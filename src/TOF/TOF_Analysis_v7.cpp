// ==========================================================================
//  TOF_Analysis.cpp — Analisi del Tempo di Volo e ricostruzione β, θ
// ==========================================================================
//
//  SCOPO:
//    A partire dai parametri di calibrazione (retta Δt₁₂ vs x e costante C)
//    e dai dati acquisiti con PMT3 in posizione Guida B (sotto la barra),
//    ricostruire evento per evento:
//      - la posizione di impatto x sulla barra
//      - il tempo di volo (TOF)
//      - la distanza percorsa l
//      - l'angolo di incidenza θ
//      - la velocità v e β = v/c
//
//  FISICA:
//    Per ogni evento con tripla coincidenza:
//
//    1) POSIZIONE DI IMPATTO sulla barra:
//       Dalla retta di calibrazione  Δt₁₂ = m·x + q   si inverte:
//         x_imp = (Δt₁₂ − q) / m
//
//    2) TEMPO DI VOLO:
//         T_meas = t₃ − (t₁ + t₂)/2
//         TOF = T_meas − C(x_imp)
//       dove C(x_imp) è la costante di offset dalla calibrazione, modellata
//       come funzione a scalini (piecewise constant).
//
//    3) DISTANZA PERCORSA (Pitagora):
//         l = √[ (x_imp − x_PMT3)² + h² ]
//       dove h = distanza verticale barra→PMT3, x_PMT3 = offset orizzontale
//       del centro di PMT3 rispetto al centro della barra.
//
//    4) ANGOLO DI INCIDENZA rispetto alla verticale:
//         θ = arctan[ (x_imp − x_PMT3) / h ]
//
//    5) VELOCITÀ e β:
//         v = l / TOF       [cm/ns]
//         β = v / c         (c = 29.9792 cm/ns)
//         1/v = TOF / l     [ns/cm]  (distribuzione più simmetrica)
//
//  UTILIZZO:
//    root -l 'TOF_Analysis.cpp("file1.xml,file2.xml,file3.xml")'
//    oppure per un singolo file:
//    root -l 'TOF_Analysis.cpp("dati_tof.xml")'
//
//  PARAMETRI DI CALIBRAZIONE:
//    I parametri m, q e i valori C(x_k) devono essere impostati nella
//    SEZIONE 1 del codice a partire dai risultati di TOF_Calibration.cpp.
//    In alternativa, possono essere caricati dal file ROOT di calibrazione
//    con la funzione LoadCalibrationFromFile().
//
//  DIPENDENZE: ROOT 6+ (TTree, TH1D, TF1, TCanvas, TGraphErrors)
//
//  AUTORI: Luca (sviluppo) + Claude (assistenza)
//  DATA:   Aprile 2026
// ==========================================================================

// ==========================================================================
//  TOF_Analysis_lead — Analisi TOF con filtro di coincidenza FPGA (FIFO)
// ==========================================================================
//
//  SCOPO:
//    Estende TOF_Analysis() con un filtraggio degli eventi DRS basato su
//    una scheda FPGA che registra segnali di due canali con relativi
//    timestamp. Solo gli eventi DRS che presentano una coincidenza nel
//    canale 1 (entro ±200 cicli di clock = ±1 µs) seguita da un segnale
//    nel canale 2 (entro 20–4000 cicli di clock) vengono analizzati.
//
//  FORMATO FILE FIFO:
//    Ogni riga contiene due colonne:
//      colonna 1: tipo evento (1 = canale 1, 2 = canale 2, 2^31 = reset buffer)
//      colonna 2: timestamp (0 … 2^30–1) oppure numero reset (2^31, 2^31+1, …)
//    I timestamp si azzerano ogni 5.36870912 s (= 2^30 × 5 ns).
//    Il primo reset valido è quello con colonna 2 = 2^31 (il secondo sarà
//    2^31+1, ecc.). I righe prima del primo reset vanno scartate.
//
//  NOME FILE FIFO:
//    Il nome del file contiene l'orario di inizio acquisizione FPGA nel
//    formato "FIFOread_YYYYMMDD-HHMMSS.txt". Questo orario è 2 ore
//    INDIETRO rispetto ai timestamp degli eventi DRS.
//
//  ALGORITMO DI FILTRAGGIO (per ogni evento DRS):
//    1. Legge l'orario DRS (campo timestamp XML: "YYYY/MM/DD HH:MM:SS.sss")
//    2. Sottrae 2 ore per ricondurlo al riferimento temporale FPGA
//    3. Calcola il tempo trascorso dall'inizio acquisizione FPGA (in secondi)
//    4. Divide per 5.36870912 → numero di buffer da saltare (troncato a intero)
//    5. Scorre i reset fino a raggiungere il buffer corrispondente
//    6. Sottrae il tempo del reset → tempo residuo in secondi
//    7. Converte in cicli di clock (÷ 5 ns) → timestamp atteso nel FIFO
//    8. Cerca nel vettore FIFO un evento di canale 1 entro ±200 clock
//    9. Se trovato, cerca un evento di canale 2 con Δt ∈ [20, 4000] clock
//   10. Solo se entrambe le condizioni sono soddisfatte, l'evento DRS è valido
//
//  UTILIZZO:
//    root -l 'TOF_Analysis.cpp("file.xml", "fifo.txt", "output.root")'
//    root -l 'TOF_Analysis.cpp("f1.xml,f2.xml", "fifo.txt", "out.root", "cal.root")'
//
//  NOTA: questa funzione esegue internamente l'intera analisi TOF_Analysis()
//  (stessi grafici, stesse quantità fisiche, stessi tagli di qualità) sugli
//  eventi filtrati. L'output ROOT ha lo stesso layout di TOF_Analysis().
// ==========================================================================


#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TROOT.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TAxis.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>


// ==========================================================================
//  SEZIONE 1: PARAMETRI CONFIGURABILI
// ==========================================================================

// --- Costanti fisiche ---
const double C_LIGHT = 29.9792458;    // Velocità della luce [cm/ns]

// --- Parametri hardware DRS4 (identici alla calibrazione) ---
const int    MAX_SAMPLES  = 1024;
const int    MAX_CHANNELS = 4;
const int    NBL_SAMPLES  = 50;
const double CFD_FRACTION = 0.20;
const double NOISE_THRESH = 5.0;

// --- Soglie di qualità per scartare eventi non fisici ---
// CLIP_V_LO: se V_min scende sotto questa soglia [mV], il segnale è saturato.
//   La DRS4 satura a circa ±500 mV. Con margine: −450 mV.
const double CLIP_V_LO = -499.0;

// CLIP_V_HI: se V_max supera questa soglia [mV], il segnale ha una componente
//   positiva anomala (oscillazione, ringing). Un impulso fisico negativo non
//   dovrebbe mai superare significativamente la baseline (~0 mV).
const double CLIP_V_HI = 100.0;

// ============================================================
//  RILEVAMENTO OSCILLAZIONE: discriminatore a isteresi (Schmitt trigger)
// ============================================================
// Contiamo il numero di crossing in SALITA sopra OSC_V_HIGH dopo il
// picco negativo, con isteresi su OSC_V_LOW per ignorare ringing fine.
//   1 crossing = singolo overshoot post-clip → ammesso (fisiologico)
//   2 crossing = doppio rimbalzo → ammesso
//   3+ crossing = oscillazione vera → scartato
const double OSC_V_HIGH     = 50.0;   // [mV] soglia alta crossing
const double OSC_V_LOW      = 25.0;   // [mV] soglia bassa isteresi
const int    OSC_NCROSS_MAX = 3;      // crossing massimi ammessi

// X_CUT_LO, X_CUT_HI: limiti fisici STRETTI della barra (-140, +140 cm).
//   La barra BC408 è lunga 280 cm e nel sistema centrato si estende da
//   -140 a +140. Un muone non può fisicamente aver attraversato la barra
//   fuori da questo intervallo, quindi un x_imp ricostruito esterno
//   indica che il modello di calibrazione (Δt12 = m·x + q) NON è
//   applicabile a quell'evento (CFD agganciato al rumore, ringing,
//   saturazione sopravvissuta al cut, ecc.).
//
//   ATTENZIONE: il taglio rigido a ±140 può tagliare anche eventi FISICI
//   ai bordi a causa della risoluzione finita su x (~3-5 cm con CFD).
//   Se serve recuperare l'efficienza ai bordi, allargare leggermente il
//   margine (es. ±145) o tenere come due cut separati: hard a ±150 e
//   "stretto" a ±140 con un flag aggiuntivo.
const double X_CUT_LO = -145.0;     // [cm]  bordo sinistro fisico
const double X_CUT_HI =  145.0;     // [cm]  bordo destro fisico

// --- Geometria del setup TOF (Guida B) ---
// h: distanza VERTICALE tra il centro della barra e il centro di PMT3 [cm]
//    Manca & Selicato usano 156 cm; il vostro setup usa ~101 cm.
//    MODIFICARE in base alla misura effettiva.
double PAR_H = 156.0;

// x_PMT3: posizione orizzontale del centro di PMT3 rispetto al centro
//         della barra [cm]. 0 = allineato al centro.
//           stimano ~139.9 cm dalla media degli x_imp, ovvero
//         circa al centro della barra nel loro riferimento (0-280).
//         Nel nostro riferimento (centro = 0): x_PMT3 ≈ 0.
//         MODIFICARE se PMT3 è disallineato.
double PAR_X_PMT3 = 0.0;

// --- Parametri di calibrazione ---
// Retta: Δt₁₂ = m · x + q   →   x = (Δt₁₂ − q) / m
// QUESTI VALORI DEVONO ESSERE AGGIORNATI DAI RISULTATI DI TOF_Calibration.cpp
double CAL_M = 0.132;     // Pendenza retta calibrazione [ns/cm]
double CAL_Q = -1.805;      // Intercetta retta calibrazione [ns]
// ============================================================
//  PARAMETRI PER IL RECUPERO SLEW RATE DEI CLIPPATI
// ============================================================
// Coerenti con TOF_Calibration_v7.cpp. Per ogni evento clippato
// si fitta linearmente il fronte di salita nella finestra
// [SLEW_U_LO, SLEW_U_HI] in valore ASSOLUTO della tensione (u = bl - v),
// si calcola lo slew rate s = du/dt e si stima l'ampiezza vera come
// A_rec = k · s, con k calibrato sui non-clippati.
//
// Selezionando soglie sotto |CLIP_V_LO| (circa 450 mV) ci si assicura
// che i campioni del fronte siano sempre sotto al plateau di clipping.

const double SLEW_U_LO = 30.0;    // [mV] limite inferiore della finestra di fit
const double SLEW_U_HI = 150.0;   // [mV] limite superiore della finestra di fit
const int    SLEW_NMIN = 4;       // Minimo numero di campioni per accettare il fit
const double SLEW_CHI2_MAX    = 1e7;  // χ²/ndf massimo accettabile
const double SLEW_SIGMA_FLOOR = 2.0;  // [mV] errore minimo per il χ²
const double SLEW_CLIP_MARGIN = 5.0;  // [mV] margine sopra |CLIP_V_LO| per la soglia CFD ricostruita

// Costanti k = A/s [ns] caricate dalla calibrazione: gK_PMT[0]=PMT1, [1]=PMT2.
// PMT3 non clippa (ampiezza moderata) e non viene calibrato → resta a 0.
static double gK_PMT[MAX_CHANNELS]      = {0.0, 0.0, 0.0, 0.0};
static double gK_PMT_err[MAX_CHANNELS]  = {0.0, 0.0, 0.0, 0.0};
static double gQ_PMT[MAX_CHANNELS]      = {0.0, 0.0, 0.0, 0.0};   // intercetta fit A = k·s + q [mV]
static double gQ_PMT_err[MAX_CHANNELS]  = {0.0, 0.0, 0.0, 0.0};
static bool   gK_calibrated             = false;


// ============================================================
//  PARAMETRI PER IL RECUPERO TIME-OVER-THRESHOLD (TOT)
// ============================================================
//
// Modello fisico: per impulso a forma fissa con coda esponenziale,
//   TOT(q) ≈ τ_eff · ln(A/q)   →   A = p0 · exp(p1 · TOT)
// con p0 ≈ q (soglia di riferimento) e p1 = 1/τ_eff.
//
// La calibrazione viene effettuata in TOF_Calibration sui non-clippati,
// e i parametri vengono caricati qui via LoadCalibrationFromFile().
// Per il recupero dei clippati, calcoliamo TOT alla soglia di
// riferimento TOT_THR_REF e applichiamo A_TOT = p0·exp(p1·TOT).

// Soglia singola di riferimento per il calcolo del TOT [mV].
// Default = 100 mV (5σ_baseline tipico, ben sotto |CLIP_V_LO|).
// Il valore esatto viene letto dal file di calibrazione tramite
// LoadCalibrationFromFile(); questo è solo un fallback.
double TOT_THR_REF = 100.0;

// Margine sopra |CLIP_V_LO| per evitare crossing nella regione
// distorta dal plateau saturato.
const double TOT_CLIP_MARGIN  = 20.0;   // [mV]

// Numero minimo di campioni tra rise e fall per accettare TOT
const int    TOT_NMIN_SAMPLES = 3;

// Cap di estrapolazione: se A_rec > TOT_AREC_MAX_FACTOR × A_max_calib,
// non si fa il recupero (estrapolazione esponenziale troppo aggressiva).
const double TOT_AREC_MAX_FACTOR = 3.0;

// Globali popolate da LoadCalibrationFromFile():
//   A_TOT(TOT) = gTOT_p0_PMT[k] · exp(gTOT_p1_PMT[k] · TOT)
// gTOT_Amax_PMT[k] è l'ampiezza massima del campione di calibrazione,
// usata come riferimento per il cap di estrapolazione.
static double gTOT_p0_PMT[MAX_CHANNELS]   = {0.0, 0.0, 0.0, 0.0};
static double gTOT_p1_PMT[MAX_CHANNELS]   = {0.0, 0.0, 0.0, 0.0};
static double gTOT_Amax_PMT[MAX_CHANNELS] = {0.0, 0.0, 0.0, 0.0};
static bool   gTOT_calibrated             = false;
// --- C(x_k): costante di offset per il TOF ---
// Coppie {posizione_x [cm], valore_C [ns]} dalle calibrazioni.
// La funzione GetC(x) restituisce il valore C del punto di calibrazione
// più vicino a x (funzione a scalini / piecewise constant).
// AGGIORNARE con i risultati del fit di C per ogni posizione.
struct CPoint {
    double x;    // Posizione di calibrazione [cm]
    double C;    // Valore di C = t₃ − (t₁+t₂)/2 misurato [ns]
};
// Modello lineare C(x) = C_lin_p0 + C_lin_p1·x, caricato da fit_params.
// Se gC_lin_loaded è true, GetC() usa questo modello al posto
// dell'interpolazione punto-a-punto (più robusto, meno parametri).
static double gC_lin_p0     = 0.0;
static double gC_lin_p0_err = 0.0;
static double gC_lin_p1     = 0.0;
static double gC_lin_p1_err = 0.0;
static bool   gC_lin_loaded = false;
std::vector<CPoint> CAL_C_POINTS = {
    {-130.0, -13.10},   // ← SOSTITUIRE con i valori misurati
    {-112.0, -13.10},
    { -84.0, -13.10},
    { -56.0, -13.10},
    { -28.0, -13.10},
    {   0.0, -13.10},
    {  28.0, -13.10},
    {  56.0, -13.10},
    {  84.0, -13.10},
    { 112.0, -13.10},
    { 130.0, -13.10}
};

// --- Parametri di binning istogrammi TOF ---
const int    NBINS_TOF   = 150;
const double TOF_LO      = -5.0;     // [ns]
const double TOF_HI      = 25.0;     // [ns]
const int    NBINS_BETA  = 200;
const double BETA_LO     = 0.0;
const double BETA_HI     = 2.5;
const int    NBINS_INVV  = 200;
const double INVV_LO     = 0.0;
const double INVV_HI     = 0.15;     // [ns/cm]
const int    NBINS_THETA = 100;
const double THETA_LO    = -1.2;     // [rad]
const double THETA_HI    = 1.2;      // [rad]
const int    NBINS_X     = 140;
const double X_LO        = -160.0;   // [cm]
const double X_HI        = 160.0;    // [cm]

// --------------------------------------------------------------------------
//  COSTANTI SPECIFICHE DEL FILTRO LEAD
// --------------------------------------------------------------------------

// Periodo di un singolo ciclo di clock FPGA [secondi]
static const double FIFO_CLK_PERIOD_S = 5.0e-9;   // 5 ns

// Durata di un buffer (2^30 cicli di clock)
static const double FIFO_BUF_DURATION_S = 5.36870912;  // = 2^30 × 5 ns

// Valore sentinella che identifica una riga di reset nel file FIFO
static const long long FIFO_RESET_SENTINEL = (1LL << 31);  // 2147483648

// Tolleranza ±ΔTs per la ricerca del canale 1 [cicli di clock]
// ±1 µs = ±1000 ns / 5 ns = ±200 cicli
static const long long FIFO_COINC_WINDOW = 200LL;

// Finestra temporale per il canale 2 dopo il canale 1 [cicli di clock]
static const long long FIFO_CH2_MIN_DT   = 20LL;
static const long long FIFO_CH2_MAX_DT   = 4000LL;

// Offset di fuso orario: i timestamp DRS sono avanti di 2 ore rispetto
// all'orologio FPGA (scritta nel nome del file FIFO)
static const int FIFO_TZ_OFFSET_S = 2 * 3600;   // 7200 s



// ==========================================================================
//  SEZIONE 2: STRUTTURE DATI (identiche alla calibrazione)
// ==========================================================================

struct ChannelData {
    int    nsamples;
    float  time[MAX_SAMPLES];
    float  voltage[MAX_SAMPLES];
    double baseline, baseline_rms, amplitude, v_min, t_min, t_cfd;
    double v_max;          // Tensione massima [mV] (per rilevare oscillazioni bipolari)
    bool   has_pulse, cfd_ok;
    bool   is_clipped;     // true se V fuori range DRS4 (saturazione)
    bool   is_oscillating; // true se troppi campioni sopra soglia (non impulso fisico)

    // ---- PASSO 7: SLEW RATE FIT (per recupero clippati) ----
    // Risultati del fit lineare del fronte di salita nella finestra
    // [SLEW_U_LO, SLEW_U_HI]. Calcolati per OGNI evento (clippato o no);
    // i non clippati servono solo per check, i clippati per la ricostruzione.
    double slew_rate;       // Pendenza |du/dt| del fronte [mV/ns]
    double slew_rate_err;   // Errore sulla pendenza
    double slew_chi2_ndf;   // χ²/ndf del fit
    int    slew_n_pts;      // Numero di campioni nella finestra
    bool   slew_ok;         // true se fit valido

// ---- Recupero CFD via slew rate ----
    // Popolati da RecoverClippedCFD(): A_rec = k · slew_rate;
    // soglia CFD ricostruita = baseline − f · A_rec; se cade sopra il
    // clipping si trova un crossing valido sul fronte.
    double amplitude_rec;     // Ampiezza ricostruita slew [mV] (0 se non recup.)
    double t_cfd_recovered;   // Tempo CFD ricostruito slew [ns]
    bool   cfd_recovered;     // true se il recupero slew ha avuto successo

    // ---- Time-Over-Threshold alla soglia di riferimento ----
    // Calcolato in AnalyzeChannel (PASSO 8). Se il segnale supera la
    // soglia TOT_THR_REF, si misurano i due crossing (rise e fall)
    // tramite interpolazione lineare e si calcola la durata.
    double tot_q;             // TOT alla soglia TOT_THR_REF [ns]
    double tot_rise;          // Tempo crossing salita [ns]
    double tot_fall;          // Tempo crossing discesa [ns]
    bool   tot_ok;            // true se TOT misurato correttamente

    // ---- Recupero CFD via TOT ----
    // Popolati da RecoverClippedCFD_TOT(): A_TOT = p0·exp(p1·TOT);
    // tempo CFD su soglia f·A_TOT, analogo al recupero slew.
    double amplitude_tot;     // Ampiezza ricostruita TOT [mV]
    double t_cfd_rec_tot;     // Tempo CFD ricostruito TOT [ns]
    bool   cfd_recovered_tot; // true se il recupero TOT è valido
};

struct EventData {
    int         serial;
    std::string timestamp;
    int         trigger_cell;
    int         board_serial;
    int         scaler[MAX_CHANNELS];
    int         nchannels;
    int         channel_ids[MAX_CHANNELS];
    ChannelData ch[MAX_CHANNELS];
};


// ==========================================================================
//  SEZIONE 3: ANALISI FORMA D'ONDA (identica alla calibrazione)
// ==========================================================================

void AnalyzeChannel(ChannelData &cd) {
    cd.has_pulse = false;  cd.cfd_ok = false;
    cd.is_clipped = false; cd.is_oscillating = false;
    cd.baseline = 0; cd.baseline_rms = 0; cd.amplitude = 0;
    cd.v_min = 0; cd.v_max = 0; cd.t_min = 0; cd.t_cfd = -999.0;

    // Inizializzazione campi slew rate (PASSO 7)
    cd.slew_rate       = 0.0;
    cd.slew_rate_err   = 0.0;
    cd.slew_chi2_ndf   = -1.0;
    cd.slew_n_pts      = 0;
    cd.slew_ok         = false;
    cd.amplitude_rec   = 0.0;
    cd.t_cfd_recovered = -999.0;
    cd.cfd_recovered   = false;

    // Inizializzazione campi TOT (PASSO 8)
    cd.tot_q             = 0.0;
    cd.tot_rise          = 0.0;
    cd.tot_fall          = 0.0;
    cd.tot_ok            = false;
    cd.amplitude_tot     = 0.0;
    cd.t_cfd_rec_tot     = -999.0;
    cd.cfd_recovered_tot = false;

    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return;
    float *t = cd.time, *v = cd.voltage;

    // Baseline
    double sum = 0, sum2 = 0;
    for (int i = 0; i < NBL_SAMPLES; i++) { sum += v[i]; sum2 += (double)v[i]*v[i]; }
    double bl = sum / NBL_SAMPLES;
    double bl_rms = sqrt(fabs(sum2/NBL_SAMPLES - bl*bl));
    cd.baseline = bl;  cd.baseline_rms = bl_rms;

    // Minimo e massimo globali
    double vmin = v[0], vmax = v[0];
    int imin = 0;
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
        if (v[i] > vmax)   vmax = v[i];
    }
    cd.v_min = vmin;  cd.v_max = vmax;  cd.t_min = t[imin];

    // ---- CONTROLLO CLIPPING ----
    // Un segnale è clippato se la tensione raggiunge i limiti della DRS4.
    // Tipico dell'evento 13/4341: V oscilla tra −500 e +500 mV.
    if (vmin < CLIP_V_LO || vmax > CLIP_V_HI) {
        cd.is_clipped = true;
    }

    // Ampiezza
    double amp = bl - vmin;
    cd.amplitude = amp;
    if (amp < NOISE_THRESH) return;
    cd.has_pulse = true;

// ---- CONTROLLO OSCILLAZIONE (Schmitt trigger) ----
    // Dopo il picco negativo, contiamo quante volte il segnale sale
    // sopra OSC_V_HIGH con isteresi su OSC_V_LOW. Un singolo overshoot
    // post-clip conta 1 crossing (fisiologico); una vera oscillazione
    // ne genera 3+.
    {
        int n_cross = 0;
        bool above = false;
        for (int i = imin + 1; i < ns; i++) {
            double vi = (double)v[i];
            if (!above) {
                if (vi > OSC_V_HIGH) { n_cross++; above = true; }
            } else {
                if (vi < OSC_V_LOW) { above = false; }
            }
        }
        if (n_cross >= OSC_NCROSS_MAX) {
            cd.is_oscillating = true;
        }
    }
// CFD forward (calcolato anche per eventi problematici — il flag
    // di qualità viene controllato nel loop principale, non qui)
    double v_thr = bl + CFD_FRACTION * (vmin - bl);
    for (int i = NBL_SAMPLES; i < imin; i++) {
        if (v[i] > v_thr && v[i+1] <= v_thr) {
            double dv = (double)v[i+1] - v[i];
            if (fabs(dv) > 1e-6) {
                cd.t_cfd = t[i] + (v_thr - v[i]) / dv * (t[i+1] - t[i]);
                cd.cfd_ok = true;
            }
            break;
        }
    }

    // ---- PASSO 7: FIT LINEARE DELLO SLEW RATE ----
    //
    // METODO (identico a TOF_Calibration_v7):
    //   1. Definiamo u(t) = baseline − v(t) (positiva nel fronte).
    //   2. Selezioniamo i campioni con SLEW_U_LO ≤ u(t) ≤ SLEW_U_HI
    //      sul fronte di salita (i ∈ [NBL_SAMPLES, imin]).
    //   3. Se ne abbiamo almeno SLEW_NMIN, fittiamo u = a + b·t con
    //      regressione lineare (errore σ_y = max(bl_rms, SLEW_SIGMA_FLOOR)).
    //   4. La pendenza b = du/dt è lo slew rate (positivo per impulsi
    //      negativi, dato che u cresce mentre v scende).
    //   5. Validazione: pendenza > 0, χ²/ndf finito e < SLEW_CHI2_MAX.
    //
    // CASI LIMITE:
    //   - Se imin è troppo vicino a NBL_SAMPLES, pochi campioni sul fronte.
    //   - Se l'ampiezza è piccola (< SLEW_U_HI), la finestra si riduce.
    //   - Per CLIPPATI, il fronte è ancora pulito sotto al plateau quindi
    //     il fit è valido — è proprio il caso che vogliamo gestire.
    {
        // Raccogli i campioni del fronte di salita nella finestra di soglie
        std::vector<double> tt, uu;
        tt.reserve(64); uu.reserve(64);

        for (int i = NBL_SAMPLES; i <= imin; i++) {
            double u_i = bl - v[i];                // tensione "positivizzata"
            if (u_i >= SLEW_U_LO && u_i <= SLEW_U_HI) {
                tt.push_back(t[i]);
                uu.push_back(u_i);
            }
        }

        cd.slew_n_pts = (int)tt.size();

        // Procedi solo se abbiamo abbastanza punti per un fit significativo
        if (cd.slew_n_pts >= SLEW_NMIN) {

            // Errore tipico per ogni campione (dominato dal rumore di
            // baseline; floor a SLEW_SIGMA_FLOOR per evitare divisioni
            // patologiche se bl_rms è anomalo)
            double sigma_y = std::max(bl_rms, SLEW_SIGMA_FLOOR);
            double w = 1.0 / (sigma_y * sigma_y);

            // Regressione lineare standard: u = a + b·t
            // Con errori uniformi (w costante), le formule chiuse sono:
            //   b = (N·Σ(t·u) − Σt·Σu) / (N·Σt² − (Σt)²)
            //   a = (Σu − b·Σt) / N
            double N = (double)cd.slew_n_pts;
            double St = 0.0, Su = 0.0, Stt = 0.0, Stu = 0.0;
            for (int i = 0; i < cd.slew_n_pts; i++) {
                St  += tt[i];
                Su  += uu[i];
                Stt += tt[i] * tt[i];
                Stu += tt[i] * uu[i];
            }
            double D = N * Stt - St * St;
            if (fabs(D) > 1e-12) {
                double slope     = (N * Stu - St * Su) / D;
                double intercept = (Su - slope * St) / N;

                // Errori sui parametri (dalla matrice di covarianza)
                double slope_err = sqrt(N / D) * sigma_y;

                // χ² con la stessa σ_y per tutti i punti
                double chi2 = 0.0;
                for (int i = 0; i < cd.slew_n_pts; i++) {
                    double resid = uu[i] - (intercept + slope * tt[i]);
                    chi2 += resid * resid * w;
                }
                int ndf = cd.slew_n_pts - 2;
                double chi2ndf = (ndf > 0) ? chi2 / ndf : -1.0;

                cd.slew_rate     = slope;
                cd.slew_rate_err = slope_err;
                cd.slew_chi2_ndf = chi2ndf;

if (slope > 0.0 && chi2ndf > 0.0 && chi2ndf < SLEW_CHI2_MAX) {
                    cd.slew_ok = true;
                }
            }
        }
    }

    // ---- PASSO 8: TIME-OVER-THRESHOLD A SOGLIA SINGOLA ----
    //
    // Misuriamo il TOT alla soglia di riferimento TOT_THR_REF (default
    // 100 mV, leggibile dal file di calibrazione). Algoritmo:
    //   1. Cerco il PRIMO crossing in salita (rise) tra NBL_SAMPLES e imin
    //   2. Cerco il PRIMO crossing in discesa (fall) DOPO imin
    //   3. Interpolazione lineare per timing sub-campione
    //   4. TOT = t_fall − t_rise
    //
    // VINCOLI DI VALIDITÀ:
    //   - q ≤ |CLIP_V_LO| − TOT_CLIP_MARGIN (la soglia non deve cadere
    //     vicino al plateau saturato, dove il fall sarebbe distorto)
    //   - amp ≥ q (segnale supera la soglia)
    //   - rise e fall trovati entrambi
    //   - durata ≥ TOT_NMIN_SAMPLES × periodo campionamento
    //
    // Per gli eventi clippati, il rise e il fall sono entrambi nelle
    // regioni del segnale fuori dal plateau (sotto soglia), quindi il
    // TOT è ben definito anche per i clippati. Questo è il punto chiave
    // che permette di usare il TOT come variabile di ricostruzione.
    {
        double q = TOT_THR_REF;
        double q_max_safe = fabs(CLIP_V_LO) - TOT_CLIP_MARGIN;

        // Stima del periodo medio di campionamento (per il check sui
        // sample minimi). Usiamo la differenza tra il primo e l'ultimo
        // sample diviso (ns-1).
        double sample_period = (ns > 1) ?
                               ((double)t[ns - 1] - (double)t[0]) / (ns - 1) :
                               0.2;  // fallback 0.2 ns = 5 GS/s
        double dt_min_tot = TOT_NMIN_SAMPLES * sample_period;

        // Procedi solo se la soglia è sicura e il segnale la raggiunge
        if (q <= q_max_safe && amp >= q) {

            // Soglia in tensione (negativa, perché impulso negativo)
            double v_thr_tot = bl - q;

            // ---- Crossing sul fronte di salita ----
            // Cerca [i, i+1] tale che v[i] > v_thr_tot ≥ v[i+1]
            // (equivalente: u[i] < q ≤ u[i+1] con u = bl - v)
            double t_rise = -1.0;
            bool   rise_found = false;
            for (int i = NBL_SAMPLES; i < imin; i++) {
                if (v[i] > v_thr_tot && v[i + 1] <= v_thr_tot) {
                    double dv = (double)v[i + 1] - v[i];
                    if (fabs(dv) > 1e-6) {
                        t_rise = t[i] + (v_thr_tot - v[i]) / dv *
                                 (t[i + 1] - t[i]);
                        rise_found = true;
                    }
                    break;
                }
            }

            // ---- Crossing sul fronte di discesa ----
            // Cerca [i, i+1] DOPO imin tale che v[i] ≤ v_thr_tot < v[i+1]
            double t_fall = -1.0;
            bool   fall_found = false;
            if (rise_found) {
                for (int i = imin; i < ns - 1; i++) {
                    if (v[i] <= v_thr_tot && v[i + 1] > v_thr_tot) {
                        double dv = (double)v[i + 1] - v[i];
                        if (fabs(dv) > 1e-6) {
                            t_fall = t[i] + (v_thr_tot - v[i]) / dv *
                                     (t[i + 1] - t[i]);
                            fall_found = true;
                        }
                        break;
                    }
                }
            }

            // ---- Validazione finale ----
            if (rise_found && fall_found) {
                double tot_val = t_fall - t_rise;
                if (tot_val >= dt_min_tot) {
                    cd.tot_rise = t_rise;
                    cd.tot_fall = t_fall;
                    cd.tot_q    = tot_val;
                    cd.tot_ok   = true;
                }
            }
        }
    }
}


// ==========================================================================
//  RECUPERO DEGLI EVENTI CLIPPATI VIA SLEW RATE
// ==========================================================================

/// RecoverClippedCFD():
///   Per un evento clippato, stima l'ampiezza vera come A_rec = k·slew_rate
///   (k calibrato globalmente per il canale), e ricalcola il tempo CFD
///   con soglia f·A_rec. Se la soglia ricostruita cade sopra il livello
///   di clipping (con margine SLEW_CLIP_MARGIN), il crossing esiste sulla
///   waveform ed è recuperabile.
///
/// PRECONDIZIONI:
///   - cd.is_clipped == true (altrimenti niente da recuperare)
///   - cd.slew_ok == true (fit dello slew rate riuscito)
///   - gK_calibrated == true e gK_PMT[ch_index] > 0 (k disponibile)
///
/// EFFETTI: popola cd.amplitude_rec, cd.t_cfd_recovered, cd.cfd_recovered.

void RecoverClippedCFD(ChannelData &cd, int ch_index) {
    if (!cd.is_clipped)                            return;
    if (!cd.slew_ok)                               return;
    if (ch_index < 0 || ch_index >= MAX_CHANNELS)  return;
    if (!gK_calibrated)                            return;
    if (gK_PMT[ch_index] <= 0.0)                   return;

// Stima ampiezza vera col modello completo: A = k·s + q
    double A_rec = gK_PMT[ch_index] * cd.slew_rate + gQ_PMT[ch_index];
    cd.amplitude_rec = A_rec;

    // Sanità: A_rec deve essere fisicamente ragionevole
    if (A_rec <= 0.0) return;

    // Soglia CFD ricostruita: deve cadere SOPRA il clipping con margine
    double v_thr = cd.baseline - CFD_FRACTION * A_rec;
    if (v_thr <= CLIP_V_LO + SLEW_CLIP_MARGIN) return;

    // Cerca il crossing sul fronte di salita
    int   ns = cd.nsamples;
    float *t = cd.time;
    float *v = cd.voltage;

    // Indice del minimo (primo campione del plateau saturato)
    double vmin = v[0];
    int    imin = 0;
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
    }

    for (int i = NBL_SAMPLES; i < imin; i++) {
        if (v[i] > v_thr && v[i + 1] <= v_thr) {
            double dv = (double)v[i + 1] - v[i];
            if (fabs(dv) > 1e-6) {
                cd.t_cfd_recovered = t[i] + (v_thr - v[i]) / dv *
                                     (t[i + 1] - t[i]);
                cd.cfd_recovered   = true;
            }
break;
        }
    }
}


// ==========================================================================
//  RECUPERO DEGLI EVENTI CLIPPATI VIA TOT
// ==========================================================================

/// RecoverClippedCFD_TOT():
///   Variante del recupero clippati basata sul Time-Over-Threshold alla
///   soglia di riferimento TOT_THR_REF (caricata dal file di calibrazione).
///
/// MODELLO: A_TOT = p0 · exp(p1 · TOT)
///   con (p0, p1) calibrati per ciascun PMT (gTOT_p0_PMT, gTOT_p1_PMT)
///   sui non-clippati e applicati per estrapolazione ai clippati.
///
/// PRECONDIZIONI:
///   - cd.is_clipped == true
///   - cd.tot_ok == true (TOT misurato e valido)
///   - gTOT_calibrated == true e gTOT_p1_PMT[ch_index] > 0
///
/// CAP DI ESTRAPOLAZIONE:
///   Se A_TOT supera TOT_AREC_MAX_FACTOR × A_max_calibrazione, la stima
///   viene rigettata (estrapolazione esponenziale troppo aggressiva).
///
/// VERIFICA SOGLIA CFD: identica a RecoverClippedCFD: la soglia f·A_TOT
///   deve cadere sopra il livello di clipping con margine SLEW_CLIP_MARGIN,
///   altrimenti il crossing non esiste sulla waveform.
///
/// EFFETTI: popola cd.amplitude_tot, cd.t_cfd_rec_tot, cd.cfd_recovered_tot.

void RecoverClippedCFD_TOT(ChannelData &cd, int ch_index) {

    // --- Sanity check delle precondizioni ---
    if (!cd.is_clipped)                            return;
    if (ch_index < 0 || ch_index >= MAX_CHANNELS)  return;
    if (!gTOT_calibrated)                          return;
    if (gTOT_p1_PMT[ch_index] <= 0.0)              return;
    if (!cd.tot_ok)                                return;

    // --- Stima dell'ampiezza ricostruita ---
    double p0 = gTOT_p0_PMT[ch_index];
    double p1 = gTOT_p1_PMT[ch_index];
    double A_tot = p0 * exp(p1 * cd.tot_q);
    cd.amplitude_tot = A_tot;

    // Sanità: ampiezza positiva
    if (A_tot <= 0.0) return;

    // Cap di estrapolazione: protegge contro TOT anomalamente lunghi
    // (afterpulse, pile-up) che farebbero esplodere l'esponenziale.
    double A_cap = TOT_AREC_MAX_FACTOR * gTOT_Amax_PMT[ch_index];
    if (A_cap > 0.0 && A_tot > A_cap) return;

    // --- Calcolo della soglia CFD ricostruita ---
    double v_thr = cd.baseline - CFD_FRACTION * A_tot;

    // Verifica CRITICA: la soglia deve stare sopra il clipping con margine
    if (v_thr <= CLIP_V_LO + SLEW_CLIP_MARGIN) return;

    // --- Ricerca del crossing sul fronte di salita ---
    int   ns = cd.nsamples;
    float *t = cd.time;
    float *v = cd.voltage;

    // Indice del minimo (primo campione del plateau saturato)
    double vmin = v[0];
    int    imin = 0;
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
    }

    for (int i = NBL_SAMPLES; i < imin; i++) {
        if (v[i] > v_thr && v[i + 1] <= v_thr) {
            double dv = (double)v[i + 1] - v[i];
            if (fabs(dv) > 1e-6) {
                cd.t_cfd_rec_tot     = t[i] + (v_thr - v[i]) / dv *
                                       (t[i + 1] - t[i]);
                cd.cfd_recovered_tot = true;
            }
            break;
        }
    }
}


// ==========================================================================
//  SEZIONE 4: PARSING XML (identico alla calibrazione)
// ==========================================================================

int ParseXML(const char* filename, std::vector<EventData> &events) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "[ERRORE] Impossibile aprire: " << filename << std::endl;
        return -1;
    }
    std::cout << "[INFO] Parsing " << filename << " ..." << std::flush;

    enum State { IDLE, IN_EVENT, IN_BOARD, IN_CHANNEL };
    State state = IDLE;
    EventData current_evt;
    int ch_idx = -1, sample_count = 0, scaler_idx = 0;
    std::string line;
    char buf[512];

    while (std::getline(infile, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        strncpy(buf, line.c_str(), sizeof(buf)-1);
        buf[sizeof(buf)-1] = '\0';

        if (state == IDLE) {
            if (line.find("<Event>") != std::string::npos) {
                current_evt = EventData();
                current_evt.nchannels = 0;
                memset(current_evt.scaler, 0, sizeof(current_evt.scaler));
                scaler_idx = 0;
                state = IN_EVENT;
            }
            continue;
        }
        if (state == IN_EVENT) {
            int val; char ts[64];
            if (sscanf(buf, " <Serial>%d</Serial>", &val) == 1)
                current_evt.serial = val;
            else if (sscanf(buf, " <Time>%63[^<]</Time>", ts) == 1)
                current_evt.timestamp = ts;
            else if (sscanf(buf, " <Board_%d>", &val) == 1) {
                current_evt.board_serial = val; state = IN_BOARD;
            }
            else if (line.find("</Event>") != std::string::npos) {
                for (int i = 0; i < current_evt.nchannels; i++)
                    AnalyzeChannel(current_evt.ch[i]);
                events.push_back(current_evt);
                state = IDLE;
            }
            continue;
        }
        if (state == IN_BOARD) {
            int val, val2;
            if (sscanf(buf, " <Trigger_Cell>%d</Trigger_Cell>", &val) == 1)
                current_evt.trigger_cell = val;
            else if (sscanf(buf, " <Scaler%d>%d</Scaler", &val, &val2) == 2) {
                if (scaler_idx < MAX_CHANNELS) current_evt.scaler[scaler_idx] = val2;
                scaler_idx++;
            }
            else if (sscanf(buf, " <CHN%d>", &val) == 1) {
                int idx = current_evt.nchannels;
                if (idx < MAX_CHANNELS) {
                    current_evt.channel_ids[idx] = val;
                    current_evt.ch[idx].nsamples = 0;
                    ch_idx = idx; current_evt.nchannels++; sample_count = 0;
                    state = IN_CHANNEL;
                }
            }
            else if (line.find("</Board_") != std::string::npos) state = IN_EVENT;
            continue;
        }
        if (state == IN_CHANNEL) {
            if (line.find("</CHN") != std::string::npos) {
                if (ch_idx >= 0) current_evt.ch[ch_idx].nsamples = sample_count;
                state = IN_BOARD; continue;
            }
            float tv, vv;
            if (ch_idx >= 0 && sample_count < MAX_SAMPLES &&
                sscanf(buf, " <Data>%f,%f</Data>", &tv, &vv) == 2) {
                current_evt.ch[ch_idx].time[sample_count] = tv;
                current_evt.ch[ch_idx].voltage[sample_count] = vv;
                sample_count++;
            }
            continue;
        }
    }
    infile.close();
    std::cout << " " << events.size() << " eventi." << std::endl;
    return (int)events.size();
}


// ==========================================================================
//  SEZIONE 5: FUNZIONE C(x) — INTERPOLAZIONE A SCALINI
// ==========================================================================

/// GetC(): restituisce il valore della costante di offset C per una data
/// posizione x sulla barra, usando un'interpolazione a scalini (piecewise
/// constant / nearest neighbor).
///
/// ALGORITMO:
///   I punti di calibrazione {x_k, C_k} definiscono una partizione della
///   barra in intervalli. Ogni intervallo è centrato su un punto x_k e
///   si estende fino ai punti medi con i vicini:
///
///     [ (x_{k-1}+x_k)/2,  (x_k+x_{k+1})/2 )  →  C = C_k
///
///   Ai bordi (x < primo punto, x > ultimo punto) si usa il valore
///   del punto estremo più vicino (estrapolazione costante).
///
///   Questa scelta è fisicamente motivata: C varia poco e in modo
///   discreto tra posizioni adiacenti, senza che abbia senso interpolare
///   linearmente tra i valori.

// double GetC(double x) {

//     int n = (int)CAL_C_POINTS.size();
//     if (n == 0) {
//         std::cerr << "[ERRORE] Nessun punto di calibrazione C definito!" << std::endl;
//         return 0.0;
//     }
//     if (n == 1) return CAL_C_POINTS[0].C;

//     // Trova il punto di calibrazione più vicino a x
//     int i_nearest = 0;
//     double min_dist = fabs(x - CAL_C_POINTS[0].x);
//     for (int i = 1; i < n; i++) {
//         double dist = fabs(x - CAL_C_POINTS[i].x);
//         if (dist < min_dist) {
//             min_dist = dist;
//             i_nearest = i;
//         }
//     }
//     return CAL_C_POINTS[i_nearest].C;
// }
double GetC(double x) {

    // Se il modello lineare C(x) = C_0 + C_1·x è stato caricato
    // dal file di calibrazione, lo usiamo (fit robusto, meno parametri,
    // evita artefatti dell'interpolazione punto-a-punto).
    if (gC_lin_loaded) {
        return gC_lin_p0 + gC_lin_p1 * x;
    }

    int n = (int)CAL_C_POINTS.size();

    // Nessun punto di calibrazione → valore di default
    // (dovrebbe essere stato impostato da LoadCalibrationFromFile o
    //  manualmente, ma se non c'è nulla restituiamo 0 con un warning).
    if (n == 0) {
        std::cerr << "[ERRORE] Nessun punto di calibrazione C definito!"
                  << std::endl;
        return 0.0;
    }
    if (n == 1) return CAL_C_POINTS[0].C;

    // ---- INTERPOLAZIONE LINEARE TRA PUNTI ADIACENTI ----
    // Strategia migliorativa rispetto al nearest-neighbor: per ogni x
    // si trova l'intervallo [x_k, x_{k+1}] che lo contiene, e si
    // restituisce
    //     C(x) = C_k + (C_{k+1} - C_k) · (x - x_k) / (x_{k+1} - x_k)
    //
    // Vantaggi:
    //   - nessuna discontinuità ai bordi degli intervalli (a differenza
    //     della funzione a scalini)
    //   - usa l'informazione di entrambi i punti vicini, pesata dalla
    //     distanza, anziché ignorare quello più lontano
    //   - ipotesi più debole possibile (assenza di salti) compatibile
    //     con i dati di calibrazione
    //
    // ESTRAPOLAZIONE: per x fuori dal range dei punti di calibrazione,
    // si usa il valore del punto estremo più vicino (clamp ai bordi).
    // Questo evita estrapolazioni lineari pericolose oltre i bordi
    // della barra; per posizioni leggermente fuori range (a causa
    // della risoluzione finita su x_imp), il valore al bordo è la
    // miglior stima ragionevole.
    //
    // PRECONDIZIONE: CAL_C_POINTS deve essere ordinato per x crescente
    // (LoadCalibrationFromFile esegue std::sort, ed è responsabilità
    // di chi modifica i punti hardcoded mantenere l'ordine).

    // Clamp ai bordi: estrapolazione costante
    if (x <= CAL_C_POINTS[0].x)         return CAL_C_POINTS[0].C;
    if (x >= CAL_C_POINTS[n - 1].x)     return CAL_C_POINTS[n - 1].C;

    // Trova l'intervallo [k, k+1] che contiene x
    // (ricerca lineare; n ~ 11 punti, costo trascurabile)
    int k = 0;
    for (int i = 0; i < n - 1; i++) {
        if (x >= CAL_C_POINTS[i].x && x <= CAL_C_POINTS[i + 1].x) {
            k = i;
            break;
        }
    }

    // Interpolazione lineare
    double x_lo = CAL_C_POINTS[k].x;
    double x_hi = CAL_C_POINTS[k + 1].x;
    double C_lo = CAL_C_POINTS[k].C;
    double C_hi = CAL_C_POINTS[k + 1].C;

    double dx = x_hi - x_lo;
    if (fabs(dx) < 1e-9) return 0.5 * (C_lo + C_hi);   // safety

    double frac = (x - x_lo) / dx;
    return C_lo + frac * (C_hi - C_lo);
}

// ==========================================================================
//  SEZIONE 6: CARICAMENTO CALIBRAZIONE DA FILE ROOT (opzionale)
// ==========================================================================

/// LoadCalibrationFromFile(): legge i parametri di calibrazione dal file
/// ROOT prodotto da TOF_Calibration.cpp.
///
/// Legge il TTree "summary" per estrarre le coppie (x_k, C_k) e
/// i parametri della retta dal fit salvato.
///
/// Restituisce true se il caricamento è riuscito, false altrimenti.

bool LoadCalibrationFromFile(const char* cal_filename) {

    TFile *fcal = TFile::Open(cal_filename, "READ");
    if (!fcal || fcal->IsZombie()) {
        std::cerr << "[WARNING] Non riesco ad aprire " << cal_filename
                  << " — uso i parametri hardcoded." << std::endl;
        return false;
    }

    // --- Leggi i punti C(x_k) dal TTree summary ---
    TTree *summary = (TTree*)fcal->Get("summary");
    if (!summary) {
        std::cerr << "[WARNING] TTree 'summary' non trovato in "
                  << cal_filename << std::endl;
        fcal->Close();
        return false;
    }

    Float_t s_x, s_C_mu;
    summary->SetBranchAddress("x",    &s_x);
    summary->SetBranchAddress("C_mu", &s_C_mu);

    CAL_C_POINTS.clear();
    int nentries = (int)summary->GetEntries();
    for (int i = 0; i < nentries; i++) {
        summary->GetEntry(i);
        CAL_C_POINTS.push_back({(double)s_x, (double)s_C_mu});
    }

    // Ordina per x crescente (necessario per GetC)
    std::sort(CAL_C_POINTS.begin(), CAL_C_POINTS.end(),
              [](const CPoint &a, const CPoint &b) { return a.x < b.x; });

    std::cout << "[INFO] Caricati " << CAL_C_POINTS.size()
              << " punti C(x) dalla calibrazione." << std::endl;
    for (auto &p : CAL_C_POINTS) {
        std::cout << "       x = " << p.x << " cm, C = " << p.C << " ns" << std::endl;
    }

// --- Leggi m, q e k_PMT[] dal TTree "fit_params" ---
    // Il TTree fit_params contiene un singolo entry con i parametri
    // del fit della retta di calibrazione (m, q) e le costanti del
    // recupero clippati via slew rate (k_PMT1, k_PMT2). Se il TTree
    // non esiste (vecchia versione di TOF_Calibration), si lascia
    // l'utente impostare i parametri a mano via SetCalibration().
    TTree *fp = (TTree*)fcal->Get("fit_params");
    if (fp && fp->GetEntries() > 0) {
        Float_t fp_m, fp_m_err, fp_q, fp_q_err;
        Float_t fp_v_eff, fp_v_eff_err;
        Float_t fp_k_PMT1, fp_k_PMT1_err;
        Float_t fp_k_PMT2, fp_k_PMT2_err;
        Float_t fp_chi2ndf;

        fp->SetBranchAddress("m",          &fp_m);
        fp->SetBranchAddress("m_err",      &fp_m_err);
        fp->SetBranchAddress("q",          &fp_q);
        fp->SetBranchAddress("q_err",      &fp_q_err);
        fp->SetBranchAddress("v_eff",      &fp_v_eff);
        fp->SetBranchAddress("v_eff_err",  &fp_v_eff_err);
        fp->SetBranchAddress("k_PMT1",     &fp_k_PMT1);
        fp->SetBranchAddress("k_PMT1_err", &fp_k_PMT1_err);
        fp->SetBranchAddress("k_PMT2",     &fp_k_PMT2);
        fp->SetBranchAddress("k_PMT2_err", &fp_k_PMT2_err);
        fp->SetBranchAddress("chi2ndf",    &fp_chi2ndf);

        fp->GetEntry(0);

        // Aggiorna le globali della retta
        CAL_M = (double)fp_m;
        CAL_Q = (double)fp_q;

        // Aggiorna le globali del recupero slew rate
        gK_PMT[0]      = (double)fp_k_PMT1;
        gK_PMT_err[0]  = (double)fp_k_PMT1_err;
        gK_PMT[1]      = (double)fp_k_PMT2;
        gK_PMT_err[1]  = (double)fp_k_PMT2_err;
        gK_PMT[2]      = 0.0;     // PMT3 non clippa, niente recupero
        gK_PMT_err[2]  = 0.0;

        // Intercetta del fit A = k·s + q (branch opzionali: se non presenti
        // nel file di calibrazione, restano a 0 → comportamento retrocompatibile)
        Float_t fp_q_PMT1 = 0, fp_q_PMT1_err = 0;
        Float_t fp_q_PMT2 = 0, fp_q_PMT2_err = 0;
        if (fp->GetBranch("q_PMT1")) {
            fp->SetBranchAddress("q_PMT1",     &fp_q_PMT1);
            fp->SetBranchAddress("q_PMT1_err", &fp_q_PMT1_err);
            fp->SetBranchAddress("q_PMT2",     &fp_q_PMT2);
            fp->SetBranchAddress("q_PMT2_err", &fp_q_PMT2_err);
            fp->GetEntry(0);  // ri-leggi l'entry con i nuovi branch
        }
        gQ_PMT[0]      = (double)fp_q_PMT1;
        gQ_PMT_err[0]  = (double)fp_q_PMT1_err;
        gQ_PMT[1]      = (double)fp_q_PMT2;
        gQ_PMT_err[1]  = (double)fp_q_PMT2_err;

        // Recupero abilitato solo se entrambi i k sono sensati
        gK_calibrated  = (gK_PMT[0] > 0.0 && gK_PMT[1] > 0.0);

        std::cout << "[INFO] Parametri retta caricati dal file:" << std::endl;
        std::cout << "       m     = " << CAL_M  << " ± " << fp_m_err
                  << " ns/cm   (χ²/ndf = " << fp_chi2ndf << ")" << std::endl;
        std::cout << "       q     = " << CAL_Q  << " ± " << fp_q_err
                  << " ns" << std::endl;
        std::cout << "       v_eff = " << fp_v_eff
                  << " ± " << fp_v_eff_err << " cm/ns" << std::endl;
        std::cout << "[INFO] Costanti slew rate caricate:" << std::endl;
        std::cout << "       k_PMT1 = " << gK_PMT[0]
                  << " ± " << gK_PMT_err[0] << " ns" << std::endl;
        std::cout << "       k_PMT2 = " << gK_PMT[1]
                  << " ± " << gK_PMT_err[1] << " ns" << std::endl;
        std::cout << "       q_PMT1 = " << gQ_PMT[0]
                  << " ± " << gQ_PMT_err[0] << " mV" << std::endl;
        std::cout << "       q_PMT2 = " << gQ_PMT[1]
                  << " ± " << gQ_PMT_err[1] << " mV" << std::endl;
std::cout << "       Recupero clippati (slew): "
                  << (gK_calibrated ? "ABILITATO" : "DISABILITATO")
                  << std::endl;

        // ---- Caricamento parametri recupero TOT ----
        // Modello: A_TOT = p0 · exp(p1 · TOT), calibrato per PMT1 e PMT2.
        // I branch sono opzionali: se non presenti (file di calibrazione
        // vecchio, senza CalibrateTOT), gTOT_calibrated resta false e
        // RecoverClippedCFD_TOT() non farà nulla — retrocompatibile.
        if (fp->GetBranch("p0_TOT_PMT1")) {
            Float_t fp_p0_TOT_PMT1, fp_p0_TOT_PMT1_err;
            Float_t fp_p1_TOT_PMT1, fp_p1_TOT_PMT1_err;
            Float_t fp_p0_TOT_PMT2, fp_p0_TOT_PMT2_err;
            Float_t fp_p1_TOT_PMT2, fp_p1_TOT_PMT2_err;
            Float_t fp_Amax_TOT_PMT1, fp_Amax_TOT_PMT2;
            Float_t fp_TOT_thr_ref;

            fp->SetBranchAddress("p0_TOT_PMT1",     &fp_p0_TOT_PMT1);
            fp->SetBranchAddress("p0_TOT_PMT1_err", &fp_p0_TOT_PMT1_err);
            fp->SetBranchAddress("p1_TOT_PMT1",     &fp_p1_TOT_PMT1);
            fp->SetBranchAddress("p1_TOT_PMT1_err", &fp_p1_TOT_PMT1_err);
            fp->SetBranchAddress("p0_TOT_PMT2",     &fp_p0_TOT_PMT2);
            fp->SetBranchAddress("p0_TOT_PMT2_err", &fp_p0_TOT_PMT2_err);
            fp->SetBranchAddress("p1_TOT_PMT2",     &fp_p1_TOT_PMT2);
            fp->SetBranchAddress("p1_TOT_PMT2_err", &fp_p1_TOT_PMT2_err);
            fp->SetBranchAddress("Amax_TOT_PMT1",   &fp_Amax_TOT_PMT1);
            fp->SetBranchAddress("Amax_TOT_PMT2",   &fp_Amax_TOT_PMT2);
            fp->SetBranchAddress("TOT_thr_ref",      &fp_TOT_thr_ref);

            fp->GetEntry(0);   // ri-leggi con i nuovi branch

            // Popola le globali TOT
            gTOT_p0_PMT[0]   = (double)fp_p0_TOT_PMT1;
            gTOT_p1_PMT[0]   = (double)fp_p1_TOT_PMT1;
            gTOT_Amax_PMT[0] = (double)fp_Amax_TOT_PMT1;
            gTOT_p0_PMT[1]   = (double)fp_p0_TOT_PMT2;
            gTOT_p1_PMT[1]   = (double)fp_p1_TOT_PMT2;
            gTOT_Amax_PMT[1] = (double)fp_Amax_TOT_PMT2;
            // PMT3 non viene recuperato via TOT
            gTOT_p0_PMT[2]   = 0.0;
            gTOT_p1_PMT[2]   = 0.0;
            gTOT_Amax_PMT[2] = 0.0;

            // Aggiorna la soglia TOT di riferimento
            if (fp_TOT_thr_ref > 0.0) {
                TOT_THR_REF = (double)fp_TOT_thr_ref;
            }

            // Calibrazione TOT valida se entrambi i PMT hanno p1 > 0
            gTOT_calibrated = (gTOT_p1_PMT[0] > 0.0 && gTOT_p1_PMT[1] > 0.0);

            std::cout << "[INFO] Parametri TOT caricati:" << std::endl;
            std::cout << "       PMT1: p0 = " << gTOT_p0_PMT[0]
                      << ", p1 = " << gTOT_p1_PMT[0]
                      << " (τ_eff = " << 1.0/gTOT_p1_PMT[0] << " ns)"
                      << ", A_max = " << gTOT_Amax_PMT[0] << " mV" << std::endl;
            std::cout << "       PMT2: p0 = " << gTOT_p0_PMT[1]
                      << ", p1 = " << gTOT_p1_PMT[1]
                      << " (τ_eff = " << 1.0/gTOT_p1_PMT[1] << " ns)"
                      << ", A_max = " << gTOT_Amax_PMT[1] << " mV" << std::endl;
            std::cout << "       Soglia TOT ref = " << TOT_THR_REF << " mV" << std::endl;
            std::cout << "       Recupero clippati (TOT): "
                      << (gTOT_calibrated ? "ABILITATO" : "DISABILITATO")
                      << std::endl;
        } else {
            std::cout << "[INFO] Branch TOT non presenti nel file di calibrazione."
                      << std::endl;
            std::cout << "       Recupero clippati (TOT): DISABILITATO" << std::endl;
            gTOT_calibrated = false;
        }
// ---- Caricamento modello lineare C(x) = C_0 + C_1·x ----
        if (fp->GetBranch("C_lin_p0")) {
            Float_t fp_Cp0, fp_Cp0_err, fp_Cp1, fp_Cp1_err, fp_Cchi2;
            fp->SetBranchAddress("C_lin_p0",      &fp_Cp0);
            fp->SetBranchAddress("C_lin_p0_err",  &fp_Cp0_err);
            fp->SetBranchAddress("C_lin_p1",      &fp_Cp1);
            fp->SetBranchAddress("C_lin_p1_err",  &fp_Cp1_err);
            fp->SetBranchAddress("C_lin_chi2ndf", &fp_Cchi2);
            fp->GetEntry(0);

            gC_lin_p0     = (double)fp_Cp0;
            gC_lin_p0_err = (double)fp_Cp0_err;
            gC_lin_p1     = (double)fp_Cp1;
            gC_lin_p1_err = (double)fp_Cp1_err;
            gC_lin_loaded = true;

            std::cout << "[INFO] Modello lineare C(x) caricato:" << std::endl;
            std::cout << "       C(x) = " << gC_lin_p0
                      << " + (" << gC_lin_p1 << ") · x  ns" << std::endl;
            std::cout << "       χ²/ndf = " << fp_Cchi2 << std::endl;
            std::cout << "       → C usata: MODELLO LINEARE" << std::endl;
        } else {
            gC_lin_loaded = false;
            std::cout << "[INFO] Branch C_lin non presenti. C usata: INTERPOLAZIONE"
                      << std::endl;
        }
    } else {
        std::cerr << "[WARNING] TTree 'fit_params' non trovato in "
                  << cal_filename << ".\n           "
                  << "Lasciate CAL_M, CAL_Q ai valori di default oppure "
                  << "usate SetCalibration(m, q) manualmente.\n           "
                  << "Il recupero dei clippati via slew rate sarà DISABILITATO."
                  << std::endl;
        gK_calibrated = false;
    }

    fcal->Close();
    return true;
}
/// SetCalibration(): imposta m, q da terminale ROOT senza ricompilare.
void SetCalibration(double m, double q) {
    CAL_M = m;
    CAL_Q = q;
    std::cout << "[INFO] Calibrazione aggiornata: m = " << m
              << " ns/cm, q = " << q << " ns" << std::endl;
    std::cout << "       v_eff = " << 2.0/fabs(m) << " cm/ns" << std::endl;
}

/// SetGeometry(): imposta h e x_PMT3 da terminale ROOT senza ricompilare.
void SetGeometry(double h, double x_pmt3 = 0.0) {
    PAR_H = h;
    PAR_X_PMT3 = x_pmt3;
    std::cout << "[INFO] Geometria aggiornata: h = " << h
              << " cm, x_PMT3 = " << x_pmt3 << " cm" << std::endl;
}


// ==========================================================================
//  FUNZIONE AUSILIARIA: parsing del nome file FIFO → orario di start FPGA
// ==========================================================================
//
//  Estrae l'orario di inizio acquisizione FPGA dal nome del file.
//  Formato atteso: "...FIFOread_YYYYMMDD-HHMMSS.txt" (il percorso può
//  contenere directory, viene estratto solo il basename).
//
//  Restituisce il numero di secondi dall'inizio della giornata (00:00:00)
//  per l'orario FPGA, oppure -1 in caso di errore.
//
//  NOTA: si usa solo la parte oraria (HH:MM:SS) perché il filtraggio
//  si basa sul tempo trascorso dall'inizio dell'acquisizione all'interno
//  della stessa sessione; la data è implicita (stessa giornata del run).

static double ParseFifoStartTime(const std::string &fifo_path) {

    // Estrai il basename (dopo l'ultima '/')
    size_t slash = fifo_path.rfind('/');
    std::string fname = (slash == std::string::npos)
                        ? fifo_path
                        : fifo_path.substr(slash + 1);

    // Cerca il pattern "FIFOread_YYYYMMDD-HHMMSS"
    // es. "FIFOread_20260507-163005.txt"
    size_t pos = fname.find("FIFOread_");
    if (pos == std::string::npos) {
        std::cerr << "[ERRORE] ParseFifoStartTime: pattern 'FIFOread_' "
                  << "non trovato in: " << fname << std::endl;
        return -1.0;
    }

    // Dopo "FIFOread_" ci sono 8 cifre di data + '-' + 6 cifre di ora
    size_t date_start = pos + 9;  // lunghezza di "FIFOread_"
    if (fname.size() < date_start + 15) {
        std::cerr << "[ERRORE] ParseFifoStartTime: nome file troppo corto: "
                  << fname << std::endl;
        return -1.0;
    }

    // Estrai i campi: YYYYMMDD-HHMMSS
    // Posizioni relative a date_start:  0-3 Y, 4-5 M, 6-7 D, 8='-', 9-10 H, 11-12 M, 13-14 S
    std::string date_part = fname.substr(date_start, 15);
    // Verifica formato: 8 cifre + '-' + 6 cifre
    bool fmt_ok = true;
    for (int i = 0; i < 15 && fmt_ok; i++) {
        if (i == 8)  { if (date_part[i] != '-') fmt_ok = false; }
        else         { if (!isdigit(date_part[i])) fmt_ok = false; }
    }
    if (!fmt_ok) {
        std::cerr << "[ERRORE] ParseFifoStartTime: formato data non valido: "
                  << date_part << std::endl;
        return -1.0;
    }

    int HH = std::stoi(date_part.substr(9, 2));
    int MM = std::stoi(date_part.substr(11, 2));
    int SS = std::stoi(date_part.substr(13, 2));

    double fpga_start_s = (double)(HH * 3600 + MM * 60 + SS);
    std::cout << "[FIFO] Orario start FPGA letto dal nome file: "
              << HH << ":" << MM << ":" << SS
              << "  (" << fpga_start_s << " s dall'inizio del giorno)"
              << std::endl;
    return fpga_start_s;
}


// ==========================================================================
//  FUNZIONE AUSILIARIA: parsing del timestamp DRS → secondi dall'inizio giorno
// ==========================================================================
//
//  Converte la stringa "YYYY/MM/DD HH:MM:SS.sss" (campo <Time> dell'XML)
//  in secondi dall'inizio della giornata, già corretta per il fuso orario
//  (sottrae FIFO_TZ_OFFSET_S = 7200 s = 2 ore).
//
//  Restituisce -1.0 in caso di errore di parsing.

static double ParseDRSTimestamp(const std::string &ts) {

    // Formato atteso: "2026/05/07 18:30:34.512"
    //                  0123456789012345678901234
    // Posizioni: YYYY=0, MM=5, DD=8, HH=11, MM=14, SS=17, .sss=19
    if (ts.size() < 19) {
        std::cerr << "[ERRORE] ParseDRSTimestamp: timestamp troppo corto: '"
                  << ts << "'" << std::endl;
        return -1.0;
    }

    int HH, MM, SS;
    double frac = 0.0;

    // sscanf con formato fisso
    int matched = sscanf(ts.c_str(),
                         "%*4d/%*2d/%*2d %2d:%2d:%lf",
                         &HH, &MM, &frac);
    if (matched < 3) {
        // Fallback: prova a leggere SS e .sss separatamente
        int SS_int = 0;
        matched = sscanf(ts.c_str(),
                         "%*4d/%*2d/%*2d %2d:%2d:%2d",
                         &HH, &MM, &SS_int);
        if (matched < 3) {
            std::cerr << "[ERRORE] ParseDRSTimestamp: parsing fallito per: '"
                      << ts << "'" << std::endl;
            return -1.0;
        }
        frac = (double)SS_int;
        // Aggiungi la parte decimale se presente
        size_t dot = ts.rfind('.');
        if (dot != std::string::npos) {
            double dec = std::stod(ts.substr(dot));
            frac += dec;
        }
    }

    double drs_s = (double)(HH * 3600 + MM * 60) + frac;

    // Correzione fuso orario: il DRS è avanti di 2 ore rispetto a FPGA
    double fpga_equiv_s = drs_s - (double)FIFO_TZ_OFFSET_S;

    return fpga_equiv_s;
}


// ==========================================================================
//  FUNZIONE AUSILIARIA: caricamento vettori FIFO
// ==========================================================================
//
//  Legge il file FIFO e popola:
//    fifo_chan  — vettore dei canali (1 o 2; i reset vengono conservati
//                 come tipo 0 per tenere traccia dei confini di buffer)
//    fifo_ts    — vettore dei timestamp (cicli di clock, già "assoluti"
//                 rispetto all'inizio del file, sommando i buffer)
//    reset_idx  — indici (nel vettore fifo_*) dove avvengono i reset
//    reset_abs_ts — timestamp assoluto (in cicli di clock dall'inizio)
//                   di ciascun reset
//
//  SCARTA le righe prima del primo reset (colonna 2 = 2^31).
//
//  Restituisce il numero di reset trovati, o -1 in caso di errore.
//
//  NOTA SULLA STRUTTURA DEI TIMESTAMP ASSOLUTI:
//    I timestamp nel file sono relativi all'interno del buffer corrente
//    (0 … 2^30−1). Per la ricerca binaria usiamo timestamp "assoluti"
//    calcolati come:
//      ts_abs = buffer_index × 2^30 + ts_relativo
//    dove buffer_index è il numero di buffer dall'inizio (0-based, scartando
//    le righe prima del primo reset). In questo modo i timestamp sono
//    strettamente crescenti e la ricerca è O(log N).

struct FifoData {
    std::vector<int>       chan;   // 1, 2  (0 = reset, non usato nella ricerca)
    std::vector<long long> ts;    // timestamp assoluti [cicli di clock]
    std::vector<size_t>    reset_pos;   // indici nel vettore dove ci sono reset
    std::vector<long long> reset_ts;    // ts assoluti dei reset (buffer × 2^30)
};

static int LoadFifo(const std::string &fifo_path, FifoData &fd) {

    std::ifstream fin(fifo_path);
    if (!fin.is_open()) {
        std::cerr << "[ERRORE] LoadFifo: impossibile aprire '" << fifo_path
                  << "'" << std::endl;
        return -1;
    }

    std::cout << "[FIFO] Lettura file: " << fifo_path << " ..." << std::flush;

    const long long BUF_SIZE  = (1LL << 30);   // 2^30 cicli di clock
    const long long RESET_VAL = (1LL << 31);   // 2^31

    fd.chan.clear();
    fd.ts.clear();
    fd.reset_pos.clear();
    fd.reset_ts.clear();

    long long col1, col2;
    bool      first_reset_seen = false;
    long long n_resets         = 0;   // numero di reset visti (0-based)
    long long n_lines          = 0;

    while (fin >> col1 >> col2) {
        n_lines++;

        // ---- Riga di reset: col1 == 2^31 ----
        if (col1 == RESET_VAL) {
            first_reset_seen = true;
            // n_resets sarà l'indice del prossimo buffer (0-based: il primo
            // reset segna la fine del buffer 0 e l'inizio del buffer 1)
            long long abs_ts = n_resets * BUF_SIZE;
            fd.reset_pos.push_back(fd.chan.size());
            fd.reset_ts.push_back(abs_ts);
            // Inseriamo nel vettore principale con chan=0 (marker di reset)
            fd.chan.push_back(0);
            fd.ts.push_back(abs_ts);
            n_resets++;
            continue;
        }

        // Scarta le righe prima del primo reset
        if (!first_reset_seen) continue;

        // ---- Riga dati normale: col1 = 1 o 2 ----
        if (col1 != 1 && col1 != 2) continue;  // ignora righe anomale

        // Timestamp assoluto: (n_resets−1) buffer già completati + ts relativo
        // Quando siamo dentro il buffer k (dopo il reset k-1), il buffer count
        // è n_resets (che conta quanti reset abbiamo visto FINORA, non quello
        // in corso). Il primo buffer dopo il reset 0 è il buffer n_resets=1,
        // quindi:
        //   ts_abs = (n_resets - 1) * BUF_SIZE + col2
        // (−1 perché n_resets viene incrementato subito dopo il reset)
        long long abs_ts = (n_resets - 1) * BUF_SIZE + col2;

        fd.chan.push_back((int)col1);
        fd.ts.push_back(abs_ts);
    }

    fin.close();

    std::cout << " " << fd.chan.size() << " voci totali"
              << " (" << fd.reset_pos.size() << " reset)" << std::endl;

    if (fd.reset_pos.empty()) {
        std::cerr << "[ERRORE] LoadFifo: nessun reset trovato nel file FIFO."
                  << std::endl;
        return -1;
    }

    return (int)fd.reset_pos.size();
}


// ==========================================================================
//  FUNZIONE AUSILIARIA: calcolo del timestamp FIFO atteso per un evento DRS
// ==========================================================================
//
//  Dato il timestamp DRS (già corretto per il fuso orario, in secondi
//  dall'inizio del giorno) e l'orario di start FPGA (stessa unità),
//  calcola il timestamp FIFO assoluto (in cicli di clock) atteso.
//
//  ALGORITMO:
//    1. elapsed = drs_s_fpga_equiv − fpga_start_s   [secondi]
//    2. n_buf   = (int)(elapsed / FIFO_BUF_DURATION_S)  [buffer da saltare, troncato]
//    3. reset_time = n_buf × FIFO_BUF_DURATION_S    [secondi]
//    4. residuo = elapsed − reset_time               [secondi entro il buffer]
//    5. ts_clk  = round(residuo / FIFO_CLK_PERIOD_S) [cicli di clock]
//    6. ts_abs  = n_buf × 2^30 + ts_clk             [cicli assoluti]
//
//  Restituisce il timestamp assoluto, o -1 in caso di errore (elapsed < 0).

static long long CalcExpectedFifoTs(double drs_s_fpga_equiv, double fpga_start_s) {

    double elapsed = drs_fpga_equiv - fpga_start_s;

    if (elapsed < 0.0) return -1LL;

    long long ts_abs = (long long)(elapsed / FIFO_CLK_PERIOD_S + 0.5);

    return ts_abs;
    
}


// ==========================================================================
//  FUNZIONE AUSILIARIA: ricerca coincidenza nel FIFO
// ==========================================================================
//
//  Dato un timestamp assoluto atteso (ts_expected), cerca nel vettore fd:
//    1. Un evento di canale 1 entro [ts_expected−FIFO_COINC_WINDOW,
//                                    ts_expected+FIFO_COINC_WINDOW]
//    2. Se trovato al timestamp ts1, cerca un evento di canale 2
//       nel range [ts1+FIFO_CH2_MIN_DT, ts1+FIFO_CH2_MAX_DT]
//
//  La ricerca usa std::lower_bound per efficienza (O(log N)).
//
//  Restituisce true se entrambe le condizioni sono soddisfatte (evento valido).

static bool FifoCoincidenceCheck(const FifoData &fd, long long ts_expected) {

    if (ts_expected < 0LL) return false;
    if (fd.ts.empty())     return false;

    // ---- Ricerca binaria nell'intervallo [ts_lo, ts_hi] per canale 1 ----
    long long ts_lo = ts_expected - FIFO_COINC_WINDOW;
    long long ts_hi = ts_expected + FIFO_COINC_WINDOW;

    // lower_bound: primo elemento con ts >= ts_lo
    auto it_lo = std::lower_bound(fd.ts.begin(), fd.ts.end(), ts_lo);
    if (it_lo == fd.ts.end()) return false;

    // Scorri gli elementi nella finestra [ts_lo, ts_hi] cercando canale 1
    long long ts1_found = -1LL;
    for (auto it = it_lo; it != fd.ts.end() && *it <= ts_hi; ++it) {
        size_t idx = (size_t)(it - fd.ts.begin());
        if (fd.chan[idx] == 1) {
            ts1_found = *it;
            break;
        }
    }

    if (ts1_found < 0LL) return false;

    // ---- Ricerca canale 2 in avanti: [ts1+CH2_MIN, ts1+CH2_MAX] ----
    long long ts2_lo = ts1_found + FIFO_CH2_MIN_DT;
    long long ts2_hi = ts1_found + FIFO_CH2_MAX_DT;

    auto it2_lo = std::lower_bound(fd.ts.begin(), fd.ts.end(), ts2_lo);
    if (it2_lo == fd.ts.end()) return false;

    for (auto it = it2_lo; it != fd.ts.end() && *it <= ts2_hi; ++it) {
        size_t idx = (size_t)(it - fd.ts.begin());
        if (fd.chan[idx] == 2) {
            return true;   // coincidenza trovata!
        }
    }

    return false;
}


// ==========================================================================
//  SEZIONE 7: FUNZIONE PRINCIPALE — ANALISI TOF
// ==========================================================================

/// TOF_Analysis(): funzione principale. Processa uno o più file XML di
/// dati TOF (acquisiti con PMT3 in Guida B) e ricostruisce per ogni
/// evento il tempo di volo, la velocità e l'angolo di incidenza.
///
/// ARGOMENTI:
///   xml_files — stringa con i percorsi dei file XML, separati da virgola.
///               Es: "run1.xml,run2.xml,run3.xml"
///   outname   — nome del file ROOT di output
///   cal_file  — (opzionale) file ROOT di calibrazione da cui caricare C(x_k)

void TOF_Analysis(const char* xml_files,
                  const char* outname  = "TOF_Analysis_output.root",
                  const char* cal_file = "") {
    gROOT->SetBatch(kTRUE); //disattiva le finestre grafiche
    std::cout << "=============================================" << std::endl;
    std::cout << "  TOF Analysis                               " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Geometria: h = " << PAR_H << " cm, x_PMT3 = "
              << PAR_X_PMT3 << " cm" << std::endl;
    std::cout << "  Calibrazione: m = " << CAL_M << " ns/cm, q = "
              << CAL_Q << " ns" << std::endl;
    std::cout << "  v_eff = " << 2.0/fabs(CAL_M) << " cm/ns" << std::endl;
    std::cout << "=============================================" << std::endl;

    // --- Caricamento calibrazione (opzionale) ---
    if (strlen(cal_file) > 0) {
        LoadCalibrationFromFile(cal_file);
    }

    // --- Parsing della lista di file XML ---
    // I file sono separati da virgola nella stringa xml_files
    std::vector<std::string> file_list;
    {
        std::string files_str(xml_files);
        std::istringstream ss(files_str);
        std::string token;
        while (std::getline(ss, token, ',')) {
            // Rimuovi spazi iniziali e finali
            size_t start = token.find_first_not_of(" \t");
            size_t end   = token.find_last_not_of(" \t");
            if (start != std::string::npos)
                file_list.push_back(token.substr(start, end - start + 1));
        }
    }

    if (file_list.empty()) {
        std::cerr << "[ERRORE] Nessun file XML specificato." << std::endl;
        return;
    }
    std::cout << "[INFO] " << file_list.size() << " file da processare." << std::endl;

    // --- Setup stile grafico ---
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetTitleSize(0.045, "xy");

    // --- Apri il file di output ---
    TFile *fout = new TFile(outname, "RECREATE");

    // --- Preparazione del TTree di output ---
    TTree *tree = new TTree("tof_data", "Dati TOF evento per evento");

    // Tempi CFD dei 3 PMT [ns]
    Float_t t_cfd[3];
    tree->Branch("t_cfd", t_cfd, "t_cfd[3]/F");

    // Ampiezze dei 3 PMT [mV]
    Float_t amp[3];
    tree->Branch("amp", amp, "amp[3]/F");

    // Differenza temporale Δt₁₂ [ns]
    Float_t dt12;
    tree->Branch("dt12", &dt12, "dt12/F");

    // Posizione di impatto ricostruita sulla barra [cm]
    Float_t x_imp;
    tree->Branch("x_imp", &x_imp, "x_imp/F");

    // T_misurato = t₃ − (t₁+t₂)/2  [ns]
    Float_t T_meas;
    tree->Branch("T_meas", &T_meas, "T_meas/F");

    // C(x_imp) usato per questo evento [ns]
    Float_t C_used;
    tree->Branch("C_used", &C_used, "C_used/F");

    // Tempo di volo = T_meas − C  [ns]
    Float_t tof;
    tree->Branch("tof", &tof, "tof/F");

    // Distanza percorsa dalla particella [cm]
    Float_t path_len;
    tree->Branch("path_len", &path_len, "path_len/F");

    // Angolo di incidenza rispetto alla verticale [rad]
    Float_t theta;
    tree->Branch("theta", &theta, "theta/F");

    // Velocità [cm/ns]
    Float_t vel;
    tree->Branch("vel", &vel, "vel/F");

    // β = v/c
    Float_t beta;
    tree->Branch("beta", &beta, "beta/F");



    // 1/v [ns/cm] — distribuzione più simmetrica di v
    Float_t inv_vel;
    tree->Branch("inv_vel", &inv_vel, "inv_vel/F");

    // Flag qualità
    Int_t good;
    tree->Branch("good", &good, "good/I");

// Indice del file sorgente (per tracciabilità)
    Int_t file_idx;
    tree->Branch("file_idx", &file_idx, "file_idx/I");

    // ============================================================
    //  BRANCH PARALLELI PER IL METODO TOT
    // ============================================================
    // Le due pipeline di analisi (slew rate e TOT) sono COMPLETAMENTE
    // INDIPENDENTI: ciascuna usa solo il proprio metodo di recupero
    // dei clippati, calcola i propri Δt₁₂, x_imp, TOF, β, ecc., e
    // applica i propri tagli. Questo evita di mescolare i due metodi
    // e permette confronti diretti a posteriori.
    //
    // Convenzione di naming: stessi nomi dei branch slew, con suffisso
    // "_tot". Per gli eventi che NON sono clippati, i due metodi
    // producono ESATTAMENTE gli stessi numeri (perché il t_cfd standard
    // viene usato in entrambi i casi). Le differenze emergono solo
    // sugli eventi clippati, dove le due ricostruzioni sono indipendenti.

    // ---- Tempi e ampiezze ricostruite via TOT ----
    Float_t t_cfd_tot[3], amp_tot[3];
    tree->Branch("t_cfd_tot", t_cfd_tot, "t_cfd_tot[3]/F");
    tree->Branch("amp_tot",   amp_tot,   "amp_tot[3]/F");

    // ---- Quantità cinematiche calcolate col metodo TOT ----
    Float_t dt12_tot, x_imp_tot;
    Float_t T_meas_tot, C_used_tot, tof_tot;
    Float_t path_len_tot, theta_tot, vel_tot, beta_tot, inv_vel_tot;
    Int_t   good_tot;
    tree->Branch("dt12_tot",     &dt12_tot,     "dt12_tot/F");
    tree->Branch("x_imp_tot",    &x_imp_tot,    "x_imp_tot/F");
    tree->Branch("T_meas_tot",   &T_meas_tot,   "T_meas_tot/F");
    tree->Branch("C_used_tot",   &C_used_tot,   "C_used_tot/F");
    tree->Branch("tof_tot",      &tof_tot,      "tof_tot/F");
    tree->Branch("path_len_tot", &path_len_tot, "path_len_tot/F");
    tree->Branch("theta_tot",    &theta_tot,    "theta_tot/F");
    tree->Branch("vel_tot",      &vel_tot,      "vel_tot/F");
    tree->Branch("beta_tot",     &beta_tot,     "beta_tot/F");
    tree->Branch("inv_vel_tot",  &inv_vel_tot,  "inv_vel_tot/F");
    tree->Branch("good_tot",     &good_tot,     "good_tot/I");

    // ---- Diagnostica: quale canale era clippato e quale metodo l'ha
    // recuperato (per analisi a posteriori senza dover rilanciare). ----
    //   clipped_chan[k]      = 1 se il canale k era clippato
    //   recov_slew[k]        = 1 se recuperato via slew rate
    //   recov_tot[k]         = 1 se recuperato via TOT
    Int_t clipped_chan[3], recov_slew[3], recov_tot[3];
    tree->Branch("clipped_chan", clipped_chan, "clipped_chan[3]/I");
    tree->Branch("recov_slew",   recov_slew,   "recov_slew[3]/I");
    tree->Branch("recov_tot",    recov_tot,    "recov_tot[3]/I");

    // --- Istogrammi ---
    TH1D *h_tof    = new TH1D("h_tof",    "Tempo di volo;TOF [ns];Conteggi",
                               NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta   = new TH1D("h_beta",   "Distribuzione #beta = v/c;#beta;Conteggi",
                               NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_norm   = new TH1D("h_beta_norm",   "Distribuzione #beta normalizzata = v/c;#beta;Conteggi",
                               NBINS_BETA, BETA_LO, BETA_HI);                           
    TH1D *h_invv   = new TH1D("h_invv",   "Distribuzione 1/v;1/v [ns/cm];Conteggi",
                               NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_theta  = new TH1D("h_theta",  "Distribuzione angolare;#theta [rad];Conteggi",
                               NBINS_THETA, THETA_LO, THETA_HI);
    TH1D *h_x_imp  = new TH1D("h_x_imp",  "Posizione di impatto;x_{imp} [cm];Conteggi",
                               NBINS_X, X_LO, X_HI);
    TH1D *h_path   = new TH1D("h_path",   "Distanza percorsa;l [cm];Conteggi",
                               100, 90, 250);
    TH1D *h_T_meas = new TH1D("h_T_meas", "T_{meas} = t_{3} - (t_{1}+t_{2})/2;"
                               "T_{meas} [ns];Conteggi",
                               150, -5, 25);

    // Istogramma 2D: TOF vs x_imp (per cercare anomalie)
    TH2D *h2_tof_x = new TH2D("h2_tof_x",
                               "TOF vs posizione;x_{imp} [cm];TOF [ns]",
                               70, X_LO, X_HI, 75, TOF_LO, TOF_HI);

// Istogramma 2D: beta vs x_imp
    TH2D *h2_beta_x = new TH2D("h2_beta_x",
                                "#beta vs posizione;x_{imp} [cm];#beta",
                                70, X_LO, X_HI, 100, BETA_LO, BETA_HI);

    // ============================================================
    //  ISTOGRAMMI PARALLELI PER IL METODO TOT
    // ============================================================
    // Stessi binning e range degli istogrammi slew, con suffisso "_tot".
    TH1D *h_tof_tot     = new TH1D("h_tof_tot",
                                   "Tempo di volo (TOT);TOF [ns];Conteggi",
                                   NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta_tot    = new TH1D("h_beta_tot",
                                   "Distribuzione #beta (TOT);#beta;Conteggi",
                                   NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_norm_tot = new TH1D("h_beta_norm_tot",
                                   "Distribuzione #beta normalizzata (TOT);"
                                   "#beta;Conteggi",
                                   NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_invv_tot    = new TH1D("h_invv_tot",
                                   "Distribuzione 1/v (TOT);1/v [ns/cm];Conteggi",
                                   NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_theta_tot   = new TH1D("h_theta_tot",
                                   "Distribuzione angolare (TOT);#theta [rad];Conteggi",
                                   NBINS_THETA, THETA_LO, THETA_HI);
    TH1D *h_x_imp_tot   = new TH1D("h_x_imp_tot",
                                   "Posizione di impatto (TOT);x_{imp} [cm];Conteggi",
                                   NBINS_X, X_LO, X_HI);
    TH1D *h_path_tot    = new TH1D("h_path_tot",
                                   "Distanza percorsa (TOT);l [cm];Conteggi",
                                   100, 90, 250);
    TH1D *h_T_meas_tot  = new TH1D("h_T_meas_tot",
                                   "T_{meas} (TOT);T_{meas} [ns];Conteggi",
                                   150, -5, 25);
    TH2D *h2_tof_x_tot  = new TH2D("h2_tof_x_tot",
                                   "TOF vs posizione (TOT);x_{imp} [cm];TOF [ns]",
                                   70, X_LO, X_HI, 75, TOF_LO, TOF_HI);
    TH2D *h2_beta_x_tot = new TH2D("h2_beta_x_tot",
                                   "#beta vs posizione (TOT);x_{imp} [cm];#beta",
                                   70, X_LO, X_HI, 100, BETA_LO, BETA_HI);


// ==================================================================
    //  CALCOLO LIMITE FISICO SUL TOF
    // ==================================================================
    // tof_min_phys = h/c è il TOF MINIMO ASSOLUTO compatibile con la
    // causalità (β ≤ 1) per qualunque traiettoria sulla barra. È usato
    // come hard cut nel TAGLIO 4. Nota: lo definisco QUI (e non come
    // const globale) perché PAR_H può essere modificato a runtime via
    // SetGeometry() prima della chiamata a TOF_Analysis().
    const double tof_min_phys = PAR_H / C_LIGHT;
    std::cout << "  TOF minimo fisico (h/c) = " << tof_min_phys
              << " ns (h = " << PAR_H << " cm, c = " << C_LIGHT
              << " cm/ns)" << std::endl;
    std::cout << "=============================================" << std::endl;

    // ==================================================================
    //  LOOP SUI FILE XML
    // ==================================================================
    int total_events    = 0;
    // Contatori per il metodo SLEW (pipeline principale)
    int total_good      = 0;
    int rej_no_cfd      = 0;
    int rej_clipped     = 0;
    int rej_oscillating = 0;
    int rej_x_out       = 0;
    int rej_tof_unphys  = 0;
    // Contatori paralleli per il metodo TOT
    int total_good_tot      = 0;
    int rej_clipped_tot     = 0;   // clippati non recuperati via TOT
    int rej_x_out_tot       = 0;   // x_imp fuori range usando t_cfd_tot
    int rej_tof_unphys_tot  = 0;   // TOF non fisico usando il metodo TOT

    for (size_t fi = 0; fi < file_list.size(); fi++) {

        std::vector<EventData> events;
        int n_parsed = ParseXML(file_list[fi].c_str(), events);
        if (n_parsed <= 0) {
            std::cerr << "[WARNING] Nessun evento in " << file_list[fi] << std::endl;
            continue;
        }

        // --- Loop sugli eventi ---
for (size_t ev = 0; ev < events.size(); ev++) {

            EventData &e = events[ev];
            total_events++;
            file_idx = (Int_t)fi;

            // ---- Init dei branch diagnostici e dei nuovi branch TOT ----
            // Questi default vengono usati per gli eventi che NON entrano
            // nel blocco di recupero/pipeline (scartati al TAGLIO 0).
            // Per gli eventi che entrano nel blocco, vengono sovrascritti.
            for (int k = 0; k < 3; k++) {
                clipped_chan[k] = 0;
                recov_slew[k]   = 0;
                recov_tot[k]    = 0;
                t_cfd_tot[k]    = -999;
                amp_tot[k]      = 0;
            }
            good_tot = 0;
            dt12_tot = -999;     x_imp_tot = -999;
            T_meas_tot = -999;   C_used_tot = -999;   tof_tot = -999;
            path_len_tot = -999; theta_tot = -999;
            vel_tot = -999;      beta_tot = -999;     inv_vel_tot = -999;

            // ============================================================
            //  TAGLIO 0: almeno 3 canali presenti con CFD valido
            // ============================================================
            if (e.nchannels < 3) { good = 0; tree->Fill(); rej_no_cfd++; continue; }

            bool all_ok = true;
            for (int k = 0; k < 3; k++) {
                t_cfd[k] = (Float_t)e.ch[k].t_cfd;
                amp[k]   = (Float_t)e.ch[k].amplitude;
                if (!e.ch[k].cfd_ok) all_ok = false;
            }

            if (!all_ok) {
                good = 0; dt12 = -999; x_imp = -999;
                T_meas = -999; C_used = -999; tof = -999;
                path_len = -999; theta = -999; vel = -999;
                beta = -999; inv_vel = -999;
                tree->Fill();
                rej_no_cfd++;
                continue;
            }
// ============================================================
// ============================================================
            //  RECUPERO CLIPPATI: SLEW e TOT IN PARALLELO
            // ============================================================
            // Eseguiamo ENTRAMBI i metodi di recupero su ogni canale
            // clippato. I risultati vanno in campi separati della struct
            // (cfd_recovered/amplitude_rec per slew, cfd_recovered_tot/
            // amplitude_tot per TOT). Le due pipeline di analisi useranno
            // ciascuna i propri risultati.
            //
            // PMT3 NON viene mai recuperato (non clippa nei nostri dati;
            // se clippa è caso anomalo). Se PMT3 è clippato, ENTRAMBE le
            // pipeline scartano l'evento.
            for (int k = 0; k < 2; k++) {
                if (e.ch[k].is_clipped) {
                    RecoverClippedCFD    (e.ch[k], k);  // metodo slew
                    RecoverClippedCFD_TOT(e.ch[k], k);  // metodo TOT
                }
            }

            // ---- Popolamento dei branch diagnostici (stato finale) ----
            for (int k = 0; k < 3; k++) {
                clipped_chan[k] = e.ch[k].is_clipped         ? 1 : 0;
                recov_slew[k]   = e.ch[k].cfd_recovered      ? 1 : 0;
                recov_tot[k]    = e.ch[k].cfd_recovered_tot  ? 1 : 0;
            }

            // ============================================================
            //  TAGLIO 2 (PRE-PIPELINE): oscillazione su qualunque canale
            // ============================================================
            // Un impulso fisico di scintillazione è un singolo picco
            // negativo breve. Se >15% dei campioni hanno |V−bl| > 20%
            // dell'ampiezza, il segnale è un'oscillazione patologica.
            // Questo taglio si applica a ENTRAMBE le pipeline (è una
            // proprietà fisica del segnale, indipendente dal metodo di
            // ricostruzione).
            bool any_osc = false;
            for (int k = 0; k < 3; k++) {
                if (e.ch[k].is_oscillating) { any_osc = true; break; }
            }
            if (any_osc) {
                // Riempi tutto a sentinella per ENTRAMBE le pipeline
                good = 0;     dt12 = -999;     x_imp = -999;
                T_meas = -999; C_used = -999; tof = -999;
                path_len = -999; theta = -999; vel = -999;
                beta = -999;  inv_vel = -999;
                good_tot = 0; dt12_tot = -999; x_imp_tot = -999;
                T_meas_tot = -999; C_used_tot = -999; tof_tot = -999;
                path_len_tot = -999; theta_tot = -999; vel_tot = -999;
                beta_tot = -999; inv_vel_tot = -999;
                for (int k = 0; k < 3; k++) {
                    t_cfd_tot[k] = -999;
                    amp_tot[k]   = 0;
                }
                tree->Fill();
                rej_oscillating++;
                continue;
            }

            // ============================================================
            //  PIPELINE PARALLELE: una lambda per metodo
            // ============================================================
            // Definiamo una lambda che prende come input i tempi e le
            // ampiezze del canale (già con eventuale recupero applicato),
            // applica TUTTI i tagli e calcola le quantità derivate.
            //
            // I parametri di output sono passati per riferimento; il
            // valore di ritorno è uno status:
            //    0 = evento "buono"
            //   -1 = clippato non recuperato
            //   -2 = x_imp fuori range
            //   -3 = TOF non fisico (ma con valori calcolati)
            //
            // PMT3 clippato è gestito internamente: la lambda controlla
            // is_clipped sul canale 2 e ritorna -1 se vero.

            auto RunPipeline = [&](
                // INPUT: tempi e ampiezze già con recupero applicato
                const Float_t  t_in[3], const Float_t  amp_in[3],
                bool           ch_recovered[2],   // se PMT1/PMT2 sono recuperati
                // OUTPUT: tutte le quantità calcolate
                Float_t &dt12_o, Float_t &x_imp_o,
                Float_t &T_meas_o, Float_t &C_used_o, Float_t &tof_o,
                Float_t &path_o, Float_t &theta_o,
                Float_t &vel_o, Float_t &beta_o, Float_t &invv_o
            ) -> int {

                // Inizializzazione output a sentinella
                dt12_o = -999; x_imp_o = -999;
                T_meas_o = -999; C_used_o = -999; tof_o = -999;
                path_o = -999; theta_o = -999;
                vel_o = -999; beta_o = -999; invv_o = -999;

                // ---- TAGLIO 1: clipping non recuperato ----
                // PMT3 mai recuperato: se clippato, scarta sempre
                if (e.ch[2].is_clipped) return -1;
                // PMT1/PMT2: se clippato, deve essere recuperato
                for (int k = 0; k < 2; k++) {
                    if (e.ch[k].is_clipped && !ch_recovered[k]) return -1;
                }

                // ---- PASSO 1: Δt₁₂ e posizione di impatto ----
                dt12_o  = t_in[0] - t_in[1];
                x_imp_o = (dt12_o - CAL_Q) / CAL_M;

                // ---- TAGLIO 3: x_imp dentro la barra ----
                if (x_imp_o < X_CUT_LO || x_imp_o > X_CUT_HI) return -2;

                // ---- PASSO 2: tempo di volo ----
                T_meas_o = t_in[2] - (t_in[0] + t_in[1]) / 2.0;
                C_used_o = (Float_t)GetC((double)x_imp_o);
                tof_o    = T_meas_o - C_used_o;

                // ---- TAGLIO 4: TOF fisicamente lecito ----
                if (tof_o < tof_min_phys) return -3;

                // ---- PASSO 3: distanza percorsa ----
                double dx = (double)x_imp_o - PAR_X_PMT3;
                path_o = (Float_t)sqrt(dx * dx + PAR_H * PAR_H);

                // ---- PASSO 4: angolo ----
                theta_o = (Float_t)atan2(dx, PAR_H);

                // ---- PASSO 5: velocità e β ----
                if (tof_o > 0.1) {
                    vel_o   = (Float_t)((double)path_o / (double)tof_o);
                    beta_o  = (Float_t)((double)vel_o / C_LIGHT);
                    invv_o  = (Float_t)((double)tof_o / (double)path_o);
                }

                return 0;
            };

            // ============================================================
            //  PIPELINE 1 — METODO SLEW RATE
            // ============================================================
            // Costruzione array di tempi/ampiezze per il metodo slew:
            // i clippati recuperati sostituiscono t_cfd standard con
            // t_cfd_recovered; i non clippati usano t_cfd direttamente.
            Float_t t_slew[3], amp_slew[3];
            bool    ch_rec_slew[2];
            for (int k = 0; k < 3; k++) {
                if (k < 2 && e.ch[k].is_clipped && e.ch[k].cfd_recovered) {
                    t_slew[k]   = (Float_t)e.ch[k].t_cfd_recovered;
                    amp_slew[k] = (Float_t)e.ch[k].amplitude_rec;
                    ch_rec_slew[k] = true;
                } else {
                    t_slew[k]   = (Float_t)e.ch[k].t_cfd;
                    amp_slew[k] = (Float_t)e.ch[k].amplitude;
                    if (k < 2) ch_rec_slew[k] = false;
                }
            }
            // Copia anche nei branch principali (slew = "metodo principale")
            for (int k = 0; k < 3; k++) {
                t_cfd[k] = t_slew[k];
                amp[k]   = amp_slew[k];
            }

            int status_slew = RunPipeline(t_slew, amp_slew, ch_rec_slew,
                                          dt12, x_imp, T_meas, C_used, tof,
                                          path_len, theta, vel, beta, inv_vel);

            // Aggiornamento contatori e flag good per slew
            switch (status_slew) {
                case  0: good = 1; total_good++;       break;
                case -1: good = 0; rej_clipped++;      break;
                case -2: good = 0; rej_x_out++;        break;
                case -3: good = 0; rej_tof_unphys++;   break;
                default: good = 0;                     break;
            }

            // ============================================================
            //  PIPELINE 2 — METODO TOT
            // ============================================================
            Float_t t_tot_arr[3], amp_tot_arr[3];
            bool    ch_rec_tot[2];
            for (int k = 0; k < 3; k++) {
                if (k < 2 && e.ch[k].is_clipped && e.ch[k].cfd_recovered_tot) {
                    t_tot_arr[k]   = (Float_t)e.ch[k].t_cfd_rec_tot;
                    amp_tot_arr[k] = (Float_t)e.ch[k].amplitude_tot;
                    ch_rec_tot[k]  = true;
                } else {
                    t_tot_arr[k]   = (Float_t)e.ch[k].t_cfd;
                    amp_tot_arr[k] = (Float_t)e.ch[k].amplitude;
                    if (k < 2) ch_rec_tot[k] = false;
                }
            }
            // Branch dedicati al metodo TOT
            for (int k = 0; k < 3; k++) {
                t_cfd_tot[k] = t_tot_arr[k];
                amp_tot[k]   = amp_tot_arr[k];
            }

            int status_tot = RunPipeline(t_tot_arr, amp_tot_arr, ch_rec_tot,
                                         dt12_tot, x_imp_tot, T_meas_tot,
                                         C_used_tot, tof_tot, path_len_tot,
                                         theta_tot, vel_tot, beta_tot,
                                         inv_vel_tot);

            switch (status_tot) {
                case  0: good_tot = 1; total_good_tot++;      break;
                case -1: good_tot = 0; rej_clipped_tot++;     break;
                case -2: good_tot = 0; rej_x_out_tot++;       break;
                case -3: good_tot = 0; rej_tof_unphys_tot++;  break;
                default: good_tot = 0;                        break;
            }

            // ============================================================
            //  RIEMPIMENTO ISTOGRAMMI
            // ============================================================
            // Per il metodo slew (istogrammi originali)
            if (good == 1) {
                h_T_meas->Fill(T_meas);
                h_x_imp->Fill(x_imp);
                if (tof > 0.1) {
                    h_tof->Fill(tof);
                    h_path->Fill(path_len);
                    h_theta->Fill(theta);
                    h2_tof_x->Fill(x_imp, tof);
                    if (beta > 0 && beta < 3.0) {
                        h_beta->Fill(beta);
                        h_beta_norm->Fill(beta);
                        h2_beta_x->Fill(x_imp, beta);
                    }
                    if (inv_vel > 0 && inv_vel < 0.2) {
                        h_invv->Fill(inv_vel);
                    }
                }
            }

            // Per il metodo TOT (istogrammi paralleli "_tot")
            if (good_tot == 1) {
                h_T_meas_tot->Fill(T_meas_tot);
                h_x_imp_tot->Fill(x_imp_tot);
                if (tof_tot > 0.1) {
                    h_tof_tot->Fill(tof_tot);
                    h_path_tot->Fill(path_len_tot);
                    h_theta_tot->Fill(theta_tot);
                    h2_tof_x_tot->Fill(x_imp_tot, tof_tot);
                    if (beta_tot > 0 && beta_tot < 3.0) {
                        h_beta_tot->Fill(beta_tot);
                        h_beta_norm_tot->Fill(beta_tot);
                        h2_beta_x_tot->Fill(x_imp_tot, beta_tot);
                    }
                    if (inv_vel_tot > 0 && inv_vel_tot < 0.2) {
                        h_invv_tot->Fill(inv_vel_tot);
                    }
                }
            }

tree->Fill();
        }
    }

    auto pct = [&](int n)->double{
        return (total_events > 0) ? 100.0 * n / total_events : 0.0;
    };

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  RIEPILOGO TAGLI DI QUALITA'                " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Eventi totali:           " << total_events << std::endl;
    std::cout << "  Scartati (no CFD):       " << rej_no_cfd
              << " (" << pct(rej_no_cfd) << "%)" << std::endl;
    std::cout << "  Scartati (oscillazione): " << rej_oscillating
              << " (" << pct(rej_oscillating) << "%)" << std::endl;
    std::cout << "  ───────── METODO SLEW ─────────────" << std::endl;
    std::cout << "  Scartati (clipped):      " << rej_clipped
              << " (" << pct(rej_clipped) << "%)" << std::endl;
    std::cout << "  Scartati (x fuori):      " << rej_x_out
              << " (" << pct(rej_x_out) << "%)" << std::endl;
    std::cout << "  Scartati (TOF<h/c):      " << rej_tof_unphys
              << " (" << pct(rej_tof_unphys) << "%)" << std::endl;
    std::cout << "  Eventi buoni (slew):     " << total_good
              << " (" << pct(total_good) << "%)" << std::endl;
    std::cout << "  ───────── METODO TOT ──────────────" << std::endl;
    std::cout << "  Scartati (clipped):      " << rej_clipped_tot
              << " (" << pct(rej_clipped_tot) << "%)" << std::endl;
    std::cout << "  Scartati (x fuori):      " << rej_x_out_tot
              << " (" << pct(rej_x_out_tot) << "%)" << std::endl;
    std::cout << "  Scartati (TOF<h/c):      " << rej_tof_unphys_tot
              << " (" << pct(rej_tof_unphys_tot) << "%)" << std::endl;
    std::cout << "  Eventi buoni (TOT):      " << total_good_tot
              << " (" << pct(total_good_tot) << "%)" << std::endl;
    std::cout << "=============================================" << std::endl;
    // ==================================================================
    //  STATISTICHE DELLA DISTRIBUZIONE β
    // ==================================================================
    // Calcolo tre statistiche complementari della distribuzione di β:
    //
    //   1) STATISTICA GREZZA (momenti empirici dell'istogramma h_beta):
    //      <β> e σ(β) calcolati come media e RMS sui contenuti dei bin.
    //      ATTENZIONE: dipendono dal range del binning [BETA_LO, BETA_HI];
    //      eventi fuori dal range non contribuiscono. Risentono delle
    //      code asimmetriche (la trasformazione β = ℓ/(c·TOF) amplifica
    //      gli errori per TOF piccoli, allungando la coda di β grandi).
    //
    //   2) FIT GAUSSIANO DEL PICCO:
    //      Fit nella regione [μ_picco − 0.3, μ_picco + 0.3] dove μ_picco
    //      è il bin di massimo. Restituisce il valore "tipico" di β e la
    //      sua dispersione naturale, robusto rispetto alle code.
    //      È la quantità più adatta per stimare la risoluzione temporale
    //      del sistema TOF se la distribuzione è approssimativamente
    //      gaussiana (controllare il χ²/ndf).
    //
    //   3) STIMA ROBUSTA (mediana e MAD = Median Absolute Deviation):
    //      Indipendente dall'assunzione di gaussianità e dalle code
    //      (50% dei dati su ciascun lato della mediana). Per una
    //      distribuzione gaussiana, σ ≈ 1.4826 · MAD.
    //
    // Le statistiche sono calcolate sull'istogramma h_beta, che contiene
    // SOLO eventi che passano TUTTI i tagli (good=1) e che hanno β nel
    // range [BETA_LO, BETA_HI]. È quindi un campione "pulito".

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  STATISTICHE DISTRIBUZIONE  β = v/c           " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Eventi nell'istogramma h_beta: "
              << (long)h_beta->GetEntries() << std::endl;

    // -------------------------------------------------------------
    //  1) Statistica grezza dell'istogramma
    // -------------------------------------------------------------
    double beta_mean_raw = h_beta->GetMean();
    double beta_rms_raw  = h_beta->GetStdDev();
    // Errore sulla media: σ/√N (assumendo eventi indipendenti)
    double beta_mean_err = (h_beta->GetEntries() > 1)
                           ? beta_rms_raw / sqrt(h_beta->GetEntries())
                           : 0.0;

    std::cout << "\n  [1] Statistica grezza (istogramma):" << std::endl;
    std::cout << "      <β>      = " << beta_mean_raw
              << " ± " << beta_mean_err << std::endl;
    std::cout << "      σ(β)     = " << beta_rms_raw << std::endl;

    // -------------------------------------------------------------
    //  2) Fit Gaussiano del picco
    // -------------------------------------------------------------
    // Strategia: trovo il bin di massimo, faccio un primo fit in una
    // finestra di ±0.3 attorno, poi (opzionalmente) raffino il fit nella
    // finestra ±2σ del primo fit per stabilità (se il picco è stretto).
    double beta_mean_fit = -1, beta_sigma_fit = -1;
    double beta_mean_fit_err = 0, beta_sigma_fit_err = 0;
    double beta_chi2_ndf = -1;

    if (h_beta->GetEntries() > 50) {  // serve un minimo di statistica
        int    bin_peak = h_beta->GetMaximumBin();
        double x_peak   = h_beta->GetBinCenter(bin_peak);

        // Range iniziale di fit: ±0.3 attorno al picco
        double fit_lo = std::max(x_peak - 0.3, BETA_LO + 0.01);
        double fit_hi = std::min(x_peak + 0.3, BETA_HI - 0.01);

        TF1 *fgaus = new TF1("fgaus_beta", "gaus", fit_lo, fit_hi);
        fgaus->SetParameters(h_beta->GetMaximum(), x_peak, 0.1);
        // Opzioni: "Q" silenzioso, "R" rispetta range, "S" salva risultato,
        // "0" non disegna sull'istogramma (lo disegnerò io dopo).
        TFitResultPtr fr = h_beta->Fit(fgaus, "QRS0");

        if (fr.Get() && fr->IsValid()) {
            beta_mean_fit      = fgaus->GetParameter(1);
            beta_sigma_fit     = fabs(fgaus->GetParameter(2));
            beta_mean_fit_err  = fgaus->GetParError(1);
            beta_sigma_fit_err = fgaus->GetParError(2);
            int ndf = fgaus->GetNDF();
            beta_chi2_ndf = (ndf > 0) ? fgaus->GetChisquare() / ndf : -1;

            std::cout << "\n  [2] Fit Gaussiano del picco "
                      << "(range [" << fit_lo << ", " << fit_hi << "]):"
                      << std::endl;
            std::cout << "      μ        = " << beta_mean_fit
                      << " ± " << beta_mean_fit_err << std::endl;
            std::cout << "      σ        = " << beta_sigma_fit
                      << " ± " << beta_sigma_fit_err << std::endl;
            std::cout << "      χ²/ndf   = " << beta_chi2_ndf
                      << "  (ndf = " << ndf << ")" << std::endl;
        } else {
            std::cout << "\n  [2] Fit Gaussiano: NON CONVERGENTE." << std::endl;
        }
    } else {
        std::cout << "\n  [2] Fit Gaussiano: SALTATO "
                  << "(troppe poche entries: "
                  << h_beta->GetEntries() << ")" << std::endl;
    }

    // -------------------------------------------------------------
    //  3) Statistica robusta: mediana e MAD
    // -------------------------------------------------------------
    // GetQuantiles() di ROOT calcola il quantile interpolato dai
    // contenuti dei bin. Per la mediana: quantile 0.5.
    // Per la MAD ho bisogno di mediana(|x − mediana|), che richiede
    // un secondo passaggio: ricostruisco i valori bin per bin pesando
    // per i conteggi e calcolo la mediana delle deviazioni assolute.
    double beta_median = -1, beta_mad = -1;

    if (h_beta->GetEntries() > 0) {
        // Mediana via quantili
        double qprob[1] = {0.5};
        double qval[1]  = {0.0};
        h_beta->GetQuantiles(1, qval, qprob);
        beta_median = qval[0];

        // MAD: costruisco un istogramma temporaneo delle deviazioni
        // assolute |β − mediana|, poi prendo la sua mediana.
        // Lo faccio bin per bin: ogni bin di h_beta contribuisce con
        // un "punto" al centro del bin con peso pari ai conteggi.
        TH1D h_dev("h_dev_beta", "|beta-median|", NBINS_BETA, 0, BETA_HI);
        for (int ib = 1; ib <= h_beta->GetNbinsX(); ++ib) {
            double xb = h_beta->GetBinCenter(ib);
            double cb = h_beta->GetBinContent(ib);
            if (cb > 0) h_dev.Fill(fabs(xb - beta_median), cb);
        }
        if (h_dev.GetEntries() > 0) {
            h_dev.GetQuantiles(1, qval, qprob);
            beta_mad = qval[0];
        }

        // Stima di σ "robusta" assumendo gaussianità: σ ≈ 1.4826 · MAD
        double sigma_from_mad = (beta_mad > 0) ? 1.4826 * beta_mad : -1;

        std::cout << "\n  [3] Statistica robusta:" << std::endl;
        std::cout << "      mediana       = " << beta_median << std::endl;
        std::cout << "      MAD           = " << beta_mad << std::endl;
        std::cout << "      1.4826·MAD    = " << sigma_from_mad
                  << "  (≈ σ se distribuzione gaussiana)" << std::endl;
    }

    std::cout << "=============================================" << std::endl;

    // ==================================================================
    //  NORMALIZZAZIONE h_beta_norm → PDF
    // ==================================================================
    // h_beta_norm viene trasformata in una funzione di densità di probabilità
    // (PDF) dividendo ogni bin per l'integrale totale calcolato "con larghezza"
    // (opzione "width" di Integral): questo corrisponde a calcolare
    //
    //     ∫ PDF(β) dβ  =  Σ_i  [contenuto_i × Δβ]  =  1
    //
    // Dopo la normalizzazione l'asse Y ha unità di [1/β], cioè la probabilità
    // per unità di β. Il vantaggio rispetto a h_beta (in conteggi assoluti) è
    // che l'istogramma normalizzato può essere confrontato direttamente con
    // distribuzioni provenienti da dataset con numero diverso di eventi acquisiti
    // (run differenti, configurazioni diverse, risultati di letteratura, ecc.).
    //
    // NOTA: la normalizzazione NON altera h_beta (usato per le statistiche e il
    // fit gaussiano precedente); agisce solo sulla copia h_beta_norm.
    if (h_beta_norm->Integral() > 0.0) {
        // Integral("width") = Σ contenuto_i * larghezza_bin_i = integrale numerico
        // dell'istogramma rispettando le larghezze dei bin (esatto per bin uniformi).
        double norm_factor = h_beta_norm->Integral("width");
        h_beta_norm->Scale(1.0 / norm_factor);
        // Aggiorna il titolo dell'asse Y per riflettere le unità della PDF
        h_beta_norm->GetYaxis()->SetTitle("PDF [1/unit. #beta]");
        std::cout << "[INFO] h_beta_norm normalizzata all'area unitaria."
                  << "  Fattore = " << norm_factor
                  << "  (entries = " << (long)h_beta_norm->GetEntries() << ")"
                  << std::endl;
} else {
        std::cout << "[WARNING] h_beta_norm vuota: normalizzazione saltata." << std::endl;
    }

    // ---- Stessa normalizzazione per h_beta_norm_tot (metodo TOT) ----
    if (h_beta_norm_tot->Integral() > 0.0) {
        double norm_factor_tot = h_beta_norm_tot->Integral("width");
        h_beta_norm_tot->Scale(1.0 / norm_factor_tot);
        h_beta_norm_tot->GetYaxis()->SetTitle("PDF [1/unit. #beta]");
        std::cout << "[INFO] h_beta_norm_tot normalizzata all'area unitaria."
                  << "  Fattore = " << norm_factor_tot
                  << "  (entries = " << (long)h_beta_norm_tot->GetEntries() << ")"
                  << std::endl;
    } else {
        std::cout << "[WARNING] h_beta_norm_tot vuota: normalizzazione saltata." << std::endl;
    }

    // ==================================================================
    //  SALVATAGGIO E CANVAS
    // ==================================================================
    fout->cd();
    tree->Write();

    // --- Canvas 1: distribuzione TOF ---
    TCanvas *c1 = new TCanvas("c_tof", "Distribuzione TOF", 800, 600);
    c1->SetGrid();
    h_tof->SetLineColor(kBlue + 1);
    h_tof->SetLineWidth(2);
    h_tof->Draw();
    c1->Write("Canvas_TOF");

// --- Canvas 2: distribuzione β con fit Gaussiano e statistiche ---
    TCanvas *c2 = new TCanvas("c_beta", "Distribuzione #beta", 900, 650);
    c2->SetGrid();
    h_beta->SetLineColor(kRed + 1);
    h_beta->SetLineWidth(2);
    h_beta->Draw();

    // Overlay del fit Gaussiano (se è stato fatto e converge)
    TF1 *fgaus_show = (TF1*)gROOT->FindObject("fgaus_beta");
    if (fgaus_show && beta_mean_fit > 0) {
        fgaus_show->SetLineColor(kBlue + 2);
        fgaus_show->SetLineWidth(2);
        fgaus_show->SetLineStyle(1);
        fgaus_show->Draw("same");
    }

    // Linea verticale a β = 1 (limite causalità)
    TLine *l_beta1 = new TLine(1.0, 0, 1.0, h_beta->GetMaximum() * 0.95);
    l_beta1->SetLineColor(kBlack);
    l_beta1->SetLineStyle(2);
    l_beta1->SetLineWidth(2);
    l_beta1->Draw("same");

    // Box con le statistiche (in alto a destra)
    // Uso TPaveText con coordinate normalizzate (NDC).
    TPaveText *pt = new TPaveText(0.60, 0.55, 0.89, 0.88, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);   // allineamento sx, centrato verticale
    pt->SetTextFont(42);
    pt->SetTextSize(0.030);
    pt->AddText(Form("Entries = %lld", (Long64_t)h_beta->GetEntries()));
    pt->AddText("");
    pt->AddText("#bf{Statistica grezza}");
    pt->AddText(Form("#LT#beta#GT = %.4f #pm %.4f",
                     beta_mean_raw, beta_mean_err));
    pt->AddText(Form("#sigma(#beta) = %.4f", beta_rms_raw));
    if (beta_mean_fit > 0) {
        pt->AddText("");
        pt->AddText("#bf{Fit gaussiano}");
        pt->AddText(Form("#mu = %.4f #pm %.4f",
                         beta_mean_fit, beta_mean_fit_err));
        pt->AddText(Form("#sigma = %.4f #pm %.4f",
                         beta_sigma_fit, beta_sigma_fit_err));
        pt->AddText(Form("#chi^{2}/ndf = %.2f", beta_chi2_ndf));
    }
    if (beta_median > 0) {
        pt->AddText("");
        pt->AddText("#bf{Robusta}");
        pt->AddText(Form("mediana = %.4f", beta_median));
        pt->AddText(Form("MAD = %.4f", beta_mad));
    }
    pt->Draw();

    c2->Write("Canvas_beta");

    // --- Canvas 2b: PDF di β normalizzata all'area unitaria ---
    // Questo canvas mostra h_beta_norm, la versione PDF di h_beta.
    // Poiché l'integrale è 1 per costruzione, la forma è direttamente
    // confrontabile con distribuzioni da altri dataset (basta sovrapporre
    // due canvas con Draw("same") o salvare il TH1D e rileggerlo).
    // La linea tratteggiata a β = 1 marca il limite di causalità.
    TCanvas *c_beta_pdf = new TCanvas("c_beta_pdf",
                                      "PDF di #beta normalizzata", 900, 650);
    c_beta_pdf->SetGrid();
    h_beta_norm->SetLineColor(kRed + 1);
    h_beta_norm->SetLineWidth(2);
    h_beta_norm->Draw("HIST");

    // Linea verticale a β = 1 (limite causalità: nessuna particella massiva
    // può viaggiare a v = c, quindi il picco fisico deve stare a β < 1).
    TLine *l_b1_pdf = new TLine(1.0, 0.0, 1.0, h_beta_norm->GetMaximum() * 0.95);
    l_b1_pdf->SetLineColor(kBlack);
    l_b1_pdf->SetLineStyle(2);   // linea tratteggiata
    l_b1_pdf->SetLineWidth(2);
    l_b1_pdf->Draw("same");

    // Etichetta "β = 1" accanto alla linea verticale
    TLatex *lat_pdf = new TLatex(1.02, h_beta_norm->GetMaximum() * 0.80, "#beta = 1");
    lat_pdf->SetTextSize(0.035);
    lat_pdf->SetTextColor(kBlack);
    lat_pdf->Draw("same");

    // Box informativo: numero di eventi originali e conferma della normalizzazione.
    // Il numero di entries è quello PRIMA della normalizzazione (Scale non cambia
    // GetEntries), quindi rappresenta la statistica reale del campione.
    TPaveText *pt_pdf = new TPaveText(0.60, 0.72, 0.89, 0.88, "NDC");
    pt_pdf->SetFillColor(0);
    pt_pdf->SetBorderSize(1);
    pt_pdf->SetTextAlign(12);
    pt_pdf->SetTextFont(42);
    pt_pdf->SetTextSize(0.030);
    pt_pdf->AddText(Form("Entries = %lld", (Long64_t)h_beta_norm->GetEntries()));
    pt_pdf->AddText("#int PDF(#beta) d#beta = 1");   // conferma normalizzazione
    pt_pdf->Draw();

    // Sovrapponi la PDF del metodo TOT sullo stesso canvas
    h_beta_norm_tot->SetLineColor(kBlue + 1);
    h_beta_norm_tot->SetLineWidth(2);
    h_beta_norm_tot->SetLineStyle(2);
    h_beta_norm_tot->Draw("HIST same");

    // Aggiorna la legenda con entrambe le curve
    TLegend *leg_pdf = new TLegend(0.15, 0.72, 0.42, 0.88);
    leg_pdf->SetTextFont(42);
    leg_pdf->SetTextSize(0.030);
    leg_pdf->AddEntry(h_beta_norm,     "Slew rate", "l");
    leg_pdf->AddEntry(h_beta_norm_tot, "TOT",       "l");
    leg_pdf->Draw();

    c_beta_pdf->Write("Canvas_beta_PDF");

    // --- Canvas 2c: PDF di β TOT separata ---
    TCanvas *c_beta_pdf_tot = new TCanvas("c_beta_pdf_tot",
                                          "PDF di #beta normalizzata (TOT)", 900, 650);
    c_beta_pdf_tot->SetGrid();
    h_beta_norm_tot->SetLineColor(kBlue + 1);
    h_beta_norm_tot->SetLineWidth(2);
    h_beta_norm_tot->SetLineStyle(1);
    h_beta_norm_tot->Draw("HIST");

    TLine *l_b1_pdf_tot = new TLine(1.0, 0.0, 1.0,
                                     h_beta_norm_tot->GetMaximum() * 0.95);
    l_b1_pdf_tot->SetLineColor(kBlack);
    l_b1_pdf_tot->SetLineStyle(2);
    l_b1_pdf_tot->SetLineWidth(2);
    l_b1_pdf_tot->Draw("same");

    TLatex *lat_pdf_tot = new TLatex(1.02,
                                      h_beta_norm_tot->GetMaximum() * 0.80,
                                      "#beta = 1");
    lat_pdf_tot->SetTextSize(0.035);
    lat_pdf_tot->SetTextColor(kBlack);
    lat_pdf_tot->Draw("same");

    TPaveText *pt_pdf_tot = new TPaveText(0.60, 0.72, 0.89, 0.88, "NDC");
    pt_pdf_tot->SetFillColor(0);
    pt_pdf_tot->SetBorderSize(1);
    pt_pdf_tot->SetTextAlign(12);
    pt_pdf_tot->SetTextFont(42);
    pt_pdf_tot->SetTextSize(0.030);
    pt_pdf_tot->AddText(Form("Entries = %lld",
                              (Long64_t)h_beta_norm_tot->GetEntries()));
    pt_pdf_tot->AddText("#int PDF(#beta) d#beta = 1");
    pt_pdf_tot->Draw();

    c_beta_pdf_tot->Write("Canvas_beta_PDF_TOT");    

    // --- Canvas 3: distribuzione 1/v ---
    TCanvas *c3 = new TCanvas("c_invv", "Distribuzione 1/v", 800, 600);
    c3->SetGrid();
    h_invv->SetLineColor(kGreen + 2);
    h_invv->SetLineWidth(2);
    h_invv->Draw();
    // Linea verticale a 1/c
    double inv_c = 1.0 / C_LIGHT;
    TLine *l_invc = new TLine(inv_c, 0, inv_c, h_invv->GetMaximum() * 0.9);
    l_invc->SetLineColor(kBlack);
    l_invc->SetLineStyle(2);
    l_invc->SetLineWidth(2);
    l_invc->Draw("same");
    c3->Write("Canvas_invv");

    // --- Canvas 4: distribuzione angolare θ ---
    TCanvas *c4 = new TCanvas("c_theta", "Distribuzione angolare", 800, 600);
    c4->SetGrid();
    h_theta->SetLineColor(kMagenta + 1);
    h_theta->SetLineWidth(2);
    h_theta->Draw();
    c4->Write("Canvas_theta");

    // --- Canvas 5: posizione di impatto ---
    TCanvas *c5 = new TCanvas("c_ximp", "Posizione di impatto", 800, 600);
    c5->SetGrid();
    h_x_imp->SetLineColor(kOrange + 1);
    h_x_imp->SetLineWidth(2);
    h_x_imp->Draw();
    c5->Write("Canvas_x_imp");

    // --- Canvas 6: TOF vs x (2D) ---
    TCanvas *c6 = new TCanvas("c_tof_vs_x", "TOF vs posizione", 900, 600);
    c6->SetGrid();
    h2_tof_x->Draw("COLZ");
    c6->Write("Canvas_TOF_vs_x");

    // --- Canvas 7: β vs x (2D) ---
    TCanvas *c7 = new TCanvas("c_beta_vs_x", "#beta vs posizione", 900, 600);
    c7->SetGrid();
    h2_beta_x->Draw("COLZ");
    // Linea orizzontale a β = 1
    TLine *l_b1 = new TLine(X_LO, 1.0, X_HI, 1.0);
    l_b1->SetLineColor(kRed);
    l_b1->SetLineStyle(2);
    l_b1->Draw("same");
    c7->Write("Canvas_beta_vs_x");

    // --- Canvas 8: pannello riassuntivo 2×3 ---
    TCanvas *c_summary = new TCanvas("c_summary", "Riassunto TOF", 1200, 800);
    c_summary->Divide(3, 2);

    c_summary->cd(1); gPad->SetGrid(); h_tof->Draw();
    c_summary->cd(2); gPad->SetGrid(); h_beta->Draw();
    c_summary->cd(3); gPad->SetGrid(); h_invv->Draw();
    c_summary->cd(4); gPad->SetGrid(); h_theta->Draw();
    c_summary->cd(5); gPad->SetGrid(); h_x_imp->Draw();
    c_summary->cd(6); gPad->SetGrid(); h_path->Draw();

    c_summary->Write("Canvas_Summary");

    // Salva tutti gli istogrammi
    h_tof->Write();
    h_beta->Write();
    h_invv->Write();
    h_theta->Write();
    h_x_imp->Write();
    h_path->Write();
    h_T_meas->Write();
    h2_tof_x->Write();
    h2_beta_x->Write();
    // Salva la PDF di β: l'asse Y è in unità di [1/β] dopo la normalizzazione,
    // quindi questo TH1D può essere letto da un'altra sessione ROOT e sovrapposto
    // a distribuzioni da altri dataset con TH1D::Draw("same") per un confronto diretto.
    h_beta_norm->Write();

    // ---- Istogrammi del metodo TOT (paralleli) ----
    h_tof_tot->Write();
    h_beta_tot->Write();
    h_invv_tot->Write();
    h_theta_tot->Write();
    h_x_imp_tot->Write();
    h_path_tot->Write();
    h_T_meas_tot->Write();
    h2_tof_x_tot->Write();
    h2_beta_x_tot->Write();
    h_beta_norm_tot->Write();

    fout->Close();

    // --- Riepilogo finale ---
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  ANALISI TOF COMPLETATA                      " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  File processati:        " << file_list.size() << std::endl;
    std::cout << "  Eventi totali:          " << total_events << std::endl;
    std::cout << "  Eventi buoni (slew):    " << total_good << std::endl;
    std::cout << "  Eventi buoni (TOT):     " << total_good_tot << std::endl;
    std::cout << "  Parametri usati:" << std::endl;
    std::cout << "    h       = " << PAR_H << " cm" << std::endl;
    std::cout << "    x_PMT3  = " << PAR_X_PMT3 << " cm" << std::endl;
    std::cout << "    m (cal) = " << CAL_M << " ns/cm" << std::endl;
    std::cout << "    q (cal) = " << CAL_Q << " ns" << std::endl;
    std::cout << "  Output: " << outname << std::endl;
    std::cout << "=============================================" << std::endl;
}


// ==========================================================================
//  FUNZIONE PRINCIPALE: TOF_Analysis_lead
// ==========================================================================
//
//  Identica a TOF_Analysis() ma con il filtro FIFO aggiuntivo.
//  Analizza solo gli eventi DRS che passano il controllo di coincidenza.
//
//  ARGOMENTI:
//    xml_files   — stringa con i percorsi dei file XML, separati da virgola
//    fifo_file   — percorso del file FIFOread (es. "FIFOread_20260507-163005.txt")
//    outname     — nome del file ROOT di output (default: "TOF_lead_output.root")
//    cal_file    — (opzionale) file ROOT di calibrazione

void TOF_Analysis_lead(const char* xml_files,
                       const char* fifo_file,
                       const char* outname  = "TOF_lead_output.root",
                       const char* cal_file = "") {

    gROOT->SetBatch(kTRUE);

    std::cout << "=============================================" << std::endl;
    std::cout << "  TOF Analysis with Lead (FIFO coincidence)  " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Geometria: h = " << PAR_H << " cm, x_PMT3 = "
              << PAR_X_PMT3 << " cm" << std::endl;
    std::cout << "  Calibrazione: m = " << CAL_M << " ns/cm, q = "
              << CAL_Q << " ns" << std::endl;
    std::cout << "  v_eff = " << 2.0/fabs(CAL_M) << " cm/ns" << std::endl;
    std::cout << "  FIFO file: " << fifo_file << std::endl;
    std::cout << "=============================================" << std::endl;

    // --- Caricamento calibrazione (opzionale) ---
    if (strlen(cal_file) > 0) {
        LoadCalibrationFromFile(cal_file);
    }

    // ---------------------------------------------------------------
    //  PASSO A: parsing orario di start FPGA dal nome del file FIFO
    // ---------------------------------------------------------------
    double fpga_start_s = ParseFifoStartTime(std::string(fifo_file));
    if (fpga_start_s < 0.0) {
        std::cerr << "[ERRORE] Impossibile determinare l'orario di start FPGA "
                  << "dal nome del file FIFO. Interruzione." << std::endl;
        return;
    }

    // ---------------------------------------------------------------
    //  PASSO B: caricamento del file FIFO in memoria
    // ---------------------------------------------------------------
    FifoData fd;
    int n_resets = LoadFifo(std::string(fifo_file), fd);
    if (n_resets < 0) {
        std::cerr << "[ERRORE] Caricamento FIFO fallito. Interruzione." << std::endl;
        return;
    }

    // ---------------------------------------------------------------
    //  PASSO C: parsing della lista di file XML
    // ---------------------------------------------------------------
    std::vector<std::string> file_list;
    {
        std::string files_str(xml_files);
        std::istringstream ss(files_str);
        std::string token;
        while (std::getline(ss, token, ',')) {
            size_t start = token.find_first_not_of(" \t");
            size_t end   = token.find_last_not_of(" \t");
            if (start != std::string::npos)
                file_list.push_back(token.substr(start, end - start + 1));
        }
    }

    if (file_list.empty()) {
        std::cerr << "[ERRORE] Nessun file XML specificato." << std::endl;
        return;
    }
    std::cout << "[INFO] " << file_list.size() << " file XML da processare." << std::endl;

    // --- Setup stile grafico ---
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetTitleSize(0.045, "xy");

    // --- Apri il file di output ---
    TFile *fout = new TFile(outname, "RECREATE");

    // ---------------------------------------------------------------
    //  PASSO D: preparazione del TTree di output (identico a TOF_Analysis)
    // ---------------------------------------------------------------
    TTree *tree = new TTree("tof_data", "Dati TOF evento per evento (lead filtered)");

    Float_t t_cfd[3];       tree->Branch("t_cfd",       t_cfd,       "t_cfd[3]/F");
    Float_t amp[3];         tree->Branch("amp",          amp,         "amp[3]/F");
    Float_t dt12;           tree->Branch("dt12",        &dt12,        "dt12/F");
    Float_t x_imp;          tree->Branch("x_imp",       &x_imp,       "x_imp/F");
    Float_t T_meas;         tree->Branch("T_meas",      &T_meas,      "T_meas/F");
    Float_t C_used;         tree->Branch("C_used",      &C_used,      "C_used/F");
    Float_t tof;            tree->Branch("tof",         &tof,         "tof/F");
    Float_t path_len;       tree->Branch("path_len",    &path_len,    "path_len/F");
    Float_t theta;          tree->Branch("theta",       &theta,       "theta/F");
    Float_t vel;            tree->Branch("vel",         &vel,         "vel/F");
    Float_t beta;           tree->Branch("beta",        &beta,        "beta/F");
    Float_t inv_vel;        tree->Branch("inv_vel",     &inv_vel,     "inv_vel/F");
    Int_t   good;           tree->Branch("good",        &good,        "good/I");
    Int_t   file_idx;       tree->Branch("file_idx",    &file_idx,    "file_idx/I");

    // Branch paralleli metodo TOT
    Float_t t_cfd_tot[3];   tree->Branch("t_cfd_tot",   t_cfd_tot,   "t_cfd_tot[3]/F");
    Float_t amp_tot[3];     tree->Branch("amp_tot",      amp_tot,     "amp_tot[3]/F");
    Float_t dt12_tot;       tree->Branch("dt12_tot",    &dt12_tot,    "dt12_tot/F");
    Float_t x_imp_tot;      tree->Branch("x_imp_tot",   &x_imp_tot,   "x_imp_tot/F");
    Float_t T_meas_tot;     tree->Branch("T_meas_tot",  &T_meas_tot,  "T_meas_tot/F");
    Float_t C_used_tot;     tree->Branch("C_used_tot",  &C_used_tot,  "C_used_tot/F");
    Float_t tof_tot;        tree->Branch("tof_tot",     &tof_tot,     "tof_tot/F");
    Float_t path_len_tot;   tree->Branch("path_len_tot",&path_len_tot,"path_len_tot/F");
    Float_t theta_tot;      tree->Branch("theta_tot",   &theta_tot,   "theta_tot/F");
    Float_t vel_tot;        tree->Branch("vel_tot",     &vel_tot,     "vel_tot/F");
    Float_t beta_tot;       tree->Branch("beta_tot",    &beta_tot,    "beta_tot/F");
    Float_t inv_vel_tot;    tree->Branch("inv_vel_tot", &inv_vel_tot, "inv_vel_tot/F");
    Int_t   good_tot;       tree->Branch("good_tot",    &good_tot,    "good_tot/I");

    // Diagnostica clipping/recupero
    Int_t clipped_chan[3], recov_slew[3], recov_tot[3];
    tree->Branch("clipped_chan", clipped_chan, "clipped_chan[3]/I");
    tree->Branch("recov_slew",   recov_slew,   "recov_slew[3]/I");
    tree->Branch("recov_tot",    recov_tot,    "recov_tot[3]/I");

    // Branch aggiuntivo: timestamp FIFO atteso (per diagnostica)
    Long64_t fifo_ts_expected;
    tree->Branch("fifo_ts_expected", &fifo_ts_expected, "fifo_ts_expected/L");
    Int_t fifo_ok;
    tree->Branch("fifo_ok", &fifo_ok, "fifo_ok/I");

    // ---------------------------------------------------------------
    //  PASSO E: istogrammi (identici a TOF_Analysis, con suffisso "_lead")
    // ---------------------------------------------------------------
    TH1D *h_tof      = new TH1D("h_tof",      "Tempo di volo (lead);TOF [ns];Conteggi",
                                 NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta     = new TH1D("h_beta",      "Distribuzione #beta (lead);#beta;Conteggi",
                                 NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_norm= new TH1D("h_beta_norm", "PDF #beta (lead);#beta;PDF",
                                 NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_invv     = new TH1D("h_invv",      "Distribuzione 1/v (lead);1/v [ns/cm];Conteggi",
                                 NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_theta    = new TH1D("h_theta",     "Distribuzione angolare (lead);#theta [rad];Conteggi",
                                 NBINS_THETA, THETA_LO, THETA_HI);
    TH1D *h_x_imp    = new TH1D("h_x_imp",     "Posizione impatto (lead);x_{imp} [cm];Conteggi",
                                 NBINS_X, X_LO, X_HI);
    TH1D *h_path     = new TH1D("h_path",      "Distanza percorsa (lead);l [cm];Conteggi",
                                 100, 90, 250);
    TH1D *h_T_meas   = new TH1D("h_T_meas",    "T_{meas} (lead);T_{meas} [ns];Conteggi",
                                 150, -5, 25);
    TH2D *h2_tof_x   = new TH2D("h2_tof_x",    "TOF vs posizione (lead);x_{imp} [cm];TOF [ns]",
                                 70, X_LO, X_HI, 75, TOF_LO, TOF_HI);
    TH2D *h2_beta_x  = new TH2D("h2_beta_x",   "#beta vs posizione (lead);x_{imp} [cm];#beta",
                                 70, X_LO, X_HI, 100, BETA_LO, BETA_HI);

    // Istogrammi TOT paralleli
    TH1D *h_tof_tot      = new TH1D("h_tof_tot",      "TOF (lead, TOT);TOF [ns];Conteggi",
                                     NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta_tot     = new TH1D("h_beta_tot",     "#beta (lead, TOT);#beta;Conteggi",
                                     NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_norm_tot= new TH1D("h_beta_norm_tot","PDF #beta (lead, TOT);#beta;PDF",
                                     NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_invv_tot     = new TH1D("h_invv_tot",     "1/v (lead, TOT);1/v [ns/cm];Conteggi",
                                     NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_theta_tot    = new TH1D("h_theta_tot",    "#theta (lead, TOT);#theta [rad];Conteggi",
                                     NBINS_THETA, THETA_LO, THETA_HI);
    TH1D *h_x_imp_tot    = new TH1D("h_x_imp_tot",   "Posizione impatto (lead, TOT);x_{imp} [cm];Conteggi",
                                     NBINS_X, X_LO, X_HI);
    TH1D *h_path_tot     = new TH1D("h_path_tot",     "Distanza percorsa (lead, TOT);l [cm];Conteggi",
                                     100, 90, 250);
    TH1D *h_T_meas_tot   = new TH1D("h_T_meas_tot",  "T_{meas} (lead, TOT);T_{meas} [ns];Conteggi",
                                     150, -5, 25);
    TH2D *h2_tof_x_tot   = new TH2D("h2_tof_x_tot",  "TOF vs posizione (lead, TOT);x_{imp} [cm];TOF [ns]",
                                     70, X_LO, X_HI, 75, TOF_LO, TOF_HI);
    TH2D *h2_beta_x_tot  = new TH2D("h2_beta_x_tot", "#beta vs posizione (lead, TOT);x_{imp} [cm];#beta",
                                     70, X_LO, X_HI, 100, BETA_LO, BETA_HI);

    // Istogrammi filtrati FIFO (lead) — metodo slew
    TH1D *h_tof_lead      = new TH1D("h_tof_lead",      "TOF (lead);TOF [ns];Conteggi",            NBINS_TOF,  TOF_LO,  TOF_HI);
    TH1D *h_beta_lead     = new TH1D("h_beta_lead",      "#beta (lead);#beta;Conteggi",              NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_norm_lead= new TH1D("h_beta_norm_lead", "PDF #beta (lead);#beta;PDF",               NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_invv_lead     = new TH1D("h_invv_lead",      "1/v (lead);1/v [ns/cm];Conteggi",         NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_theta_lead    = new TH1D("h_theta_lead",     "#theta (lead);#theta [rad];Conteggi",      NBINS_THETA,THETA_LO,THETA_HI);
    TH1D *h_x_imp_lead    = new TH1D("h_x_imp_lead",     "x_{imp} (lead);x_{imp} [cm];Conteggi",    NBINS_X,    X_LO,    X_HI);
    TH1D *h_path_lead     = new TH1D("h_path_lead",      "Distanza (lead);l [cm];Conteggi",         100, 90, 250);
    TH1D *h_T_meas_lead   = new TH1D("h_T_meas_lead",    "T_{meas} (lead);T_{meas} [ns];Conteggi",  150, -5, 25);
    TH2D *h2_tof_x_lead   = new TH2D("h2_tof_x_lead",   "TOF vs x (lead);x_{imp} [cm];TOF [ns]",   70, X_LO, X_HI, 75, TOF_LO, TOF_HI);
    TH2D *h2_beta_x_lead  = new TH2D("h2_beta_x_lead",  "#beta vs x (lead);x_{imp} [cm];#beta",    70, X_LO, X_HI, 100, BETA_LO, BETA_HI);

    // Istogrammi filtrati FIFO (lead) — metodo TOT
    TH1D *h_tof_tot_lead      = new TH1D("h_tof_tot_lead",      "TOF (lead,TOT);TOF [ns];Conteggi",          NBINS_TOF,  TOF_LO,  TOF_HI);
    TH1D *h_beta_tot_lead     = new TH1D("h_beta_tot_lead",     "#beta (lead,TOT);#beta;Conteggi",            NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_norm_tot_lead= new TH1D("h_beta_norm_tot_lead","PDF #beta (lead,TOT);#beta;PDF",             NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_invv_tot_lead     = new TH1D("h_invv_tot_lead",     "1/v (lead,TOT);1/v [ns/cm];Conteggi",       NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_theta_tot_lead    = new TH1D("h_theta_tot_lead",    "#theta (lead,TOT);#theta [rad];Conteggi",    NBINS_THETA,THETA_LO,THETA_HI);
    TH1D *h_x_imp_tot_lead    = new TH1D("h_x_imp_tot_lead",   "x_{imp} (lead,TOT);x_{imp} [cm];Conteggi",  NBINS_X,    X_LO,    X_HI);
    TH1D *h_path_tot_lead     = new TH1D("h_path_tot_lead",     "Distanza (lead,TOT);l [cm];Conteggi",       100, 90, 250);
    TH1D *h_T_meas_tot_lead   = new TH1D("h_T_meas_tot_lead",  "T_{meas} (lead,TOT);T_{meas} [ns];Conteggi",150, -5, 25);
    TH2D *h2_tof_x_tot_lead   = new TH2D("h2_tof_x_tot_lead",  "TOF vs x (lead,TOT);x_{imp} [cm];TOF [ns]", 70, X_LO, X_HI, 75, TOF_LO, TOF_HI);
    TH2D *h2_beta_x_tot_lead  = new TH2D("h2_beta_x_tot_lead", "#beta vs x (lead,TOT);x_{imp} [cm];#beta",  70, X_LO, X_HI, 100, BETA_LO, BETA_HI);
    
    // Limite fisico TOF
    const double tof_min_phys = PAR_H / C_LIGHT;
    std::cout << "  TOF minimo fisico (h/c) = " << tof_min_phys
              << " ns" << std::endl;
    std::cout << "  Finestra coincidenza:   ± " << FIFO_COINC_WINDOW
              << " clock (± " << FIFO_COINC_WINDOW * 5.0e-3 << " µs)" << std::endl;
    std::cout << "  Finestra canale 2:      [" << FIFO_CH2_MIN_DT
              << ", " << FIFO_CH2_MAX_DT << "] clock = ["
              << FIFO_CH2_MIN_DT * 5.0e-3 << ", " << FIFO_CH2_MAX_DT * 5.0e-3
              << "] µs" << std::endl;
    std::cout << "=============================================" << std::endl;

    // ---------------------------------------------------------------
    //  PASSO F: loop sui file XML
    // ---------------------------------------------------------------
    int total_events    = 0;
    int rej_fifo        = 0;    // scartati dal filtro FIFO
    int rej_fifo_notime = 0;    // timestamp FIFO non calcolabile (evento prima dell'FPGA)
    // Contatori identici a TOF_Analysis
    int total_good      = 0;
    int rej_no_cfd      = 0;
    int rej_clipped     = 0;
    int rej_oscillating = 0;
    int rej_x_out       = 0;
    int rej_tof_unphys  = 0;
    int total_good_tot      = 0;
    int rej_clipped_tot     = 0;
    int rej_x_out_tot       = 0;
    int rej_tof_unphys_tot  = 0;
    

    for (size_t fi = 0; fi < file_list.size(); fi++) {

        std::vector<EventData> events;
        int n_parsed = ParseXML(file_list[fi].c_str(), events);
        if (n_parsed <= 0) {
            std::cerr << "[WARNING] Nessun evento in " << file_list[fi] << std::endl;
            continue;
        }

        // --- Loop sugli eventi DRS ---
        for (size_t ev = 0; ev < events.size(); ev++) {

            EventData &e = events[ev];
            total_events++;
            file_idx = (Int_t)fi;

            // Inizializzazione branch diagnostici
            for (int k = 0; k < 3; k++) {
                clipped_chan[k] = 0;
                recov_slew[k]   = 0;
                recov_tot[k]    = 0;
                t_cfd_tot[k]    = -999;
                amp_tot[k]      = 0;
            }
            good_tot = 0;
            dt12_tot = -999;     x_imp_tot = -999;
            T_meas_tot = -999;   C_used_tot = -999;   tof_tot = -999;
            path_len_tot = -999; theta_tot = -999;
            vel_tot = -999;      beta_tot = -999;     inv_vel_tot = -999;
            fifo_ts_expected = -1LL;
            fifo_ok = 0;

            // ============================================================
            //  FILTRO FIFO — eseguito prima di qualsiasi altro taglio
            // ============================================================
            //
            // Se l'evento non ha un timestamp DRS valido o non supera la
            // coincidenza, viene scartato subito (non va nel TTree).
            // Questo rende l'output identico a TOF_Analysis() ma solo
            // sugli eventi che passano il filtro hardware.

            // Calcola il timestamp FIFO atteso per questo evento DRS
            // Calcola il timestamp FIFO atteso
            double drs_fpga_equiv = ParseDRSTimestamp(e.timestamp);
            if (drs_fpga_equiv < 0.0) {
            fifo_ok = 0;
            fifo_ts_expected = -1LL;
            rej_fifo_notime++;
            } else {
             long long ts_exp = CalcExpectedFifoTs(drs_fpga_equiv, fpga_start_s);
             fifo_ts_expected = ts_exp;
            if (ts_exp < 0LL) {
                fifo_ok = 0;
                rej_fifo_notime++;
            } else if (!FifoCoincidenceCheck(fd, ts_exp)) {
                fifo_ok = 0;
                rej_fifo++;
            } else {
                fifo_ok = 1;
            }
            }

            static int debug_count = 0;
            if (debug_count < 10) {
                std::cout << "[DEBUG] Event " << e.serial
                        << "  ts_DRS='" << e.timestamp << "'"
                        << "  drs_fpga_equiv=" << drs_fpga_equiv
                        << "  fpga_start_s=" << fpga_start_s
                        << "  ts_expected=" << fifo_ts_expected
                        << "  fifo_ok=" << fifo_ok
                        << std::endl;
                debug_count++;
            }
                // NON ci sono più continue qui: si prosegue sempre con la pipeline

            // ============================================================
            //  DA QUI IN POI: identico a TOF_Analysis()
            // ============================================================

            // TAGLIO 0: almeno 3 canali con CFD valido
            if (e.nchannels < 3) {
                good = 0; tree->Fill(); rej_no_cfd++; continue;
            }

            bool all_ok = true;
            for (int k = 0; k < 3; k++) {
                t_cfd[k] = (Float_t)e.ch[k].t_cfd;
                amp[k]   = (Float_t)e.ch[k].amplitude;
                if (!e.ch[k].cfd_ok) all_ok = false;
            }

            if (!all_ok) {
                good = 0; dt12 = -999; x_imp = -999;
                T_meas = -999; C_used = -999; tof = -999;
                path_len = -999; theta = -999; vel = -999;
                beta = -999; inv_vel = -999;
                tree->Fill();
                rej_no_cfd++;
                continue;
            }

            // Recupero clippati: SLEW e TOT in parallelo
            for (int k = 0; k < 2; k++) {
                if (e.ch[k].is_clipped) {
                    RecoverClippedCFD    (e.ch[k], k);
                    RecoverClippedCFD_TOT(e.ch[k], k);
                }
            }

            // Branch diagnostici
            for (int k = 0; k < 3; k++) {
                clipped_chan[k] = e.ch[k].is_clipped         ? 1 : 0;
                recov_slew[k]   = e.ch[k].cfd_recovered      ? 1 : 0;
                recov_tot[k]    = e.ch[k].cfd_recovered_tot  ? 1 : 0;
            }

            // TAGLIO 2: oscillazione su qualunque canale
            bool any_osc = false;
            for (int k = 0; k < 3; k++) {
                if (e.ch[k].is_oscillating) { any_osc = true; break; }
            }
            if (any_osc) {
                good = 0;     dt12 = -999;     x_imp = -999;
                T_meas = -999; C_used = -999; tof = -999;
                path_len = -999; theta = -999; vel = -999;
                beta = -999;  inv_vel = -999;
                good_tot = 0; dt12_tot = -999; x_imp_tot = -999;
                T_meas_tot = -999; C_used_tot = -999; tof_tot = -999;
                path_len_tot = -999; theta_tot = -999; vel_tot = -999;
                beta_tot = -999; inv_vel_tot = -999;
                for (int k = 0; k < 3; k++) { t_cfd_tot[k] = -999; amp_tot[k] = 0; }
                tree->Fill();
                rej_oscillating++;
                continue;
            }

            // Lambda pipeline (identica a TOF_Analysis)
            auto RunPipeline = [&](
                const Float_t  t_in[3], const Float_t  amp_in[3],
                bool           ch_recovered[2],
                Float_t &dt12_o, Float_t &x_imp_o,
                Float_t &T_meas_o, Float_t &C_used_o, Float_t &tof_o,
                Float_t &path_o, Float_t &theta_o,
                Float_t &vel_o, Float_t &beta_o, Float_t &invv_o
            ) -> int {

                dt12_o = -999; x_imp_o = -999;
                T_meas_o = -999; C_used_o = -999; tof_o = -999;
                path_o = -999; theta_o = -999;
                vel_o = -999; beta_o = -999; invv_o = -999;

                if (e.ch[2].is_clipped) return -1;
                for (int k = 0; k < 2; k++) {
                    if (e.ch[k].is_clipped && !ch_recovered[k]) return -1;
                }

                dt12_o  = t_in[0] - t_in[1];
                x_imp_o = (dt12_o - CAL_Q) / CAL_M;

                if (x_imp_o < X_CUT_LO || x_imp_o > X_CUT_HI) return -2;

                T_meas_o = t_in[2] - (t_in[0] + t_in[1]) / 2.0;
                C_used_o = (Float_t)GetC((double)x_imp_o);
                tof_o    = T_meas_o - C_used_o;

                if (tof_o < tof_min_phys) return -3;

                double dx = (double)x_imp_o - PAR_X_PMT3;
                path_o = (Float_t)sqrt(dx * dx + PAR_H * PAR_H);
                theta_o = (Float_t)atan2(dx, PAR_H);

                if (tof_o > 0.1) {
                    vel_o   = (Float_t)((double)path_o / (double)tof_o);
                    beta_o  = (Float_t)((double)vel_o / C_LIGHT);
                    invv_o  = (Float_t)((double)tof_o / (double)path_o);
                }

                return 0;
            };

            // Pipeline 1 — metodo SLEW
            Float_t t_slew[3], amp_slew[3];
            bool    ch_rec_slew[2];
            for (int k = 0; k < 3; k++) {
                if (k < 2 && e.ch[k].is_clipped && e.ch[k].cfd_recovered) {
                    t_slew[k]      = (Float_t)e.ch[k].t_cfd_recovered;
                    amp_slew[k]    = (Float_t)e.ch[k].amplitude_rec;
                    ch_rec_slew[k] = true;
                } else {
                    t_slew[k]   = (Float_t)e.ch[k].t_cfd;
                    amp_slew[k] = (Float_t)e.ch[k].amplitude;
                    if (k < 2) ch_rec_slew[k] = false;
                }
            }
            for (int k = 0; k < 3; k++) { t_cfd[k] = t_slew[k]; amp[k] = amp_slew[k]; }

            int status_slew = RunPipeline(t_slew, amp_slew, ch_rec_slew,
                                          dt12, x_imp, T_meas, C_used, tof,
                                          path_len, theta, vel, beta, inv_vel);
            switch (status_slew) {
                case  0: good = 1; total_good++;       break;
                case -1: good = 0; rej_clipped++;      break;
                case -2: good = 0; rej_x_out++;        break;
                case -3: good = 0; rej_tof_unphys++;   break;
                default: good = 0;                     break;
            }

            // Pipeline 2 — metodo TOT
            Float_t t_tot_arr[3], amp_tot_arr[3];
            bool    ch_rec_tot_arr[2];
            for (int k = 0; k < 3; k++) {
                if (k < 2 && e.ch[k].is_clipped && e.ch[k].cfd_recovered_tot) {
                    t_tot_arr[k]      = (Float_t)e.ch[k].t_cfd_rec_tot;
                    amp_tot_arr[k]    = (Float_t)e.ch[k].amplitude_tot;
                    ch_rec_tot_arr[k] = true;
                } else {
                    t_tot_arr[k]   = (Float_t)e.ch[k].t_cfd;
                    amp_tot_arr[k] = (Float_t)e.ch[k].amplitude;
                    if (k < 2) ch_rec_tot_arr[k] = false;
                }
            }
            for (int k = 0; k < 3; k++) { t_cfd_tot[k] = t_tot_arr[k]; amp_tot[k] = amp_tot_arr[k]; }

            int status_tot = RunPipeline(t_tot_arr, amp_tot_arr, ch_rec_tot_arr,
                                         dt12_tot, x_imp_tot, T_meas_tot,
                                         C_used_tot, tof_tot, path_len_tot,
                                         theta_tot, vel_tot, beta_tot, inv_vel_tot);
            switch (status_tot) {
                case  0: good_tot = 1; total_good_tot++;      break;
                case -1: good_tot = 0; rej_clipped_tot++;     break;
                case -2: good_tot = 0; rej_x_out_tot++;       break;
                case -3: good_tot = 0; rej_tof_unphys_tot++;  break;
                default: good_tot = 0;                        break;
            }

            // Riempimento istogrammi — identico a TOF_Analysis
            if (good == 1) {
                h_T_meas->Fill(T_meas);
                h_x_imp->Fill(x_imp);
                if (tof > 0.1) {
                    h_tof->Fill(tof);
                    h_path->Fill(path_len);
                    h_theta->Fill(theta);
                    h2_tof_x->Fill(x_imp, tof);
                    if (beta > 0 && beta < 3.0) {
                        h_beta->Fill(beta);
                        h_beta_norm->Fill(beta);
                        h2_beta_x->Fill(x_imp, beta);
                    }
                    if (inv_vel > 0 && inv_vel < 0.2) h_invv->Fill(inv_vel);
                }
                // Aggiunta: riempimento lead slew
                if (fifo_ok == 1) {
                    h_T_meas_lead->Fill(T_meas);
                    h_x_imp_lead->Fill(x_imp);
                    if (tof > 0.1) {
                        h_tof_lead->Fill(tof);
                        h_path_lead->Fill(path_len);
                        h_theta_lead->Fill(theta);
                        h2_tof_x_lead->Fill(x_imp, tof);
                        if (beta > 0 && beta < 3.0) {
                            h_beta_lead->Fill(beta);
                            h_beta_norm_lead->Fill(beta);
                            h2_beta_x_lead->Fill(x_imp, beta);
                         }
                        if (inv_vel > 0 && inv_vel < 0.2) h_invv_lead->Fill(inv_vel);
                        }
                    }
            }
            if (good_tot == 1) {
                h_T_meas_tot->Fill(T_meas_tot);
                h_x_imp_tot->Fill(x_imp_tot);
                if (tof_tot > 0.1) {
                    h_tof_tot->Fill(tof_tot);
                    h_path_tot->Fill(path_len_tot);
                    h_theta_tot->Fill(theta_tot);
                    h2_tof_x_tot->Fill(x_imp_tot, tof_tot);
                    if (beta_tot > 0 && beta_tot < 3.0) {
                        h_beta_tot->Fill(beta_tot);
                        h_beta_norm_tot->Fill(beta_tot);
                        h2_beta_x_tot->Fill(x_imp_tot, beta_tot);
                    }
                    if (inv_vel_tot > 0 && inv_vel_tot < 0.2) h_invv_tot->Fill(inv_vel_tot);
                }
                 // Aggiunta: riempimento lead TOT
                if (fifo_ok == 1) {
                h_T_meas_tot_lead->Fill(T_meas_tot);
                h_x_imp_tot_lead->Fill(x_imp_tot);
                if (tof_tot > 0.1) {
                    h_tof_tot_lead->Fill(tof_tot);
                    h_path_tot_lead->Fill(path_len_tot);
                    h_theta_tot_lead->Fill(theta_tot);
                    h2_tof_x_tot_lead->Fill(x_imp_tot, tof_tot);
                    if (beta_tot > 0 && beta_tot < 3.0) {
                        h_beta_tot_lead->Fill(beta_tot);
                        h_beta_norm_tot_lead->Fill(beta_tot);
                        h2_beta_x_tot_lead->Fill(x_imp_tot, beta_tot);
            }
            if (inv_vel_tot > 0 && inv_vel_tot < 0.2) h_invv_tot_lead->Fill(inv_vel_tot);
        }
    }
            }

            tree->Fill();

        }  // fine loop eventi
    }  // fine loop file

    // ---------------------------------------------------------------
    //  RIEPILOGO TAGLI
    // ---------------------------------------------------------------
    auto pct = [&](int n) -> double {
        return (total_events > 0) ? 100.0 * n / total_events : 0.0;
    };

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  RIEPILOGO TAGLI DI QUALITA' (LEAD)         " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Eventi totali letti:       " << total_events << std::endl;
    std::cout << "  Scartati (fifo notime):    " << rej_fifo_notime
              << " (" << pct(rej_fifo_notime) << "%)" << std::endl;
    std::cout << "  Scartati (FIFO no coinc.): " << rej_fifo
              << " (" << pct(rej_fifo) << "%)" << std::endl;
    int total_after_fifo = total_events - rej_fifo - rej_fifo_notime;
    std::cout << "  Dopo filtro FIFO:          " << total_after_fifo << std::endl;
    std::cout << "  Scartati (no CFD):         " << rej_no_cfd
              << " (" << pct(rej_no_cfd) << "%)" << std::endl;
    std::cout << "  Scartati (oscillazione):   " << rej_oscillating
              << " (" << pct(rej_oscillating) << "%)" << std::endl;
    std::cout << "  ───────── METODO SLEW ─────────────" << std::endl;
    std::cout << "  Scartati (clipped):        " << rej_clipped
              << " (" << pct(rej_clipped) << "%)" << std::endl;
    std::cout << "  Scartati (x fuori):        " << rej_x_out
              << " (" << pct(rej_x_out) << "%)" << std::endl;
    std::cout << "  Scartati (TOF<h/c):        " << rej_tof_unphys
              << " (" << pct(rej_tof_unphys) << "%)" << std::endl;
    std::cout << "  Eventi buoni (slew):       " << total_good
              << " (" << pct(total_good) << "%)" << std::endl;
    std::cout << "  ───────── METODO TOT ──────────────" << std::endl;
    std::cout << "  Scartati (clipped):        " << rej_clipped_tot
              << " (" << pct(rej_clipped_tot) << "%)" << std::endl;
    std::cout << "  Scartati (x fuori):        " << rej_x_out_tot
              << " (" << pct(rej_x_out_tot) << "%)" << std::endl;
    std::cout << "  Scartati (TOF<h/c):        " << rej_tof_unphys_tot
              << " (" << pct(rej_tof_unphys_tot) << "%)" << std::endl;
    std::cout << "  Eventi buoni (TOT):        " << total_good_tot
              << " (" << pct(total_good_tot) << "%)" << std::endl;
    std::cout << "=============================================" << std::endl;

    // ---------------------------------------------------------------
    //  STATISTICHE β — identiche a TOF_Analysis
    // ---------------------------------------------------------------
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  STATISTICHE DISTRIBUZIONE  β = v/c (LEAD)  " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Eventi nell'istogramma h_beta: "
              << (long)h_beta->GetEntries() << std::endl;

    double beta_mean_raw = h_beta->GetMean();
    double beta_rms_raw  = h_beta->GetStdDev();
    double beta_mean_err = (h_beta->GetEntries() > 1)
                           ? beta_rms_raw / sqrt(h_beta->GetEntries()) : 0.0;
    std::cout << "\n  [1] Statistica grezza:" << std::endl;
    std::cout << "      <β>    = " << beta_mean_raw << " ± " << beta_mean_err << std::endl;
    std::cout << "      σ(β)   = " << beta_rms_raw << std::endl;

    double beta_mean_fit = -1, beta_sigma_fit = -1;
    double beta_mean_fit_err = 0, beta_sigma_fit_err = 0;
    double beta_chi2_ndf = -1;

    if (h_beta->GetEntries() > 50) {
        int    bin_peak = h_beta->GetMaximumBin();
        double x_peak   = h_beta->GetBinCenter(bin_peak);
        double fit_lo = std::max(x_peak - 0.3, BETA_LO + 0.01);
        double fit_hi = std::min(x_peak + 0.3, BETA_HI - 0.01);
        TF1 *fgaus = new TF1("fgaus_beta_lead", "gaus", fit_lo, fit_hi);
        fgaus->SetParameters(h_beta->GetMaximum(), x_peak, 0.1);
        TFitResultPtr fr = h_beta->Fit(fgaus, "QRS0");
        if (fr.Get() && fr->IsValid()) {
            beta_mean_fit      = fgaus->GetParameter(1);
            beta_sigma_fit     = fabs(fgaus->GetParameter(2));
            beta_mean_fit_err  = fgaus->GetParError(1);
            beta_sigma_fit_err = fgaus->GetParError(2);
            int ndf = fgaus->GetNDF();
            beta_chi2_ndf = (ndf > 0) ? fgaus->GetChisquare() / ndf : -1;
            std::cout << "\n  [2] Fit Gaussiano del picco:" << std::endl;
            std::cout << "      μ       = " << beta_mean_fit << " ± " << beta_mean_fit_err << std::endl;
            std::cout << "      σ       = " << beta_sigma_fit << " ± " << beta_sigma_fit_err << std::endl;
            std::cout << "      χ²/ndf  = " << beta_chi2_ndf << std::endl;
        } else {
            std::cout << "\n  [2] Fit Gaussiano: NON CONVERGENTE." << std::endl;
        }
    } else {
        std::cout << "\n  [2] Fit Gaussiano: SALTATO (poche entries)." << std::endl;
    }

    double beta_median = -1, beta_mad = -1;
    if (h_beta->GetEntries() > 0) {
        double qprob[1] = {0.5};
        double qval[1]  = {0.0};
        h_beta->GetQuantiles(1, qval, qprob);
        beta_median = qval[0];
        TH1D h_dev("h_dev_beta_lead", "|beta-median|", NBINS_BETA, 0, BETA_HI);
        for (int ib = 1; ib <= h_beta->GetNbinsX(); ++ib) {
            double xb = h_beta->GetBinCenter(ib);
            double cb = h_beta->GetBinContent(ib);
            if (cb > 0) h_dev.Fill(fabs(xb - beta_median), cb);
        }
        if (h_dev.GetEntries() > 0) {
            h_dev.GetQuantiles(1, qval, qprob);
            beta_mad = qval[0];
        }
        double sigma_from_mad = (beta_mad > 0) ? 1.4826 * beta_mad : -1;
        std::cout << "\n  [3] Statistica robusta:" << std::endl;
        std::cout << "      mediana    = " << beta_median << std::endl;
        std::cout << "      MAD        = " << beta_mad << std::endl;
        std::cout << "      1.4826·MAD = " << sigma_from_mad << std::endl;
    }
    std::cout << "=============================================" << std::endl;

    // Normalizzazione PDF β
    if (h_beta_norm->Integral() > 0.0) {
        double nf = h_beta_norm->Integral("width");
        h_beta_norm->Scale(1.0 / nf);
        h_beta_norm->GetYaxis()->SetTitle("PDF [1/unit. #beta]");
        std::cout << "[INFO] h_beta_norm normalizzata (fattore = " << nf << ")" << std::endl;
    }
    if (h_beta_norm_tot->Integral() > 0.0) {
        double nf = h_beta_norm_tot->Integral("width");
        h_beta_norm_tot->Scale(1.0 / nf);
        h_beta_norm_tot->GetYaxis()->SetTitle("PDF [1/unit. #beta]");
        std::cout << "[INFO] h_beta_norm_tot normalizzata (fattore = " << nf << ")" << std::endl;
    }
    if (h_beta_norm_lead->Integral() > 0.0) {
        double nf = h_beta_norm_lead->Integral("width");
        h_beta_norm_lead->Scale(1.0 / nf);
        h_beta_norm_lead->GetYaxis()->SetTitle("PDF [1/unit. #beta]");
    }
    if (h_beta_norm_tot_lead->Integral() > 0.0) {
        double nf = h_beta_norm_tot_lead->Integral("width");
        h_beta_norm_tot_lead->Scale(1.0 / nf);
        h_beta_norm_tot_lead->GetYaxis()->SetTitle("PDF [1/unit. #beta]");
    }

    // ---------------------------------------------------------------
    //  SALVATAGGIO E CANVAS — identici a TOF_Analysis
    // ---------------------------------------------------------------
    fout->cd();
    tree->Write();

    // Canvas 1: TOF
    TCanvas *c1 = new TCanvas("c_tof", "Distribuzione TOF (lead)", 800, 600);
    c1->SetGrid();
    h_tof->SetLineColor(kBlue + 1); h_tof->SetLineWidth(2); h_tof->Draw();
    c1->Write("Canvas_TOF");

    // Canvas 2: β con fit gaussiano e statistiche
    TCanvas *c2 = new TCanvas("c_beta", "Distribuzione #beta (lead)", 900, 650);
    c2->SetGrid();
    h_beta->SetLineColor(kRed + 1); h_beta->SetLineWidth(2); h_beta->Draw();

    TF1 *fgaus_show = (TF1*)gROOT->FindObject("fgaus_beta_lead");
    if (fgaus_show && beta_mean_fit > 0) {
        fgaus_show->SetLineColor(kBlue + 2); fgaus_show->SetLineWidth(2);
        fgaus_show->SetLineStyle(1); fgaus_show->Draw("same");
    }
    TLine *l_beta1 = new TLine(1.0, 0, 1.0, h_beta->GetMaximum() * 0.95);
    l_beta1->SetLineColor(kBlack); l_beta1->SetLineStyle(2); l_beta1->SetLineWidth(2);
    l_beta1->Draw("same");

    TPaveText *pt = new TPaveText(0.60, 0.55, 0.89, 0.88, "NDC");
    pt->SetFillColor(0); pt->SetBorderSize(1); pt->SetTextAlign(12);
    pt->SetTextFont(42);  pt->SetTextSize(0.030);
    pt->AddText(Form("Entries = %lld", (Long64_t)h_beta->GetEntries()));
    pt->AddText("");
    pt->AddText("#bf{Statistica grezza}");
    pt->AddText(Form("#LT#beta#GT = %.4f #pm %.4f", beta_mean_raw, beta_mean_err));
    pt->AddText(Form("#sigma(#beta) = %.4f", beta_rms_raw));
    if (beta_mean_fit > 0) {
        pt->AddText(""); pt->AddText("#bf{Fit gaussiano}");
        pt->AddText(Form("#mu = %.4f #pm %.4f", beta_mean_fit, beta_mean_fit_err));
        pt->AddText(Form("#sigma = %.4f #pm %.4f", beta_sigma_fit, beta_sigma_fit_err));
        pt->AddText(Form("#chi^{2}/ndf = %.2f", beta_chi2_ndf));
    }
    if (beta_median > 0) {
        pt->AddText(""); pt->AddText("#bf{Robusta}");
        pt->AddText(Form("mediana = %.4f", beta_median));
        pt->AddText(Form("MAD = %.4f", beta_mad));
    }
    pt->Draw();
    c2->Write("Canvas_beta");

    // Canvas 2b: PDF β normalizzata (slew + TOT sovrapposti)
    TCanvas *c_beta_pdf = new TCanvas("c_beta_pdf", "PDF #beta normalizzata (lead)", 900, 650);
    c_beta_pdf->SetGrid();
    h_beta_norm->SetLineColor(kRed + 1); h_beta_norm->SetLineWidth(2);
    h_beta_norm->Draw("HIST");
    TLine *l_b1_pdf = new TLine(1.0, 0.0, 1.0, h_beta_norm->GetMaximum() * 0.95);
    l_b1_pdf->SetLineColor(kBlack); l_b1_pdf->SetLineStyle(2); l_b1_pdf->SetLineWidth(2);
    l_b1_pdf->Draw("same");
    TLatex *lat_pdf = new TLatex(1.02, h_beta_norm->GetMaximum() * 0.80, "#beta = 1");
    lat_pdf->SetTextSize(0.035); lat_pdf->Draw("same");
    h_beta_norm_tot->SetLineColor(kBlue + 1); h_beta_norm_tot->SetLineWidth(2);
    h_beta_norm_tot->SetLineStyle(2); h_beta_norm_tot->Draw("HIST same");
    TLegend *leg_pdf = new TLegend(0.15, 0.72, 0.42, 0.88);
    leg_pdf->SetTextFont(42); leg_pdf->SetTextSize(0.030);
    leg_pdf->AddEntry(h_beta_norm,     "Slew rate", "l");
    leg_pdf->AddEntry(h_beta_norm_tot, "TOT",       "l");
    leg_pdf->Draw();
    TPaveText *pt_pdf = new TPaveText(0.60, 0.72, 0.89, 0.88, "NDC");
    pt_pdf->SetFillColor(0); pt_pdf->SetBorderSize(1); pt_pdf->SetTextAlign(12);
    pt_pdf->SetTextFont(42); pt_pdf->SetTextSize(0.030);
    pt_pdf->AddText(Form("Entries = %lld", (Long64_t)h_beta_norm->GetEntries()));
    pt_pdf->AddText("#int PDF(#beta) d#beta = 1");
    pt_pdf->Draw();
    c_beta_pdf->Write("Canvas_beta_PDF");

    // Canvas 2c: PDF β TOT separata
    TCanvas *c_beta_pdf_tot = new TCanvas("c_beta_pdf_tot",
                                          "PDF #beta normalizzata TOT (lead)", 900, 650);
    c_beta_pdf_tot->SetGrid();
    h_beta_norm_tot->SetLineStyle(1); h_beta_norm_tot->SetLineColor(kBlue + 1);
    h_beta_norm_tot->SetLineWidth(2); h_beta_norm_tot->Draw("HIST");
    TLine *l_b1_pdf_tot = new TLine(1.0, 0.0, 1.0, h_beta_norm_tot->GetMaximum() * 0.95);
    l_b1_pdf_tot->SetLineColor(kBlack); l_b1_pdf_tot->SetLineStyle(2);
    l_b1_pdf_tot->SetLineWidth(2); l_b1_pdf_tot->Draw("same");
    TLatex *lat_pdf_tot = new TLatex(1.02, h_beta_norm_tot->GetMaximum() * 0.80, "#beta = 1");
    lat_pdf_tot->SetTextSize(0.035); lat_pdf_tot->Draw("same");
    c_beta_pdf_tot->Write("Canvas_beta_PDF_TOT");

    // Canvas 3: 1/v
    TCanvas *c3 = new TCanvas("c_invv", "Distribuzione 1/v (lead)", 800, 600);
    c3->SetGrid();
    h_invv->SetLineColor(kGreen + 2); h_invv->SetLineWidth(2); h_invv->Draw();
    double inv_c = 1.0 / C_LIGHT;
    TLine *l_invc = new TLine(inv_c, 0, inv_c, h_invv->GetMaximum() * 0.9);
    l_invc->SetLineColor(kBlack); l_invc->SetLineStyle(2); l_invc->SetLineWidth(2);
    l_invc->Draw("same");
    c3->Write("Canvas_invv");

    // Canvas 4: θ
    TCanvas *c4 = new TCanvas("c_theta", "Distribuzione angolare (lead)", 800, 600);
    c4->SetGrid();
    h_theta->SetLineColor(kMagenta + 1); h_theta->SetLineWidth(2); h_theta->Draw();
    c4->Write("Canvas_theta");

    // Canvas 5: x_imp
    TCanvas *c5 = new TCanvas("c_ximp", "Posizione di impatto (lead)", 800, 600);
    c5->SetGrid();
    h_x_imp->SetLineColor(kOrange + 1); h_x_imp->SetLineWidth(2); h_x_imp->Draw();
    c5->Write("Canvas_x_imp");

    // Canvas 6: TOF vs x (2D)
    TCanvas *c6 = new TCanvas("c_tof_vs_x", "TOF vs posizione (lead)", 900, 600);
    c6->SetGrid(); h2_tof_x->Draw("COLZ"); c6->Write("Canvas_TOF_vs_x");

    // Canvas 7: β vs x (2D)
    TCanvas *c7 = new TCanvas("c_beta_vs_x", "#beta vs posizione (lead)", 900, 600);
    c7->SetGrid(); h2_beta_x->Draw("COLZ");
    TLine *l_b1 = new TLine(X_LO, 1.0, X_HI, 1.0);
    l_b1->SetLineColor(kRed); l_b1->SetLineStyle(2); l_b1->Draw("same");
    c7->Write("Canvas_beta_vs_x");

    // Canvas 8: pannello riassuntivo 2×3
    TCanvas *c_summary = new TCanvas("c_summary", "Riassunto TOF (lead)", 1200, 800);
    c_summary->Divide(3, 2);
    c_summary->cd(1); gPad->SetGrid(); h_tof->Draw();
    c_summary->cd(2); gPad->SetGrid(); h_beta->Draw();
    c_summary->cd(3); gPad->SetGrid(); h_invv->Draw();
    c_summary->cd(4); gPad->SetGrid(); h_theta->Draw();
    c_summary->cd(5); gPad->SetGrid(); h_x_imp->Draw();
    c_summary->cd(6); gPad->SetGrid(); h_path->Draw();
    c_summary->Write("Canvas_Summary");

    TCanvas *c_summary_lead = new TCanvas("c_summary_lead", "Riassunto TOF (lead)", 1200, 800);
    c_summary_lead->Divide(3, 2);
    c_summary_lead->cd(1); gPad->SetGrid(); h_tof_lead->Draw();
    c_summary_lead->cd(2); gPad->SetGrid(); h_beta_lead->Draw();
    c_summary_lead->cd(3); gPad->SetGrid(); h_invv_lead->Draw();
    c_summary_lead->cd(4); gPad->SetGrid(); h_theta_lead->Draw();
    c_summary_lead->cd(5); gPad->SetGrid(); h_x_imp_lead->Draw();
    c_summary_lead->cd(6); gPad->SetGrid(); h_path_lead->Draw();
    c_summary_lead->Write("Canvas_Summary_lead");

    // Canvas confronto: beta senza e con filtro FIFO sovrapposti
    TCanvas *c_compare = new TCanvas("c_compare", "Confronto #beta: tutti vs lead", 900, 650);
    c_compare->SetGrid();
    h_beta_norm->SetLineColor(kRed + 1);  h_beta_norm->SetLineWidth(2);
    h_beta_norm->Draw("HIST");
    h_beta_norm_lead->SetLineColor(kBlue + 1); h_beta_norm_lead->SetLineWidth(2);
    h_beta_norm_lead->SetLineStyle(2);
    h_beta_norm_lead->Draw("HIST same");
    TLegend *leg_cmp = new TLegend(0.15, 0.72, 0.50, 0.88);
    leg_cmp->SetTextFont(42); leg_cmp->SetTextSize(0.030);
    leg_cmp->AddEntry(h_beta_norm,      "Tutti gli eventi (slew)", "l");
    leg_cmp->AddEntry(h_beta_norm_lead, "Filtrati FIFO (slew)",    "l");
    leg_cmp->Draw();
    c_compare->Write("Canvas_beta_confronto");
    
    // Salva tutti gli istogrammi
    h_tof->Write();      h_beta->Write();     h_beta_norm->Write();
    h_invv->Write();     h_theta->Write();    h_x_imp->Write();
    h_path->Write();     h_T_meas->Write();   h2_tof_x->Write();
    h2_beta_x->Write();
    h_tof_tot->Write();  h_beta_tot->Write(); h_beta_norm_tot->Write();
    h_invv_tot->Write(); h_theta_tot->Write();h_x_imp_tot->Write();
    h_path_tot->Write(); h_T_meas_tot->Write();h2_tof_x_tot->Write();
    h2_beta_x_tot->Write();

    h_tof_lead->Write();      h_beta_lead->Write();     h_beta_norm_lead->Write();
    h_invv_lead->Write();     h_theta_lead->Write();    h_x_imp_lead->Write();
    h_path_lead->Write();     h_T_meas_lead->Write();   h2_tof_x_lead->Write();
    h2_beta_x_lead->Write();
    h_tof_tot_lead->Write();      h_beta_tot_lead->Write();     h_beta_norm_tot_lead->Write();
    h_invv_tot_lead->Write();     h_theta_tot_lead->Write();    h_x_imp_tot_lead->Write();
    h_path_tot_lead->Write();     h_T_meas_tot_lead->Write();   h2_tof_x_tot_lead->Write();
    h2_beta_x_tot_lead->Write();

    fout->Close();

    // Riepilogo finale
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  ANALISI TOF LEAD COMPLETATA                 " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  File XML processati:    " << file_list.size() << std::endl;
    std::cout << "  File FIFO:             " << fifo_file << std::endl;
    std::cout << "  Reset FPGA trovati:    " << n_resets << std::endl;
    std::cout << "  Eventi totali DRS:     " << total_events << std::endl;
    std::cout << "  Filtrati dal FIFO:     " << rej_fifo
              << " (" << pct(rej_fifo) << "%)" << std::endl;
    std::cout << "  Eventi buoni (slew):   " << total_good << std::endl;
    std::cout << "  Eventi buoni (TOT):    " << total_good_tot << std::endl;
    std::cout << "  Eventi buoni lead (slew): " << (long)h_tof_lead->GetEntries() << std::endl;
    std::cout << "  Eventi buoni lead (TOT):  " << (long)h_tof_tot_lead->GetEntries() << std::endl;
    std::cout << "  Output: " << outname << std::endl;
    std::cout << "=============================================" << std::endl;
}
