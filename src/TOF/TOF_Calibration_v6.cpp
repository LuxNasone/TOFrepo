// ==========================================================================
//  TOF_Calibration.cpp — Calibrazione dell'apparato Time-of-Flight
// ==========================================================================
//
//  SCOPO:
//    Determinare la velocità effettiva della luce nella barra (v_eff) e la
//    costante di offset C necessaria per calcolare il tempo di volo.
//
//  FISICA:
//    La barra scintillatrice (BICRON BC408, 280 cm) è letta da PMT1 e PMT2
//    ai due estremi. Un terzo scintillatore mobile (PMT3) è posizionato in
//    contatto con la barra (Guida A) in diverse posizioni note x_k.
//
//    Definendo x la coordinata lungo la barra dal centro (x ∈ [-140, 140] cm),
//    i tempi di arrivo del segnale ai tre PMT sono:
//
//      t₁ = t₀ + d₁(x)/v + δ₁       (d₁ = distanza scintillazione → PMT1)
//      t₂ = t₀ + d₂(x)/v + δ₂       (d₂ = distanza scintillazione → PMT2)
//      t₃ = t₀ + TOF      + δ₃       (TOF ≈ 0 in configurazione Guida A)
//
//    dove t₀ è l'istante della scintillazione, δᵢ i ritardi strumentali.
//
//    1) RETTA DI CALIBRAZIONE: Δt₁₂ = t₁ - t₂ = (2/v)·x + (δ₁ − δ₂) + cost
//       → il fit lineare Δt₁₂ vs x dà:   pendenza m = ±2/v   →   v = 2/|m|
//       → l'intercetta q = offset relativo PMT1-PMT2
//
//    2) COSTANTE C:  C = t₃ − (t₁+t₂)/2
//       Nella configurazione di calibrazione (PMT3 a contatto, TOF ≈ 0):
//         C = δ₃ − (δ₁+δ₂)/2 − L/(2v)
//       Questa costante è indipendente da x (entro le incertezze).
//       In fase di misura TOF:   TOF = t₃ − (t₁+t₂)/2 − C
//
//  UTILIZZO:
//    root -l 'TOF_Calibration.cpp("/percorso/cartella/dati/")'
//
//    I file XML devono avere nomi nella forma "N130.xml", "X0.xml", ecc.
//    dove N indica posizioni negative e X positive (in cm dal centro barra).
//
//  OUTPUT:
//    - File ROOT con TTree per ogni posizione (dati evento per evento)
//    - Istogrammi fittati delle distribuzioni Δt₁₂ e C per ogni posizione
//    - Canvas con retta di calibrazione Δt₁₂ vs x (+ pannello residui)
//    - Canvas con C vs x (verifica uniformità)
//    - TTree riassuntivo con i risultati dei fit
//
//  DIPENDENZE: ROOT 6+ (usa TTree, TH1D, TF1, TCanvas, TGraphErrors)
//
// 
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
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
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
//  SEZIONE 1: COSTANTI E PARAMETRI CONFIGURABILI
// ==========================================================================

// --- Parametri hardware DRS4 ---
const int MAX_SAMPLES  = 1024;   // Celle del chip DRS4 (fisso, NON modificare)
const int MAX_CHANNELS = 4;     // Canali della DRS4 Eval Board

// --- Parametri di analisi della forma d'onda ---
const int    NBL_SAMPLES   = 50;    // Campioni per il calcolo della baseline (~10 ns a 5 GSPS)
const double CFD_FRACTION  = 0.2;   // Frazione CFD (0.5 = soglia al 50% dell'ampiezza)
const double NOISE_THRESH  = 5.0;   // Soglia minima ampiezza [mV] per dichiarare un impulso
// --- Soglie di qualità per scartare eventi non fisici ---
//
// Due flag indipendenti, ciascuno con la propria firma fisica:
//
//   is_clipped     ← V_min < CLIP_V_LO (saturazione del bottom della DRS4)
//   is_oscillating ← V_max > OSC_VMAX_THRESH (escursione positiva anomala)
//
// Un impulso fisico di scintillazione è MONOPOLARE NEGATIVO: V_max è
// vicina alla baseline (pochi mV positivi al massimo, anche con rumore
// e overshoot fisiologico). Una V_max significativamente positiva è
// la firma più diretta di un evento patologico (afterpulse, EMI,
// breakdown, ringing).

// CLIP_V_LO: tensione minima [mV] sotto cui il segnale è saturato.
//   La DRS4 satura tipicamente a circa −500 mV. Soglia con margine.
//   ATTENZIONE: aggiornare se il range della DRS4 è impostato diversamente.
const double CLIP_V_LO       = -450.0;

// OSC_VMAX_THRESH: tensione massima [mV] sopra cui il segnale è
//   considerato oscillante. Soglia conservativa: ~10 σ del rumore
//   baseline (tipicamente 2-5 mV RMS) e ~20× il V_max di un impulso
//   fisico pulito (inclusi rumore e overshoot fisiologico).
const double OSC_VMAX_THRESH =   100.0;

// Flag globali per abilitare/disabilitare i tagli (utile per analisi
// comparative). Default: tutti abilitati.
const bool ENABLE_CFD_CUT  = true;   // scarta se nessun tempo utilizzabile
                                     //  (né CFD standard né CFD ricostruito)
const bool ENABLE_CLIP_CUT = false;  // FALSE = il clipping NON causa scarto
                                     //   diretto; gli eventi clippati vengono
                                     //   recuperati via slew rate quando
                                     //   possibile, e scartati altrimenti
                                     //   tramite la verifica di all_t_ok.
                                     //  TRUE  = comportamento "puristico":
                                     //   scarta tutti gli eventi clippati,
                                     //   indipendentemente dal recupero
                                     //   (utile per confronti di sistematica).
const bool ENABLE_OSC_CUT  = true;   // scarta se V_max > OSC_VMAX_THRESH
// ============================================================
//  PARAMETRI PER IL RECUPERO DEI SEGNALI CLIPPATI
//  (calibrazione "ampiezza ↔ slew rate" e CFD ricostruito)
// ============================================================
//
// Strategia: per gli impulsi non clippati, ampiezza A e slew rate s
// del fronte di salita sono legate dalla forma normalizzata g(t):
//
//      v(t) = A · g(t − t₀)   ⇒   dv/dt = A · g'(t − t₀)
//
// Misurando s in una finestra fissa di tensione (lontana sia dal
// rumore di baseline sia dal livello di clipping), s è proporzionale
// ad A: A = k · s, con k caratteristica del PMT (e della sua catena).
//
// Una volta determinato k da eventi NON clippati, lo applichiamo agli
// eventi clippati per stimare A_rec = k · s e quindi il tempo CFD
// usando la frazione CFD_FRACTION sull'ampiezza ricostruita.

// --- Finestra di tensione per il fit lineare del fronte di salita ---
// Espressa come u = baseline − v (cioè ampiezza istantanea dalla baseline,
// sempre positiva). Limiti scelti per:
//   - SLEW_U_LO sopra il rumore baseline (~3-5 mV RMS) e sopra il piede
//     non lineare del leading edge (rampa iniziale curva).
//   - SLEW_U_HI ben sotto |CLIP_V_LO| (= 450 mV), così la finestra è
//     SEMPRE non clippata anche per impulsi enormi.
const double SLEW_U_LO     = 30.0;   // [mV] limite inferiore della finestra
const double SLEW_U_HI     = 250.0;  // [mV] limite superiore della finestra

// Numero minimo di campioni nella finestra per accettare il fit lineare.
// A 5 GS/s (0.2 ns/sample) e rise time ~3 ns, la finestra contiene
// tipicamente 8-12 campioni per impulsi grandi, di più per impulsi medi.
// 5 è una soglia di sicurezza per garantire stabilità del fit.
const int    SLEW_NMIN     = 5;

// Soglia di qualità sul χ²/ndf del fit lineare. Eventi con χ²/ndf > soglia
// vengono scartati dalla calibrazione e dalla ricostruzione. Il χ² è
// calcolato usando come errore tipico del campione max(bl_rms, SLEW_SIGMA_FLOOR).
const double SLEW_CHI2_MAX = 3.0;

// Floor per σ_y nel calcolo del χ²: evita che eventi con baseline RMS
// anomalamente bassa (sotto 2 mV) abbiano χ² gonfiati artificialmente
// dalla pura curvatura fisica del fronte di salita.
const double SLEW_SIGMA_FLOOR = 2.0; // [mV]

// Margine sopra |CLIP_V_LO| sotto cui la soglia CFD ricostruita è
// considerata "troppo vicina al clipping" e l'evento viene scartato.
// Es: se |CLIP_V_LO| = 450 mV e margine = 20 mV, accettiamo solo
// soglie CFD nella regione v_thr > -430 mV (cioè u_thr < 430 mV).
const double SLEW_CLIP_MARGIN = 20.0;  // [mV]

// --- Range di posizioni usate per la calibrazione di k ---
// Solo dataset con |x_k| <= CALIB_K_X_ABS_MAX contribuiscono al fit
// di k. Giustificazione: a queste posizioni, entrambi i PMT vedono
// abbondanti eventi NON clippati su un range di ampiezze sufficiente
// per determinare la pendenza A = k·s. Posizioni più estreme (±112,
// ±130) sono dominate da clipping su un PMT e quindi non utili per
// la calibrazione (anche se sono i target dell'applicazione).
const double CALIB_K_X_ABS_MAX = 84.0;  // [cm]


// ============================================================
//  COSTANTI DI CALIBRAZIONE k_PMT (popolate da CalibrateSlewK)
// ============================================================
// Indici 0,1,2 = PMT1, PMT2, PMT3. PMT3 è quello mobile sopra la barra
// in configurazione di calibrazione e tipicamente NON è clippato
// (vede solo i fotoni "diretti" dello scintillatore di trigger),
// quindi non viene calibrato per il recupero — il suo k resta a zero.
//
// Le variabili sono "static" a scope di file: persistono tra chiamate
// di funzioni dentro questo .cpp ma non sono visibili da altri file.
static double gK_PMT[MAX_CHANNELS]      = {0.0, 0.0, 0.0, 0.0};   // k [mV / (mV/ns)] = [ns]
static double gK_PMT_err[MAX_CHANNELS]  = {0.0, 0.0, 0.0, 0.0};   // errore di k
static bool   gK_calibrated             = false;                  // flag stato globale

// --- Parametri di binning istogrammi ---
const int    NBINS_DT  = 120;       // Numero di bin per gli istogrammi Δt
const double DT_LO     = -30.0;     // Limite inferiore Δt [ns]
const double DT_HI     = 30.0;      // Limite superiore Δt [ns]
const int    NBINS_C   = 120;       // Numero di bin per gli istogrammi di C
const double C_LO      = -30.0;     // Limite inferiore C [ns]
const double C_HI      = 0.0;       // Limite superiore C [ns]

// --- Errore sulla posizione x ---
// Metà dell'estensione del PMT3 nella direzione parallela alla barra.
//  3 cm  Modificare se PMT3 ha dimensioni diverse.
const double DX_POSITION = 3.0;     // Incertezza sulla posizione [cm]
// ============================================================
//  PARAMETRI ERRORE SULLA POSIZIONE (modello posizione-dipendente)
// ============================================================
//
// L'errore totale su x_k è la somma in quadratura di due termini:
//
//   σ_x² = σ_parallax²(x_k) + σ_tape²(x_k)
//
// 1) σ_parallax(x_k) — viene letto dal file di output del Monte Carlo
//    di accettanza (TOF_AcceptanceMC.cpp::ParallaxScan), che fornisce
//    σ stocastica della posizione di scintillazione nella barra dato
//    x_PMT3 letto col metro. Include sia la dimensione finita di PMT3
//    (~ℓ/√12) che la parallasse geometrica dovuta agli angoli ~cos^α(θ).
//
// 2) σ_tape(x_k) — errore di lettura del metro a nastro, modellato
//    in due modi alternativi (selezionabili via TAPE_MODEL_INCREMENTAL):
//
//    ► MODALITÀ A "lettura singola" (default, conservativa):
//      Ogni posizione ≠ 0 ha errore di lettura costante σ_read.
//      La posizione 0 è il riferimento → errore nullo.
//
//    ► MODALITÀ B "letture incrementali":
//      Si aggiungono σ_read² in quadratura per ogni "passo" del
//      metro percorso dal riferimento (es. ogni 28 cm).
//      Es. a x_k = 84 cm: σ_tape = σ_read · √3.

// Errore di lettura del metro a nastro per ogni "evento" di lettura.
// Tipicamente 0.1 cm = 1 mm per un metro di buona qualità.
const double SIGMA_TAPE_READ = 0.1;     // [cm]

// Posizione di riferimento del metro (errore nullo a questa posizione)
const double TAPE_REFERENCE_X = 0.0;    // [cm]

// Modalità di accumulo dell'errore del metro:
//   false → modalità A "lettura singola" (1 mm per qualunque x ≠ 0)
//   true  → modalità B "letture incrementali" (1 mm × √n_step)
const bool TAPE_MODEL_INCREMENTAL = false;

// Lunghezza di un singolo passo di posizionamento del metro [cm].
// Usata SOLO in modalità B per calcolare il numero di passi da 0 a x_k.
// Per il vostro setup le posizioni sono multipli di 28 cm.
const double TAPE_STEP_LENGTH = 28.0;   // [cm]


// ============================================================
//  FILE DI PARALLASSE (output del MC TOF_AcceptanceMC.cpp)
// ============================================================
//
// Path al file ROOT prodotto da ParallaxScan(...) per la
// configurazione di calibrazione (Guida A, PMT3 sopra). Contiene:
//   - TGraph "g_sigma_x" con σ_parallax vs x_PMT3
//   - TGraph "g_bias_x"   con bias vs x_PMT3
//   - TTree  "parallax"
//
// Se il file non è disponibile (stringa vuota o file mancante),
// il codice usa il valore costante PARALLAX_FALLBACK come errore,
// emettendo un warning. Il default è il valore osservato al centro
// della distribuzione σ_x dal MC in Guida A (~3.95 cm).

const char*  PARALLAX_FILE      = "";        // path file ROOT MC (vuoto = disabilitato)
const double PARALLAX_FALLBACK  = 3.95;      // [cm] valore costante se file non disponibile
// --- Lista posizioni di calibrazione ---
// Formato: {"nome_file_senza_estensione", posizione_in_cm_dal_centro}
// Modificare in base alle posizioni effettivamente acquisite.
struct PositionInfo {
    std::string name;    // Nome del file XML (senza .xml)
    double      x_cm;    // Posizione dal centro della barra [cm]
};

// Posizioni predefinite (modificare se necessario)
std::vector<PositionInfo> DEFAULT_POSITIONS = {
    {"N130", -130.0},
    {"N112", -112.0},
    {"N84",   -84.0},
    {"N56",   -56.0},
    {"N28",   -28.0},
    {"X0",      0.0},
    {"X28",    28.0},
    {"X56",    56.0},
    {"X84",    84.0},
    {"X112",  112.0},
    {"X130",  130.0}
};


// ==========================================================================
//  SEZIONE 2: STRUTTURE DATI
// ==========================================================================

/// Dati di un singolo canale di un singolo evento.
/// Contiene la forma d'onda grezza e le grandezze estratte dall'analisi.
/// Dati di un singolo canale di un singolo evento.
/// Contiene la forma d'onda grezza e le grandezze estratte dall'analisi.
struct ChannelData {
    int    nsamples;                  // Numero di campioni letti (tipicamente 1024)
    float  time[MAX_SAMPLES];        // Tempi dei campioni [ns], calibrati dalla DRS4
    float  voltage[MAX_SAMPLES];     // Tensione dei campioni [mV]

    // Grandezze estratte da AnalyzeChannel():
    double baseline;      // Media tensione nei primi NBL_SAMPLES [mV]
    double baseline_rms;  // Deviazione standard della baseline [mV]
    double amplitude;     // Ampiezza: baseline - V_min [mV] (positiva)
    double v_min;         // Tensione minima (picco negativo) [mV]
    double v_max;         // Tensione massima [mV] — usato per rilevare bipolarità
    double t_min;         // Tempo del campione con V minima [ns]
    double t_cfd;         // Tempo CFD interpolato [ns] (-999 se non valido)

    // Flag di qualità:
    bool   has_pulse;       // true se amplitude > NOISE_THRESH
    bool   cfd_ok;          // true se il CFD ha trovato un crossing valido
    bool   is_clipped;      // true se V_min<CLIP_V_LO o V_max>CLIP_V_HI (saturazione)
    bool   is_oscillating;  // true se troppi campioni "grandi" (impulso non fisico)

    // ---- Grandezze del fit lineare dello slew rate (PASSO 7) ----
    // Calcolate per OGNI evento (clippato o no) per permettere sia la
    // calibrazione di k (sui non clippati) sia il recupero dei clippati.
    double slew_rate;       // pendenza b del fit lineare u(t) = a + b·t [mV/ns]
                            //   u = baseline − voltage (positivo per impulsi negativi)
    double slew_rate_err;   // errore di b dal fit (formula OLS standard) [mV/ns]
    double slew_chi2_ndf;   // χ²/ndf del fit (con σ_y = max(bl_rms, SLEW_SIGMA_FLOOR))
    int    slew_n_pts;      // numero di campioni nella finestra [SLEW_U_LO, SLEW_U_HI]
    bool   slew_ok;         // true se: n_pts ≥ SLEW_NMIN, slope > 0,
                            //          χ²/ndf ∈ (0, SLEW_CHI2_MAX)

    // ---- Quantità di recupero per eventi clippati (PASSO 8) ----
    // Popolate solo da RecoverClippedCFD() durante ProcessDataset(),
    // dopo che gK_PMT[] è stato calibrato.
    double amplitude_rec;   // ampiezza ricostruita: A_rec = k_PMT · slew_rate [mV]
    double t_cfd_recovered; // tempo CFD calcolato sulla soglia f·A_rec [ns]
    bool   cfd_recovered;   // true se il CFD ricostruito è valido e usabile
};

/// Dati completi di un evento DRS4 (header + canali).
struct EventData {
    int         serial;                    // Numero progressivo (parte da 1)
    std::string timestamp;                 // Data/ora "YYYY/MM/DD HH:MM:SS.mmm"
    int         trigger_cell;              // Cella dove l'onda domino si è fermata
    int         board_serial;              // Numero seriale della scheda
    int         scaler[MAX_CHANNELS];      // Rate scaler per canale [Hz]
    int         nchannels;                 // Numero canali presenti (1–4)
    int         channel_ids[MAX_CHANNELS]; // Numeri dei canali (1-based: CHN1=1, ecc.)
    ChannelData ch[MAX_CHANNELS];          // Dati dei canali (0-based)
};

/// Risultati del fit di una distribuzione (Δt₁₂ o C) per una singola posizione.
/// Usato per assemblare i punti della retta di calibrazione.
struct FitResult {
    double center;       // Centro della distribuzione dal fit [ns]
    double center_err;   // Errore sul centro [ns]
    double width;        // Larghezza (FWHM per Breit-Wigner) [ns]
    double width_err;    // Errore sulla larghezza [ns]
    double chi2_ndf;     // Chi²/ndf del fit (diagnostica)
    int    nentries;     // Numero di eventi nell'istogramma
    bool   fit_ok;       // true se il fit è convergito con risultati sensati
};


// ==========================================================================
//  SEZIONE 3: ANALISI DELLA FORMA D'ONDA
// ==========================================================================

/// AnalyzeChannel(): estrae da un canale baseline, ampiezza, tempo CFD,
/// e i flag di qualità is_clipped e is_oscillating.
///
/// FLAG DI QUALITÀ:
///   - is_clipped: il segnale ha raggiunto il range della DRS4. La V_min
///     non è il vero picco ma il valore di saturazione, e la CFD calcolata
///     su un fronte di discesa "tagliato" è bias-ata. Anche un V_max grande
///     (segnale bipolare) indica problemi (afterpulse, EMI).
///   - is_oscillating: il segnale ha troppi campioni "grandi" rispetto a un
///     impulso fisico. Tipico di scariche nel PMT, interferenze EMI, o
///     saturazione dell'elettronica che genera un'oscillazione.
///
/// L'algoritmo CFD viene eseguito comunque (per conservare diagnostica nel
/// TTree), ma il loop principale di ProcessDataset() scarterà gli eventi
/// con flag positivi.

void AnalyzeChannel(ChannelData &cd) {

    // --- Inizializzazione a valori "non calcolato" ---
    cd.has_pulse      = false;
    cd.cfd_ok         = false;
    cd.is_clipped     = false;
    cd.is_oscillating = false;
    cd.baseline       = 0.0;
    cd.baseline_rms   = 0.0;
    cd.amplitude      = 0.0;
    cd.v_min          = 0.0;
    cd.v_max          = 0.0;
    cd.t_min          = 0.0;
    cd.t_cfd          = -999.0;

    // Inizializzazione delle grandezze di slew rate e di recupero clipping
    cd.slew_rate       = 0.0;
    cd.slew_rate_err   = 0.0;
    cd.slew_chi2_ndf   = -1.0;     // -1 = "non calcolato" (qualunque valore < 0)
    cd.slew_n_pts      = 0;
    cd.slew_ok         = false;
    cd.amplitude_rec   = 0.0;
    cd.t_cfd_recovered = -999.0;
    cd.cfd_recovered   = false;

    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return;  // Troppo pochi campioni

    float *t = cd.time;
    float *v = cd.voltage;

    // ---- PASSO 1: BASELINE ----
    // Media e deviazione standard dei primi NBL_SAMPLES campioni.
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < NBL_SAMPLES; i++) {
        sum  += v[i];
        sum2 += (double)v[i] * v[i];
    }
    double bl     = sum / NBL_SAMPLES;
    double bl_rms = sqrt(fabs(sum2 / NBL_SAMPLES - bl * bl));
    cd.baseline     = bl;
    cd.baseline_rms = bl_rms;

    // ---- PASSO 2: RICERCA DEL MINIMO E MASSIMO GLOBALI ----
    // Il minimo è il picco negativo dell'impulso. Il massimo serve per
    // rilevare segnali bipolari (oscillanti) che non sono impulsi fisici.
    double vmin = v[0], vmax = v[0];
    int    imin = 0;
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
        if (v[i] > vmax)   vmax = v[i];
    }
    cd.v_min = vmin;
    cd.v_max = vmax;
    cd.t_min = t[imin];

// ---- PASSO 3: RILEVAMENTO CLIPPING (saturazione DRS4) ----
    // V_min vicino al limite inferiore del range significa che il segnale
    // è saturato: il vero picco è stato "tagliato" dalla saturazione e la
    // CFD calcolata sul fronte di discesa è sistematicamente bias-ata.
    if (vmin < CLIP_V_LO) {
        cd.is_clipped = true;
    }

    // ---- PASSO 4: AMPIEZZA E SOGLIA DI RUMORE ----
    double amp = bl - vmin;
    cd.amplitude = amp;

    if (amp < NOISE_THRESH) return;  // Sotto soglia: nessun impulso reale
    cd.has_pulse = true;

// ---- PASSO 5: RILEVAMENTO OSCILLAZIONE (escursione positiva anomala) ----
    // Un impulso di scintillazione del PMT è monopolare negativo per
    // costruzione. V_max sostanzialmente positiva (≥ 50 mV) non ha
    // spiegazione fisica e indica un evento patologico: afterpulse,
    // EMI captata sui cavi, breakdown nel tubo, oscillazioni elettroniche.
    // L'oscillazione patologica può manifestarsi anche senza saturare
    // la DRS4 (V_min può essere benissimo > CLIP_V_LO), quindi questo
    // criterio è indipendente dal flag is_clipped.
    if (vmax > OSC_VMAX_THRESH) {
        cd.is_oscillating = true;
    }

    // ---- PASSO 6: CFD (Constant Fraction Discriminator) ----
    // Soglia CFD: una frazione dell'ampiezza sotto la baseline.
    //   v_thr = bl + CFD_FRACTION * (vmin - bl) = bl - CFD_FRACTION * amp
    //
    // Ricerca IN AVANTI dalla fine della baseline al minimo, primo intervallo
    // [i, i+1] dove la tensione attraversa la soglia verso il basso. Trovato
    // il crossing, interpolazione lineare sub-campione.
    //
    // La ricerca in avanti (e non backward dal minimo) garantisce di
    // trovare il crossing sul vero leading edge, evitando falsi match
    // su strutture parassite.
    double v_thr = bl + CFD_FRACTION * (vmin - bl);

    for (int i = NBL_SAMPLES; i < imin; i++) {
        if (v[i] > v_thr && v[i + 1] <= v_thr) {
            double dv = (double)v[i + 1] - v[i];
            if (fabs(dv) > 1e-6) {
                cd.t_cfd  = t[i] + (v_thr - v[i]) / dv * (t[i + 1] - t[i]);
                cd.cfd_ok = true;
            }
            break;
        }
    }

    // ---- PASSO 7: FIT LINEARE DELLO SLEW RATE ----
    // Obiettivo: misurare la pendenza b del fronte di salita in una
    // finestra di tensione fissa, lontana sia dal rumore baseline che
    // dal livello di clipping.
    //
    // Definizione: u(t) = baseline − v(t)  [sempre positiva durante il
    // fronte di salita di un impulso negativo].
    //
    // Algoritmo:
    //   1. Iteriamo sui campioni dalla fine della baseline al picco
    //      (i = NBL_SAMPLES … imin) — NON oltre il picco, dove inizia
    //      il fronte di discesa o (se clippato) il plateau saturato.
    //   2. Selezioniamo i campioni con SLEW_U_LO ≤ u_i ≤ SLEW_U_HI.
    //   3. Se ne abbiamo almeno SLEW_NMIN, fittiamo u = a + b·t con
    //      regressione lineare ordinaria (OLS), formula chiusa.
    //   4. Calcoliamo χ² assumendo errore σ_y = max(bl_rms, SLEW_SIGMA_FLOOR)
    //      uguale per tutti i campioni (il floor evita χ² gonfiato per
    //      eventi con bl_rms anomalamente piccolo).
    //   5. Validiamo: pendenza positiva e χ²/ndf < SLEW_CHI2_MAX.
    //
    // Questo passo viene eseguito SEMPRE (clippato o non clippato),
    // perché lo slew rate dei NON clippati serve per calibrare k,
    // e quello dei CLIPPATI serve per il recupero in PASSO 8.
    {
        // Raccolta campioni: scorriamo solo il fronte di salita.
        // Limite superiore: imin (campione del minimo, cioè del picco).
        // Per eventi clippati il "picco" è il primo campione del
        // plateau saturato — ma nel codice attuale imin è comunque
        // dentro il plateau, e i campioni precedenti coprono il fronte.
        std::vector<double> t_pts;
        std::vector<double> u_pts;
        t_pts.reserve(32);   // pre-allocazione: tipicamente <30 campioni
        u_pts.reserve(32);

        for (int i = NBL_SAMPLES; i <= imin; i++) {
            double u = bl - (double)v[i];     // ampiezza istantanea dalla baseline
            if (u >= SLEW_U_LO && u <= SLEW_U_HI) {
                t_pts.push_back((double)t[i]);
                u_pts.push_back(u);
            }
        }

        cd.slew_n_pts = (int)t_pts.size();

        // Procediamo al fit solo se abbiamo abbastanza campioni
        if (cd.slew_n_pts >= SLEW_NMIN) {

            // Somme per la regressione lineare OLS (formule chiuse).
            // Riferimento: Bevington & Robinson, "Data Reduction and
            // Error Analysis", cap. 6, eq. (6.13)-(6.16).
            int    N      = cd.slew_n_pts;
            double sum_t  = 0.0;
            double sum_u  = 0.0;
            double sum_tt = 0.0;   // Σ t_i²
            double sum_tu = 0.0;   // Σ t_i · u_i

            for (int j = 0; j < N; j++) {
                sum_t  += t_pts[j];
                sum_u  += u_pts[j];
                sum_tt += t_pts[j] * t_pts[j];
                sum_tu += t_pts[j] * u_pts[j];
            }

            // Determinante del sistema normale.
            // D = N·Σt² − (Σt)² = N · Var(t) · N (a meno di costanti).
            // Se D è degenero (tutti i t_i uguali, impossibile in pratica),
            // saltiamo il fit.
            double D = (double)N * sum_tt - sum_t * sum_t;

            if (fabs(D) > 1e-12) {

                // Coefficienti del fit: u = a + b·t
                double slope     = ((double)N * sum_tu - sum_t * sum_u) / D;
                double intercept = (sum_u - slope * sum_t) / (double)N;

                // Errore tipico del campione (con floor per stabilità)
                double sigma_y = std::max(bl_rms, SLEW_SIGMA_FLOOR);

                // Errore della pendenza dalla formula OLS:
                //   σ_b² = σ_y² · N / D
                double slope_err = sigma_y * sqrt((double)N / D);

                // Calcolo χ² rispetto al modello lineare
                double chi2 = 0.0;
                for (int j = 0; j < N; j++) {
                    double resid = u_pts[j] - (intercept + slope * t_pts[j]);
                    chi2 += (resid * resid) / (sigma_y * sigma_y);
                }
                int    ndf     = N - 2;     // 2 parametri liberi (a, b)
                double chi2ndf = (ndf > 0) ? chi2 / (double)ndf : -1.0;

                // Salvataggio nei campi della struct
                cd.slew_rate     = slope;
                cd.slew_rate_err = slope_err;
                cd.slew_chi2_ndf = chi2ndf;

                // Validazione: pendenza fisica (positiva) e χ²/ndf accettabile
                if (slope > 0.0 && chi2ndf > 0.0 && chi2ndf < SLEW_CHI2_MAX) {
                    cd.slew_ok = true;
                }
            }
        }
    }
}
// ==========================================================================
//  SEZIONE 3b: RECUPERO DEGLI EVENTI CLIPPATI VIA SLEW RATE
// ==========================================================================

/// RecoverClippedCFD():
///   Per un canale clippato di un dato PMT, prova a ricostruire il tempo CFD
///   usando l'ampiezza stimata da A_rec = k_PMT · slew_rate.
///
/// PRECONDIZIONI:
///   - cd.is_clipped == true (altrimenti non c'è nulla da recuperare;
///                            usa cd.t_cfd standard).
///   - cd.slew_ok == true   (il fit dello slew rate ha buona qualità).
///   - gK_calibrated == true e gK_PMT[ch_index] > 0 (k disponibile).
///
/// COMPORTAMENTO:
///   - Calcola A_rec = k · slew_rate.
///   - Calcola la soglia CFD: v_thr = baseline − CFD_FRACTION · A_rec.
///   - Verifica che v_thr stia ABBONDANTEMENTE sopra il livello di clipping
///     (margine SLEW_CLIP_MARGIN). Se cade nella regione saturata, non c'è
///     un crossing fisico misurabile sulla forma d'onda → NON si recupera.
///   - Cerca il crossing sul fronte di salita (campioni in cui la tensione
///     scende attraverso v_thr). Usa interpolazione lineare sub-campione.
///   - Salva il risultato in cd.t_cfd_recovered e setta cd.cfd_recovered.
///
/// PARAMETRI:
///   cd       — riferimento al ChannelData da recuperare
///   ch_index — indice del PMT (0=PMT1, 1=PMT2, 2=PMT3) per leggere gK_PMT[]
///
/// NOTE IMPORTANTI:
///   - PMT3 (ch_index=2) di norma NON è clippato in fase di calibrazione,
///     quindi gK_PMT[2] resta a 0 e questa funzione ritorna senza fare nulla.
///   - In caso di mancata calibrazione (gK_calibrated == false), la funzione
///     non fa nulla; cd.cfd_recovered resta false e l'evento sarà scartato
///     dai tagli di qualità a valle (perché is_clipped è true).

void RecoverClippedCFD(ChannelData &cd, int ch_index) {

    // --- Sanity check delle precondizioni ---
    if (!cd.is_clipped)                       return;  // evento non clippato: niente da fare
    if (!cd.slew_ok)                          return;  // slew rate non affidabile
    if (!gK_calibrated)                       return;  // calibrazione non fatta
    if (ch_index < 0 || ch_index >= MAX_CHANNELS) return; // indice fuori range
    if (gK_PMT[ch_index] <= 0.0)              return;  // k non disponibile per questo canale

    // --- Stima dell'ampiezza ricostruita ---
    // A_rec = k · s, con k in [ns] e s in [mV/ns] → A_rec in [mV].
    double A_rec = gK_PMT[ch_index] * cd.slew_rate;
    cd.amplitude_rec = A_rec;

    // Sanità fisica: A_rec deve essere maggiore dell'ampiezza apparente
    // (clippata). Se non lo è, qualcosa non torna — non recuperare.
    if (A_rec <= cd.amplitude) return;

    // --- Calcolo della soglia CFD ricostruita ---
    // v_thr [mV] = baseline − CFD_FRACTION · A_rec.
    // Per impulsi negativi, v_thr < baseline.
    double v_thr = cd.baseline - CFD_FRACTION * A_rec;

    // --- Verifica che la soglia stia SOPRA il livello di clipping ---
    // Se v_thr ≤ CLIP_V_LO + margine, la soglia cade nella regione
    // saturata (o troppo vicina al limite per essere affidabile) e
    // un crossing fisico non esiste sulla forma d'onda.
    if (v_thr <= CLIP_V_LO + SLEW_CLIP_MARGIN) return;

    // --- Ricerca del crossing sul fronte di salita ---
    // Stessa logica del CFD originale (in AnalyzeChannel PASSO 6):
    // ricerca IN AVANTI dalla fine della baseline al picco.
    // Usiamo i puntatori float interni alla struct.
    int   ns = cd.nsamples;
    float *t = cd.time;
    float *v = cd.voltage;

    // Trova l'indice del minimo (= primo campione del plateau saturato
    // se clippato). Lo ricalcoliamo qui per pulizia, anche se è già
    // disponibile come cd.t_min — preferiamo usare gli indici grezzi.
    double vmin = v[0];
    int    imin = 0;
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
    }

    for (int i = NBL_SAMPLES; i < imin; i++) {
        // Cerca il primo intervallo [i, i+1] dove v scende attraverso v_thr
        if (v[i] > v_thr && v[i + 1] <= v_thr) {
            double dv = (double)v[i + 1] - v[i];
            if (fabs(dv) > 1e-6) {
                // Interpolazione lineare sub-campione
                cd.t_cfd_recovered = t[i] + (v_thr - v[i]) / dv * (t[i + 1] - t[i]);
                cd.cfd_recovered   = true;
            }
            break;
        }
    }
}
// ==========================================================================
//  SEZIONE 4: PARSING XML DRS4
// ==========================================================================

/// ParseXML(): legge un file XML della DRS4 e popola il vettore di eventi.
///
/// Il parser usa una macchina a stati finiti con 4 stati:
///   IDLE       → attesa di un <Event>
///   IN_EVENT   → lettura header (Serial, Time)
///   IN_BOARD   → lettura dati scheda (Trigger_Cell, Scaler, canali)
///   IN_CHANNEL → lettura campioni <Data>tempo,tensione</Data>
///
/// Restituisce il numero di eventi letti, o -1 in caso di errore.

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
    int ch_idx       = -1;   // Indice canale corrente in ch[] (0-based)
    int sample_count = 0;    // Contatore campioni nel canale
    int scaler_idx   = 0;    // Contatore progressivo scaler
    std::string line;
    char buf[512];

    while (std::getline(infile, line)) {

        // Rimuovi '\r' se line ending Windows
        if (!line.empty() && line.back() == '\r') line.pop_back();

        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';

        // --- STATO IDLE: attesa di <Event> ---
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

        // --- STATO IN_EVENT: lettura header ---
        if (state == IN_EVENT) {
            int val;
            char ts[64];

            if (sscanf(buf, " <Serial>%d</Serial>", &val) == 1) {
                current_evt.serial = val;
            }
            else if (sscanf(buf, " <Time>%63[^<]</Time>", ts) == 1) {
                current_evt.timestamp = ts;
            }
            else if (sscanf(buf, " <Board_%d>", &val) == 1) {
                current_evt.board_serial = val;
                state = IN_BOARD;
            }
            else if (line.find("</Event>") != std::string::npos) {
                // Evento completo: analizza tutti i canali
                for (int i = 0; i < current_evt.nchannels; i++) {
                    AnalyzeChannel(current_evt.ch[i]);
                }
                events.push_back(current_evt);
                state = IDLE;
            }
            continue;
        }

        // --- STATO IN_BOARD: lettura dati scheda ---
        if (state == IN_BOARD) {
            int val, val2;

            if (sscanf(buf, " <Trigger_Cell>%d</Trigger_Cell>", &val) == 1) {
                current_evt.trigger_cell = val;
            }
            else if (sscanf(buf, " <Scaler%d>%d</Scaler", &val, &val2) == 2) {
                if (scaler_idx < MAX_CHANNELS) {
                    current_evt.scaler[scaler_idx] = val2;
                }
                scaler_idx++;
            }
            else if (sscanf(buf, " <CHN%d>", &val) == 1) {
                int idx = current_evt.nchannels;
                if (idx < MAX_CHANNELS) {
                    current_evt.channel_ids[idx] = val;
                    current_evt.ch[idx].nsamples = 0;
                    ch_idx = idx;
                    current_evt.nchannels++;
                    sample_count = 0;
                    state = IN_CHANNEL;
                }
            }
            else if (line.find("</Board_") != std::string::npos) {
                state = IN_EVENT;
            }
            continue;
        }

        // --- STATO IN_CHANNEL: lettura campioni ---
        if (state == IN_CHANNEL) {
            if (line.find("</CHN") != std::string::npos) {
                if (ch_idx >= 0) {
                    current_evt.ch[ch_idx].nsamples = sample_count;
                }
                state = IN_BOARD;
                continue;
            }

            float tv, vv;
            if (ch_idx >= 0 && sample_count < MAX_SAMPLES &&
                sscanf(buf, " <Data>%f,%f</Data>", &tv, &vv) == 2) {
                current_evt.ch[ch_idx].time[sample_count]    = tv;
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
//  SEZIONE 4b: CALIBRAZIONE k_PMT (ampiezza ↔ slew rate)
// ==========================================================================

/// CalibrateSlewK():
///   Determina la costante di proporzionalità k tra ampiezza e slew rate
///   per ciascun PMT (PMT1 e PMT2), usando i dataset di calibrazione con
///   |x_k| ≤ CALIB_K_X_ABS_MAX.
///
/// METODO:
///   1. Per ogni dataset selezionato, parsa l'XML e per ogni evento estrae
///      i ChannelData già analizzati da AnalyzeChannel (incluso slew rate).
///   2. Per PMT1 e PMT2 separatamente, raccoglie le coppie (s, A) di
///      eventi che soddisfano TUTTI i criteri:
///        - has_pulse == true   (impulso fisico, sopra soglia di rumore)
///        - is_clipped == false (NO clipping: A misurata è quella vera)
///        - slew_ok == true     (fit slew rate riuscito con buona qualità)
///        - is_oscillating == false (esclude eventi patologici)
///   3. Costruisce un TGraph (s, A) per ciascun PMT e fitta:
///        - una volta con intercetta libera (controllo di sistematica),
///        - una volta con A = k · s (forzato per origine, k = parametro).
///   4. Salva k1, k2 nelle variabili globali gK_PMT[] e gK_PMT_err[].
///   5. Salva i grafici nel file ROOT di output per ispezione visiva.
///
/// PARAMETRI:
///   folder    — cartella contenente i file XML
///   positions — lista delle posizioni note (DEFAULT_POSITIONS o altra)
///   fout      — file ROOT aperto in scrittura per salvare i canvas
///
/// EFFETTO COLLATERALE:
///   gK_PMT[0..1], gK_PMT_err[0..1], gK_calibrated vengono aggiornate.
///
/// NOTA: questa funzione PARSA gli stessi XML che verranno parsati di nuovo
///   dal loop principale di TOF_Calibration. Il doppio parsing è inefficiente
///   ma mantiene il codice modulare. Costo tipico: ~30 secondi extra totali.

void CalibrateSlewK(const char* folder,
                    const std::vector<PositionInfo> &positions,
                    TFile* fout) {

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  CALIBRAZIONE SLEW RATE → AMPIEZZA           " << std::endl;
    std::cout << "  Range posizioni: |x| ≤ "
              << CALIB_K_X_ABS_MAX << " cm                     " << std::endl;
    std::cout << "=============================================" << std::endl;

    // ---- TGraph (s, A) per ciascun PMT ----
    // Nota: TGraph riempito uno punto alla volta con SetPoint(i, x, y).
    // Useremo l'indice progressivo n_pts_PMTk per riempirlo.
    TGraph *g_PMT1 = new TGraph();   g_PMT1->SetName("g_calib_slew_PMT1");
    TGraph *g_PMT2 = new TGraph();   g_PMT2->SetName("g_calib_slew_PMT2");
    g_PMT1->SetTitle("Calibrazione PMT1: A vs slew rate;"
                     "Slew rate s [mV/ns];Ampiezza A [mV]");
    g_PMT2->SetTitle("Calibrazione PMT2: A vs slew rate;"
                     "Slew rate s [mV/ns];Ampiezza A [mV]");

    int n_pts_PMT1 = 0;
    int n_pts_PMT2 = 0;

    // Contatori diagnostici
    int n_files_used      = 0;
    int n_events_total    = 0;
    int n_events_kept_p1  = 0;
    int n_events_kept_p2  = 0;

    // ---- Loop sulle posizioni di calibrazione ----
    for (size_t ip = 0; ip < positions.size(); ip++) {

        const PositionInfo &pos = positions[ip];

        // Selezione: solo dataset con |x| ≤ CALIB_K_X_ABS_MAX
        if (fabs(pos.x_cm) > CALIB_K_X_ABS_MAX) continue;

        // Risoluzione del path al file XML (stessa logica di TOF_Calibration)
        std::string xml_path = std::string(folder) + "/" + pos.name + ".xml";
        std::ifstream test_file(xml_path.c_str());
        if (!test_file.good()) {
            xml_path = std::string(folder) + pos.name + ".xml";
            test_file.open(xml_path.c_str());
            if (!test_file.good()) {
                std::cerr << "[WARNING] CalibrateSlewK: file non trovato per "
                          << pos.name << " — salto." << std::endl;
                continue;
            }
        }
        test_file.close();

        // Parsing dell'XML (chiama già AnalyzeChannel per ogni canale)
        std::vector<EventData> events;
        int n_parsed = ParseXML(xml_path.c_str(), events);
        if (n_parsed <= 0) continue;

        n_files_used++;
        n_events_total += (int)events.size();

        // ---- Loop sugli eventi: estrai (s, A) per PMT1 e PMT2 ----
        for (size_t ev = 0; ev < events.size(); ev++) {
            const EventData &e = events[ev];
            if (e.nchannels < 2) continue;

            // Iteriamo sui due PMT della barra (indici 0 e 1)
            for (int k = 0; k < 2; k++) {
                const ChannelData &cd = e.ch[k];

                // Criteri di selezione per la calibrazione di k:
                // - impulso reale (sopra soglia di rumore)
                // - NON clippato (l'ampiezza misurata è quella vera)
                // - slew rate fit OK (qualità del fit lineare)
                // - non oscillante (esclude patologie)
                if (!cd.has_pulse)        continue;
                if ( cd.is_clipped)       continue;
                if (!cd.slew_ok)          continue;
                if ( cd.is_oscillating)   continue;

                // Aggiungi al TGraph del PMT corrispondente
                if (k == 0) {
                    g_PMT1->SetPoint(n_pts_PMT1, cd.slew_rate, cd.amplitude);
                    n_pts_PMT1++;
                    n_events_kept_p1++;
                } else {
                    g_PMT2->SetPoint(n_pts_PMT2, cd.slew_rate, cd.amplitude);
                    n_pts_PMT2++;
                    n_events_kept_p2++;
                }
            }
        }
    }

    std::cout << "[INFO] CalibrateSlewK: "
              << n_files_used << " dataset, "
              << n_events_total << " eventi totali" << std::endl;
    std::cout << "       PMT1: " << n_events_kept_p1
              << " punti utili per il fit" << std::endl;
    std::cout << "       PMT2: " << n_events_kept_p2
              << " punti utili per il fit" << std::endl;

    // ---- Verifica statistica minima ----
    if (n_pts_PMT1 < 100 || n_pts_PMT2 < 100) {
        std::cerr << "[ERRORE] CalibrateSlewK: troppi pochi punti per il fit. "
                  << "Calibrazione FALLITA, recupero clippati DISABILITATO."
                  << std::endl;
        gK_calibrated = false;
        return;
    }

    // ---- Stima del range degli slew rate per inizializzare il fit ----
    // Useremo questi limiti anche per il range del TF1.
    double s_min_p1, s_max_p1, A_min_p1, A_max_p1;
    double s_min_p2, s_max_p2, A_min_p2, A_max_p2;
    g_PMT1->ComputeRange(s_min_p1, A_min_p1, s_max_p1, A_max_p1);
    g_PMT2->ComputeRange(s_min_p2, A_min_p2, s_max_p2, A_max_p2);

    // ---- Fit con intercetta libera (sistematica): A = k·s + q ----
    // Se q viene compatibile con 0 entro 2σ, possiamo usare il fit vincolato.
    // Altrimenti, la presenza di q ≠ 0 è un'indicazione di un offset di
    // baseline o di un piede non lineare residuo nella finestra di fit.
    auto FitOne = [&](TGraph* gr, const char* name,
                      double s_min, double s_max,
                      double &k_out, double &k_err_out,
                      double &q_out, double &q_err_out,
                      double &chi2ndf_out) {

        // Fit lineare libero: A = q + k·s (range [s_min, s_max])
        TF1 *f_lin = new TF1(Form("%s_lin", name),
                             "[0] + [1]*x", s_min, s_max);
        f_lin->SetParName(0, "q");
        f_lin->SetParName(1, "k");
        f_lin->SetParameters(0.0, A_max_p1 / std::max(s_max, 1.0));
        gr->Fit(f_lin, "Q R N");      // N = non disegnare (lo facciamo dopo)

        q_out     = f_lin->GetParameter(0);
        q_err_out = f_lin->GetParError(0);
        k_out     = f_lin->GetParameter(1);
        k_err_out = f_lin->GetParError(1);
        chi2ndf_out = (f_lin->GetNDF() > 0) ?
                      f_lin->GetChisquare() / (double)f_lin->GetNDF() : -1.0;

        delete f_lin;
    };

    auto FitOneOrigin = [&](TGraph* gr, const char* name,
                            double s_min, double s_max,
                            double &k_out, double &k_err_out,
                            double &chi2ndf_out) {

        // Fit lineare attraverso l'origine: A = k·s (range [s_min, s_max])
        TF1 *f_orig = new TF1(Form("%s_orig", name),
                              "[0]*x", s_min, s_max);
        f_orig->SetParName(0, "k");
        f_orig->SetParameter(0, A_max_p1 / std::max(s_max, 1.0));
        gr->Fit(f_orig, "Q R N +");   // + = aggiungi (mantiene anche il libero)

        k_out     = f_orig->GetParameter(0);
        k_err_out = f_orig->GetParError(0);
        chi2ndf_out = (f_orig->GetNDF() > 0) ?
                      f_orig->GetChisquare() / (double)f_orig->GetNDF() : -1.0;

        // Salva il TF1 nella TGraph (così viene scritto nel file)
        gr->GetListOfFunctions()->Add(f_orig);
    };

    // ---- Esecuzione fit ----
    double k1, k1_err, q1, q1_err, chi2_lin1;
    double k2, k2_err, q2, q2_err, chi2_lin2;
    FitOne(g_PMT1, "fit_PMT1", s_min_p1, s_max_p1,
           k1, k1_err, q1, q1_err, chi2_lin1);
    FitOne(g_PMT2, "fit_PMT2", s_min_p2, s_max_p2,
           k2, k2_err, q2, q2_err, chi2_lin2);

    double k1_orig, k1_orig_err, chi2_orig1;
    double k2_orig, k2_orig_err, chi2_orig2;
    FitOneOrigin(g_PMT1, "fit_PMT1", s_min_p1, s_max_p1,
                 k1_orig, k1_orig_err, chi2_orig1);
    FitOneOrigin(g_PMT2, "fit_PMT2", s_min_p2, s_max_p2,
                 k2_orig, k2_orig_err, chi2_orig2);

    // ---- Salva nelle globali (usiamo il fit vincolato per origine) ----
    gK_PMT[0]      = k1_orig;
    gK_PMT_err[0]  = k1_orig_err;
    gK_PMT[1]      = k2_orig;
    gK_PMT_err[1]  = k2_orig_err;
    gK_calibrated  = true;

    // ---- Stampa diagnostica ----
    std::cout << "\n  RISULTATI CALIBRAZIONE k_PMT:" << std::endl;
    std::cout << "  ----------------------------------------" << std::endl;
    std::cout << "  PMT1 (fit libero):    A = ("
              << Form("%.3f ± %.3f", k1, k1_err) << ") · s + ("
              << Form("%.2f ± %.2f", q1, q1_err) << ")  χ²/ndf="
              << Form("%.2f", chi2_lin1) << std::endl;
    std::cout << "  PMT1 (fit per origine): k = "
              << Form("%.3f ± %.3f", k1_orig, k1_orig_err)
              << "  χ²/ndf=" << Form("%.2f", chi2_orig1) << std::endl;
    std::cout << "  PMT2 (fit libero):    A = ("
              << Form("%.3f ± %.3f", k2, k2_err) << ") · s + ("
              << Form("%.2f ± %.2f", q2, q2_err) << ")  χ²/ndf="
              << Form("%.2f", chi2_lin2) << std::endl;
    std::cout << "  PMT2 (fit per origine): k = "
              << Form("%.3f ± %.3f", k2_orig, k2_orig_err)
              << "  χ²/ndf=" << Form("%.2f", chi2_orig2) << std::endl;

    // Test di sistematica sull'intercetta
    if (fabs(q1) > 2.0 * q1_err) {
        std::cout << "  [WARNING] PMT1: intercetta q = " << Form("%.2f", q1)
                  << " mV non compatibile con 0 entro 2σ → possibile offset"
                  << " di baseline o piede non lineare residuo." << std::endl;
    }
    if (fabs(q2) > 2.0 * q2_err) {
        std::cout << "  [WARNING] PMT2: intercetta q = " << Form("%.2f", q2)
                  << " mV non compatibile con 0 entro 2σ → possibile offset"
                  << " di baseline o piede non lineare residuo." << std::endl;
    }

    // ---- Costruzione canvas di output ----
    TCanvas *c_calib_k = new TCanvas("c_calib_k",
                                      "Calibrazione k_PMT: A vs slew rate",
                                      1200, 500);
    c_calib_k->Divide(2, 1);

    c_calib_k->cd(1);
    gPad->SetGrid(1, 1);
    g_PMT1->SetMarkerStyle(20);
    g_PMT1->SetMarkerSize(0.4);
    g_PMT1->SetMarkerColor(kBlue + 1);
    g_PMT1->Draw("AP");

    TPaveText *pt1 = new TPaveText(0.15, 0.65, 0.55, 0.88, "NDC");
    pt1->SetFillColorAlpha(kWhite, 0.85);
    pt1->SetTextFont(42);
    pt1->SetTextSize(0.04);
    pt1->SetTextAlign(12);
    pt1->AddText("PMT1: A = k_{1} #upoint s");
    pt1->AddText(Form("k_{1} = %.3f #pm %.3f ns", k1_orig, k1_orig_err));
    pt1->AddText(Form("#chi^{2}/ndf = %.2f", chi2_orig1));
    pt1->AddText(Form("N punti = %d", n_pts_PMT1));
    pt1->Draw();

    c_calib_k->cd(2);
    gPad->SetGrid(1, 1);
    g_PMT2->SetMarkerStyle(20);
    g_PMT2->SetMarkerSize(0.4);
    g_PMT2->SetMarkerColor(kRed + 1);
    g_PMT2->Draw("AP");

    TPaveText *pt2 = new TPaveText(0.15, 0.65, 0.55, 0.88, "NDC");
    pt2->SetFillColorAlpha(kWhite, 0.85);
    pt2->SetTextFont(42);
    pt2->SetTextSize(0.04);
    pt2->SetTextAlign(12);
    pt2->AddText("PMT2: A = k_{2} #upoint s");
    pt2->AddText(Form("k_{2} = %.3f #pm %.3f ns", k2_orig, k2_orig_err));
    pt2->AddText(Form("#chi^{2}/ndf = %.2f", chi2_orig2));
    pt2->AddText(Form("N punti = %d", n_pts_PMT2));
    pt2->Draw();

    c_calib_k->Update();
    if (fout) {
        fout->cd();
        c_calib_k->Write("Calibration_k_slew_to_amplitude");
        g_PMT1->Write();
        g_PMT2->Write();
    }

    std::cout << "  [DONE] Calibrazione k completata." << std::endl;
    std::cout << "=============================================\n" << std::endl;
}
// ==========================================================================
//  SEZIONE 5: FIT DI UNA DISTRIBUZIONE TEMPORALE
// ==========================================================================

/// EstimateHistFWHM(): stima la FWHM direttamente dall'istogramma.
///
/// Metodo robusto contro outlier: parte dal bin con il massimo e scende
/// a destra e sinistra fino a trovare i punti a metà altezza, interpolando
/// linearmente tra bin adiacenti. Insensibile alle code sparse.

double EstimateHistFWHM(TH1D* h) {

    int max_bin = h->GetMaximumBin();
    double half_max = h->GetBinContent(max_bin) / 2.0;
    int nbins = h->GetNbinsX();

    // Bordo sinistro: scendi dal picco verso sinistra
    double x_left = h->GetBinCenter(1);
    for (int i = max_bin; i > 1; i--) {
        if (h->GetBinContent(i) < half_max) {
            double y1 = h->GetBinContent(i);
            double y2 = h->GetBinContent(i + 1);
            double x1 = h->GetBinCenter(i);
            double x2 = h->GetBinCenter(i + 1);
            if (fabs(y2 - y1) > 1e-6)
                x_left = x1 + (half_max - y1) / (y2 - y1) * (x2 - x1);
            else
                x_left = (x1 + x2) / 2.0;
            break;
        }
    }

    // Bordo destro: scendi dal picco verso destra
    double x_right = h->GetBinCenter(nbins);
    for (int i = max_bin; i < nbins; i++) {
        if (h->GetBinContent(i) < half_max) {
            double y1 = h->GetBinContent(i - 1);
            double y2 = h->GetBinContent(i);
            double x1 = h->GetBinCenter(i - 1);
            double x2 = h->GetBinCenter(i);
            if (fabs(y2 - y1) > 1e-6)
                x_right = x1 + (half_max - y1) / (y2 - y1) * (x2 - x1);
            else
                x_right = (x1 + x2) / 2.0;
            break;
        }
    }

    double fwhm = x_right - x_left;
    if (fwhm < h->GetBinWidth(1)) fwhm = 2.0 * h->GetStdDev();
    return fwhm;
}


/// FitLorentzian(): fitta un istogramma con una Lorentziana (Cauchy)
/// usando una parametrizzazione CUSTOM con parametri fisici intuitivi
/// e log-likelihood Poissoniana sull'INTERO range dell'istogramma.
///
/// FORMA FUNZIONALE:
///
///                       (Γ/2)²
///   f(x) = A · ─────────────────────
///              (x − x₀)² + (Γ/2)²
///
///   [0] = A      → altezza del picco (conteggi al centro, f(x₀) = A)
///   [1] = x₀     → posizione del centro [ns]
///   [2] = Γ      → FWHM (Full Width at Half Maximum) [ns]
///
/// PERCHÈ QUESTA PARAMETRIZZAZIONE e non "breitwigner" di ROOT:
///   La "breitwigner" di ROOT usa p0 come normalizzazione dell'AREA,
///   non come altezza del picco. Al centro: BW(x₀) = p0 · 2/(πΓ),
///   quindi p0 ha significato non intuitivo e l'inizializzazione è
///   difficile. Con la formula custom, p0 = A = altezza del picco →
///   l'inizializzazione è banale: A ≈ conteggi nel bin più alto.
///
/// PERCHÈ RANGE COMPLETO e non ristretto:
///   La Lorentziana ha code ∼ 1/x² che decadono lentamente. Gli eventi
///   sparsi nelle code (anche quelli "outlier" da CFD errato) sono
///   COERENTI con una distribuzione Lorentziana. Il fit di log-likelihood
///   Poissoniana gestisce correttamente i bin con 0, 1, 2 conteggi
///   (a differenza del χ² che necessita bin ben popolati). Restringere
///   il range causerebbe problemi di normalizzazione e perdita di
///   informazione sulle code.
///
/// INIZIALIZZAZIONE:
///   A   ← altezza del bin più alto (esatto per definizione)
///   x₀  ← posizione del bin più alto
///   Γ   ← FWHM misurata direttamente dall'istogramma (robusta)

FitResult FitLorentzian(TH1D* h, const char* fit_name = "bw_fit") {

    FitResult res;
    res.fit_ok = false;
    res.nentries = (int)h->GetEntries();

    // Se l'istogramma ha troppi pochi eventi, non fittare
    if (res.nentries < 30) {
        std::cerr << "[WARNING] Troppi pochi eventi (" << res.nentries
                  << ") per il fit di " << h->GetName() << std::endl;
        res.center    = h->GetMean();
        res.center_err = h->GetMeanError();
        res.width     = h->GetStdDev();
        res.width_err = h->GetStdDevError();
        res.chi2_ndf  = -1.0;
        return res;
    }

    // ---- Stime iniziali dai dati ----
    int    imax     = h->GetMaximumBin();
    double peak_hgt = h->GetBinContent(imax);           // altezza picco [conteggi]
    double peak_pos = h->GetBinCenter(imax);             // posizione picco [ns]
    double hist_fwhm = EstimateHistFWHM(h);              // FWHM robusta [ns]

    // Fallback: se la FWHM dall'istogramma non è sensata, usa il 5% del range
    double xlo = h->GetXaxis()->GetXmin();
    double xhi = h->GetXaxis()->GetXmax();
    if (hist_fwhm <= 0 || hist_fwhm > (xhi - xlo)) {
        hist_fwhm = 0.05 * (xhi - xlo);
    }

    // ---- Costruzione della TF1 con formula Lorentziana custom ----
    // Range di fit: INTERO istogramma (la log-likelihood gestisce le code)
    TF1 *fL = new TF1(fit_name,
        "[0]*([2]/2.)*([2]/2.)/((x-[1])*(x-[1]) + ([2]/2.)*([2]/2.))",
        xlo, xhi);

    fL->SetParName(0, "A_peak");
    fL->SetParName(1, "x0");
    fL->SetParName(2, "Gamma");

    // ---- Inizializzazione parametri ----
    // Diretta e intuitiva grazie alla parametrizzazione custom:
    //   A = altezza picco (f(x₀) = A per definizione)
    //   x₀ = posizione del massimo dell'istogramma
    //   Γ = FWHM misurata dall'istogramma
    fL->SetParameters(peak_hgt, peak_pos, hist_fwhm);

    // Limiti morbidi per evitare regioni non fisiche
    fL->SetParLimits(0, 0.1, 100.0 * peak_hgt);          // A positivo
    fL->SetParLimits(1, xlo, xhi);                        // x₀ nel range
    fL->SetParLimits(2, 1e-3 * (xhi - xlo), xhi - xlo);  // Γ positivo

    // ---- Esecuzione del fit ----
    // "L" = log-likelihood Poissoniana binnata (gestisce bin con 0-1-2 conteggi)
    // "S" = salva TFitResultPtr
    // "R" = usa il range della TF1
    // "Q" = quiet
    int fit_status = h->Fit(fL, "L S R Q");

    // ---- Estrazione risultati ----
    res.center     = fL->GetParameter(1);         // x₀ = centro distribuzione
    res.center_err = fL->GetParError(1);
    res.width      = fL->GetParameter(2);          // Γ = FWHM
    res.width_err  = fL->GetParError(2);
    res.chi2_ndf   = (fL->GetNDF() > 0) ?
                     fL->GetChisquare() / fL->GetNDF() : -1.0;
    res.fit_ok     = (fit_status == 0);

    // Stampa diagnostica
    std::cout << "  [FIT] " << h->GetName()
              << ": x0 = " << res.center << " ± " << res.center_err
              << " ns, Γ = " << res.width << " ± " << res.width_err
              << " ns, χ²/ndf = " << res.chi2_ndf << std::endl;

    return res;
}

/// FitGaussianCore(): fitta un istogramma con una Gaussiana ristretta al
/// "core" della distribuzione, definito come l'intervallo che contiene
/// il 90% (default) degli eventi attorno alla mediana.
///
/// MOTIVAZIONE FISICA:
///   Dopo aver applicato i tagli di qualità a livello di forma d'onda
///   (clipping, oscillazione), la distribuzione delle differenze temporali
///   (Δt₁₂, Δt₁₃, Δt₂₃, C) ha un core dominato dal jitter Gaussiano del
///   timing — somma di molti contributi indipendenti (jitter elettronico,
///   fluttuazioni del numero di fotoelettroni, rumore del PMT). Per il
///   teorema del limite centrale, il core è effettivamente Gaussiano.
///
///   Le code residue contengono ancora alcuni effetti sistematici:
///     - time-walk residuo (rise time variabile con la posizione)
///     - eventi con CFD su un fronte distorto da pile-up minore
///     - fluttuazioni di calibrazione della DRS4
///   Questi effetti deformano la distribuzione asimmetricamente e gonfiano
///   la deviazione standard apparente. Restringendo il fit al 90% centrale,
///   si ottiene una σ "pulita" che riflette la VERA risoluzione temporale
///   dell'apparato per la quantità misurata.
///
/// FORMA FUNZIONALE:
///                           ⎛   (x − μ)²  ⎞
///   f(x) = A · exp ⎜ − ───────── ⎟
///                           ⎝     2σ²      ⎠
///
///   [0] = A   → altezza del picco (conteggi al centro, f(μ) = A)
///   [1] = μ   → posizione del centro [ns]   ← VALORE MEDIO
///   [2] = σ   → deviazione standard [ns]    ← RISOLUZIONE TEMPORALE
///
///   Per riferimento: FWHM = 2√(2 ln 2) σ ≈ 2.355 σ
///
/// RANGE DEL FIT:
///   I quantili 5% e 95% dell'istogramma definiscono il core 90%.
///   ROOT li calcola con TH1::GetQuantiles(). Il range del fit è quindi
///   [Q05, Q95] e tutti i bin fuori da questo intervallo sono ignorati.
///
/// INIZIALIZZAZIONE PARAMETRI:
///   A   ← altezza del bin più alto (esatto per definizione)
///   μ   ← posizione del bin più alto
///   σ   ← FWHM_istogramma / 2.355  (FWHM stimata robustamente da
///         EstimateHistFWHM() — già definita altrove nel codice)
///
/// PARAMETRI:
///   h              — istogramma da fittare
///   core_fraction  — frazione di eventi nel core (default 0.90 → quantili 5%–95%)
///   fit_name       — nome interno della TF1 (per evitare conflitti)

FitResult FitGaussianCore(TH1D* h,
                          double core_fraction = 0.90,
                          const char* fit_name = "gaus_fit") {

    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h->GetEntries();

    // ---- Controllo statistica minima ----
    // Sotto i 30 eventi il fit non è affidabile: ritorna le statistiche
    // grezze dell'istogramma come fallback.
    if (res.nentries < 30) {
        std::cerr << "[WARNING] Troppi pochi eventi (" << res.nentries
                  << ") per il fit di " << h->GetName() << std::endl;
        res.center     = h->GetMean();
        res.center_err = h->GetMeanError();
        res.width      = h->GetStdDev();
        res.width_err  = h->GetStdDevError();
        res.chi2_ndf   = -1.0;
        return res;
    }

    // ---- Calcolo dei quantili che definiscono il core ----
    // Per core_fraction = 0.90: quantili 5% e 95%.
    // GetQuantiles() restituisce i valori di x ai quantili richiesti.
    double q_lo = 0.5 - core_fraction / 2.0;  // 0.05
    double q_hi = 0.5 + core_fraction / 2.0;  // 0.95
    const int nq = 2;
    double xq[nq] = { q_lo, q_hi };
    double yq[nq] = { 0.0, 0.0 };
    h->GetQuantiles(nq, yq, xq);
    double x_lo = yq[0];   // Estremo inferiore del core
    double x_hi = yq[1];   // Estremo superiore del core

    // ---- Sanity check sul range ----
    // Se il range è inconsistente (es. istogramma con tutti gli eventi
    // in un solo bin), torna al range completo dell'istogramma.
    if (x_hi - x_lo < 5.0 * h->GetBinWidth(1)) {
        std::cerr << "[WARNING] Range core troppo piccolo per "
                  << h->GetName() << ", uso range completo." << std::endl;
        x_lo = h->GetXaxis()->GetXmin();
        x_hi = h->GetXaxis()->GetXmax();
    }

    // ---- Stime iniziali dai dati ----
    int    imax       = h->GetMaximumBin();
    double peak_hgt   = h->GetBinContent(imax);              // altezza picco [conteggi]
    double peak_pos   = h->GetBinCenter(imax);                // posizione picco [ns]
    double hist_fwhm  = EstimateHistFWHM(h);                  // FWHM robusta [ns]
    double sigma_init = hist_fwhm / 2.355;                    // σ ≈ FWHM / 2.355

    // Fallback se la stima è inconsistente
    if (sigma_init <= 0 || sigma_init > 0.5 * (x_hi - x_lo)) {
        sigma_init = 0.25 * (x_hi - x_lo);
    }

    // ---- Costruzione della TF1 Gaussiana ----
    // Nota: usiamo la formula esplicita (non "gaus" di ROOT) per avere
    // controllo completo sui nomi dei parametri e sul significato fisico
    // di [0] (altezza del picco, non normalizzazione).
    TF1 *fG = new TF1(fit_name,
        "[0]*TMath::Exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",
        x_lo, x_hi);

    fG->SetParName(0, "A_peak");
    fG->SetParName(1, "mu");
    fG->SetParName(2, "sigma");

    fG->SetParameters(peak_hgt, peak_pos, sigma_init);

    // Limiti morbidi per evitare regioni non fisiche
    fG->SetParLimits(0, 0.1, 100.0 * peak_hgt);            // A positivo
    fG->SetParLimits(1, x_lo, x_hi);                        // μ nel range
    fG->SetParLimits(2, h->GetBinWidth(1), x_hi - x_lo);    // σ tra 1 bin e range

    // ---- Esecuzione del fit ----
    // "L" = log-likelihood Poissoniana binnata
    //       (gestisce correttamente bin con pochi conteggi anche se nel core
    //        sono quasi sempre ben popolati, è comunque la scelta più rigorosa)
    // "S" = salva TFitResultPtr per uso futuro
    // "R" = usa il range della TF1 (cioè [x_lo, x_hi])
    // "Q" = quiet (sopprimi output di Minuit)
    int fit_status = h->Fit(fG, "L S R Q");

    // ---- Estrazione risultati ----
    // ATTENZIONE: width = σ (NON FWHM). La risoluzione temporale è σ direttamente.
    // Per ottenere FWHM, moltiplicare per 2.355.
    res.center     = fG->GetParameter(1);             // μ = centro distribuzione
    res.center_err = fG->GetParError(1);
    res.width      = fG->GetParameter(2);              // σ = risoluzione temporale
    res.width_err  = fG->GetParError(2);
    res.chi2_ndf   = (fG->GetNDF() > 0) ?
                     fG->GetChisquare() / fG->GetNDF() : -1.0;
    res.fit_ok     = (fit_status == 0);

    // Stampa diagnostica
    std::cout << "  [GAUSS-FIT] " << h->GetName()
              << ": μ = " << res.center << " ± " << res.center_err
              << " ns, σ = " << res.width << " ± " << res.width_err
              << " ns (range [" << x_lo << ", " << x_hi << "])"
              << ", χ²/ndf = " << res.chi2_ndf << std::endl;

    return res;
}
// ==========================================================================
//  SEZIONE 5b: COMPONENTI DELL'ERRORE SULLA POSIZIONE
// ==========================================================================

/// GetParallaxSigmaFromFile():
///   Restituisce σ_parallax [cm] alla posizione x [cm], leggendo il
///   TGraph "g_sigma_x" dal file ROOT prodotto da ParallaxScan().
///   Usa interpolazione lineare se x non è esattamente sulla griglia.
///
///   Comportamento di fallback:
///     - Se file vuoto o non apribile      → ritorna PARALLAX_FALLBACK
///     - Se TGraph "g_sigma_x" non presente → ritorna PARALLAX_FALLBACK
///     - Se x è fuori dal range della griglia → estrapolazione lineare
///       (TGraph::Eval estrapola linearmente; emette un warning)
///
///   IMPORTANTE: per efficienza, il file ROOT viene aperto UNA SOLA VOLTA
///   e mantenuto aperto in una variabile statica. Le chiamate successive
///   riutilizzano la stessa istanza. Il file viene chiuso automaticamente
///   alla fine del programma.

double GetParallaxSigmaFromFile(double x_k, const char* parallax_file) {

    // Variabili statiche: persistono tra chiamate successive
    static TFile  *s_file  = nullptr;     // istanza del file
    static TGraph *s_graph = nullptr;     // TGraph caricato
    static std::string s_loaded_path = "";  // path attualmente caricato
    static bool   s_warned = false;        // per stampare il warning una volta sola

    // Path richiesto come stringa (evita problemi con strcmp)
    std::string requested_path = (parallax_file ? parallax_file : "");

    // Caso 1: nessun file fornito → fallback con warning una sola volta
    if (requested_path.empty()) {
        if (!s_warned) {
            std::cerr << "[WARNING] PARALLAX_FILE non specificato, "
                      << "uso valore costante " << PARALLAX_FALLBACK
                      << " cm per σ_parallax." << std::endl;
            s_warned = true;
        }
        return PARALLAX_FALLBACK;
    }

    // Caso 2: il path è cambiato dall'ultima chiamata → ricarica
    if (requested_path != s_loaded_path) {
        // Pulisci l'istanza precedente (se esisteva)
        if (s_file) { s_file->Close(); delete s_file; s_file = nullptr; s_graph = nullptr; }

        // Apri il nuovo file
        s_file = TFile::Open(parallax_file, "READ");
        if (!s_file || s_file->IsZombie()) {
            std::cerr << "[WARNING] Impossibile aprire file parallasse '"
                      << requested_path << "', uso fallback "
                      << PARALLAX_FALLBACK << " cm." << std::endl;
            if (s_file) { delete s_file; s_file = nullptr; }
            s_loaded_path = requested_path;  // memorizza per non ritentare
            return PARALLAX_FALLBACK;
        }

        // Estrai il TGraph σ_parallax(x)
        s_graph = (TGraph*)s_file->Get("g_sigma_x");
        if (!s_graph) {
            std::cerr << "[WARNING] TGraph 'g_sigma_x' non trovato in '"
                      << requested_path << "', uso fallback "
                      << PARALLAX_FALLBACK << " cm." << std::endl;
            s_file->Close();
            delete s_file;
            s_file = nullptr;
            s_loaded_path = requested_path;
            return PARALLAX_FALLBACK;
        }

        s_loaded_path = requested_path;
        std::cout << "[INFO] Caricato file di parallasse: "
                  << requested_path << " (" << s_graph->GetN() << " punti)"
                  << std::endl;
    }

    // Caso 3: file e grafico OK → interpola linearmente
    // TGraph::Eval() restituisce il valore del grafico interpolato
    // linearmente tra i due punti più vicini. Estrapola fuori range.
    if (s_graph) {
        return s_graph->Eval(x_k);
    }

    // Fallback finale (non dovrebbe mai accadere)
    return PARALLAX_FALLBACK;
}


/// GetTapeSigma():
///   Restituisce σ_tape [cm] alla posizione x_k, secondo il modello
///   selezionato da TAPE_MODEL_INCREMENTAL.

double GetTapeSigma(double x_k) {

    double dx_from_ref = fabs(x_k - TAPE_REFERENCE_X);

    // Posizione di riferimento: errore zero per definizione
    if (dx_from_ref < 1e-6) return 0.0;

    if (!TAPE_MODEL_INCREMENTAL) {
        // Modalità A "lettura singola": ogni posizione ≠ 0 ha un solo
        // evento di lettura → σ costante.
        return SIGMA_TAPE_READ;
    } else {
        // Modalità B "letture incrementali": σ_read² accumula in
        // quadratura per ogni passo di TAPE_STEP_LENGTH cm percorso
        // dal riferimento.
        double n_step_d = dx_from_ref / TAPE_STEP_LENGTH;
        // Si arrotonda per eccesso (anche un passo parziale conta come uno
        // step di posizionamento). Es. da 0 a 30 cm: 30/28 = 1.07 → 2 passi.
        // Modificare in round/floor a seconda della convenzione adottata.
        int n_step = (int)ceil(n_step_d - 1e-6);
        if (n_step < 1) n_step = 1;
        return SIGMA_TAPE_READ * sqrt((double)n_step);
    }
}


/// ComputeSigmaX():
///   Funzione PRINCIPALE: restituisce l'errore totale σ_x [cm] per
///   la posizione x_k, combinando in quadratura parallasse e metro.
///
///     σ_x²(x_k) = σ_parallax²(x_k) + σ_tape²(x_k)
///
///   parallax_file: path al file ROOT del MC. Stringa vuota = disabilitato.

double ComputeSigmaX(double x_k, const char* parallax_file) {
    double sigma_par  = GetParallaxSigmaFromFile(x_k, parallax_file);
    double sigma_tape = GetTapeSigma(x_k);
    return sqrt(sigma_par * sigma_par + sigma_tape * sigma_tape);
}
// ==========================================================================
//  SEZIONE 6: ELABORAZIONE DI UN SINGOLO DATASET (una posizione x_k)
// ==========================================================================

/// ProcessDataset(): per una data posizione x_k, legge il file XML,
/// calcola le differenze temporali evento per evento, riempie un TTree,
/// crea istogrammi di Δt₁₂ e C, e fitta le distribuzioni.
///
/// Restituisce una coppia {FitResult per Δt₁₂, FitResult per C}.
///
/// MAPPATURA CANALI:
///   CH1 (idx 0) → PMT1 (estremo sinistro, x = -L/2)
///   CH2 (idx 1) → PMT2 (estremo destro, x = +L/2)
///   CH3 (idx 2) → PMT3 (mobile, posizione x_k)
///   CH4 (idx 3) → trigger NIM (non usato per il timing)
///
/// N.B.: la mappatura si adatta automaticamente al numero di canali presenti
/// nel file, ma richiede che PMT1, PMT2, PMT3 siano nei primi 3 canali.

struct DatasetResults {
    FitResult dt12;     // Risultato fit Δt₁₂ = t₁ − t₂
    FitResult dt13;     // Risultato fit Δt₁₃ = t₁ − t₃
    FitResult dt23;     // Risultato fit Δt₂₃ = t₂ − t₃
    FitResult C_fit;    // Risultato fit C = t₃ − (t₁+t₂)/2
    int       n_good;   // Numero di eventi con CFD valido su tutti e 3 i canali
};

DatasetResults ProcessDataset(const char* xml_path,
                              const std::string &tree_name,
                              TFile* outfile) {

    DatasetResults result;
    result.n_good = 0;

    // ---- Parsing del file XML ----
    std::vector<EventData> events;
    int n_parsed = ParseXML(xml_path, events);
    if (n_parsed <= 0) {
        std::cerr << "[ERRORE] Nessun evento letto da " << xml_path << std::endl;
        return result;
    }

    // ---- Preparazione del TTree ----
    // Salviamo per ogni evento le grandezze di interesse.
    // Questo permette analisi successive senza dover ri-parsare l'XML.
    outfile->cd();
    TTree* tree = new TTree(tree_name.c_str(),
                            Form("Dati TOF posizione %s", tree_name.c_str()));

    // Branch per i tempi CFD dei 3 PMT [ns]
    Float_t t_cfd[3];
    tree->Branch("t_cfd", t_cfd, "t_cfd[3]/F");

    // Branch per le ampiezze dei 3 PMT [mV]
    Float_t amp[3];
    tree->Branch("amp", amp, "amp[3]/F");

    // Branch per le baseline dei 3 PMT [mV]
    Float_t bl[3];
    tree->Branch("bl", bl, "bl[3]/F");

    // Branch per le differenze temporali [ns]
    Float_t dt12, dt13, dt23;
    tree->Branch("dt12", &dt12, "dt12/F");   // t₁ − t₂
    tree->Branch("dt13", &dt13, "dt13/F");   // t₁ − t₃
    tree->Branch("dt23", &dt23, "dt23/F");   // t₂ − t₃

    // Branch per la costante C = t₃ − (t₁+t₂)/2  [ns]
    Float_t C_val;
    tree->Branch("C", &C_val, "C/F");

    // Flag di qualità (1 = buono, 0 = scartato)
    Int_t good_event;
    tree->Branch("good", &good_event, "good/I");

    // ---- Istogrammi ----
    // Creiamo gli istogrammi qui, li riempiamo nel loop sugli eventi,
    // e poi li fittiamo alla fine. I nomi includono la posizione per
    // evitare conflitti nel file ROOT.
    TH1D *h_dt12 = new TH1D(Form("h_dt12_%s", tree_name.c_str()),
                             Form("#Delta t_{12} @ %s;#Delta t_{12} = t_{1} - t_{2} [ns];Conteggi",
                                  tree_name.c_str()),
                             NBINS_DT, DT_LO, DT_HI);

    TH1D *h_dt13 = new TH1D(Form("h_dt13_%s", tree_name.c_str()),
                             Form("#Delta t_{13} @ %s;#Delta t_{13} = t_{1} - t_{3} [ns];Conteggi",
                                  tree_name.c_str()),
                             NBINS_DT, DT_LO, DT_HI);

    TH1D *h_dt23 = new TH1D(Form("h_dt23_%s", tree_name.c_str()),
                             Form("#Delta t_{23} @ %s;#Delta t_{23} = t_{2} - t_{3} [ns];Conteggi",
                                  tree_name.c_str()),
                             NBINS_DT, DT_LO, DT_HI);

    TH1D *h_C = new TH1D(Form("h_C_%s", tree_name.c_str()),
                          Form("C = t_{3} - (t_{1}+t_{2})/2 @ %s;"
                               "C = t_{3} - #frac{t_{1}+t_{2}}{2} [ns];Conteggi",
                               tree_name.c_str()),
                          NBINS_C, C_LO, C_HI);

    // Branch flag clipping/oscillazione (per diagnostica nel TTree)
    Int_t clip_flag, osc_flag, recovered_flag;
    tree->Branch("clip",      &clip_flag,      "clip/I");        // 1 se almeno un canale clippato
    tree->Branch("osc",       &osc_flag,       "osc/I");         // 1 se almeno un canale oscillante
    tree->Branch("recovered", &recovered_flag, "recovered/I");   // 1 se evento recuperato via slew

    // Branch ampiezze ricostruite (per analisi a posteriori)
    Float_t amp_rec[3];
    tree->Branch("amp_rec", amp_rec, "amp_rec[3]/F");   // 0 se non recuperato

    // Contatori delle ragioni di scarto (per diagnostica statistica)
    int n_no_cfd      = 0;   // CFD fallito su qualche canale
    int n_clipped     = 0;   // clipping su qualche canale
    int n_oscillating = 0;   // oscillazione su qualche canale
    // ---- Loop sugli eventi ----
    for (size_t ev = 0; ev < events.size(); ev++) {

        EventData &e = events[ev];

        // Verifica che ci siano almeno 3 canali
        if (e.nchannels < 3) {
            good_event = 0;
            for (int k = 0; k < 3; k++) { t_cfd[k] = -999; amp[k] = 0; bl[k] = 0; }
            dt12 = -999; dt13 = -999; dt23 = -999; C_val = -999;
            tree->Fill();
            continue;
        }

        // ---- Tentativo di RECUPERO degli eventi clippati ----
        // Per PMT1 e PMT2 (indici 0,1) che sono clippati, prova a ricostruire
        // il tempo CFD usando il metodo dello slew rate (k_PMT calibrato).
        // PMT3 (indice 2) NON viene recuperato (è il PMT mobile, tipicamente
        // non clippato; il suo k non è stato calibrato).
        // Modifica il riferimento all'EventData per scrivere nei suoi campi.
        EventData &e_mod = events[ev];
        for (int k = 0; k < 2; k++) {
            if (e_mod.ch[k].is_clipped) {
                RecoverClippedCFD(e_mod.ch[k], k);
            }
        }

        // ---- Estrazione dei tempi e dei flag di qualità ----
        // Per ogni canale, decidiamo quale tempo usare:
        //   - Se NON clippato: t_cfd standard (cd.t_cfd) — richiede cfd_ok
        //   - Se CLIPPATO + recuperato: t_cfd_recovered (PMT1, PMT2 soltanto)
        //   - Altrimenti: -999 → canale non utilizzabile per questo evento
        //
        // t_used[k] è il tempo "migliore" per il canale k.
        Float_t t_used[3];
        bool    used_ok[3];
        bool    ch_recovered[3] = {false, false, false};

        bool any_clipped     = false;
        bool any_oscillating = false;
        bool any_recovered   = false;

        for (int k = 0; k < 3; k++) {
            const ChannelData &cd = e_mod.ch[k];

            t_cfd[k]   = (Float_t)cd.t_cfd;
            amp[k]     = (Float_t)cd.amplitude;
            bl[k]      = (Float_t)cd.baseline;
            amp_rec[k] = (Float_t)cd.amplitude_rec;   // 0 se non recuperato

            if (cd.is_clipped)     any_clipped     = true;
            if (cd.is_oscillating) any_oscillating = true;

            // Selezione del tempo "migliore" per questo canale
            if (!cd.is_clipped) {
                // Canale non clippato: usa CFD standard
                t_used[k]  = (Float_t)cd.t_cfd;
                used_ok[k] = cd.cfd_ok;
            } else if (cd.cfd_recovered) {
                // Canale clippato ma recuperato via slew rate
                t_used[k]       = (Float_t)cd.t_cfd_recovered;
                used_ok[k]      = true;
                ch_recovered[k] = true;
                any_recovered   = true;
            } else {
                // Canale clippato e non recuperabile (slew_ok=false, soglia
                // CFD nel clip, o k non calibrato): tempo non usabile
                t_used[k]  = -999.0f;
                used_ok[k] = false;
            }
        }

        bool all_t_ok = used_ok[0] && used_ok[1] && used_ok[2];

        // Flag aggregati per il TTree
        clip_flag      = any_clipped     ? 1 : 0;
        osc_flag       = any_oscillating ? 1 : 0;
        recovered_flag = any_recovered   ? 1 : 0;

        // ---- Decisione di qualità (NUOVA logica) ----
        // L'evento è "buono" se:
        //   1. Tutti e tre i canali hanno un tempo utilizzabile (CFD standard
        //      OPPURE CFD ricostruito).
        //   2. Nessun canale ha rilevato oscillazione patologica.
        //
        // Notare che un canale "clippato" non implica più automaticamente
        // lo scarto: se il recupero è andato a buon fine, l'evento sopravvive.
        bool good_quality = true;
        if (ENABLE_CFD_CUT && !all_t_ok)      good_quality = false;
        if (ENABLE_OSC_CUT && any_oscillating) good_quality = false;
        // ENABLE_CLIP_CUT non viene più applicato come scarto duro:
        // il clipping è gestito tramite il recupero. Lo manteniamo nei
        // contatori per diagnostica, e per uso "puristico" (se l'utente
        // setta ENABLE_CLIP_CUT = true E vuole il vecchio comportamento,
        // basta scartare anche gli eventi recuperati con: any_clipped).
        // Per ora consideriamo: se ENABLE_CLIP_CUT è true E l'evento ha
        // canali clippati NON recuperati su PMT1/PMT2, allora già fallisce
        // per all_t_ok. Comportamento equivalente in pratica.

        // Aggiornamento contatori delle ragioni di scarto (informativi)
        if (!all_t_ok)        n_no_cfd++;
        if (any_clipped)      n_clipped++;
        if (any_oscillating)  n_oscillating++;

        if (good_quality) {
            // Calcolo delle differenze temporali e di C usando t_used[]
            // (= t_cfd o t_cfd_recovered, a seconda del canale).
            dt12  = t_used[0] - t_used[1];                    // Δt₁₂ = t₁ − t₂
            dt13  = t_used[0] - t_used[2];                    // Δt₁₃ = t₁ − t₃
            dt23  = t_used[1] - t_used[2];                    // Δt₂₃ = t₂ − t₃
            C_val = t_used[2] - (t_used[0] + t_used[1]) / 2.0; // C = t₃ − (t₁+t₂)/2

            good_event = 1;
            result.n_good++;

            // Riempi gli istogrammi solo con eventi buoni
            h_dt12->Fill(dt12);
            h_dt13->Fill(dt13);
            h_dt23->Fill(dt23);
            h_C->Fill(C_val);

        } else {
            // Evento scartato: valori sentinella nel TTree
            dt12 = -999; dt13 = -999; dt23 = -999; C_val = -999;
            good_event = 0;
        }

        tree->Fill();
    }

    // ---- Stampa riepilogativa con breakdown delle ragioni di scarto ----
    // Aggiunto un contatore di eventi recuperati via slew rate.
    int n_recovered = 0;
    for (size_t ev = 0; ev < events.size(); ev++) {
        const EventData &e = events[ev];
        if (e.nchannels < 3) continue;
        bool rec = false;
        for (int k = 0; k < 2; k++) {
            if (e.ch[k].is_clipped && e.ch[k].cfd_recovered) rec = true;
        }
        if (rec) n_recovered++;
    }

    std::cout << "[INFO] " << tree_name << ": "
              << result.n_good << "/" << events.size()
              << " eventi buoni"
              << " (di cui " << n_recovered << " recuperati via slew rate)"
              << std::endl;
    std::cout << "       Categorie: "
              << n_no_cfd      << " no-CFD/non-utilizzabili, "
              << n_clipped     << " con clipping, "
              << n_oscillating << " con oscillazione "
              << "(possono sovrapporsi)" << std::endl;
// ---- Fit delle distribuzioni con Gaussiana sul 90% centrale ----
    // FitGaussianCore() esclude le code 5% per parte e fitta una Gaussiana
    // pura. La σ restituita è la risoluzione temporale del core.
    result.dt12  = FitGaussianCore(h_dt12, 0.90, Form("fit_dt12_%s", tree_name.c_str()));
    result.dt13  = FitGaussianCore(h_dt13, 0.90, Form("fit_dt13_%s", tree_name.c_str()));
    result.dt23  = FitGaussianCore(h_dt23, 0.90, Form("fit_dt23_%s", tree_name.c_str()));
    result.C_fit = FitGaussianCore(h_C,    0.90, Form("fit_C_%s",    tree_name.c_str()));
    // ---- Salvataggio ----
    outfile->cd();
    tree->Write();
    h_dt12->Write();
    h_dt13->Write();
    h_dt23->Write();
    h_C->Write();

    return result;
}


// ==========================================================================
//  SEZIONE 7: FUNZIONE PRINCIPALE DI CALIBRAZIONE
// ==========================================================================

/// TOF_Calibration(): funzione principale. Processa tutti i dataset,
/// assembla la retta di calibrazione, e genera l'output grafico.
///
/// ARGOMENTI:
///   folder  — percorso della cartella contenente i file XML
///             (es. "/home/lux_n/TOF_DRS/")
///   outname — nome del file ROOT di output
///             (default: "TOF_Calibration_output.root")
///
/// UTILIZZO:
///   root -l 'TOF_Calibration.cpp("/percorso/dati/")'
///   oppure dalla console ROOT:
///     .L TOF_Calibration.cpp
///     TOF_Calibration("/percorso/dati/")

void TOF_Calibration(const char* folder,
                     const char* outname = "TOF_Calibration_output.root",
                     const char* parallax_file = nullptr) {

    std::cout << "=============================================" << std::endl;
    std::cout << "  TOF Calibration                            " << std::endl;
    std::cout << "  Folder: " << folder                         << std::endl;
    std::cout << "  Output: " << outname                        << std::endl;
    std::cout << "=============================================" << std::endl;

    // Setup stile grafico
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);          // Mostra parametri fit nei canvas
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetTitleSize(0.045, "xy");
    gStyle->SetTitleOffset(1.0, "y");
// --- Risoluzione del file di parallasse ---
    // Priorità: argomento esplicito > costante PARALLAX_FILE > vuoto
    const char* effective_parallax_file = parallax_file
                                         ? parallax_file
                                         : PARALLAX_FILE;

    std::cout << "  Parallax file: "
              << (effective_parallax_file && strlen(effective_parallax_file) > 0
                  ? effective_parallax_file
                  : "(non specificato — uso fallback costante)")
              << std::endl;
    std::cout << "  Tape model:    "
              << (TAPE_MODEL_INCREMENTAL
                  ? Form("incrementale (%.1f mm × √n_step, passo %.0f cm)",
                         SIGMA_TAPE_READ * 10.0, TAPE_STEP_LENGTH)
                  : Form("lettura singola (%.1f mm per x ≠ 0)",
                         SIGMA_TAPE_READ * 10.0))
              << std::endl;
    std::cout << "=============================================" << std::endl;
    // Apri il file di output ROOT
    TFile* fout = new TFile(outname, "RECREATE");
    if (!fout->IsOpen()) {
        std::cerr << "[ERRORE] Impossibile creare " << outname << std::endl;
        return;
    }

    // ---- Vettori per raccogliere i risultati posizione per posizione ----
    std::vector<double> x_pos;        // Posizioni x_k [cm]
    std::vector<double> dx_pos;       // Errori su x [cm]
    std::vector<double> dt12_mean;    // Centri dei fit Δt₁₂ [ns]
    std::vector<double> dt12_err;     // Errori sui centri [ns]
    std::vector<double> dt12_width;   // FWHM Δt₁₂ [ns] (= risoluzione temporale)
    std::vector<double> dt12_werr;    // Errori su FWHM [ns]
    std::vector<double> C_mean;       // Centri dei fit C [ns]
    std::vector<double> C_err;        // Errori sui centri C [ns]
    std::vector<std::string> labels;  // Etichette posizioni (per asse X)

    // ---- TTree riassuntivo con tutti i risultati dei fit ----
    fout->cd();
    TTree* summary = new TTree("summary", "Risultati calibrazione per posizione");
    Float_t s_x, s_dx;
    Float_t s_dt12_mu, s_dt12_mu_err, s_dt12_gamma, s_dt12_gamma_err;
    Float_t s_dt13_mu, s_dt13_mu_err, s_dt13_gamma, s_dt13_gamma_err;
    Float_t s_dt23_mu, s_dt23_mu_err, s_dt23_gamma, s_dt23_gamma_err;
    Float_t s_C_mu, s_C_mu_err, s_C_gamma, s_C_gamma_err;
    Int_t   s_n_good;
    Char_t  s_name[32];

    summary->Branch("x",            &s_x,              "x/F");
    summary->Branch("dx",           &s_dx,             "dx/F");
    summary->Branch("dt12_mu",      &s_dt12_mu,        "dt12_mu/F");
    summary->Branch("dt12_mu_err",  &s_dt12_mu_err,    "dt12_mu_err/F");
    summary->Branch("dt12_gamma",   &s_dt12_gamma,     "dt12_gamma/F");
    summary->Branch("dt12_gamma_err", &s_dt12_gamma_err, "dt12_gamma_err/F");
    summary->Branch("dt13_mu",      &s_dt13_mu,        "dt13_mu/F");
    summary->Branch("dt13_mu_err",  &s_dt13_mu_err,    "dt13_mu_err/F");
    summary->Branch("dt13_gamma",   &s_dt13_gamma,     "dt13_gamma/F");
    summary->Branch("dt13_gamma_err", &s_dt13_gamma_err, "dt13_gamma_err/F");
    summary->Branch("dt23_mu",      &s_dt23_mu,        "dt23_mu/F");
    summary->Branch("dt23_mu_err",  &s_dt23_mu_err,    "dt23_mu_err/F");
    summary->Branch("dt23_gamma",   &s_dt23_gamma,     "dt23_gamma/F");
    summary->Branch("dt23_gamma_err", &s_dt23_gamma_err, "dt23_gamma_err/F");
    summary->Branch("C_mu",         &s_C_mu,           "C_mu/F");
    summary->Branch("C_mu_err",     &s_C_mu_err,       "C_mu_err/F");
    summary->Branch("C_gamma",      &s_C_gamma,        "C_gamma/F");
    summary->Branch("C_gamma_err",  &s_C_gamma_err,    "C_gamma_err/F");
    summary->Branch("n_good",       &s_n_good,         "n_good/I");
    summary->Branch("name",         s_name,            "name[32]/C");

// ==================================================================
    //  CALIBRAZIONE k_PMT (slew rate ↔ ampiezza)
    // ==================================================================
    // Eseguita PRIMA del loop principale: serve per recuperare gli eventi
    // clippati durante il processing. Usa solo i dataset con |x| ≤
    // CALIB_K_X_ABS_MAX (per default 84 cm) per garantire abbondanza di
    // eventi NON clippati su entrambi i PMT.
    //
    // Effetto: popola le variabili globali gK_PMT[0..1] e setta
    // gK_calibrated = true. Se la calibrazione fallisce (statistica
    // insufficiente), gK_calibrated resta false e nessun evento clippato
    // verrà recuperato (comportamento equivalente al codice originale,
    // tranne che per il fatto che è ENABLE_CLIP_CUT a controllare lo scarto).
    CalibrateSlewK(folder, DEFAULT_POSITIONS, fout);


    // ==================================================================
    //  LOOP SULLE POSIZIONI
    // ==================================================================
    for (size_t ip = 0; ip < DEFAULT_POSITIONS.size(); ip++) {

        const PositionInfo &pos = DEFAULT_POSITIONS[ip];

        // Costruisci il percorso completo del file XML
        std::string xml_path = std::string(folder) + "/" + pos.name + ".xml";

        // Verifica che il file esista prima di provare a parsarlo
        std::ifstream test_file(xml_path.c_str());
        if (!test_file.good()) {
            // Prova senza lo slash finale nel folder
            xml_path = std::string(folder) + pos.name + ".xml";
            test_file.open(xml_path.c_str());
            if (!test_file.good()) {
                std::cerr << "[WARNING] File non trovato: " << pos.name
                          << ".xml — salto questa posizione." << std::endl;
                continue;
            }
        }
        test_file.close();

        // Processa il dataset
        DatasetResults dres = ProcessDataset(xml_path.c_str(), pos.name, fout);

        // Controlla che il fit Δt₁₂ sia riuscito
        if (!dres.dt12.fit_ok && dres.n_good < 30) {
            std::cerr << "[WARNING] Fit fallito per " << pos.name
                      << " — salto questa posizione per la calibrazione." << std::endl;
            continue;
        }

        // Raccogli i risultati per la retta di calibrazione
        x_pos.push_back(pos.x_cm);
// Errore su x_k posizione-dipendente (parallasse + metro a nastro)
        double sx_k = ComputeSigmaX(pos.x_cm, effective_parallax_file);
        dx_pos.push_back(sx_k);
        dt12_mean.push_back(dres.dt12.center);
        dt12_err.push_back(dres.dt12.center_err);
        dt12_width.push_back(dres.dt12.width);
        dt12_werr.push_back(dres.dt12.width_err);
        C_mean.push_back(dres.C_fit.center);
        C_err.push_back(dres.C_fit.center_err);
        labels.push_back(pos.name);

        // Riempi il TTree riassuntivo
        s_x  = pos.x_cm;
        s_dx = sx_k;   // riusa lo stesso valore appena calcolato
        s_dt12_mu        = dres.dt12.center;
        s_dt12_mu_err    = dres.dt12.center_err;
        s_dt12_gamma     = dres.dt12.width;
        s_dt12_gamma_err = dres.dt12.width_err;
        s_dt13_mu        = dres.dt13.center;
        s_dt13_mu_err    = dres.dt13.center_err;
        s_dt13_gamma     = dres.dt13.width;
        s_dt13_gamma_err = dres.dt13.width_err;
        s_dt23_mu        = dres.dt23.center;
        s_dt23_mu_err    = dres.dt23.center_err;
        s_dt23_gamma     = dres.dt23.width;
        s_dt23_gamma_err = dres.dt23.width_err;
        s_C_mu           = dres.C_fit.center;
        s_C_mu_err       = dres.C_fit.center_err;
        s_C_gamma        = dres.C_fit.width;
        s_C_gamma_err    = dres.C_fit.width_err;
        s_n_good         = dres.n_good;
        strncpy(s_name, pos.name.c_str(), 31);
        s_name[31] = '\0';
        summary->Fill();

        std::cout << "[RESULT] " << pos.name
                  << ": x = " << pos.x_cm << " cm"
                  << ", Δt₁₂ = " << dres.dt12.center
                  << " ± " << dres.dt12.center_err << " ns"
                  << ", FWHM = " << dres.dt12.width << " ns"
                  << ", C = " << dres.C_fit.center
                  << " ± " << dres.C_fit.center_err << " ns"
                  << std::endl;
    }

    // Salva il TTree riassuntivo
    fout->cd();
    summary->Write();

    // ==================================================================
    //  VERIFICA: abbastanza punti per la retta di calibrazione?
    // ==================================================================
    int npts = (int)x_pos.size();
    if (npts < 3) {
        std::cerr << "[ERRORE] Solo " << npts << " posizioni valide — "
                  << "impossibile fare il fit lineare." << std::endl;
        fout->Close();
        return;
    }

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  Calibrazione con " << npts << " posizioni     " << std::endl;
    std::cout << "=============================================" << std::endl;


    // ==================================================================
    //  CANVAS 1: RETTA DI CALIBRAZIONE Δt₁₂ vs x
    // ==================================================================
    // Plottiamo Δt₁₂ sull'asse Y e la posizione x sull'asse X,
    // perché x è la variabile controllata (indipendente) e Δt₁₂ è misurata.
    // Il fit è: Δt₁₂ = m·x + q   con   m = 2/v_eff  (o -2/v_eff)

    TCanvas *c_cal = new TCanvas("c_calibration",
                                  "Calibrazione #Delta t_{12} vs x",
                                  900, 700);

    // Pannello superiore: dati + fit
    c_cal->Divide(1, 2);

    c_cal->cd(1);
    gPad->SetPad(0.0, 0.35, 1.0, 1.0);
    gPad->SetBottomMargin(0.02);
    gPad->SetTopMargin(0.08);
    gPad->SetGrid(1, 1);

    TGraphErrors *g_cal = new TGraphErrors(npts,
                                           x_pos.data(),      // asse X: posizione [cm]
                                           dt12_mean.data(),   // asse Y: Δt₁₂ [ns]
                                           dx_pos.data(),      // errore X: incertezza posizione
                                           dt12_err.data());   // errore Y: errore dal fit

    g_cal->SetMarkerStyle(20);
    g_cal->SetMarkerSize(1.0);
    g_cal->SetMarkerColor(kBlue + 1);
    g_cal->SetLineColor(kBlue + 1);
    g_cal->SetTitle(";Posizione x [cm];#Delta t_{12} = t_{1} - t_{2} [ns]");
    g_cal->GetXaxis()->SetLabelSize(0);  // Etichette X nascoste (le mostra il pannello residui)

    // Fit lineare: Δt₁₂ = m·x + q
    TF1 *f_lin = new TF1("f_lin", "[0] + [1]*x", -150, 150);
    f_lin->SetParNames("q (offset)", "m (2/v)");
    g_cal->Fit(f_lin, "S R Q");

    g_cal->Draw("AP");

    // Box con parametri di fit e velocità derivata
    double m_fit   = f_lin->GetParameter(1);
    double m_err   = f_lin->GetParError(1);
    double q_fit   = f_lin->GetParameter(0);
    double q_err   = f_lin->GetParError(0);
    double chi2ndf = (f_lin->GetNDF() > 0) ? f_lin->GetChisquare() / f_lin->GetNDF() : -1;

    // Velocità: v = 2/|m|, con propagazione errore: σ_v = 2·σ_m / m²
    double v_eff     = 2.0 / fabs(m_fit);
    double v_eff_err = 2.0 * m_err / (m_fit * m_fit);

    TPaveText *pt_cal = new TPaveText(0.15, 0.62, 0.55, 0.92, "NDC");
    pt_cal->SetFillColorAlpha(kWhite, 0.85);
    pt_cal->SetLineColor(kBlack);
    pt_cal->SetTextFont(42);
    pt_cal->SetTextSize(0.045);
    pt_cal->SetTextAlign(12);
    pt_cal->AddText(Form("Fit: #Delta t_{12} = m #upoint x + q"));
    pt_cal->AddText(Form("m = %.5f #pm %.5f ns/cm", m_fit, m_err));
    pt_cal->AddText(Form("q = %.3f #pm %.3f ns", q_fit, q_err));
    pt_cal->AddText(Form("#chi^{2}/ndf = %.2f / %d", f_lin->GetChisquare(), f_lin->GetNDF()));
    pt_cal->AddText(Form("v_{eff} = %.2f #pm %.2f cm/ns", v_eff, v_eff_err));
    pt_cal->Draw();

    // Pannello inferiore: residui del fit
    c_cal->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.35);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.25);
    gPad->SetGrid(1, 1);

    // Calcolo dei residui: res_i = Δt₁₂_misurato − Δt₁₂_fittato
    std::vector<double> residui(npts);
    for (int i = 0; i < npts; i++) {
        residui[i] = dt12_mean[i] - f_lin->Eval(x_pos[i]);
    }

    TGraphErrors *g_res = new TGraphErrors(npts,
                                           x_pos.data(),
                                           residui.data(),
                                           dx_pos.data(),
                                           dt12_err.data());
    g_res->SetMarkerStyle(20);
    g_res->SetMarkerSize(1.0);
    g_res->SetMarkerColor(kBlue + 1);
    g_res->SetLineColor(kBlue + 1);
    g_res->SetTitle(";Posizione x [cm];Residui [ns]");
    g_res->GetXaxis()->SetTitleSize(0.08);
    g_res->GetXaxis()->SetLabelSize(0.07);
    g_res->GetYaxis()->SetTitleSize(0.08);
    g_res->GetYaxis()->SetLabelSize(0.07);
    g_res->GetYaxis()->SetTitleOffset(0.5);
    g_res->Draw("AP");

    // Linea orizzontale a zero
    double xmin_r = *std::min_element(x_pos.begin(), x_pos.end()) - 10;
    double xmax_r = *std::max_element(x_pos.begin(), x_pos.end()) + 10;
    TLine *l_zero = new TLine(xmin_r, 0, xmax_r, 0);
    l_zero->SetLineColor(kRed);
    l_zero->SetLineStyle(2);
    l_zero->Draw("same");

    c_cal->Update();

    // Salva il canvas
    fout->cd();
    c_cal->Write("Calibration_Dt12_vs_x");


    // ==================================================================
    //  CANVAS 2: RISOLUZIONE TEMPORALE FWHM(Δt₁₂) vs x
    // ==================================================================
    TCanvas *c_res = new TCanvas("c_resolution",
                                  "Risoluzione temporale #Gamma(#Delta t_{12}) vs x",
                                  800, 500);
    c_res->SetGrid(1, 1);

    TGraphErrors *g_width = new TGraphErrors(npts,
                                             x_pos.data(),
                                             dt12_width.data(),
                                             dx_pos.data(),
                                             dt12_werr.data());
    g_width->SetMarkerStyle(21);
    g_width->SetMarkerSize(1.0);
    g_width->SetMarkerColor(kRed + 1);
    g_width->SetLineColor(kRed + 1);
    g_width->SetTitle("Risoluzione temporale vs posizione;"
                      "Posizione x [cm];"
                      "#sigma(#Delta t_{12}) [ns]");
    g_width->Draw("AP");
    c_res->Update();

    fout->cd();
    c_res->Write("Resolution_FWHM_vs_x");


    // ==================================================================
    //  CANVAS 3: COSTANTE C vs x (verifica uniformità)
    // ==================================================================
    TCanvas *c_C = new TCanvas("c_C_vs_x",
                                "Costante C = t_{3} - (t_{1}+t_{2})/2 vs x",
                                800, 500);
    c_C->SetGrid(1, 1);

    TGraphErrors *g_C = new TGraphErrors(npts,
                                         x_pos.data(),
                                         C_mean.data(),
                                         dx_pos.data(),
                                         C_err.data());
    g_C->SetMarkerStyle(22);
    g_C->SetMarkerSize(1.2);
    g_C->SetMarkerColor(kGreen + 2);
    g_C->SetLineColor(kGreen + 2);
    g_C->SetTitle("Costante di offset C vs posizione;"
                  "Posizione x [cm];"
                  "C = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]");

    // Fit con una costante per verificare l'uniformità
    TF1 *f_const = new TF1("f_const", "[0]", -150, 150);
    f_const->SetParName(0, "C");
    g_C->Fit(f_const, "S R Q");
    g_C->Draw("AP");

    // Box con il valore di C
    double C_offset     = f_const->GetParameter(0);
    double C_offset_err = f_const->GetParError(0);
    double C_chi2ndf    = (f_const->GetNDF() > 0) ?
                          f_const->GetChisquare() / f_const->GetNDF() : -1;

    TPaveText *pt_C = new TPaveText(0.15, 0.72, 0.55, 0.92, "NDC");
    pt_C->SetFillColorAlpha(kWhite, 0.85);
    pt_C->SetLineColor(kBlack);
    pt_C->SetTextFont(42);
    pt_C->SetTextSize(0.045);
    pt_C->SetTextAlign(12);
    pt_C->AddText(Form("C = %.3f #pm %.3f ns", C_offset, C_offset_err));
    pt_C->AddText(Form("#chi^{2}/ndf = %.2f / %d",
                       f_const->GetChisquare(), f_const->GetNDF()));
    pt_C->Draw();

    c_C->Update();
    fout->cd();
    c_C->Write("C_vs_x");


    // ==================================================================
    //  CANVAS 4: GALLERIA ISTOGRAMMI Δt₁₂ (tutti i dataset)
    // ==================================================================
    {
        // Determina il layout della griglia
        int ncols = (npts <= 4) ? npts : (npts <= 9 ? 3 : 4);
        int nrows = (npts + ncols - 1) / ncols;

        TCanvas *c_gallery = new TCanvas("c_dt12_gallery",
                                          "Distribuzioni #Delta t_{12}",
                                          350 * ncols, 300 * nrows);
        c_gallery->Divide(ncols, nrows);

        for (int i = 0; i < npts; i++) {
            c_gallery->cd(i + 1);
            gPad->SetGrid(1, 1);
            // Recupera l'istogramma dal file
            TH1D *h = (TH1D*)fout->Get(Form("h_dt12_%s", labels[i].c_str()));
            if (h) {
                h->SetLineColor(kBlue + 1);
                h->SetFillColorAlpha(kBlue - 9, 0.3);
                h->Draw();
            }
        }
        c_gallery->Update();
        fout->cd();
        c_gallery->Write("Gallery_Dt12");
    }


    // ==================================================================
    //  CANVAS 5: GALLERIA ISTOGRAMMI C (tutti i dataset)
    // ==================================================================
    {
        int ncols = (npts <= 4) ? npts : (npts <= 9 ? 3 : 4);
        int nrows = (npts + ncols - 1) / ncols;

        TCanvas *c_gal_C = new TCanvas("c_C_gallery",
                                        "Distribuzioni C",
                                        350 * ncols, 300 * nrows);
        c_gal_C->Divide(ncols, nrows);

        for (int i = 0; i < npts; i++) {
            c_gal_C->cd(i + 1);
            gPad->SetGrid(1, 1);
            TH1D *h = (TH1D*)fout->Get(Form("h_C_%s", labels[i].c_str()));
            if (h) {
                h->SetLineColor(kGreen + 2);
                h->SetFillColorAlpha(kGreen - 9, 0.3);
                h->Draw();
            }
        }
        c_gal_C->Update();
        fout->cd();
        c_gal_C->Write("Gallery_C");
    }


    // ==================================================================
    //  RIEPILOGO FINALE
    // ==================================================================
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  RISULTATI CALIBRAZIONE                      " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Velocità effettiva: v_eff = "
              << Form("%.2f ± %.2f", v_eff, v_eff_err) << " cm/ns" << std::endl;
    std::cout << "  Offset Δt₁₂:       q     = "
              << Form("%.3f ± %.3f", q_fit, q_err) << " ns" << std::endl;
    std::cout << "  Costante C:        C     = "
              << Form("%.3f ± %.3f", C_offset, C_offset_err) << " ns" << std::endl;
    std::cout << "  χ²/ndf (retta):           = "
              << Form("%.2f", chi2ndf) << std::endl;
    std::cout << "  χ²/ndf (C = cost):        = "
              << Form("%.2f", C_chi2ndf) << std::endl;
    std::cout << "\n  Output salvato in: " << outname << std::endl;
    std::cout << "  TTree 'summary' con tutti i parametri dei fit." << std::endl;
    std::cout << "=============================================" << std::endl;

    // Chiudi il file di output
    fout->Close();

    std::cout << "\n[DONE] Calibrazione completata." << std::endl;
}
