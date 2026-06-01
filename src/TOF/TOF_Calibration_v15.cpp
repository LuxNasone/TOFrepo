// ==========================================================================
//  TOF_Calibration_v13.cpp — Calibrazione dell'apparato Time-of-Flight
// ==========================================================================
//
//  SCOPO:
//    Determinare la velocità effettiva della luce nella barra (v_eff) e la
//    costante di offset C necessaria per calcolare il tempo di volo.
//
//  FISICA:
//    La barra scintillatrice (BC408, 280 cm) è letta da PMT1 e PMT2 ai due
//    estremi. Un terzo scintillatore mobile (PMT3) è posizionato in contatto
//    con la barra (Guida A) in diverse posizioni note x_k.
//
//    Definendo x la coordinata lungo la barra dal centro (x ∈ [-140,140] cm):
//      t1 = t0 + d1(x)/v + delta1
//      t2 = t0 + d2(x)/v + delta2
//      t3 = t0 + TOF     + delta3        (TOF ~ 0 in Guida A)
//
//    1) RETTA DI CALIBRAZIONE: Dt12 = t1-t2 = (2/v)*x + cost
//       fit lineare Dt12 vs x -> pendenza m = 2/v -> v_eff = 2/|m|
//    2) COSTANTE C: C = t3 - (t1+t2)/2  (~ costante in x entro le incertezze)
//       In fase di misura TOF:  TOF = t3 - (t1+t2)/2 - C(x)
//
// ==========================================================================
//  VERSIONE 13 — PARADIGMA TOT-POLINOMIALE UNICO
// ==========================================================================
//
//  Rispetto alla v12 questa versione adotta UN SOLO modello ufficiale per la
//  ricostruzione dei segnali clippati (saturati dalla DRS4):
//
//      A(TOT) = a0 + a1*TOT + a2*TOT^2 + a3*TOT^3 ,   con a0 = q_thr = 50 mV
//
//  E' stato RIMOSSO completamente:
//    - il recupero clippati via "slew rate" (CalibrateSlewK, RecoverClippedCFD,
//      tutte le costanti SLEW_*, le globali gK/gQ, i campi slew_* di ChannelData,
//      il fit lineare del fronte in AnalyzeChannel, il canvas slew);
//    - il vecchio modello esponenziale A = q + p0*(exp(p1*TOT)-1);
//    - la calibrazione PURE_TOT separata (CalibrateTOT_Pure): ora la modalita'
//      diagnostica PURE_TOT usa ESATTAMENTE gli stessi coefficienti polinomiali
//      della calibrazione principale;
//    - il flag USE_POLY_TOT (il polinomio cubico e' ora l'unico modello);
//    - le costanti legacy OSC_VMAX_THRESH / OSC_USE_LEGACY (il rilevamento di
//      oscillazione e' solo lo Schmitt trigger con isteresi).
//
//  La calibrazione produce TRE metodi diagnostici, che CONDIVIDONO lo stesso
//  modello A(TOT) e differiscono solo per QUALI segnali vengono ricostruiti:
//    - noclip     : solo segnali NON clippati, nessun recupero  (riferimento)
//    - hybrid_tot : CFD classico sui non clippati + ricostruzione TOT sui soli
//                   clippati                                    (METODO FISICO)
//    - pure_tot   : ricostruzione TOT applicata a TUTTI i segnali (diagnostica
//                   di simmetria: verifica se il TOT introduce bias in C(x))
//
//  Il TTree "fit_params" del file ROOT di output e' il CONTRATTO con la macro
//  di analisi TOF_Analysis_v9.cpp: contiene m, q, v_eff, il modello quadratico
//  C(x) = p0+p1*x+p2*x^2 (principale), il modello costante <C> (confronto),
//  il vecchio C(x) lineare (legacy), e i coefficienti polinomiali A(TOT) per
//  PMT1/2/3, tutti relativi al metodo principale hybrid_tot.
//
//  UTILIZZO:
//    root -l 'TOF_Calibration_v13.cpp("/percorso/cartella/dati/")'
//    I file XML hanno nomi tipo "N130.xml", "X0.xml", ... (N = posizioni
//    negative, X = positive, in cm dal centro barra).
//
//  DIPENDENZE: ROOT 6+
//  AUTORI: Luca (sviluppo) + Claude (assistenza)
//  DATA:   Maggio 2026
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
#include <TH2.h>
#include <TProfile.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>

// --- Header aggiuntivi per la cache di parsing (Punto 1) ---
#include <TSystem.h>     // gSystem->mkdir() per creare la cartella di cache
#include <TParameter.h>  // marcatore di versione del formato della cache
#include <sys/stat.h>    // stat() per leggere il tempo di modifica dei file
#include <TList.h>       // gestione lista funzioni dell'istogramma (Punto 2)
#include <TMultiGraph.h> // overlay di piu' TGraph con auto-scaling assi (Punto 4)
// ==========================================================================
//  SEZIONE 1: COSTANTI E PARAMETRI CONFIGURABILI
// ==========================================================================

// --- Parametri hardware DRS4 ---
const int MAX_SAMPLES  = 1024;   // Celle del chip DRS4 (fisso, NON modificare)
const int MAX_CHANNELS = 4;      // Canali della DRS4 Eval Board

// --- Parametri di analisi della forma d'onda ---
const int    NBL_SAMPLES   = 50;     // Campioni per il calcolo della baseline (~10 ns a 5 GS/s)
// FRAZIONE CFD: definizione operativa dell'istante di arrivo del segnale.
// DEVE essere identica in TOF_Calibration_v13.cpp e TOF_Analysis_v9.cpp:
// la retta Dt12(x) e la funzione C(x) dipendono da questa scelta, quindi una
// frazione diversa fra calibrazione e analisi introdurrebbe un bias sistematico.
const double CFD_FRACTION  = 0.15;   // Frazione CFD (0.15 = soglia al 15% dell'ampiezza)
const double NOISE_THRESH  = 5.0;    // Soglia minima ampiezza [mV] per dichiarare un impulso

// CLIP_V_LO: tensione minima [mV] sotto cui il segnale e' considerato saturato.
//   La DRS4 satura tipicamente attorno a -500 mV: soglia posta con margine.
const double CLIP_V_LO     = -499.0;

// ============================================================
//  RILEVAMENTO OSCILLAZIONE: discriminatore a isteresi (Schmitt trigger)
// ============================================================
//  Un impulso fisico di scintillazione e' MONOPOLARE NEGATIVO. Contiamo quante
//  volte il segnale ATTRAVERSA IN SALITA la soglia OSC_V_HIGH dopo il picco
//  negativo: un singolo overshoot post-clip conta 1, un'oscillazione vera con
//  N cicli conta N. L'isteresi (un nuovo crossing si conta solo dopo che il
//  segnale e' ridisceso sotto OSC_V_LOW) sopprime il ringing fine.
//  NOTA: in v13 questo e' l'UNICO criterio di oscillazione (rimosso il vecchio
//  criterio "v_max > soglia" che era impreciso).
const double OSC_V_HIGH      = 50.0;   // [mV] soglia alta di crossing
const double OSC_V_LOW       = 25.0;   // [mV] soglia bassa di isteresi
const int    OSC_NCROSS_MAX  = 3;      // >= 3 crossing -> oscillazione vera -> scarto

// Flag globali per abilitare/disabilitare i tagli di qualita'.
const bool ENABLE_CFD_CUT  = true;   // scarta se manca un tempo utilizzabile su un canale
const bool ENABLE_OSC_CUT  = true;   // scarta se un canale e' oscillante

// ============================================================
//  PARAMETRI TIME-OVER-THRESHOLD (TOT) — recupero clippati
// ============================================================
//  Per impulsi a forma fissa A*g(t), la durata sopra una soglia q (TOT) e' una
//  funzione monotona crescente di A. Calibrando A vs TOT su segnali NON
//  clippati (dove l'ampiezza vera e' nota) si ottiene A(TOT); applicando A(TOT)
//  ai clippati se ne stima l'ampiezza e quindi il tempo CFD.

// Soglie [mV] a cui si calcola il TOT per ogni evento (multi-soglia per studi
// di sistematica). Tutte >= ~10 sigma_baseline e <= |CLIP_V_LO| - TOT_CLIP_MARGIN.
const int    NTOT_THR        = 5;
const double TOT_THR_MV[NTOT_THR] = { 50.0, 100.0, 150.0, 200.0, 300.0 };

// Indice (0-based) della soglia di RIFERIMENTO usata per:
//   - il fit di calibrazione A(TOT);
//   - il recupero ampiezza dei clippati (hybrid_tot);
//   - la ricostruzione PURE_TOT.
// In v13 e' 0 -> 50 mV: una sola soglia, una sola calibrazione, usata da tutti
// i metodi (vedi guida sez. 3.1). Corrisponde a -50 mV sul segnale negativo.
const int    TOT_CALIB_THR_IDX = 0;   // -> TOT_THR_MV[0] = 50 mV

// Margine sopra |CLIP_V_LO| per dichiarare una soglia "troppo vicina al clip":
// se q > |CLIP_V_LO| - TOT_CLIP_MARGIN il TOT a quella soglia non viene
// calcolato (il fall crossing cadrebbe nel plateau saturato).
const double TOT_CLIP_MARGIN = 5.0;  // [mV]

// Numero minimo di campioni tra rise e fall per accettare il TOT.
const int    TOT_NMIN_SAMPLES = 3;

// Cap di estrapolazione: A_rec accettata solo se <= TOT_AREC_MAX_FACTOR * Amax,
// dove Amax e' l'ampiezza massima del campione di calibrazione per quel PMT.
// Protegge contro TOT anomalamente lunghi (pile-up, deriva di baseline).
const double TOT_AREC_MAX_FACTOR = 3.0;

// ============================================================
//  MODELLO A(TOT): POLINOMIO CUBICO (unico modello ufficiale)
// ============================================================
//      A(TOT) = a0 + a1*TOT + a2*TOT^2 + a3*TOT^3
//  con a0 = q_thr = TOT_THR_MV[TOT_CALIB_THR_IDX] FISSO nel fit.
//  Significato fisico del vincolo a0 = q_thr: un segnale che resta sopra
//  soglia per tempo nullo (TOT -> 0) e' appena tangente alla soglia, quindi la
//  sua ampiezza e' ~ q_thr.
const int    TOT_POLY_DEGREE  = 3;                    // grado del polinomio (cubico)
const int    TOT_POLY_NPAR    = TOT_POLY_DEGREE + 1;  // numero di coefficienti = 4
const int    TOT_POLY_MAX_DEG = 5;                    // dimensione array statici (NON superare)
const bool   FIX_A0_TO_QTHR   = true;                 // se true, a0 viene fissato a q_thr nel fit

// ============================================================
//  FIT POLINOMIALE A PROFILO MEDIANO
// ============================================================
//  Invece di fittare il polinomio sulla NUVOLA grezza (TOT, A) — che ha enorme
//  dispersione intrinseca — lo si fitta sul PROFILO MEDIANO: i punti puliti
//  vengono divisi in bin EQUIPOPOLATI di TOT e per ogni bin si estrae la
//  mediana di TOT, la mediana di A e l'errore della mediana (MAD*1.4826/sqrt(N)).
const int    TOT_PROFILE_NBINS    = 30;   // numero di bin equipopolati del profilo
const int    TOT_PROFILE_MIN_N    = 5;    // eventi minimi per dichiarare valido un bin
const double TOT_PROFILE_AERR_MIN = 1.0;  // [mV] floor sull'errore della mediana di A

// ============================================================
//  ISTOGRAMMI DIAGNOSTICI dA = A_rec(TOT) - A_vera
// ============================================================
//  Popolati SOLO su eventi NON clippati (dove A_vera e' ben definita): servono
//  a quantificare bias ed errore sistematico del recupero TOT.
const int    NBINS_DA   = 200;       // bin per gli istogrammi 1D di dA
const double DA_LO      = -100.0;    // [mV] limite inferiore dA
const double DA_HI      =  100.0;    // [mV] limite superiore dA
const int    NBINS_DA2_A  = 50;      // bin asse X (A_vera) dei TH2
const double DA2_A_LO     = 0.0;     // [mV]
const double DA2_A_HI     = 500.0;   // [mV]
const int    NBINS_DA2_DA = 100;     // bin asse Y (dA) dei TH2
const double DA2_DA_LO    = -100.0;  // [mV]
const double DA2_DA_HI    = 100.0;   // [mV]
// Asse X del TH2 dA vs TOT
const int    NBINS_DA2_TOT = 60;     // bin asse X (TOT) [ns]
const double DA2_TOT_LO    = 0.0;
const double DA2_TOT_HI    = 60.0;

// --- Parametri di binning istogrammi temporali ---
// Range ADATTIVO per gli istogrammi Dt12/Dt13/Dt23 (e C).
// Il range viene calcolato per ogni posizione x_k a partire dalla mediana
// e dalla MAD (Median Absolute Deviation) dei dati raccolti nel primo pass.
// Questo garantisce che i bin siano concentrati attorno al picco centrale,
// indipendentemente dalla posizione del PMT3 sulla barra.
const int    NBINS_DT         = 50;   // bin per gli istogrammi Dt (range adattivo)
const double DT_RANGE_NSIGMA  = 5.0;   // semi-intervallo in unita' di sigma robusto (MAD)
const double DT_RANGE_MIN     = 3.0;   // semi-intervallo minimo [ns]
const int    NBINS_C   = 50;       // bin per gli istogrammi di C
const double C_LO      = -16.0;     // [ns]
const double C_HI      =   -8.0;     // [ns]

// --- Errore sulla posizione x ---
const double DX_POSITION = 3.0;     // [cm] meta' estensione PMT3 lungo la barra (fallback)

// ============================================================
//  PARAMETRI ERRORE SULLA POSIZIONE (modello posizione-dipendente)
// ============================================================
//  sigma_x^2(x_k) = sigma_parallax^2(x_k) + sigma_tape^2(x_k)
//  - sigma_parallax letto da file ROOT del Monte Carlo di accettanza;
//  - sigma_tape errore di lettura del metro, modello "lettura singola" oppure
//    "letture incrementali" (sigma_read^2 accumula in quadratura per ogni passo).
const double SIGMA_TAPE_READ      = 0.1;   // [cm] errore di lettura del metro
const double TAPE_REFERENCE_X     = 0.0;   // [cm] posizione di riferimento (errore nullo)
const bool   TAPE_MODEL_INCREMENTAL = true;// false = lettura singola, true = incrementale
const double TAPE_STEP_LENGTH     = 28.0;  // [cm] lunghezza di un passo (solo modello incrementale)

// File ROOT del MC di parallasse. Stringa vuota = usa il fallback costante.
const char*  PARALLAX_FILE      = "";      // path file ROOT MC (vuoto = disabilitato)
const double PARALLAX_FALLBACK  = 3.95;    // [cm] valore costante se file non disponibile

// --- Lista posizioni di calibrazione ---
struct PositionInfo {
    std::string name;    // Nome del file XML (senza .xml)
    double      x_cm;    // Posizione dal centro della barra [cm]
};

std::vector<PositionInfo> DEFAULT_POSITIONS = {
    {"N130", -130.0},
    {"N112", -112.0},
    {"N98",   -98.0},
    {"N84",   -84.0},
    {"N70",   -70.0},
    {"N56",   -56.0},
    {"N42",   -42.0},
    {"N28",   -28.0},
    {"N14",   -14.0},
    {"X0",      0.0},
    {"X14",    14.0},
    {"X28",    28.0},
    {"X42",    42.0},
    {"X56",    56.0},
    {"X70",    70.0},
    {"X84",    84.0},
    {"X98",    98.0},
    {"X112",  112.0},
    {"X130",  130.0}
};

// Range di posizioni usate per CALIBRARE A(TOT).
// Solo dataset con |x_k| <= TOT_CALIB_X_ABS_MAX contribuiscono al fit del
// polinomio: a queste posizioni entrambi i PMT vedono abbondanti eventi NON
// clippati su un range di ampiezze sufficiente. (In v12 questa costante si
// chiamava CALIB_K_X_ABS_MAX e riguardava lo slew-rate, ora rimosso.)
const double TOT_CALIB_X_ABS_MAX = 56.0;  // [cm]


// ==========================================================================
//  PARAMETRI DELLA CACHE DI PARSING (XML -> ROOT)
// ==========================================================================
//  Il parsing testuale degli XML DRS4 e' la fase piu' lenta del programma.
//  Per evitarlo a ogni run, le FORME D'ONDA GREZZE di ogni dataset vengono
//  salvate in un file ROOT di cache (uno per posizione) dentro la sottocartella
//  PARSE_CACHE_SUBDIR, accanto agli XML. Alle run successive, se la cache esiste
//  ed e' piu' recente dell'XML, gli eventi vengono caricati da li'.
//
//  SCELTA DI PROGETTO: nella cache si salvano SOLO le forme d'onda (time[],
//  voltage[]) e l'header, NON le grandezze calcolate da AnalyzeChannel. Al
//  caricamento AnalyzeChannel viene RIESEGUITO. Cosi' la cache resta valida
//  anche se si cambiano CFD_FRACTION, le soglie TOT o l'algoritmo di analisi:
//  cio' che e' congelato e' solo l'input grezzo, non la fisica derivata.
const bool  USE_PARSE_CACHE    = true;   // false = disattiva del tutto la cache
const bool  FORCE_REPARSE      = false;  // true  = ignora la cache e riparsa (poi la riscrive)
const char* PARSE_CACHE_SUBDIR = "parsed_cache";  // sottocartella della cache
const int   PARSE_CACHE_FORMAT = 1;      // versione del formato (cambiala per invalidare)


// ==========================================================================
//  RANGE LINEARE RISTRETTO PER LA VELOCITA' (Punto 3)
// ==========================================================================
//  Finestra di posizioni in cui il rise time di ENTRAMBI i PMT e' ben
//  approssimato da una relazione lineare in x (vedi rise_range.png). Fuori da
//  questa finestra, attenuazione/dispersione della luce e clipping curvano
//  Dt12(x) e gonfiano il time-walk. Restringere il fit a questo intervallo da'
//  una stima di v_eff meno contaminata dalle non-linearita' di bordo; il
//  confronto con la stima globale e' una diagnostica della sistematica.
const double LINRANGE_X_LO = -84.0;   // [cm] estremo inferiore (incluso)
const double LINRANGE_X_HI =  70.0;   // [cm] estremo superiore (incluso)

// ==========================================================================
//  GLOBALI DI CALIBRAZIONE TOT (popolate da CalibrateTOT)
// ==========================================================================
//  Modello valutato: A(TOT) = sum_k gTOT_poly_PMT[ch][k] * TOT^k  (k = 0..3).
//  In v13 i tre metodi (noclip / hybrid_tot / pure_tot) condividono questi
//  stessi coefficienti: non esistono piu' insiemi separati per il PURE_TOT.
//  Gli array sono dimensionati a MAX_CHANNELS: gli indici 0,1,2 = PMT1,2,3.
//  PMT3 (indice 2) viene calibrato anche lui perche' entra nel calcolo di
//  C = t3 - (t1+t2)/2 e una sua eventuale ricostruzione modifica il TOF.
static double gTOT_poly_PMT       [MAX_CHANNELS][TOT_POLY_MAX_DEG + 1] = {{0.0}};
static double gTOT_poly_err_PMT   [MAX_CHANNELS][TOT_POLY_MAX_DEG + 1] = {{0.0}};
static double gTOT_poly_chi2ndf_PMT[MAX_CHANNELS]                      = {0.0};
static double gTOT_Amax_PMT       [MAX_CHANNELS]                       = {0.0}; // A massima nel calib (cap)
static bool   gTOT_calibrated     [MAX_CHANNELS]                       = {false};// flag per canale
static int    gTOT_poly_degree    = TOT_POLY_DEGREE;  // grado effettivo (per coerenza con l'analisi)

// ============================================================
//  ISTOGRAMMI DIAGNOSTICI AGGREGATI dA (su tutte le posizioni)
// ============================================================
//  Vivono a scope di file perche' vengono riempiti dentro ProcessDataset (una
//  posizione alla volta) ma aggregati su TUTTE le posizioni. Creati in
//  TOF_Calibration prima del loop, popolati in ProcessDataset, scritti alla fine.
static TH1D* g_h_dA       [MAX_CHANNELS] = {nullptr, nullptr, nullptr, nullptr};
static TH2F* g_h_dA_vs_A  [MAX_CHANNELS] = {nullptr, nullptr, nullptr, nullptr};
static TH2F* g_h_dA_vs_TOT[MAX_CHANNELS] = {nullptr, nullptr, nullptr, nullptr};


// ==========================================================================
//  SEZIONE 2: STRUTTURE DATI
// ==========================================================================

/// Dati di un singolo canale di un singolo evento.
/// Contiene la forma d'onda grezza e le grandezze estratte dall'analisi.
struct ChannelData {
    int    nsamples;                  // Numero di campioni letti (tipicamente 1024)
    float  time[MAX_SAMPLES];         // Tempi dei campioni [ns], calibrati dalla DRS4
    float  voltage[MAX_SAMPLES];      // Tensione dei campioni [mV]

    // ---- Grandezze estratte da AnalyzeChannel() ----
    double baseline;      // Media tensione nei primi NBL_SAMPLES [mV]
    double baseline_rms;  // Deviazione standard della baseline [mV]
    double amplitude;     // Ampiezza: baseline - V_min [mV] (positiva)
    double v_min;         // Tensione minima (picco negativo) [mV]
    double v_max;         // Tensione massima [mV] (solo diagnostica)
    double t_min;         // Tempo del campione con V minima [ns]
    double t_cfd;         // Tempo CFD interpolato [ns] (-999 se non valido)

    // ---- Flag di qualita' ----
    bool   has_pulse;       // true se amplitude > NOISE_THRESH
    bool   cfd_ok;          // true se il CFD ha trovato un crossing valido
    bool   is_clipped;      // true se V_min < CLIP_V_LO (saturazione DRS4)
    bool   is_oscillating;  // true se n_pos_crossings >= OSC_NCROSS_MAX

    // ---- Diagnostica del rilevamento oscillazione ----
    int    n_pos_crossings; // crossing in salita di OSC_V_HIGH dopo il picco (con isteresi)

    // ---- Time-Over-Threshold multi-soglia ----
    // tot_q[k]    = durata del segnale sopra la soglia TOT_THR_MV[k] [ns]
    // tot_rise/fall[k] = istanti dei due crossing [ns]
    // tot_ok[k]   = true se entrambi i crossing sono validi e separati di
    //               almeno TOT_NMIN_SAMPLES campioni.
    double tot_q   [NTOT_THR];
    double tot_rise[NTOT_THR];
    double tot_fall[NTOT_THR];
    bool   tot_ok  [NTOT_THR];

    // ---- Recupero clippati via TOT polinomiale (metodo principale hybrid_tot) ----
    // Popolati SOLO sui canali clippati da RecoverClippedCFD_TOT().
    double amplitude_tot;     // ampiezza ricostruita A(TOT) [mV]
    double t_cfd_rec_tot;     // tempo CFD ricalcolato sulla soglia f*A(TOT) [ns]
    bool   cfd_recovered_tot; // true se il recupero del clippato e' riuscito

    // ---- Metodo diagnostico PURE_TOT (ricostruzione TOT su TUTTI i canali) ----
    // Popolati da ComputePureTOTTime() per OGNI canale con tot_ok valido,
    // a prescindere dal clipping. Usano lo STESSO polinomio A(TOT).
    double amplitude_pure_tot; // ampiezza ricostruita via TOT [mV]
    double t_cfd_pure_tot;     // tempo CFD su soglia f*A_PURE_TOT [ns]
    bool   cfd_pure_tot_ok;    // true se la ricostruzione PURE_TOT e' valida
    // ---- Rise Time T90-T10 (diagnostica dipendenza posizione) ----
    // Tempo impiegato dal segnale per passare dal 10% al 90% dell'ampiezza,
    // misurato sul fronte di salita (leading edge). Calcolato da AnalyzeChannel.
    // Per segnali negativi: T10 e' l'istante in cui V scende sotto bl-0.10*A,
    // T90 l'istante in cui scende sotto bl-0.90*A. Rise time = T90 - T10 > 0.
    // Valido solo se entrambi i crossing sono trovati sull'onda digitalizzata.
    double rise_time;     // T90 - T10 [ns] (-999 se non valido)
    bool   rise_time_ok;  // true se entrambi i crossing sono stati trovati
    // ---- Crossing a soglie fisse 10%, 30% (per correzione Pietro&Rick) ----
    // t_10, t_30 sono gli istanti in cui il segnale negativo attraversa
    // sul fronte di salita le soglie:
    //   V_10 = baseline - 0.10 * amp     (vicino alla baseline, alto SNR)
    //   V_30 = baseline - 0.30 * amp     (frazione tipica del CFT classico)
    // La differenza dT = t_30 - t_10 e' un rise time "effettivo" sulla
    // porzione di fronte rilevante per il timing CFD, piu' stabile del
    // canonico T_90 - T_10 (che include anche la regione vicino al picco
    // dove ringing e overshoot rendono T_90 rumoroso).
    // Sono i mattoni elementari della correzione del time-walk residuo
    // dovuto al rise time variabile (vedi Pietro&Rick eq. 18-19).
    double t_10;       // istante crossing al 10% sul leading edge [ns]
    double t_30;       // istante crossing al 30% sul leading edge [ns]
    bool   t10_ok;     // true se t_10 e' stato trovato
    bool   t30_ok;     // true se t_30 e' stato trovato
};

/// Dati completi di un evento DRS4 (header + canali).
struct EventData {
    int         serial;                    // Numero progressivo (parte da 1)
    std::string timestamp;                 // Data/ora "YYYY/MM/DD HH:MM:SS.mmm"
    int         trigger_cell;              // Cella dove l'onda domino si e' fermata
    int         board_serial;              // Numero seriale della scheda
    int         scaler[MAX_CHANNELS];      // Rate scaler per canale [Hz]
    int         nchannels;                 // Numero canali presenti (1-4)
    int         channel_ids[MAX_CHANNELS]; // Numeri dei canali (1-based)
    ChannelData ch[MAX_CHANNELS];          // Dati dei canali (0-based)
};

/// Risultato del fit di una distribuzione (centro + larghezza).
struct FitResult {
    double center;       // Centro della distribuzione dal fit [ns]
    double center_err;   // Errore sul centro [ns]
    double width;        // Larghezza (sigma per la gaussiana di core) [ns]
    double width_err;    // Errore sulla larghezza [ns]
    double chi2_ndf;     // chi2/ndf del fit (diagnostica)
    int    nentries;     // Numero di eventi nell'istogramma
    bool   fit_ok;       // true se il fit e' convergito con risultati sensati
};

/// Risultati dei fit di un SINGOLO metodo di ricostruzione, per una posizione.
struct MethodResults {
    FitResult dt12;    // fit Dt12 = t1 - t2
    FitResult dt13;    // fit Dt13 = t1 - t3
    FitResult dt23;    // fit Dt23 = t2 - t3
    FitResult C;       // fit C = t3 - (t1+t2)/2
    int       n_good;  // eventi buoni per questo metodo
};

/// Risultati completi di una posizione x_k: i tre metodi affiancati.
struct DatasetResults {
    double      x;     // posizione [cm]
    double      dx;    // incertezza sulla posizione [cm]
    std::string name;  // nome del dataset (= nome file XML senza estensione)

    MethodResults noclip;      // solo segnali non clippati (riferimento)
    MethodResults hybrid_tot;  // CFD + recupero TOT sui clippati (METODO FISICO)
    MethodResults pure_tot;    // ricostruzione TOT su tutti i canali (diagnostica)
    // METODO "corrected": parte da hybrid_tot e sottrae il bias del rise
    // time variabile alla Pietro&Rick. L'evento e' valido per questo metodo
    // SOLO SE: (i) e' valido per hybrid_tot, (ii) NESSUN canale e' clippato
    // (perche' su un clippato amp e' sottostimata -> dT = t30 - t10 non
    // rappresenta il rise time fisico), (iii) t_10 e t_30 sono entrambi
    // disponibili su tutti e 3 i PMT.
    MethodResults corrected;   // hybrid_tot con correzione P&R applicata
};


// ==========================================================================
//  SEZIONE 3: ANALISI DELLA FORMA D'ONDA
// ==========================================================================

/// AnalyzeChannel(): estrae da un canale baseline, ampiezza, tempo CFD, i flag
/// di qualita' (is_clipped, is_oscillating) e il Time-Over-Threshold multi-soglia.
///
/// FLUSSO (v13, senza slew-rate):
///   1. Baseline (media + RMS dei primi NBL_SAMPLES campioni).
///   2. Ricerca di v_min (picco negativo) e v_max.
///   3. Flag is_clipped tramite CLIP_V_LO.
///   4. Ampiezza = baseline - v_min, flag has_pulse tramite NOISE_THRESH.
///   5. Rilevamento oscillazione: Schmitt trigger con isteresi (OSC_V_HIGH,
///      OSC_V_LOW, OSC_NCROSS_MAX). Unico criterio (rimosso il vecchio "v_max").
///   6. Tempo CFD standard (per i segnali non clippati).
///   7. Time-Over-Threshold per tutte le soglie in TOT_THR_MV[].
///
/// Il CFD e il TOT vengono calcolati SEMPRE (anche sui clippati e oscillanti)
/// per conservare diagnostica nel TTree; sara' ProcessDataset() a decidere quali
/// eventi/canali usare in base ai flag.
void AnalyzeChannel(ChannelData &cd) {

    // --- Inizializzazione a valori "non calcolato" ---
    cd.has_pulse       = false;
    cd.cfd_ok          = false;
    cd.is_clipped      = false;
    cd.is_oscillating  = false;
    cd.baseline        = 0.0;
    cd.baseline_rms    = 0.0;
    cd.amplitude       = 0.0;
    cd.v_min           = 0.0;
    cd.v_max           = 0.0;
    cd.t_min           = 0.0;
    cd.t_cfd           = -999.0;
    cd.n_pos_crossings = 0;

    // Inizializzazione campi TOT
    for (int kth = 0; kth < NTOT_THR; kth++) {
        cd.tot_q   [kth] = 0.0;
        cd.tot_rise[kth] = 0.0;
        cd.tot_fall[kth] = 0.0;
        cd.tot_ok  [kth] = false;
    }
    // Inizializzazione campi di recupero TOT (hybrid_tot)
    cd.amplitude_tot     = 0.0;
    cd.t_cfd_rec_tot     = -999.0;
    cd.cfd_recovered_tot = false;
    // Inizializzazione campi PURE_TOT (diagnostica)
    cd.amplitude_pure_tot = 0.0;
    cd.t_cfd_pure_tot     = -999.0;
    cd.cfd_pure_tot_ok    = false;

// Inizializzazione rise time e crossing per la correzione Pietro&Rick
    cd.rise_time    = -999.0;
    cd.rise_time_ok = false;
    cd.t_10         = -999.0;
    cd.t_30         = -999.0;
    cd.t10_ok       = false;
    cd.t30_ok       = false;
    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return;  // Troppo pochi campioni: nulla da analizzare

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
    // V_min vicino al limite inferiore del range -> il vero picco e' "tagliato".
    if (vmin < CLIP_V_LO) cd.is_clipped = true;

    // ---- PASSO 4: AMPIEZZA E SOGLIA DI RUMORE ----
    double amp = bl - vmin;
    cd.amplitude = amp;
    if (amp < NOISE_THRESH) return;  // Sotto soglia: nessun impulso reale
    cd.has_pulse = true;

    // ---- PASSO 5: RILEVAMENTO OSCILLAZIONE (Schmitt trigger con isteresi) ----
    //   Stato "basso"  (above=false): si attende una salita sopra OSC_V_HIGH;
    //   trovata -> si conta UN crossing e si passa allo stato "alto".
    //   Stato "alto"   (above=true):  si attende una discesa sotto OSC_V_LOW;
    //   trovata -> si riarma lo Schmitt trigger SENZA contare nulla.
    //   La banda morta [OSC_V_LOW, OSC_V_HIGH] assorbe il ringing fine.
    //   Si analizza solo la coda DOPO il picco (i > imin): il fronte di salita
    //   di un impulso negativo va verso valori sempre piu' negativi e non
    //   genera crossing positivi.
    {
        int  n_cross = 0;
        bool above   = false;   // false = regione bassa, true = regione alta
        for (int i = imin + 1; i < ns; i++) {
            double vi = (double)v[i];
            if (!above) {
                if (vi > OSC_V_HIGH) { n_cross++; above = true; }
            } else {
                if (vi < OSC_V_LOW)  { above = false; }
            }
        }
        cd.n_pos_crossings = n_cross;
        if (n_cross >= OSC_NCROSS_MAX) cd.is_oscillating = true;
    }

    // ---- PASSO 6: CFD (Constant Fraction Discriminator) ----
    // Soglia CFD: una frazione CFD_FRACTION dell'ampiezza sotto la baseline.
    //   v_thr = baseline - CFD_FRACTION * amp
    // Ricerca IN AVANTI dalla fine della baseline al minimo: primo intervallo
    // [i, i+1] dove la tensione attraversa la soglia verso il basso, con
    // interpolazione lineare sub-campione.
    double v_thr = bl + CFD_FRACTION * (vmin - bl);   // = bl - CFD_FRACTION*amp
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
// ---- PASSO 6b: RISE TIME T90-T10 + CROSSING T10/T30 (correzione P&R) ----
    //
    // Calcoliamo nello stesso passaggio TRE istanti di crossing sul fronte di
    // salita (leading edge) del segnale negativo:
    //
    //   t_10  -> attraversamento del 10% dell'ampiezza (V_10 = bl - 0.10*amp)
    //   t_30  -> attraversamento del 30% dell'ampiezza (V_30 = bl - 0.30*amp)
    //   t_90  -> attraversamento del 90% dell'ampiezza (V_90 = bl - 0.90*amp)
    //
    // Il rise time canonico T90-T10 (cd.rise_time) misura la durata del
    // fronte di salita ed e' usato come diagnostica.
    //
    // Le grandezze t_10 e t_30 vengono SALVATE separatamente perche' sono
    // gli ingredienti elementari della correzione di Pietro&Rick:
    //
    //     dT_n = t_30,n - t_10,n         (rise time "effettivo" 10->30%)
    //     tof_bias = dT_3 - (dT_1+dT_2)/2
    //
    // Sottrarre dT_n al tempo CFD di ciascun canale equivale ad "ancorare"
    // l'istante di arrivo al crossing del 10% (zona di alto SNR, poco
    // sensibile a oscillazioni della baseline) per tutti i canali,
    // eliminando una grossa parte del time-walk residuo causato dal fatto
    // che il rise time dipende dalla posizione x sulla barra
    // (effetto di dispersione geometrica/cromatica della luce in BC408).
    //
    // ATTENZIONE AI CLIPPATI: se il segnale e' clippato l'ampiezza
    // misurata amp = bl - V_min e' SOTTOSTIMATA (il vero picco e' tagliato
    // dal plateau della DRS4). I livelli V_10, V_30, V_90 calcolati da
    // questa amp NON corrispondono alle frazioni della vera ampiezza, e i
    // crossing trovati su di essi sono biased. Sara' la logica di
    // ProcessDataset a richiedere is_clipped = false su tutti e 3 i canali
    // prima di applicare la correzione P&R.
    {
        double v_10 = bl - 0.10 * amp;   // 10% (vicino alla baseline)
        double v_30 = bl - 0.30 * amp;   // 30% (tipica soglia CFT)
        double v_90 = bl - 0.90 * amp;   // 90% (vicino al picco)

        // --- Crossing del livello 10% sul fronte di salita ---
        // Si esce al primo intervallo [i, i+1] dove il segnale, scendendo,
        // attraversa v_10. L'interpolazione lineare sub-campione e'
        // identica a quella usata per il CFD e per i TOT (coerenza
        // metodologica del timing).
        double t_10_loc = -1.0;
        bool   found_10 = false;
        for (int i = NBL_SAMPLES; i < imin; i++) {
            if (v[i] > v_10 && v[i + 1] <= v_10) {
                double dv = (double)v[i + 1] - v[i];
                if (fabs(dv) > 1e-6) {
                    t_10_loc = t[i] + (v_10 - v[i]) / dv * (t[i + 1] - t[i]);
                    found_10 = true;
                }
                break;
            }
        }

        // --- Crossing del livello 30% sul fronte di salita ---
        double t_30_loc = -1.0;
        bool   found_30 = false;
        for (int i = NBL_SAMPLES; i < imin; i++) {
            if (v[i] > v_30 && v[i + 1] <= v_30) {
                double dv = (double)v[i + 1] - v[i];
                if (fabs(dv) > 1e-6) {
                    t_30_loc = t[i] + (v_30 - v[i]) / dv * (t[i + 1] - t[i]);
                    found_30 = true;
                }
                break;
            }
        }

        // --- Crossing del livello 90% sul fronte di salita (per rise time) ---
        double t_90_loc = -1.0;
        bool   found_90 = false;
        for (int i = NBL_SAMPLES; i < imin; i++) {
            if (v[i] > v_90 && v[i + 1] <= v_90) {
                double dv = (double)v[i + 1] - v[i];
                if (fabs(dv) > 1e-6) {
                    t_90_loc = t[i] + (v_90 - v[i]) / dv * (t[i + 1] - t[i]);
                    found_90 = true;
                }
                break;
            }
        }

        // Salvataggio dei crossing 10% e 30% (con consistenza t_30 > t_10).
        if (found_10) {
            cd.t_10   = t_10_loc;
            cd.t10_ok = true;
        }
        if (found_30 && found_10 && t_30_loc > t_10_loc) {
            cd.t_30   = t_30_loc;
            cd.t30_ok = true;
        }

        // Rise time canonico T90-T10 (logica originale invariata)
        if (found_10 && found_90 && t_90_loc > t_10_loc) {
            cd.rise_time    = t_90_loc - t_10_loc;
            cd.rise_time_ok = true;
        }
    }
    // ---- PASSO 7: TIME-OVER-THRESHOLD MULTI-SOGLIA ----
    // Per ciascuna soglia q in TOT_THR_MV[]:
    //   - rise = primo crossing sul fronte di salita (tra NBL_SAMPLES e imin);
    //   - fall = primo crossing sul fronte di discesa (tra imin e ns-1);
    //   - TOT(q) = t_fall - t_rise (con interpolazione lineare sub-campione).
    // Variabile di lavoro: u(t) = baseline - v(t), positiva durante l'impulso.
    {
        const double q_max_safe = fabs(CLIP_V_LO) - TOT_CLIP_MARGIN;
        double sample_period = (ns > 1)
                             ? ((double)t[ns - 1] - (double)t[0]) / (ns - 1)
                             : 0.2;   // fallback 0.2 ns = 5 GS/s
        double dt_min_tot = TOT_NMIN_SAMPLES * sample_period;

        for (int kth = 0; kth < NTOT_THR; kth++) {
            double q = TOT_THR_MV[kth];
            if (q > q_max_safe) continue;   // soglia troppo vicina al clip
            if (amp < q)        continue;   // il segnale non raggiunge la soglia

            double v_thr_tot = bl - q;

            // --- crossing sul fronte di salita ---
            double t_rise = -1.0; bool rise_found = false;
            for (int i = NBL_SAMPLES; i < imin; i++) {
                if (v[i] > v_thr_tot && v[i + 1] <= v_thr_tot) {
                    double dv = (double)v[i + 1] - v[i];
                    if (fabs(dv) > 1e-6) {
                        t_rise = t[i] + (v_thr_tot - v[i]) / dv * (t[i + 1] - t[i]);
                        rise_found = true;
                    }
                    break;
                }
            }
            if (!rise_found) continue;

            // --- crossing sul fronte di discesa ---
            double t_fall = -1.0; bool fall_found = false;
            for (int i = imin; i < ns - 1; i++) {
                if (v[i] <= v_thr_tot && v[i + 1] > v_thr_tot) {
                    double dv = (double)v[i + 1] - v[i];
                    if (fabs(dv) > 1e-6) {
                        t_fall = t[i] + (v_thr_tot - v[i]) / dv * (t[i + 1] - t[i]);
                        fall_found = true;
                    }
                    break;
                }
            }
            if (!fall_found) continue;

            double tot_val = t_fall - t_rise;
            if (tot_val < dt_min_tot) continue;   // troppo stretto

            cd.tot_rise[kth] = t_rise;
            cd.tot_fall[kth] = t_fall;
            cd.tot_q   [kth] = tot_val;
            cd.tot_ok  [kth] = true;
        }
    }
}


// ==========================================================================
//  SEZIONE 4: MODELLO A(TOT) E RECUPERO DEI SEGNALI CLIPPATI
// ==========================================================================
//
//  In v13 esiste UN SOLO modello A(TOT): il polinomio cubico vincolato.
//  Tutte le funzioni di ricostruzione (recupero clippati hybrid_tot e modalita'
//  diagnostica PURE_TOT) chiamano la stessa EvalAmplitudeFromTOT(), garantendo
//  che calibrazione e ricostruzione siano coerenti per costruzione.

/// EvalAmplitudeFromTOT(): valuta A(TOT) in mV usando il polinomio cubico
/// calibrato per il canale ch_index, con lo schema di Horner.
///
/// INPUT:
///   tot      : durata sopra soglia (TOT) misurata sull'evento [ns]
///   ch_index : indice canale (0=PMT1, 1=PMT2, 2=PMT3)
///
/// RITORNA: ampiezza ricostruita [mV], oppure -1.0 se la calibrazione non e'
///   disponibile per quel canale o se il TOT non e' fisico (tot <= 0).
///
/// SCHEMA DI HORNER:
///   A = a0 + TOT*(a1 + TOT*(a2 + TOT*a3))
///   Numericamente stabile e veloce. Il coefficiente gTOT_poly_PMT[ch][0]
///   coincide con q_thr entro precisione numerica se FIX_A0_TO_QTHR == true.
double EvalAmplitudeFromTOT(double tot, int ch_index) {
    if (ch_index < 0 || ch_index >= MAX_CHANNELS) return -1.0;
    if (!gTOT_calibrated[ch_index])               return -1.0;
    if (tot <= 0.0)                               return -1.0;

    const double *coeff = gTOT_poly_PMT[ch_index];
    double A = coeff[gTOT_poly_degree];
    for (int k = gTOT_poly_degree - 1; k >= 0; --k) {
        A = A * tot + coeff[k];
    }
    return A;
}

/// TOTCalibrationAvailable(): true se la calibrazione polinomiale A(TOT) e'
/// disponibile e utilizzabile per il canale richiesto.
///
/// Condizioni minime:
///   - 0 <= ch_index < MAX_CHANNELS
///   - gTOT_calibrated[ch_index] == true
///   - gTOT_poly_degree == TOT_POLY_DEGREE (grado atteso)
///   - gTOT_Amax_PMT[ch_index] > 0 (esiste un cap di estrapolazione valido)
bool TOTCalibrationAvailable(int ch_index) {
    if (ch_index < 0 || ch_index >= MAX_CHANNELS) return false;
    if (!gTOT_calibrated[ch_index])               return false;
    if (gTOT_poly_degree != TOT_POLY_DEGREE)      return false;
    if (gTOT_Amax_PMT[ch_index] <= 0.0)           return false;
    return true;
}

/// RecoverClippedCFD_TOT(): unica funzione di recupero dei segnali clippati.
///   Per un canale CLIPPATO stima l'ampiezza vera dal TOT alla soglia di
///   riferimento e ricalcola il tempo CFD sulla soglia f*A_rec.
///
/// ALGORITMO (guida sez. 11.3):
///   1. Se il canale non e' clippato, non fa nulla (si usa il CFD standard).
///   2. Se il canale e' oscillante, fallisce.
///   3. Se il TOT alla soglia di riferimento non e' valido, fallisce.
///   4. A_rec = EvalAmplitudeFromTOT(tot_ref, ch_index); deve essere > 0.
///   5. Cap di estrapolazione: A_rec <= TOT_AREC_MAX_FACTOR * Amax.
///   6. Soglia CFD: V_cfd = baseline - CFD_FRACTION * A_rec.
///   7. V_cfd deve cadere SOPRA il plateau saturato:
///        V_cfd > CLIP_V_LO + TOT_CLIP_MARGIN
///      altrimenti il crossing non esiste sulla forma d'onda -> fallisce.
///   8. Ricerca del crossing sul fronte di salita con interpolazione lineare.
///
/// EFFETTO: popola cd.amplitude_tot, cd.t_cfd_rec_tot, cd.cfd_recovered_tot.
void RecoverClippedCFD_TOT(ChannelData &cd, int ch_index) {

    // --- Precondizioni ---
    if (!cd.is_clipped)                   return;  // niente da recuperare
    if (cd.is_oscillating)                return;  // segnale patologico
    if (ch_index < 0 || ch_index >= MAX_CHANNELS) return;
    if (!TOTCalibrationAvailable(ch_index))       return;

    const int kref = TOT_CALIB_THR_IDX;
    if (kref < 0 || kref >= NTOT_THR) return;
    if (!cd.tot_ok[kref])             return;   // TOT alla soglia ref non misurato

    double tot_ref = cd.tot_q[kref];

    // --- Stima dell'ampiezza ricostruita (polinomio cubico) ---
    double A_rec = EvalAmplitudeFromTOT(tot_ref, ch_index);
    if (A_rec <= 0.0) return;
    cd.amplitude_tot = A_rec;

    // --- Cap di estrapolazione (anti-outlier) ---
    double A_cap = TOT_AREC_MAX_FACTOR * gTOT_Amax_PMT[ch_index];
    if (A_cap > 0.0 && A_rec > A_cap) return;

    // --- Soglia CFD ricostruita ---
    double v_thr = cd.baseline - CFD_FRACTION * A_rec;

    // VERIFICA CRITICA: la soglia CFD deve cadere SOPRA il livello di clipping
    // (con margine), altrimenti il crossing e' dentro il plateau saturo e il
    // tempo sarebbe sbagliato per costruzione. In quel caso NON si recupera.
    if (v_thr <= CLIP_V_LO + TOT_CLIP_MARGIN) return;

    // --- Ricerca del crossing sul fronte di salita ---
    int   ns = cd.nsamples;
    float *t = cd.time;
    float *v = cd.voltage;

    double vmin = v[0];
    int    imin = 0;
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
    }

    for (int i = NBL_SAMPLES; i < imin; i++) {
        if (v[i] > v_thr && v[i + 1] <= v_thr) {
            double dv = (double)v[i + 1] - v[i];
            if (fabs(dv) > 1e-6) {
                cd.t_cfd_rec_tot     = t[i] + (v_thr - v[i]) / dv * (t[i + 1] - t[i]);
                cd.cfd_recovered_tot = true;
            }
            break;
        }
    }
}

/// ComputePureTOTTime(): modalita' DIAGNOSTICA PURE_TOT.
///   Per OGNI canale (clippato o no), calcola un tempo CFD basato
///   ESCLUSIVAMENTE sulla ricostruzione TOT dell'ampiezza, usando lo STESSO
///   polinomio A(TOT) del metodo principale.
///
/// RAZIONALE: il metodo hybrid_tot tratta in modo diverso i canali clippati e
///   non clippati; questa asimmetria, antisimmetrica in x, puo' generare una
///   pendenza spuria in C(x). Applicando il TOT a TUTTI i canali in modo
///   simmetrico, qualunque bias di ricostruzione si cancella in larga parte in
///   C = t3 - (t1+t2)/2. Confrontando C_pure_tot(x) con C_noclip(x) si misura
///   quanto il TOT e' "neutro".
///
/// ALGORITMO (guida sez. 12.2):
///   1. has_pulse == true
///   2. is_oscillating == false
///   3. tot_ok[TOT_CALIB_THR_IDX] == true
///   4. A_pure = EvalAmplitudeFromTOT(tot_ref, ch_index)
///   5. cap di estrapolazione con TOT_AREC_MAX_FACTOR
///   6. soglia CFD V_cfd = baseline - CFD_FRACTION * A_pure
///   7. verifica che V_cfd sia nel tratto digitalizzato (sopra il plateau,
///      sempre vero per i non clippati, controllato per i clippati)
///   8. ricerca del crossing CFD sul fronte di salita
///
/// EFFETTO: popola cd.amplitude_pure_tot, cd.t_cfd_pure_tot, cd.cfd_pure_tot_ok.
void ComputePureTOTTime(ChannelData &cd, int ch_index) {

    // --- Precondizioni ---
    if (!cd.has_pulse)                    return;
    if (cd.is_oscillating)                return;
    if (ch_index < 0 || ch_index >= MAX_CHANNELS) return;
    if (!TOTCalibrationAvailable(ch_index))       return;

const int kref = TOT_CALIB_THR_IDX;
    if (kref < 0 || kref >= NTOT_THR) return;

    // --- Gestione del caso "TOT non disponibile" ---
    // Per i canali NON clippati senza TOT valido (tipicamente ampiezza piu'
    // piccola della soglia di riferimento, oppure TOT troppo corto perche'
    // il segnale appena la supera): facciamo FALLBACK su t_cfd standard.
    // Questa scelta garantisce che pure_tot abbia entries >= noclip:
    //  - canali non clippati con TOT valido  -> tempo ricostruito via TOT
    //  - canali non clippati senza TOT       -> tempo CFD standard (= noclip)
    //  - canali clippati con TOT valido      -> tempo ricostruito via TOT
    //  - canali clippati senza TOT           -> fail (non ricostruibili)
    // Sul confronto C_pure vs C_noclip questo fallback non altera le
    // posizioni in cui pure_tot opera "puramente": dove TOT funziona su tutti
    // i 3 canali la ricostruzione e' simmetrica come prima; dove non funziona
    // su qualche canale non clippato, pure_tot coincide col CFD su quel canale.
    if (!cd.tot_ok[kref]) {
        if (cd.is_clipped) return;       // clippato senza TOT: non ricostruibile
        if (!cd.cfd_ok)    return;       // non clippato ma CFD pure fallito
        cd.amplitude_pure_tot = cd.amplitude;
        cd.t_cfd_pure_tot     = cd.t_cfd;
        cd.cfd_pure_tot_ok    = true;
        return;
    }

    double tot_ref = cd.tot_q[kref];

    // --- Stima dell'ampiezza ricostruita (stesso polinomio del metodo principale) ---
    double A_pure = EvalAmplitudeFromTOT(tot_ref, ch_index);
    if (A_pure <= 0.0) return;
    cd.amplitude_pure_tot = A_pure;

    // --- Cap di estrapolazione (anti-outlier) ---
    double A_cap = TOT_AREC_MAX_FACTOR * gTOT_Amax_PMT[ch_index];
    if (A_cap > 0.0 && A_pure > A_cap) return;

    // --- Soglia CFD ricostruita ---
    double v_thr = cd.baseline - CFD_FRACTION * A_pure;

    // Per i CLIPPATI: la soglia CFD deve cadere SOPRA il plateau di clipping.
    // Per i NON clippati il check e' automaticamente soddisfatto (lo applichiamo
    // comunque per uniformita').
    if (cd.is_clipped) {
        if (v_thr <= CLIP_V_LO + TOT_CLIP_MARGIN) return;
    }

    // --- Ricerca del crossing sul fronte di salita ---
    int   ns = cd.nsamples;
    float *t = cd.time;
    float *v = cd.voltage;

    double vmin = v[0];
    int    imin = 0;
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
    }

    for (int i = NBL_SAMPLES; i < imin; i++) {
        if (v[i] > v_thr && v[i + 1] <= v_thr) {
            double dv = (double)v[i + 1] - v[i];
            if (fabs(dv) > 1e-6) {
                cd.t_cfd_pure_tot  = t[i] + (v_thr - v[i]) / dv * (t[i + 1] - t[i]);
                cd.cfd_pure_tot_ok = true;
            }
            break;
        }
    }
}


// ==========================================================================
//  SEZIONE 5: PARSING XML DRS4
// ==========================================================================

/// ParseXML(): legge un file XML della DRS4 e popola il vettore di eventi.
///
/// Il parser usa una macchina a stati finiti con 4 stati:
///   IDLE       -> attesa di un <Event>
///   IN_EVENT   -> lettura header (Serial, Time)
///   IN_BOARD   -> lettura dati scheda (Trigger_Cell, Scaler, canali)
///   IN_CHANNEL -> lettura campioni <Data>tempo,tensione</Data>
///
/// Per ogni evento completo chiama AnalyzeChannel() su tutti i canali presenti.
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
//  SEZIONE 5-bis: CACHE DI PARSING (XML -> ROOT) — vedi Punto 1
// ==========================================================================
//  Tre funzioni:
//    - GetFileMTime()     : tempo di modifica di un file (per validare la cache)
//    - WriteDatasetCache(): scrive le forme d'onda grezze in un file ROOT
//    - ReadDatasetCache() : rilegge le forme d'onda e RIESEGUE AnalyzeChannel
//    - LoadOrParseDataset(): orchestratore (cache se valida, altrimenti parsing)

/// GetFileMTime(): ritorna il tempo di ultima modifica del file [secondi epoch],
/// oppure -1 se il file non esiste / non e' accessibile.
static long GetFileMTime(const char* path) {
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    return (long)st.st_mtime;
}

/// WriteDatasetCache(): salva le forme d'onda grezze di tutti gli eventi in un
/// file ROOT di cache. NON salva le grandezze calcolate da AnalyzeChannel: solo
/// l'header e gli array time[]/voltage[] per canale.
///
/// Layout del TTree "events_cache" (un'entry per evento):
///   serial, trigger_cell, board_serial, nchannels        (Int_t)
///   channel_ids[MAX_CHANNELS], scaler[MAX_CHANNELS]       (Int_t)
///   nsamp[MAX_CHANNELS]                                   (Int_t)
///   wf_t[MAX_CHANNELS][MAX_SAMPLES]                        (Float_t) tempi [ns]
///   wf_v[MAX_CHANNELS][MAX_SAMPLES]                        (Float_t) tensioni [mV]
///
/// RITORNA: true se la scrittura e' riuscita.
/// NOTA: salva/ripristina gDirectory per non interferire con il file ROOT
///       principale (fout), che resta aperto durante tutta la calibrazione.
bool WriteDatasetCache(const char* cache_path,
                       const std::vector<EventData>& events) {

    TDirectory* save = gDirectory;   // memorizza la directory ROOT corrente (fout)

    TFile* fc = new TFile(cache_path, "RECREATE");
    if (!fc || !fc->IsOpen()) {
        if (save) save->cd();
        delete fc;
        return false;
    }
    fc->cd();

    // Marcatore di versione del formato (per invalidazione futura).
    TParameter<int> fmt("PARSE_CACHE_FORMAT", PARSE_CACHE_FORMAT);
    fmt.Write();

    TTree* t = new TTree("events_cache", "Cache forme d'onda DRS4 parsate");

    Int_t   serial, trigger_cell, board_serial, nchannels;
    Int_t   channel_ids[MAX_CHANNELS], scaler[MAX_CHANNELS], nsamp[MAX_CHANNELS];
    // Array statici: ~32 KB totali, evitano di allocare/liberare a ogni run.
    static Float_t wf_t[MAX_CHANNELS][MAX_SAMPLES];
    static Float_t wf_v[MAX_CHANNELS][MAX_SAMPLES];

    t->Branch("serial",       &serial,       "serial/I");
    t->Branch("trigger_cell", &trigger_cell, "trigger_cell/I");
    t->Branch("board_serial", &board_serial, "board_serial/I");
    t->Branch("nchannels",    &nchannels,    "nchannels/I");
    t->Branch("channel_ids",  channel_ids,   Form("channel_ids[%d]/I", MAX_CHANNELS));
    t->Branch("scaler",       scaler,        Form("scaler[%d]/I",      MAX_CHANNELS));
    t->Branch("nsamp",        nsamp,         Form("nsamp[%d]/I",       MAX_CHANNELS));
    t->Branch("wf_t",         wf_t,          Form("wf_t[%d][%d]/F", MAX_CHANNELS, MAX_SAMPLES));
    t->Branch("wf_v",         wf_v,          Form("wf_v[%d][%d]/F", MAX_CHANNELS, MAX_SAMPLES));

    for (const auto& e : events) {
        serial       = e.serial;
        trigger_cell = e.trigger_cell;
        board_serial = e.board_serial;
        nchannels    = e.nchannels;
        for (int k = 0; k < MAX_CHANNELS; k++) {
            // Per i canali non presenti si scrivono zeri (verranno ignorati).
            channel_ids[k] = (k < e.nchannels) ? e.channel_ids[k]  : 0;
            scaler[k]      = (k < e.nchannels) ? e.scaler[k]        : 0;
            int ns         = (k < e.nchannels) ? e.ch[k].nsamples   : 0;
            if (ns > MAX_SAMPLES) ns = MAX_SAMPLES;
            nsamp[k] = ns;
            for (int i = 0; i < MAX_SAMPLES; i++) {
                wf_t[k][i] = (i < ns) ? e.ch[k].time[i]    : 0.0f;
                wf_v[k][i] = (i < ns) ? e.ch[k].voltage[i] : 0.0f;
            }
        }
        t->Fill();
    }

    fc->cd();
    t->Write();
    fc->Close();
    delete fc;

    if (save) save->cd();   // ripristina la directory ROOT principale (fout)
    return true;
}

/// ReadDatasetCache(): rilegge le forme d'onda da un file di cache e RICOSTRUISCE
/// il vettore di EventData. Per ogni evento RIESEGUE AnalyzeChannel sui canali,
/// in modo che le grandezze derivate siano sempre coerenti con il codice
/// corrente (la cache contiene solo l'input grezzo).
///
/// RITORNA: numero di eventi caricati, oppure -1 se la cache non e' valida
///   (file assente, formato incompatibile, TTree mancante).
int ReadDatasetCache(const char* cache_path, std::vector<EventData>& events) {

    TDirectory* save = gDirectory;

    TFile* fc = new TFile(cache_path, "READ");
    if (!fc || !fc->IsOpen()) {
        if (save) save->cd();
        delete fc;
        return -1;
    }

    // Verifica della versione del formato.
    TParameter<int>* fmt = (TParameter<int>*)fc->Get("PARSE_CACHE_FORMAT");
    if (!fmt || fmt->GetVal() != PARSE_CACHE_FORMAT) {
        std::cerr << "[CACHE] Formato incompatibile in " << cache_path
                  << " -> riparsing." << std::endl;
        fc->Close(); delete fc;
        if (save) save->cd();
        return -1;
    }

    TTree* t = (TTree*)fc->Get("events_cache");
    if (!t) {
        fc->Close(); delete fc;
        if (save) save->cd();
        return -1;
    }

    Int_t   serial, trigger_cell, board_serial, nchannels;
    Int_t   channel_ids[MAX_CHANNELS], scaler[MAX_CHANNELS], nsamp[MAX_CHANNELS];
    static Float_t wf_t[MAX_CHANNELS][MAX_SAMPLES];
    static Float_t wf_v[MAX_CHANNELS][MAX_SAMPLES];

    t->SetBranchAddress("serial",       &serial);
    t->SetBranchAddress("trigger_cell", &trigger_cell);
    t->SetBranchAddress("board_serial", &board_serial);
    t->SetBranchAddress("nchannels",    &nchannels);
    t->SetBranchAddress("channel_ids",  channel_ids);
    t->SetBranchAddress("scaler",       scaler);
    t->SetBranchAddress("nsamp",        nsamp);
    t->SetBranchAddress("wf_t",         wf_t);
    t->SetBranchAddress("wf_v",         wf_v);

    Long64_t nentries = t->GetEntries();
    events.clear();
    events.reserve((size_t)nentries);

    for (Long64_t ie = 0; ie < nentries; ie++) {
        t->GetEntry(ie);

        EventData e = EventData();   // value-init: tutti i campi a zero
        e.serial       = serial;
        e.trigger_cell = trigger_cell;
        e.board_serial = board_serial;
        e.nchannels    = nchannels;
        for (int k = 0; k < MAX_CHANNELS; k++) {
            e.channel_ids[k] = channel_ids[k];
            e.scaler[k]      = scaler[k];
            int ns = nsamp[k];
            if (ns > MAX_SAMPLES) ns = MAX_SAMPLES;
            e.ch[k].nsamples = ns;
            for (int i = 0; i < ns; i++) {
                e.ch[k].time[i]    = wf_t[k][i];
                e.ch[k].voltage[i] = wf_v[k][i];
            }
        }

        // Rieseguo l'analisi della forma d'onda: la cache NON la salva, cosi'
        // restano coerenti CFD/TOT/rise time con il codice corrente.
        for (int k = 0; k < e.nchannels; k++) AnalyzeChannel(e.ch[k]);

        events.push_back(e);
    }

    fc->Close();
    delete fc;
    if (save) save->cd();
    return (int)events.size();
}

/// LoadOrParseDataset(): punto d'ingresso unico per ottenere gli eventi di un
/// dataset. Sostituisce le chiamate dirette a ParseXML in CalibrateTOT e
/// ProcessDataset.
///
/// LOGICA:
///   1. La cartella della cache e' "<dir XML>/PARSE_CACHE_SUBDIR/<name>.root".
///   2. Se USE_PARSE_CACHE && !FORCE_REPARSE e la cache esiste ed e' piu'
///      recente (o coeva) dell'XML -> carica da cache.
///   3. Altrimenti -> ParseXML (lento) e, se la cache e' attiva, la (ri)scrive.
///
/// INPUT:
///   xml_path : path completo del file XML (gia' risolto dal chiamante)
///   name     : nome del dataset senza estensione (es. "N112"), per il file cache
///   events   : vettore di output
/// RITORNA: numero di eventi, oppure <= 0 in caso di errore di parsing.
int LoadOrParseDataset(const char* xml_path, const std::string& name,
                       std::vector<EventData>& events) {

    // Directory dei dati = directory del file XML.
    std::string xmls(xml_path);
    size_t slash = xmls.find_last_of("/\\");
    std::string dir = (slash == std::string::npos)
                    ? std::string(".")
                    : xmls.substr(0, slash);
    std::string cache_dir  = dir + "/" + PARSE_CACHE_SUBDIR;
    std::string cache_path = cache_dir + "/" + name + ".root";

    // ---- Tentativo di caricamento da cache ----
    if (USE_PARSE_CACHE && !FORCE_REPARSE) {
        long t_xml   = GetFileMTime(xml_path);
        long t_cache = GetFileMTime(cache_path.c_str());
        // Cache valida se esiste (t_cache >= 0) ed e' >= tempo XML.
        if (t_cache >= 0 && t_cache >= t_xml) {
            std::cout << "[CACHE] Carico " << name << " da "
                      << cache_path << " ..." << std::flush;
            int nc = ReadDatasetCache(cache_path.c_str(), events);
            if (nc > 0) {
                std::cout << " " << nc << " eventi (da cache)." << std::endl;
                return nc;
            }
            std::cout << " cache non valida, riparsing." << std::endl;
        }
    }

    // ---- Parsing dell'XML (lento) ----
    int n = ParseXML(xml_path, events);
    if (n <= 0) return n;

    // ---- Scrittura della cache per le run successive ----
    if (USE_PARSE_CACHE) {
        gSystem->mkdir(cache_dir.c_str(), kTRUE);   // crea la cartella (ricorsivo)
        if (WriteDatasetCache(cache_path.c_str(), events))
            std::cout << "[CACHE] Salvato " << cache_path << std::endl;
        else
            std::cerr << "[CACHE] Impossibile scrivere " << cache_path
                      << " (parsing comunque OK)." << std::endl;
    }
    return n;
}

// ==========================================================================
//  SEZIONE 6: CALIBRAZIONE DEL MODELLO A(TOT)
// ==========================================================================

/// CalibrateTOT(): calibrazione UNICA del modello A(TOT) per PMT1, PMT2, PMT3.
///
///   Determina i coefficienti del polinomio cubico
///       A(TOT) = a0 + a1*TOT + a2*TOT^2 + a3*TOT^3
///   con a0 = q_thr fissato (se FIX_A0_TO_QTHR == true), usando SOLO segnali
///   NON clippati dei dataset centrali (|x_k| <= TOT_CALIB_X_ABS_MAX).
///
/// PERCHE' PMT3: anche PMT3 viene calibrato, perche' entra in
///   C = t3 - (t1+t2)/2 e una sua ricostruzione diversa modifica il TOF.
///
/// METODO (guida sez. 9):
///   1. Parsing degli XML dei dataset centrali.
///   2. Per ogni canale (0,1,2) e ogni evento, selezione:
///        has_pulse && !is_clipped && !is_oscillating &&
///        tot_ok[TOT_CALIB_THR_IDX] && amplitude > q_thr
///      Le coppie (TOT, A) finiscono in un TGraph per canale.
///   3. Cleaning robusto (mediana + MAD su TOT, taglio anti-clip su A).
///   4. Profilo MEDIANO in bin equipopolati di TOT.
///   5. Fit polinomiale cubico del profilo (TGraphErrors -> chi2 pesato), con
///      a0 fissato a q_thr.
///   6. Check di monotonia (warning non bloccante).
///   7. Salvataggio in gTOT_poly_PMT[], gTOT_poly_err_PMT[],
///      gTOT_poly_chi2ndf_PMT[], gTOT_Amax_PMT[], gTOT_calibrated[].
///   8. Canvas "Calibration_TOT_to_amplitude" con 3 pad (PMT1/2/3).
///
/// EFFETTI: popola le globali gTOT_* e scrive un canvas nel file ROOT.
void CalibrateTOT(const char* folder,
                  const std::vector<PositionInfo> &positions,
                  TFile* fout) {

    const int    kref  = TOT_CALIB_THR_IDX;
    const double q_thr = TOT_THR_MV[kref];

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  CALIBRAZIONE A(TOT)  (polinomio cubico)     " << std::endl;
    std::cout << "  Soglia di riferimento: q_thr = " << q_thr << " mV" << std::endl;
    std::cout << "  Range posizioni: |x| <= " << TOT_CALIB_X_ABS_MAX << " cm" << std::endl;
    std::cout << "  Canali calibrati: PMT1, PMT2, PMT3           " << std::endl;
    std::cout << "=============================================" << std::endl;

    // Un TGraph (TOT, A) per ciascuno dei 3 PMT.
    TGraph *g_PMT[3];
    int     n_pts[3] = {0, 0, 0};
    double  A_max[3] = {0.0, 0.0, 0.0};
    const int   colors[3] = {kBlue + 1, kRed + 1, kGreen + 2};
    const char* pmt_tag[3] = {"PMT1", "PMT2", "PMT3"};
    for (int k = 0; k < 3; k++) {
        g_PMT[k] = new TGraph();
        g_PMT[k]->SetName(Form("g_calib_TOT_%s", pmt_tag[k]));
        g_PMT[k]->SetTitle(Form("Calibrazione %s: A vs TOT(q=%.0f mV);"
                                "TOT [ns];Ampiezza A [mV]", pmt_tag[k], q_thr));
    }

    int n_files_used = 0, n_events_total = 0;

    // ---- Loop sui dataset centrali: raccolta delle coppie (TOT, A) ----
    for (const auto &pos : positions) {

        if (fabs(pos.x_cm) > TOT_CALIB_X_ABS_MAX) continue;  // solo dataset centrali

        // Risoluzione del path al file XML (con o senza slash finale).
        std::string xml_path = std::string(folder) + "/" + pos.name + ".xml";
        std::ifstream test_file(xml_path.c_str());
        if (!test_file.good()) {
            xml_path = std::string(folder) + pos.name + ".xml";
            test_file.open(xml_path.c_str());
            if (!test_file.good()) {
                std::cerr << "[WARNING] CalibrateTOT: file non trovato per "
                          << pos.name << " — salto." << std::endl;
                continue;
            }
        }
        test_file.close();

        std::vector<EventData> events;
        // Cache di parsing: la prima run scrive i .root; le successive caricano.
        int n_parsed = LoadOrParseDataset(xml_path.c_str(), pos.name, events);
        if (n_parsed <= 0) continue;
        n_files_used++;
        n_events_total += (int)events.size();

        for (size_t ev = 0; ev < events.size(); ev++) {
            const EventData &e = events[ev];
            int nch = std::min(e.nchannels, 3);
            for (int k = 0; k < nch; k++) {
                const ChannelData &cd = e.ch[k];
                // Selezione: solo segnali NON clippati con TOT valido e
                // ampiezza sopra la soglia (per i clippati A non e' nota).
                if (!cd.has_pulse)        continue;
                if ( cd.is_clipped)       continue;
                if ( cd.is_oscillating)   continue;
                if (!cd.tot_ok[kref])     continue;
                if ( cd.amplitude <= q_thr) continue;

                g_PMT[k]->SetPoint(n_pts[k], cd.tot_q[kref], cd.amplitude);
                n_pts[k]++;
                if (cd.amplitude > A_max[k]) A_max[k] = cd.amplitude;
            }
        }
    }

    std::cout << "[INFO] CalibrateTOT: " << n_files_used << " dataset, "
              << n_events_total << " eventi totali" << std::endl;
    for (int k = 0; k < 3; k++) {
        std::cout << "       " << pmt_tag[k] << ": " << n_pts[k]
                  << " punti utili (A_max = " << Form("%.0f", A_max[k]) << " mV)"
                  << std::endl;
    }

    // ----------------------------------------------------------------------
    //  Lambda FitOneTOTPoly: fit polinomiale cubico sul PROFILO MEDIANO.
    //
    //  La nuvola grezza (TOT, A) ha enorme dispersione intrinseca: a parita' di
    //  TOT, A varia di ~50-100 mV per fluttuazioni di forma d'onda, statistica
    //  fotoelettronica del PMT, rumore di baseline. Il fit del profilo mediano
    //  estrae <A | TOT> con errori statistici corretti.
    //
    //  Passi: (1) cleaning robusto (MAD su TOT, taglio anti-clip su A),
    //         (2) profilo mediano in bin equipopolati di TOT,
    //         (3) fit polinomiale del profilo (TGraphErrors -> chi2 pesato),
    //             con a0 fissato a q_thr se FIX_A0_TO_QTHR.
    //  Il profilo e il fit vengono AGGIUNTI al TGraph originale come oggetti
    //  decorativi (visibili nel canvas finale).
    // ----------------------------------------------------------------------
    auto FitOneTOTPoly = [&](TGraph* gr, const char* name,
                             double coeffs_out[TOT_POLY_MAX_DEG + 1],
                             double cerrs_out [TOT_POLY_MAX_DEG + 1],
                             double &chi2ndf_out) -> bool {

        for (int k = 0; k <= TOT_POLY_MAX_DEG; k++) {
            coeffs_out[k] = 0.0;
            cerrs_out [k] = 0.0;
        }
        chi2ndf_out = -1.0;

        const int N = gr->GetN();
        if (N < 100) {
            std::cerr << "  [WARNING] " << name << ": troppi pochi punti ("
                      << N << ") per il fit polinomiale grado "
                      << TOT_POLY_DEGREE << "." << std::endl;
            return false;
        }

        // ---- Step 1: estrazione punti dal TGraph ----
        std::vector<double> tot_v(N), A_v(N);
        for (int i = 0; i < N; i++) {
            double xx, yy;
            gr->GetPoint(i, xx, yy);
            tot_v[i] = xx;
            A_v[i]   = yy;
        }

        // ---- Step 2: cleaning robusto su TOT (mediana + MAD) ----
        std::vector<double> tot_sorted = tot_v;
        std::sort(tot_sorted.begin(), tot_sorted.end());
        double tot_median = tot_sorted[N / 2];

        std::vector<double> abs_dev(N);
        for (int i = 0; i < N; i++) abs_dev[i] = fabs(tot_v[i] - tot_median);
        std::sort(abs_dev.begin(), abs_dev.end());
        double MAD = abs_dev[N / 2];

        const double N_MAD     = 5.0;                       // finestra 5 sigma_MAD
        double sigma_MAD       = MAD * 1.4826;
        double tot_window      = std::max(N_MAD * sigma_MAD, 3.0);
        double tot_lo_cut      = tot_median - tot_window;
        double tot_hi_cut      = tot_median + tot_window;

        // Tagli su A: escludo eventi vicini al clip (sommita' spianata, TOT
        // sovrastimato) e quelli appena sopra soglia (rumorosi).
        const double A_MAX_CLEAN   = fabs(CLIP_V_LO) - 150.0;  // ~350 mV
        const double A_MIN_CLEAN   = q_thr + 10.0;             // soglia + 10 mV
        const double TOT_MIN_CLEAN = 1.0;

        std::vector<double> tot_clean, A_clean;
        tot_clean.reserve(N);
        A_clean.reserve(N);
        for (int i = 0; i < N; i++) {
            if (tot_v[i] < tot_lo_cut || tot_v[i] > tot_hi_cut) continue;
            if (tot_v[i] < TOT_MIN_CLEAN)                       continue;
            if (A_v[i] < A_MIN_CLEAN || A_v[i] > A_MAX_CLEAN)   continue;
            tot_clean.push_back(tot_v[i]);
            A_clean.push_back(A_v[i]);
        }
        int n_clean = (int)tot_clean.size();

        std::cout << "  [" << name << "] cleaning: " << N << " -> " << n_clean
                  << " punti (TOT median = " << Form("%.1f", tot_median)
                  << " ns, MAD = " << Form("%.1f", MAD) << " ns)" << std::endl;

        if (n_clean < 100) {
            std::cerr << "  [WARNING] " << name
                      << ": troppi pochi punti dopo cleaning (" << n_clean
                      << ")." << std::endl;
            return false;
        }

        // ---- Step 3: profilo mediano in bin equipopolati di TOT ----
        std::vector<int> sort_idx(n_clean);
        for (int i = 0; i < n_clean; i++) sort_idx[i] = i;
        std::sort(sort_idx.begin(), sort_idx.end(),
                  [&tot_clean](int a, int b) { return tot_clean[a] < tot_clean[b]; });

        // Adatto il numero di bin se la statistica e' bassa.
        int nbins_eff = TOT_PROFILE_NBINS;
        while (nbins_eff > 5 && n_clean / nbins_eff < TOT_PROFILE_MIN_N) nbins_eff--;
        int npb_target = n_clean / nbins_eff;

        std::vector<double> prof_tot, prof_A, prof_errA;
        prof_tot.reserve(nbins_eff);
        prof_A.reserve(nbins_eff);
        prof_errA.reserve(nbins_eff);

        for (int b = 0; b < nbins_eff; b++) {
            int i_start = b * npb_target;
            int i_end   = (b == nbins_eff - 1) ? n_clean : (b + 1) * npb_target;
            int n_in_bin = i_end - i_start;
            if (n_in_bin < TOT_PROFILE_MIN_N) continue;

            std::vector<double> bin_tot(n_in_bin), bin_A(n_in_bin);
            for (int i = 0; i < n_in_bin; i++) {
                int idx_global = sort_idx[i_start + i];
                bin_tot[i] = tot_clean[idx_global];
                bin_A  [i] = A_clean  [idx_global];
            }
            std::sort(bin_tot.begin(), bin_tot.end());
            std::sort(bin_A.begin(),   bin_A.end());
            double tot_med_bin = bin_tot[n_in_bin / 2];
            double A_med_bin   = bin_A  [n_in_bin / 2];

            // Errore sulla mediana di A = MAD(A)*1.4826/sqrt(N_bin), con floor.
            std::vector<double> A_dev(n_in_bin);
            for (int i = 0; i < n_in_bin; i++) A_dev[i] = fabs(bin_A[i] - A_med_bin);
            std::sort(A_dev.begin(), A_dev.end());
            double A_MAD = A_dev[n_in_bin / 2];
            double A_err = A_MAD * 1.4826 / sqrt((double)n_in_bin);
            if (A_err < TOT_PROFILE_AERR_MIN) A_err = TOT_PROFILE_AERR_MIN;

            prof_tot.push_back (tot_med_bin);
            prof_A.push_back   (A_med_bin);
            prof_errA.push_back(A_err);
        }

        int n_prof = (int)prof_tot.size();
        std::cout << "  [" << name << "] profilo: " << n_prof
                  << " bin equipopolati (~" << npb_target << " eventi/bin)"
                  << std::endl;

        if (n_prof < TOT_POLY_DEGREE + 2) {
            std::cerr << "  [WARNING] " << name << ": profilo con troppi pochi bin ("
                      << n_prof << ")." << std::endl;
            return false;
        }

        // ---- Step 4: fit polinomiale cubico del profilo ----
        TGraphErrors *gr_prof = new TGraphErrors(n_prof,
                                                 prof_tot.data(),
                                                 prof_A.data(),
                                                 nullptr,            // errori su TOT trascurabili
                                                 prof_errA.data());
        gr_prof->SetName(Form("%s_profile", gr->GetName()));
        gr_prof->SetTitle(Form("%s_profile;TOT [ns];A_{mediano} [mV]", gr->GetName()));

        double tot_min_p, tot_max_p, A_min_p, A_max_p;
        gr_prof->ComputeRange(tot_min_p, A_min_p, tot_max_p, A_max_p);

        TString polname = Form("pol%d", TOT_POLY_DEGREE);
        TF1 *f_poly = new TF1(Form("%s_poly", name), polname.Data(), 0.0, tot_max_p);

        // Vincolo a0 = q_thr: A(TOT=0) = q_thr (asintoto fisico).
        if (FIX_A0_TO_QTHR) f_poly->FixParameter(0, q_thr);

        // Fit: "Q" quiet, "R" usa il range del TF1, "N" non aggiunge il TF1 al
        // graph (lo agganciamo manualmente a gr per la visualizzazione).
        // chi2 pesato automatico perche' gr_prof e' un TGraphErrors.
        gr_prof->Fit(f_poly, "Q R N");

        for (int k = 0; k <= TOT_POLY_DEGREE; k++) {
            coeffs_out[k] = f_poly->GetParameter(k);
            cerrs_out [k] = f_poly->GetParError(k);
        }
        chi2ndf_out = (f_poly->GetNDF() > 0)
                    ? f_poly->GetChisquare() / (double)f_poly->GetNDF() : -1.0;

        // ---- Step 5: check di monotonia (warning non bloccante) ----
        // Per i clippati il polinomio estrapola fino a ~1.5*TOT_max: una
        // non-monotonia in quel range distorcerebbe la ricostruzione di A.
        auto evalDeriv = [&](double T) -> double {
            double d = 0.0;
            for (int k = TOT_POLY_DEGREE; k >= 1; k--) d = d * T + k * coeffs_out[k];
            return d;
        };
        const double T_check_lo = tot_min_p;
        const double T_check_hi = 1.5 * tot_max_p;
        const int    NCHECK     = 100;
        int    n_neg = 0;
        double T_first_neg = -1.0;
        for (int i = 0; i < NCHECK; i++) {
            double T = T_check_lo + (T_check_hi - T_check_lo) * i / (double)(NCHECK - 1);
            if (evalDeriv(T) <= 0.0) {
                n_neg++;
                if (T_first_neg < 0.0) T_first_neg = T;
            }
        }
        if (n_neg > 0) {
            std::cerr << "  [WARNING] " << name << ": polinomio NON monotono in "
                      << n_neg << "/" << NCHECK << " punti del range ["
                      << Form("%.1f", T_check_lo) << ", " << Form("%.1f", T_check_hi)
                      << "] ns (primo a TOT = " << Form("%.2f", T_first_neg)
                      << " ns)." << std::endl;
        }

        // ---- Step 6: visualizzazione (decora il TGraph originale) ----
        TF1 *f_show = new TF1(Form("%s_show", name), polname.Data(), 0, 1.2 * tot_max_p);
        for (int k = 0; k <= TOT_POLY_DEGREE; k++) f_show->SetParameter(k, coeffs_out[k]);
        f_show->SetLineColor(kMagenta + 1);
        f_show->SetLineWidth(3);
        gr->GetListOfFunctions()->Add(f_show);

        gr_prof->SetMarkerStyle(20);
        gr_prof->SetMarkerSize(1.0);
        gr_prof->SetMarkerColor(kBlack);
        gr_prof->SetLineColor(kBlack);
        gr_prof->SetLineWidth(2);
        gr->GetListOfFunctions()->Add(gr_prof);  // di proprieta' della TList di gr

        delete f_poly;
        return true;
    };

    // ---- Esecuzione del fit per i 3 PMT ----
    double coeffs[3][TOT_POLY_MAX_DEG + 1] = {{0.0}};
    double cerrs [3][TOT_POLY_MAX_DEG + 1] = {{0.0}};
    double chi2  [3] = {-1.0, -1.0, -1.0};
    bool   ok    [3] = {false, false, false};

    for (int k = 0; k < 3; k++) {
        ok[k] = FitOneTOTPoly(g_PMT[k], Form("fit_TOT_%s", pmt_tag[k]),
                              coeffs[k], cerrs[k], chi2[k]);
    }

    // ---- Salvataggio nelle globali ----
    gTOT_poly_degree = TOT_POLY_DEGREE;
    for (int k = 0; k < 3; k++) {
        for (int j = 0; j <= TOT_POLY_MAX_DEG; j++) {
            gTOT_poly_PMT[k][j]     = ok[k] ? coeffs[k][j] : 0.0;
            gTOT_poly_err_PMT[k][j] = ok[k] ? cerrs [k][j] : 0.0;
        }
        gTOT_poly_chi2ndf_PMT[k] = chi2[k];
        gTOT_Amax_PMT[k]         = A_max[k];
        gTOT_calibrated[k]       = ok[k];
    }

    // ---- Stampa diagnostica ----
    std::cout << "\n  RISULTATI CALIBRAZIONE A(TOT) (polinomio grado "
              << TOT_POLY_DEGREE << "):" << std::endl;
    std::cout << "  ----------------------------------------" << std::endl;
    for (int k = 0; k < 3; k++) {
        std::cout << "  " << pmt_tag[k] << ": ";
        if (!ok[k]) { std::cout << "FIT FALLITO" << std::endl; continue; }
        std::cout << "chi2/ndf = " << Form("%.2f", chi2[k]) << std::endl;
        for (int j = 0; j <= TOT_POLY_DEGREE; j++) {
            std::cout << "        a" << j << " = "
                      << Form("%+.4e +/- %.4e", coeffs[k][j], cerrs[k][j])
                      << (j > 0 ? Form("  [mV/ns^%d]", j) : "  [mV]")
                      << std::endl;
        }
    }

    // ---- Canvas di output: 3 pad (PMT1/2/3) ----
    TCanvas *c_calib_tot = new TCanvas("c_calib_tot",
                                       "Calibrazione A(TOT)", 1500, 500);
    c_calib_tot->Divide(3, 1);

    for (int k = 0; k < 3; k++) {
        c_calib_tot->cd(k + 1);
        gPad->SetGrid(1, 1);
        g_PMT[k]->SetMarkerStyle(20);
        g_PMT[k]->SetMarkerSize(0.4);
        g_PMT[k]->SetMarkerColor(colors[k]);
        g_PMT[k]->Draw("AP");   // disegna nuvola + profilo + fit (oggetti agganciati)

        // Linea orizzontale A = q_thr (asintoto fisico richiesto dal vincolo a0).
        double xr_lo, yr_lo, xr_hi, yr_hi;
        g_PMT[k]->ComputeRange(xr_lo, yr_lo, xr_hi, yr_hi);
        TLine *l_q = new TLine(xr_lo, q_thr, xr_hi, q_thr);
        l_q->SetLineColor(kGray + 2);
        l_q->SetLineStyle(2);
        l_q->Draw("same");

        TPaveText *pt = new TPaveText(0.13, 0.55, 0.62, 0.89, "NDC");
        pt->SetFillColorAlpha(kWhite, 0.85);
        pt->SetTextFont(42);
        pt->SetTextSize(0.035);
        pt->SetTextAlign(12);
        pt->AddText(Form("%s: A = #sum_{k=0}^{%d} a_{k} #upoint TOT^{k}",
                         pmt_tag[k], TOT_POLY_DEGREE));
        if (ok[k]) {
            for (int j = 0; j <= TOT_POLY_DEGREE; j++) {
                pt->AddText(Form("a_{%d} = %+.3e #pm %.1e", j, coeffs[k][j], cerrs[k][j]));
            }
            pt->AddText(Form("#chi^{2}/ndf = %.2f", chi2[k]));
            pt->AddText(Form("a_{0} fissato a q_{thr} = %.0f mV", q_thr));
        } else {
            pt->AddText("Fit fallito");
        }
        pt->AddText(Form("N punti = %d   A_{max} = %.0f mV", n_pts[k], A_max[k]));
        pt->Draw();
    }

    c_calib_tot->Update();
    if (fout) {
        fout->cd();
        c_calib_tot->Write("Calibration_TOT_to_amplitude");
    }

    std::cout << "  [DONE] Calibrazione A(TOT) completata." << std::endl;
    std::cout << "=============================================\n" << std::endl;
}


// ==========================================================================
//  SEZIONE 7: FIT DELLE DISTRIBUZIONI TEMPORALI
// ==========================================================================

/// EstimateHistFWHM(): stima la FWHM direttamente dall'istogramma.
/// Metodo robusto contro outlier: parte dal bin con il massimo e scende a
/// destra e a sinistra fino ai punti a meta' altezza, interpolando linearmente.
double EstimateHistFWHM(TH1D* h) {

    int    max_bin  = h->GetMaximumBin();
    double half_max = h->GetBinContent(max_bin) / 2.0;
    int    nbins    = h->GetNbinsX();

    // Bordo sinistro
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

    // Bordo destro
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

/// FitLorentzian(): fitta un istogramma con una Lorentziana (Cauchy) usando
/// una parametrizzazione custom con parametri fisici intuitivi e
/// log-likelihood Poissoniana sull'intero range.
///
///                       (Gamma/2)^2
///   f(x) = A * ---------------------------
///              (x - x0)^2 + (Gamma/2)^2
///
///   [0] = A      altezza del picco (f(x0) = A)
///   [1] = x0     posizione del centro [ns]
///   [2] = Gamma  FWHM [ns]
///
/// Utile per le code sparse: la Lorentziana ha code ~ 1/x^2 e la
/// log-likelihood Poissoniana gestisce correttamente i bin poco popolati.
/// In v13 NON e' usata nel flusso principale (vedi FitGaussianCore), ma e'
/// mantenuta come strumento di confronto/diagnostica.
FitResult FitLorentzian(TH1D* h, const char* fit_name = "bw_fit") {

    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h->GetEntries();

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

    int    imax      = h->GetMaximumBin();
    double peak_hgt  = h->GetBinContent(imax);
    double peak_pos  = h->GetBinCenter(imax);
    double hist_fwhm = EstimateHistFWHM(h);

    double xlo = h->GetXaxis()->GetXmin();
    double xhi = h->GetXaxis()->GetXmax();
    if (hist_fwhm <= 0 || hist_fwhm > (xhi - xlo)) hist_fwhm = 0.05 * (xhi - xlo);

    TF1 *fL = new TF1(fit_name,
        "[0]*([2]/2.)*([2]/2.)/((x-[1])*(x-[1]) + ([2]/2.)*([2]/2.))",
        xlo, xhi);
    fL->SetParName(0, "A_peak");
    fL->SetParName(1, "x0");
    fL->SetParName(2, "Gamma");
    fL->SetParameters(peak_hgt, peak_pos, hist_fwhm);
    fL->SetParLimits(0, 0.1, 100.0 * peak_hgt);
    fL->SetParLimits(1, xlo, xhi);
    fL->SetParLimits(2, 1e-3 * (xhi - xlo), xhi - xlo);

    // "L" log-likelihood Poissoniana, "S" salva TFitResultPtr, "R" range, "Q" quiet.
    int fit_status = h->Fit(fL, "L S R Q");

    res.center     = fL->GetParameter(1);
    res.center_err = fL->GetParError(1);
    res.width      = fL->GetParameter(2);
    res.width_err  = fL->GetParError(2);
    res.chi2_ndf   = (fL->GetNDF() > 0) ? fL->GetChisquare() / fL->GetNDF() : -1.0;
    res.fit_ok     = (fit_status == 0);
    return res;
}

/// FitGaussianCore(): fitta un istogramma con una Gaussiana ristretta al
/// "core" della distribuzione, definito come l'intervallo che contiene il 90%
/// (default) degli eventi attorno alla mediana.
///
/// MOTIVAZIONE: dopo i tagli di qualita' a livello di forma d'onda, il core
/// delle differenze temporali e' dominato dal jitter Gaussiano del timing
/// (teorema del limite centrale). Le code residue contengono time-walk
/// residuo, ricostruzioni imperfette, fluttuazioni di calibrazione DRS4, e
/// gonfiano la deviazione standard apparente. Restringendo il fit al 90%
/// centrale si ottiene una sigma "pulita" = vera risoluzione temporale.
///
///   [0] = A   altezza del picco
///   [1] = mu  centro [ns]   <- valore restituito come "center"
///   [2] = sigma            <- risoluzione temporale, restituita come "width"
///   (FWHM = 2.355 * sigma)
FitResult FitGaussianCore(TH1D* h,
                          double core_fraction = 0.90,
                          const char* fit_name = "gaus_fit") {

    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h->GetEntries();

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

    // Quantili che definiscono il core (per 0.90: 5% e 95%).
    double q_lo = 0.5 - core_fraction / 2.0;
    double q_hi = 0.5 + core_fraction / 2.0;
    const int nq = 2;
    double xq[nq] = { q_lo, q_hi };
    double yq[nq] = { 0.0, 0.0 };
    h->GetQuantiles(nq, yq, xq);
    double x_lo = yq[0];
    double x_hi = yq[1];

    if (x_hi - x_lo < 5.0 * h->GetBinWidth(1)) {
        std::cerr << "[WARNING] Range core troppo piccolo per "
                  << h->GetName() << ", uso range completo." << std::endl;
        x_lo = h->GetXaxis()->GetXmin();
        x_hi = h->GetXaxis()->GetXmax();
    }

    int    imax       = h->GetMaximumBin();
    double peak_hgt   = h->GetBinContent(imax);
    double peak_pos   = h->GetBinCenter(imax);
    double hist_fwhm  = EstimateHistFWHM(h);
    double sigma_init = hist_fwhm / 2.355;
    if (sigma_init <= 0 || sigma_init > 0.5 * (x_hi - x_lo))
        sigma_init = 0.25 * (x_hi - x_lo);

    TF1 *fG = new TF1(fit_name,
        "[0]*TMath::Exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",
        x_lo, x_hi);
    fG->SetParName(0, "A_peak");
    fG->SetParName(1, "mu");
    fG->SetParName(2, "sigma");
    fG->SetParameters(peak_hgt, peak_pos, sigma_init);
    fG->SetParLimits(0, 0.1, 100.0 * peak_hgt);
    fG->SetParLimits(1, x_lo, x_hi);
    fG->SetParLimits(2, h->GetBinWidth(1), x_hi - x_lo);

// "L" log-likelihood Poissoniana, "S" salva TFitResultPtr, "R" range, "Q" quiet.
    int fit_status = h->Fit(fG, "L S R Q");

    // --- FIT ITERATIVO: refit ristretto a [mu - 2 sigma, mu + 2 sigma] ---
    // Su distribuzioni con code asimmetriche, il fit Poisson-log-likelihood
    // sull'intero core (90%) puo' spostare il centro fittato verso la coda
    // piu' popolata, allontanandolo dal massimo dell'istogramma. Una seconda
    // iterazione, centrata sul mu del primo fit e con finestra simmetrica
    // +/-2 sigma, stabilizza la stima del picco rimuovendo gran parte
    // dell'effetto delle code. Se la finestra ristretta risulta troppo
    // piccola (meno di 5 bin), si mantiene il primo fit.
    if (fit_status == 0 && fG->GetParameter(2) > 0) {
        double mu_1   = fG->GetParameter(1);
        double sig_1  = fabs(fG->GetParameter(2));
        double fit2_lo = std::max(mu_1 - 2.0 * sig_1, x_lo);
        double fit2_hi = std::min(mu_1 + 2.0 * sig_1, x_hi);
        if (fit2_hi - fit2_lo > 5.0 * h->GetBinWidth(1)) {
            fG->SetRange(fit2_lo, fit2_hi);
            fG->SetParLimits(1, fit2_lo, fit2_hi);
            // Re-inizializza l'ampiezza al massimo locale (puo' essere
            // diversa dal massimo globale se le code "alzano" il bin max).
            int    imax_loc = h->GetMaximumBin();   // resta il bin del massimo
            double hmax_loc = h->GetBinContent(imax_loc);
            fG->SetParameter(0, hmax_loc);
            fG->SetParameter(1, h->GetBinCenter(imax_loc));
            fG->SetParameter(2, sig_1);
            fit_status = h->Fit(fG, "L S R Q");
        }
    }

res.center     = fG->GetParameter(1);
    res.center_err = fG->GetParError(1);
    res.width      = fG->GetParameter(2);
    res.width_err  = fG->GetParError(2);
    res.chi2_ndf   = (fG->GetNDF() > 0) ? fG->GetChisquare() / fG->GetNDF() : -1.0;
    res.fit_ok     = (fit_status == 0);
    return res;
}

// ==========================================================================
//  SEZIONE 7-bis: FIT DEL RISE TIME CON CODA DESTRA (EMG) — vedi Punto 2
// ==========================================================================
//
//  Le distribuzioni di rise time T90-T10 hanno una CODA A DESTRA (fotoni in
//  ritardo per cammini ottici lunghi e riflessioni nella barra). Una gaussiana
//  simmetrica e' un modello sbagliato: per seguire la coda sposta il centro e
//  gonfia la sigma. Usiamo invece la Exponentially-Modified Gaussian (EMG):
//
//      f(x) = A * exp( 0.5*(sigma/tau)^2 - (x-mu)/tau )
//               * erfc( (sigma/tau - (x-mu)/sigma) / sqrt(2) )
//
//  Parametri (par[]):
//    par[0] = A     scala (NON l'altezza del picco)
//    par[1] = mu    centro della componente gaussiana [ns]
//    par[2] = sigma jitter gaussiano del timing [ns]
//    par[3] = tau   costante della coda esponenziale destra [ns]
//
//  Relazioni utili: media = mu + tau ; la moda (picco) e' in (mu, mu+tau).

/// RiseTimeEMG(): funzione EMG per TF1, con ramo asintotico numericamente
/// stabile. Sul lato sinistro (x << mu) il prodotto exp(...)*erfc(...) darebbe
/// overflow*underflow; lo si sostituisce con l'identita' esatta
///     exp(arg_exp)*erfc(z) = exp(-u^2/2) / (z*sqrt(pi))   (per z grande)
/// dove u = (x-mu)/sigma. La' l'EMG coincide con la gaussiana pura.
double RiseTimeEMG(double* xx, double* par) {
    double x   = xx[0];
    double A   = par[0];
    double mu  = par[1];
    double sig = par[2];
    double tau = par[3];
    if (sig <= 0.0 || tau <= 0.0) return 0.0;

    double r = sig / tau;
    double u = (x - mu) / sig;
    double z = (r - u) / TMath::Sqrt2();    // argomento di erfc

    if (z > 6.0) {
        // Ramo asintotico stabile (lontano dalla coda -> gaussiana pura).
        double g = TMath::Exp(-0.5 * u * u);
        return A * g / (z * TMath::Sqrt(TMath::Pi()));
    }
    double arg_exp = 0.5 * r * r - (x - mu) / tau;
    return A * TMath::Exp(arg_exp) * TMath::Erfc(z);
}

/// FitRiseTimeEMG(): fitta un istogramma di rise time con l'EMG e ritorna un
/// FitResult dove:
///   center = MODA della curva (= picco), via TF1::GetMaximumX. E' il valore
///            "tipico" del rise time, non distorto dalla coda destra.
///   width  = sigma gaussiano (jitter di timing).
///   center_err = errore su mu (proxy: a primo ordine d(moda)/d(mu) ~ 1).
///
/// PARAMETRI:
///   fit_name   : nome univoco della TF1.
///   line_color : colore della curva ( >=0 per disegnarla colorata, -1 = default).
///   attach     : se true la TF1 resta agganciata all'istogramma e viene
///                disegnata; se false si usa l'opzione "N" (solo numeri).
///
/// FALLBACK: con < 50 eventi, o se l'EMG non converge, ripiega su
///   FitGaussianCore (con nome distinto per non sovrascrivere la TF1 EMG).
FitResult FitRiseTimeEMG(TH1D* h,
                         const char* fit_name = "emg_fit",
                         int line_color = -1,
                         bool attach = true) {

    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h->GetEntries();

    // Con pochi eventi l'EMG a 4 parametri e' instabile: fallback gaussiano.
    if (res.nentries < 50) {
        return FitGaussianCore(h, 0.90, fit_name);
    }

    // ---- Stima iniziale dei parametri dalla statistica dell'istogramma ----
    int    imax     = h->GetMaximumBin();
    double peak_hgt = h->GetBinContent(imax);
    double peak_pos = h->GetBinCenter(imax);
    double rms      = h->GetStdDev();
    double xlo      = h->GetXaxis()->GetXmin();
    double xhi      = h->GetXaxis()->GetXmax();
    double bw       = h->GetBinWidth(1);

    double sig0 = std::max(0.4 * rms, 2.0 * bw);   // gaussiana stretta
    double tau0 = std::max(0.5 * rms, 2.0 * bw);   // coda destra moderata
    double mu0  = peak_pos - 0.5 * tau0;           // mu e' a sinistra del picco

    TF1* fEMG = new TF1(fit_name, RiseTimeEMG, xlo, xhi, 4);
    fEMG->SetParNames("A_scale", "mu", "sigma", "tau");
    // A e' una scala, non l'altezza: la inizializzo perche' f(picco) ~ peak_hgt.
    fEMG->SetParameters(1.0, mu0, sig0, tau0);
    double shape_at_peak = fEMG->Eval(peak_pos);
    double A0 = (shape_at_peak > 1e-12) ? peak_hgt / shape_at_peak : peak_hgt;
    fEMG->SetParameter(0, A0);

    // Limiti fisici.
    fEMG->SetParLimits(0, 1e-3 * peak_hgt, 1e4 * peak_hgt);
    fEMG->SetParLimits(1, xlo, xhi);
    fEMG->SetParLimits(2, bw, 0.7 * (xhi - xlo));
    fEMG->SetParLimits(3, 0.05 * bw, xhi - xlo);

    if (line_color >= 0) {
        fEMG->SetLineColor(line_color);
        fEMG->SetLineWidth(2);
    }

    // "L" log-likelihood Poissoniana (gestisce i bin poco popolati della coda),
    // "S" salva, "R" range della TF1, "Q" silenzioso. "N" = niente disegno.
    TString opt = attach ? "L S R Q" : "L S R Q N";
    int fit_status = h->Fit(fEMG, opt.Data());

    bool sane = (fit_status == 0)
             && (fEMG->GetParameter(2) > 0.0)
             && (fEMG->GetParameter(3) > 0.0);

    if (!sane) {
        std::cerr << "[WARNING] EMG non convergente per " << h->GetName()
                  << " -> fallback FitGaussianCore." << std::endl;
        // Stacca e cancella la TF1 EMG fallita per non disegnare una curva
        // spuria, poi ripiega sulla gaussiana di core (nome distinto).
        if (attach && h->GetListOfFunctions())
            h->GetListOfFunctions()->Remove(fEMG);
        delete fEMG;
        return FitGaussianCore(h, 0.90, Form("%s_gfb", fit_name));
    }

    // ---- Valore centrale = MODA (picco) della curva EMG ----
    double mode = fEMG->GetMaximumX(xlo, xhi);

    res.center     = mode;
    res.center_err = fEMG->GetParError(1);          // proxy: errore su mu
    res.width      = fabs(fEMG->GetParameter(2));   // sigma gaussiano
    res.width_err  = fEMG->GetParError(2);
    res.chi2_ndf   = (fEMG->GetNDF() > 0)
                   ? fEMG->GetChisquare() / fEMG->GetNDF() : -1.0;
    res.fit_ok     = true;
    return res;
}


// ==========================================================================
//  SEZIONE 8: COMPONENTI DELL'ERRORE SULLA POSIZIONE
// ==========================================================================

/// GetParallaxSigmaFromFile(): restituisce sigma_parallax [cm] alla posizione
/// x_k leggendo il TGraph "g_sigma_x" dal file ROOT del Monte Carlo di
/// accettanza. Il file viene aperto una sola volta (variabili statiche).
/// Se il file non e' disponibile o non contiene il grafico -> PARALLAX_FALLBACK.
double GetParallaxSigmaFromFile(double x_k, const char* parallax_file) {

    static TFile  *s_file  = nullptr;
    static TGraph *s_graph = nullptr;
    static std::string s_loaded_path = "";
    static bool   s_warned = false;

    std::string requested_path = (parallax_file ? parallax_file : "");

    if (requested_path.empty()) {
        if (!s_warned) {
            std::cerr << "[WARNING] PARALLAX_FILE non specificato, uso valore "
                      << "costante " << PARALLAX_FALLBACK << " cm per sigma_parallax."
                      << std::endl;
            s_warned = true;
        }
        return PARALLAX_FALLBACK;
    }

    if (requested_path != s_loaded_path) {
        if (s_file) { s_file->Close(); delete s_file; s_file = nullptr; s_graph = nullptr; }

        s_file = TFile::Open(parallax_file, "READ");
        if (!s_file || s_file->IsZombie()) {
            std::cerr << "[WARNING] Impossibile aprire file parallasse '"
                      << requested_path << "', uso fallback "
                      << PARALLAX_FALLBACK << " cm." << std::endl;
            if (s_file) { delete s_file; s_file = nullptr; }
            s_loaded_path = requested_path;
            return PARALLAX_FALLBACK;
        }

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
        std::cout << "[INFO] Caricato file di parallasse: " << requested_path
                  << " (" << s_graph->GetN() << " punti)" << std::endl;
    }

    if (s_graph) return s_graph->Eval(x_k);   // interpolazione lineare
    return PARALLAX_FALLBACK;
}

/// GetTapeSigma(): restituisce sigma_tape [cm] alla posizione x_k secondo il
/// modello selezionato da TAPE_MODEL_INCREMENTAL.
double GetTapeSigma(double x_k) {

    double dx_from_ref = fabs(x_k - TAPE_REFERENCE_X);
    if (dx_from_ref < 1e-6) return 0.0;   // posizione di riferimento: errore nullo

    if (!TAPE_MODEL_INCREMENTAL) {
        // Modello "lettura singola": ogni posizione != 0 ha errore costante.
        return SIGMA_TAPE_READ;
    } else {
        // Modello "letture incrementali": sigma_read^2 accumula in quadratura
        // per ogni passo di TAPE_STEP_LENGTH cm percorso dal riferimento.
        double n_step_d = dx_from_ref / TAPE_STEP_LENGTH;
        int    n_step   = (int)ceil(n_step_d - 1e-6);
        if (n_step < 1) n_step = 1;
        return SIGMA_TAPE_READ * sqrt((double)n_step);
    }
}

/// ComputeSigmaX(): errore totale sulla posizione, somma in quadratura di
/// parallasse e lettura del metro:  sigma_x^2 = sigma_parallax^2 + sigma_tape^2.
double ComputeSigmaX(double x_k, const char* parallax_file) {
    double sigma_par  = GetParallaxSigmaFromFile(x_k, parallax_file);
    double sigma_tape = GetTapeSigma(x_k);
    return sqrt(sigma_par * sigma_par + sigma_tape * sigma_tape);
}


// ==========================================================================
//  SEZIONE 9: ELABORAZIONE DI UN SINGOLO DATASET (una posizione x_k)
// ==========================================================================
//
//  ProcessDataset() legge il file XML di una posizione, e per ogni evento
//  calcola le differenze temporali con i TRE metodi (noclip, hybrid_tot,
//  pure_tot), riempie un TTree diagnostico e 12 istogrammi (Dt12/Dt13/Dt23/C
//  per ciascun metodo), che vengono poi fittati con FitGaussianCore.
//
//  MAPPATURA CANALI:
//    CH (idx 0) -> PMT1 (estremo sinistro, x = -L/2)
//    CH (idx 1) -> PMT2 (estremo destro, x = +L/2)
//    CH (idx 2) -> PMT3 (mobile, posizione x_k)
//    CH (idx 3) -> trigger NIM (non usato per il timing)
//
//  DEFINIZIONE DEI TRE METODI (condividono lo stesso polinomio A(TOT)):
//    noclip     : tempo = t_cfd; un canale e' valido solo se has_pulse &&
//                 cfd_ok && !is_clipped. I clippati NON vengono recuperati.
//                 E' il riferimento "pulito" privo di ricostruzione.
//    hybrid_tot : se il canale NON e' clippato -> tempo = t_cfd (valido se
//                 cfd_ok); se e' clippato -> tempo = t_cfd_rec_tot (valido se
//                 cfd_recovered_tot). E' il METODO FISICO ufficiale.
//    pure_tot   : tempo = t_cfd_pure_tot per OGNI canale (valido se
//                 cfd_pure_tot_ok). Ricostruzione TOT applicata simmetricamente
//                 a tutti i canali: diagnostica di simmetria.
//
//  TAGLI A LIVELLO DI EVENTO:
//    - servono almeno 3 canali e has_pulse su PMT1,PMT2,PMT3;
//    - se ENABLE_OSC_CUT e un canale tra i 3 e' oscillante, l'evento viene
//      scartato per TUTTI i metodi (un'oscillazione e' patologica a monte
//      della ricostruzione, indipendentemente dal metodo).
//  TAGLIO A LIVELLO DI METODO:
//    - se ENABLE_CFD_CUT, un metodo accetta l'evento solo se ha un tempo
//      valido su TUTTI e 3 i canali.

DatasetResults ProcessDataset(const char* xml_path,
                              const PositionInfo &pos,
                              TFile* outfile) {

    DatasetResults result;
    result.x    = pos.x_cm;
    result.name = pos.name;
    // dx viene riempito dal chiamante (TOF_Calibration) con ComputeSigmaX.
    result.dx   = DX_POSITION;
    result.noclip.n_good     = 0;
    result.hybrid_tot.n_good = 0;
    result.pure_tot.n_good   = 0;
    result.corrected.n_good  = 0;

    const std::string tname = pos.name;

// ---- Caricamento eventi: da cache .root se disponibile, altrimenti XML ----
    // I dataset centrali (|x| <= 56) sono gia' stati messi in cache da
    // CalibrateTOT entro questa stessa run, quindi qui vengono ricaricati e non
    // riparsati: il doppio parsing dei centrali sparisce.
    std::vector<EventData> events;
    int n_parsed = LoadOrParseDataset(xml_path, pos.name, events);
    if (n_parsed <= 0) {
        std::cerr << "[ERRORE] Nessun evento letto da " << xml_path << std::endl;
        return result;
    }

    // ---- Preparazione del TTree diagnostico ----
    outfile->cd();
    TTree* tree = new TTree(tname.c_str(),
                            Form("Dati TOF posizione %s", tname.c_str()));

    // Grandezze grezze di forma d'onda per i 3 PMT.
    Float_t t_cfd[3], amp[3], bl[3], tot_ref[3];
    Int_t   clipped[3], osc[3], n_cross[3];
    tree->Branch("t_cfd",   t_cfd,   "t_cfd[3]/F");    // tempo CFD standard [ns]
    tree->Branch("amp",     amp,     "amp[3]/F");      // ampiezza [mV]
    tree->Branch("bl",      bl,      "bl[3]/F");       // baseline [mV]
    tree->Branch("tot_ref", tot_ref, "tot_ref[3]/F"); // TOT alla soglia di rif. [ns] (0 se non valido)
    tree->Branch("clipped", clipped, "clipped[3]/I"); // 1 se canale clippato
    tree->Branch("osc",     osc,     "osc[3]/I");     // 1 se canale oscillante
    tree->Branch("n_cross", n_cross, "n_cross[3]/I"); // crossing Schmitt (diagnostica)
    // Rise time T90-T10 per i 3 PMT [ns].
    Float_t rt[3];
    tree->Branch("rise_time", rt, "rise_time[3]/F");

    // Ricostruzione TOT sui clippati (metodo hybrid_tot).
    Float_t amp_tot[3], t_rec_tot[3];
    Int_t   rec_tot_ok[3];
    tree->Branch("amp_tot",    amp_tot,    "amp_tot[3]/F");      // A ricostruita [mV]
    tree->Branch("t_rec_tot",  t_rec_tot,  "t_rec_tot[3]/F");    // tempo CFD ricostruito [ns]
    tree->Branch("rec_tot_ok", rec_tot_ok, "rec_tot_ok[3]/I");   // 1 se recupero clippato riuscito

    // Ricostruzione TOT su tutti i canali (metodo pure_tot).
    Float_t amp_pure[3], t_pure[3];
    Int_t   pure_ok[3];
    tree->Branch("amp_pure", amp_pure, "amp_pure[3]/F");   // A ricostruita PURE_TOT [mV]
    tree->Branch("t_pure",   t_pure,   "t_pure[3]/F");     // tempo CFD PURE_TOT [ns]
    tree->Branch("pure_ok",  pure_ok,  "pure_ok[3]/I");    // 1 se PURE_TOT valido per il canale

    // Differenze temporali e C per i 3 metodi + flag di evento buono.
    Float_t dt12_noclip, dt13_noclip, dt23_noclip, C_noclip;
    Float_t dt12_hybrid, dt13_hybrid, dt23_hybrid, C_hybrid;
    Float_t dt12_pure,   dt13_pure,   dt23_pure,   C_pure;
    Int_t   good_noclip, good_hybrid, good_pure;
    tree->Branch("dt12_noclip", &dt12_noclip, "dt12_noclip/F");
    tree->Branch("dt13_noclip", &dt13_noclip, "dt13_noclip/F");
    tree->Branch("dt23_noclip", &dt23_noclip, "dt23_noclip/F");
    tree->Branch("C_noclip",    &C_noclip,    "C_noclip/F");
    tree->Branch("good_noclip", &good_noclip, "good_noclip/I");
    tree->Branch("dt12_hybrid", &dt12_hybrid, "dt12_hybrid/F");
    tree->Branch("dt13_hybrid", &dt13_hybrid, "dt13_hybrid/F");
    tree->Branch("dt23_hybrid", &dt23_hybrid, "dt23_hybrid/F");
    tree->Branch("C_hybrid",    &C_hybrid,    "C_hybrid/F");
    tree->Branch("good_hybrid", &good_hybrid, "good_hybrid/I");
    tree->Branch("dt12_pure",   &dt12_pure,   "dt12_pure/F");
    tree->Branch("dt13_pure",   &dt13_pure,   "dt13_pure/F");
    tree->Branch("dt23_pure",   &dt23_pure,   "dt23_pure/F");
    tree->Branch("C_pure",      &C_pure,      "C_pure/F");
    tree->Branch("good_pure",   &good_pure,   "good_pure/I");
    // --- Crossing t_10, t_30 per PMT 1/2/3 (servono per la correzione P&R) ---
    // Salvati nel TTree per consentire ri-analisi offline della correzione
    // (es. con un coefficiente diverso, o un'altra coppia di soglie).
    Float_t t10_arr[3], t30_arr[3];
    Int_t   t10_ok_arr[3], t30_ok_arr[3];
    tree->Branch("t_10",     t10_arr,    "t_10[3]/F");      // [ns]
    tree->Branch("t_30",     t30_arr,    "t_30[3]/F");      // [ns]
    tree->Branch("t10_ok",   t10_ok_arr, "t10_ok[3]/I");
    tree->Branch("t30_ok",   t30_ok_arr, "t30_ok[3]/I");

    // --- Metodo "corrected": Dt_ij e C corretti per il bias P&R ---
    // L'evento e' good_corr=1 solo se tutti i 3 canali sono non clippati e
    // hanno t_10 e t_30 validi (oltre a essere good_hybrid).
    Float_t dt12_corr, dt13_corr, dt23_corr, C_corr;
    Float_t tof_bias_evt;   // bias P&R: dT_3 - (dT_1+dT_2)/2 [ns]
    Int_t   good_corr;
    tree->Branch("dt12_corr",   &dt12_corr,    "dt12_corr/F");
    tree->Branch("dt13_corr",   &dt13_corr,    "dt13_corr/F");
    tree->Branch("dt23_corr",   &dt23_corr,    "dt23_corr/F");
    tree->Branch("C_corr",      &C_corr,       "C_corr/F");
    tree->Branch("tof_bias",    &tof_bias_evt, "tof_bias/F");
    tree->Branch("good_corr",   &good_corr,    "good_corr/I");

// ---- Vettori per la raccolta dei valori (PASS 1) ----
    // Il range degli istogrammi Dt12/Dt13/Dt23/C e' calcolato DOPO il primo
    // loop sugli eventi a partire dalla mediana e dalla MAD dei dati raccolti.
    // Questo evita il range fisso [-30,30] che spreca bin: gli istogrammi
    // avranno ~150 bin concentrati attorno al picco, con una risoluzione
    // per bin ottimale per il fit gaussiano.
    // PASS 1: i valori vengono raccolti in vettori.
    // PASS 2 (dopo il loop): range adattivo -> creazione istogrammi -> fill.
    std::vector<float> v_dt12_nc, v_dt13_nc, v_dt23_nc, v_C_nc;
    std::vector<float> v_dt12_hy, v_dt13_hy, v_dt23_hy, v_C_hy;
    std::vector<float> v_dt12_pu, v_dt13_pu, v_dt23_pu, v_C_pu;
    // Vettori per i valori "corrected" (range adattivo nel PASS 2)
    std::vector<float> v_dt12_co, v_dt13_co, v_dt23_co, v_C_co;
    // Pre-alloco memoria: il dataset tipico ha ~1000 eventi buoni.
    const size_t est_size = events.size();
    v_dt12_nc.reserve(est_size); v_dt13_nc.reserve(est_size);
    v_dt23_nc.reserve(est_size); v_C_nc.reserve(est_size);
    v_dt12_hy.reserve(est_size); v_dt13_hy.reserve(est_size);
    v_dt23_hy.reserve(est_size); v_C_hy.reserve(est_size);
    v_dt12_pu.reserve(est_size); v_dt13_pu.reserve(est_size);
    v_dt23_pu.reserve(est_size); v_C_pu.reserve(est_size);
    v_dt12_co.reserve(est_size); v_dt13_co.reserve(est_size);
    v_dt23_co.reserve(est_size); v_C_co.reserve(est_size);
// ---- Vettori per la raccolta del rise time per PMT (solo non clippati) ----
    std::vector<float> v_rt_pmt[3];
    for (int k = 0; k < 3; k++) v_rt_pmt[k].reserve(est_size);

    // Contatori diagnostici delle ragioni di scarto.
    int n_too_few_ch = 0;   // meno di 3 canali, o has_pulse mancante
    int n_osc_event  = 0;   // evento scartato per oscillazione
    int n_clipped_ev = 0;   // eventi con almeno un canale clippato

    const int kref = TOT_CALIB_THR_IDX;

    // ====================================================================
    //  LOOP SUGLI EVENTI
    // ====================================================================
    for (size_t ev = 0; ev < events.size(); ev++) {

        EventData &e = events[ev];

        // Default "non valido" per tutti i branch (sovrascritti se l'evento e' usabile).
        for (int k = 0; k < 3; k++) {
            t_cfd[k] = -999.0f; amp[k] = 0.0f; bl[k] = 0.0f; tot_ref[k] = 0.0f;
            clipped[k] = 0; osc[k] = 0; n_cross[k] = 0;
            amp_tot[k] = 0.0f; t_rec_tot[k] = -999.0f; rec_tot_ok[k] = 0;
            amp_pure[k] = 0.0f; t_pure[k] = -999.0f; pure_ok[k] = 0;
        }
        dt12_noclip = dt13_noclip = dt23_noclip = C_noclip = -999.0f;
        dt12_hybrid = dt13_hybrid = dt23_hybrid = C_hybrid = -999.0f;
        dt12_pure   = dt13_pure   = dt23_pure   = C_pure   = -999.0f;
        good_noclip = good_hybrid = good_pure = 0;
        // Default "non valido" per i branch della correzione P&R.
        for (int k = 0; k < 3; k++) {
            t10_arr[k]    = -999.0f;
            t30_arr[k]    = -999.0f;
            t10_ok_arr[k] = 0;
            t30_ok_arr[k] = 0;
        }
        dt12_corr = dt13_corr = dt23_corr = C_corr = -999.0f;
        tof_bias_evt = -999.0f;
        good_corr = 0;

        // --- Taglio: servono almeno 3 canali ---
        if (e.nchannels < 3) {
            n_too_few_ch++;
            tree->Fill();
            continue;
        }

        // --- Esecuzione della ricostruzione TOT ---
        // hybrid_tot: recupero solo dei canali clippati (PMT1,2,3).
        // pure_tot:   ricostruzione TOT applicata a TUTTI i 3 canali.
        for (int k = 0; k < 3; k++) {
            if (e.ch[k].is_clipped) RecoverClippedCFD_TOT(e.ch[k], k);
            ComputePureTOTTime(e.ch[k], k);
        }

        // --- Riempimento dei branch grezzi + costruzione dei tempi per metodo ---
        bool any_osc      = false;
        bool any_clipped  = false;
        bool has_pulse_all = true;

        // t_use[metodo][canale] e validita'.
        Float_t t_nc[3], t_hy[3], t_pu[3];
        bool    ok_nc[3], ok_hy[3], ok_pu[3];

        for (int k = 0; k < 3; k++) {
            const ChannelData &cd = e.ch[k];

            t_cfd[k]   = (Float_t)cd.t_cfd;
            amp[k]     = (Float_t)cd.amplitude;
            bl[k]      = (Float_t)cd.baseline;
            tot_ref[k] = cd.tot_ok[kref] ? (Float_t)cd.tot_q[kref] : 0.0f;
            clipped[k] = cd.is_clipped     ? 1 : 0;
            osc[k]     = cd.is_oscillating ? 1 : 0;
            n_cross[k] = cd.n_pos_crossings;
            rt[k]      = cd.rise_time_ok ? (Float_t)cd.rise_time : -999.0f;
            amp_tot[k]    = (Float_t)cd.amplitude_tot;
            t_rec_tot[k]  = (Float_t)cd.t_cfd_rec_tot;
            rec_tot_ok[k] = cd.cfd_recovered_tot ? 1 : 0;

            amp_pure[k]   = (Float_t)cd.amplitude_pure_tot;
            t_pure[k]     = (Float_t)cd.t_cfd_pure_tot;
            pure_ok[k]    = cd.cfd_pure_tot_ok ? 1 : 0;
            // Popola gli array del TTree per i crossing 10% e 30%
            t10_arr[k]    = cd.t10_ok ? (Float_t)cd.t_10 : -999.0f;
            t30_arr[k]    = cd.t30_ok ? (Float_t)cd.t_30 : -999.0f;
            t10_ok_arr[k] = cd.t10_ok ? 1 : 0;
            t30_ok_arr[k] = cd.t30_ok ? 1 : 0;

            if (cd.is_oscillating) any_osc     = true;
            if (cd.is_clipped)     any_clipped = true;
            if (!cd.has_pulse)     has_pulse_all = false;

            // --- Metodo noclip: solo segnali non clippati con CFD valido ---
            if (cd.has_pulse && cd.cfd_ok && !cd.is_clipped) {
                t_nc[k]  = (Float_t)cd.t_cfd;
                ok_nc[k] = true;
            } else {
                t_nc[k]  = -999.0f;
                ok_nc[k] = false;
            }

            // --- Metodo hybrid_tot: CFD sui non clippati, recupero TOT sui clippati ---
            if (!cd.is_clipped) {
                t_hy[k]  = (Float_t)cd.t_cfd;
                ok_hy[k] = (cd.has_pulse && cd.cfd_ok);
            } else if (cd.cfd_recovered_tot) {
                t_hy[k]  = (Float_t)cd.t_cfd_rec_tot;
                ok_hy[k] = true;
            } else {
                t_hy[k]  = -999.0f;
                ok_hy[k] = false;
            }

            // --- Metodo pure_tot: ricostruzione TOT su tutti i canali ---
            if (cd.cfd_pure_tot_ok) {
                t_pu[k]  = (Float_t)cd.t_cfd_pure_tot;
                ok_pu[k] = true;
            } else {
                t_pu[k]  = -999.0f;
                ok_pu[k] = false;
            }
        }
// ---- Raccolta rise time per istogrammi (solo segnali puliti) ----
        for (int k = 0; k < 3; k++) {
            const ChannelData &cd_rt = e.ch[k];
            if (cd_rt.has_pulse && !cd_rt.is_clipped && !cd_rt.is_oscillating
                && cd_rt.rise_time_ok) {
                v_rt_pmt[k].push_back((float)cd_rt.rise_time);
            }
        }
        if (any_clipped) n_clipped_ev++;

        // --- Taglio di evento: has_pulse su tutti e 3 i canali ---
        if (!has_pulse_all) {
            n_too_few_ch++;
            tree->Fill();
            continue;
        }

        // --- Taglio di evento: oscillazione (vale per TUTTI i metodi) ---
        // Un'oscillazione e' una patologia della forma d'onda a monte della
        // ricostruzione: se ENABLE_OSC_CUT, l'evento e' inutilizzabile da
        // qualunque metodo.
        if (ENABLE_OSC_CUT && any_osc) {
            n_osc_event++;
            tree->Fill();
            continue;
        }

        // --- Validita' "tutti i canali" per ciascun metodo ---
        // Un metodo accetta l'evento solo se ha un tempo valido su TUTTI e 3 i
        // canali: una differenza temporale con un canale a -999 non avrebbe
        // senso fisico. Il flag ENABLE_CFD_CUT resta come interruttore
        // concettuale ma in pratica la richiesta dei 3 tempi e' sempre
        // necessaria per costruire Dt12/Dt13/Dt23/C.
        bool all_nc = ok_nc[0] && ok_nc[1] && ok_nc[2];
        bool all_hy = ok_hy[0] && ok_hy[1] && ok_hy[2];
        bool all_pu = ok_pu[0] && ok_pu[1] && ok_pu[2];

// --- Metodo NOCLIP ---
        if (all_nc) {
            dt12_noclip = t_nc[0] - t_nc[1];
            dt13_noclip = t_nc[0] - t_nc[2];
            dt23_noclip = t_nc[1] - t_nc[2];
            C_noclip    = t_nc[2] - (t_nc[0] + t_nc[1]) / 2.0f;
            good_noclip = 1;
            result.noclip.n_good++;
            v_dt12_nc.push_back(dt12_noclip);
            v_dt13_nc.push_back(dt13_noclip);
            v_dt23_nc.push_back(dt23_noclip);
            v_C_nc   .push_back(C_noclip);
        }

        // --- Metodo HYBRID_TOT (metodo fisico ufficiale) ---
        if (all_hy) {
            dt12_hybrid = t_hy[0] - t_hy[1];
            dt13_hybrid = t_hy[0] - t_hy[2];
            dt23_hybrid = t_hy[1] - t_hy[2];
            C_hybrid    = t_hy[2] - (t_hy[0] + t_hy[1]) / 2.0f;
            good_hybrid = 1;
            result.hybrid_tot.n_good++;
            v_dt12_hy.push_back(dt12_hybrid);
            v_dt13_hy.push_back(dt13_hybrid);
            v_dt23_hy.push_back(dt23_hybrid);
            v_C_hy   .push_back(C_hybrid);
        }

        // --- Metodo PURE_TOT (diagnostica di simmetria) ---
        if (all_pu) {
            dt12_pure = t_pu[0] - t_pu[1];
            dt13_pure = t_pu[0] - t_pu[2];
            dt23_pure = t_pu[1] - t_pu[2];
            C_pure    = t_pu[2] - (t_pu[0] + t_pu[1]) / 2.0f;
            good_pure = 1;
            result.pure_tot.n_good++;
            v_dt12_pu.push_back(dt12_pure);
            v_dt13_pu.push_back(dt13_pure);
            v_dt23_pu.push_back(dt23_pure);
            v_C_pu   .push_back(C_pure);
        }
// ----------------------------------------------------------------
        //  METODO CORRECTED (correzione Pietro&Rick del rise time variabile)
        // ----------------------------------------------------------------
        // Algoritmo (eq. 18-19 del documento P&R):
        //
        //   Per ciascun canale n, dT_n = t_30,n - t_10,n e' una stima del
        //   rise time effettivo nella regione del fronte rilevante per il
        //   timing. Sottraendolo al tempo CFD del canale "ancoriamo" il
        //   timing al crossing del 10% (zona di alto SNR sopra la baseline),
        //   rimuovendo la dipendenza dal rise time variabile e quindi dalla
        //   posizione x del cosmico sulla barra.
        //
        //   t_corr,n = t_cfd,n - dT_n
        //
        //   Dt_ij^corr = Dt_ij^hyb - (dT_i - dT_j)
        //   C^corr     = C^hyb    - [dT_3 - (dT_1+dT_2)/2]
        //
        // CRITERI DI ACCETTAZIONE (piu' stringenti di hybrid_tot):
        //   (1) l'evento e' gia' valido per hybrid_tot (all_hy);
        //   (2) NESSUN canale e' clippato (su un clippato amp e' sottostimata
        //       quindi le soglie 10%/30% non corrispondono alle frazioni
        //       vere e dT_n e' biased);
        //   (3) t_10 e t_30 sono entrambi disponibili su PMT1, PMT2 e PMT3.
        //
        // La statistica sara' ridotta rispetto a hybrid_tot, soprattutto
        // alle posizioni estreme (|x| ~ 130 cm) dove un PMT vede ampiezze
        // grandi e clippa spesso. Questo e' il prezzo della pulizia
        // metodologica: la correzione viene applicata SOLO dove le sue
        // ipotesi sono soddisfatte.
        bool no_clip_any =  !e.ch[0].is_clipped
                         && !e.ch[1].is_clipped
                         && !e.ch[2].is_clipped;
        bool t1030_all_ok =  e.ch[0].t10_ok && e.ch[0].t30_ok
                          && e.ch[1].t10_ok && e.ch[1].t30_ok
                          && e.ch[2].t10_ok && e.ch[2].t30_ok;

        if (all_hy && no_clip_any && t1030_all_ok) {
            // dT_n = t_30,n - t_10,n per i 3 canali
            double dT1 = e.ch[0].t_30 - e.ch[0].t_10;
            double dT2 = e.ch[1].t_30 - e.ch[1].t_10;
            double dT3 = e.ch[2].t_30 - e.ch[2].t_10;

            // Bias P&R del TOF (= tof_bias dell'eq. 18)
            double tof_bias_v = dT3 - 0.5 * (dT1 + dT2);

            // Tempi corretti = t_cfd - dT_n (matematicamente equivalente
            // a usare t_10 al posto di t_cfd se il modello esponenziale
            // fosse perfetto; qui restiamo data-driven e usiamo i tempi
            // CFD esistenti sottraendo dT)
            // double t1_c = t_hy[0] - dT1;  // (NB: non servono espliciti,
            // double t2_c = t_hy[1] - dT2;  //  uso direttamente l'algebra
            // double t3_c = t_hy[2] - dT3;  //  sui Dt_ij gia' calcolati)

            dt12_corr = dt12_hybrid - (Float_t)(dT1 - dT2);
            dt13_corr = dt13_hybrid - (Float_t)(dT1 - dT3);
            dt23_corr = dt23_hybrid - (Float_t)(dT2 - dT3);
            C_corr    = C_hybrid    - (Float_t)tof_bias_v;
            tof_bias_evt = (Float_t)tof_bias_v;
            good_corr = 1;
            result.corrected.n_good++;

            v_dt12_co.push_back(dt12_corr);
            v_dt13_co.push_back(dt13_corr);
            v_dt23_co.push_back(dt23_corr);
            v_C_co   .push_back(C_corr);
        }
        // --- Diagnostica dA = A_rec(TOT) - A_vera (solo canali NON clippati) ---
        // Sui non clippati A_vera = amplitude e' nota: confrontandola con
        // A_rec = EvalAmplitudeFromTOT(TOT) si misura bias ed errore del
        // recupero TOT. Riempie gli istogrammi globali aggregati su tutte le
        // posizioni (creati in TOF_Calibration).
        for (int k = 0; k < 3; k++) {
            const ChannelData &cd = e.ch[k];
            if (cd.is_clipped)    continue;
            if (!cd.has_pulse)    continue;
            if (!cd.tot_ok[kref]) continue;
            double A_rec = EvalAmplitudeFromTOT(cd.tot_q[kref], k);
            if (A_rec < 0.0)      continue;
            double dA = A_rec - cd.amplitude;
            if (g_h_dA[k])        g_h_dA[k]->Fill(dA);
            if (g_h_dA_vs_A[k])   g_h_dA_vs_A[k]->Fill(cd.amplitude, dA);
            if (g_h_dA_vs_TOT[k]) g_h_dA_vs_TOT[k]->Fill(cd.tot_q[kref], dA);
        }

        tree->Fill();
    }
// ====================================================================
    //  PASS 2: RANGE ADATTIVI, CREAZIONE E RIEMPIMENTO DEGLI ISTOGRAMMI
    // ====================================================================
    // Per ciascuno dei 12 vettori (4 variabili × 3 metodi), calcola la
    // mediana e la MAD (Median Absolute Deviation), determina il range
    // ottimale centrato sulla mediana, crea l'istogramma e lo riempie.
    // Il fattore 1.4826 converte la MAD in sigma-equivalente per una
    // distribuzione gaussiana: sigma_est = 1.4826 * MAD.
    // Il range e': [mediana - nsigma*sigma_est, mediana + nsigma*sigma_est]
    // con un floor di DT_RANGE_MIN per distribuzioni molto strette.

    // Lambda per calcolare il range adattivo da un vettore di valori.
    auto computeRange = [](std::vector<float> &vals,
                           double nsigma, double range_min,
                           double &lo, double &hi) {
        if (vals.empty()) { lo = -10.0; hi = 10.0; return; }
        std::sort(vals.begin(), vals.end());
        int n = (int)vals.size();
        double median = (n % 2 == 0)
                      ? 0.5 * ((double)vals[n/2 - 1] + (double)vals[n/2])
                      : (double)vals[n/2];
        // MAD = mediana delle deviazioni assolute dalla mediana
        std::vector<double> devs(n);
        for (int i = 0; i < n; i++) devs[i] = fabs((double)vals[i] - median);
        std::sort(devs.begin(), devs.end());
        double mad = (n % 2 == 0)
                   ? 0.5 * (devs[n/2 - 1] + devs[n/2])
                   : devs[n/2];
        double sigma_est  = 1.4826 * mad;   // sigma gaussiano equivalente
        double half_range = std::max(nsigma * sigma_est, range_min);
        lo = median - half_range;
        hi = median + half_range;
    };

    // Lambda per creare un istogramma con range adattivo e riempirlo.
    auto makeAndFillH = [&](const char* var, const char* method,
                            std::vector<float> &vals,
                            const char* axis) -> TH1D* {
        double lo, hi;
        computeRange(vals, DT_RANGE_NSIGMA, DT_RANGE_MIN, lo, hi);
        TH1D* h = new TH1D(
            Form("h_%s_%s_%s", var, method, tname.c_str()),
            Form("%s (%s) @ %s;%s;Conteggi", var, method, tname.c_str(), axis),
            NBINS_DT, lo, hi);
        for (size_t i = 0; i < vals.size(); i++) h->Fill(vals[i]);
        return h;
    };

    // Creazione e riempimento degli istogrammi con range adattivo.
    // noclip
    TH1D* h_dt12_nc = makeAndFillH("dt12", "noclip", v_dt12_nc, "#Delta t_{12} [ns]");
    TH1D* h_dt13_nc = makeAndFillH("dt13", "noclip", v_dt13_nc, "#Delta t_{13} [ns]");
    TH1D* h_dt23_nc = makeAndFillH("dt23", "noclip", v_dt23_nc, "#Delta t_{23} [ns]");
    TH1D* h_C_nc    = makeAndFillH("C",    "noclip", v_C_nc,    "C [ns]");
    // hybrid_tot
    TH1D* h_dt12_hy = makeAndFillH("dt12", "hybrid", v_dt12_hy, "#Delta t_{12} [ns]");
    TH1D* h_dt13_hy = makeAndFillH("dt13", "hybrid", v_dt13_hy, "#Delta t_{13} [ns]");
    TH1D* h_dt23_hy = makeAndFillH("dt23", "hybrid", v_dt23_hy, "#Delta t_{23} [ns]");
    TH1D* h_C_hy    = makeAndFillH("C",    "hybrid", v_C_hy,    "C [ns]");
    // pure_tot
    TH1D* h_dt12_pu = makeAndFillH("dt12", "pure",   v_dt12_pu, "#Delta t_{12} [ns]");
    TH1D* h_dt13_pu = makeAndFillH("dt13", "pure",   v_dt13_pu, "#Delta t_{13} [ns]");
    TH1D* h_dt23_pu = makeAndFillH("dt23", "pure",   v_dt23_pu, "#Delta t_{23} [ns]");
    TH1D* h_C_pu    = makeAndFillH("C",    "pure",   v_C_pu,    "C [ns]");
    // corrected (correzione P&R sul rise time variabile)
    TH1D* h_dt12_co = makeAndFillH("dt12", "corr", v_dt12_co, "#Delta t_{12} [ns]");
    TH1D* h_dt13_co = makeAndFillH("dt13", "corr", v_dt13_co, "#Delta t_{13} [ns]");
    TH1D* h_dt23_co = makeAndFillH("dt23", "corr", v_dt23_co, "#Delta t_{23} [ns]");
    TH1D* h_C_co    = makeAndFillH("C",    "corr", v_C_co,    "C [ns]");
    // ---- Stampa riepilogativa ----
    std::cout << "[INFO] " << tname << ": eventi totali = " << events.size()
              << std::endl;
    std::cout << "       buoni:  noclip = " << result.noclip.n_good
              << ",  hybrid_tot = " << result.hybrid_tot.n_good
              << ",  pure_tot = " << result.pure_tot.n_good << std::endl;
    std::cout << "       scarti: " << n_too_few_ch << " senza 3 canali/impulso, "
              << n_osc_event << " per oscillazione, "
              << n_clipped_ev << " con clipping (recuperabili)" << std::endl;

    // ---- Fit delle 12 distribuzioni con Gaussiana sul 90% centrale ----
    result.noclip.dt12 = FitGaussianCore(h_dt12_nc, 0.90, Form("fit_dt12_noclip_%s", tname.c_str()));
    result.noclip.dt13 = FitGaussianCore(h_dt13_nc, 0.90, Form("fit_dt13_noclip_%s", tname.c_str()));
    result.noclip.dt23 = FitGaussianCore(h_dt23_nc, 0.90, Form("fit_dt23_noclip_%s", tname.c_str()));
    result.noclip.C    = FitGaussianCore(h_C_nc,    0.90, Form("fit_C_noclip_%s",    tname.c_str()));

    result.hybrid_tot.dt12 = FitGaussianCore(h_dt12_hy, 0.90, Form("fit_dt12_hybrid_%s", tname.c_str()));
    result.hybrid_tot.dt13 = FitGaussianCore(h_dt13_hy, 0.90, Form("fit_dt13_hybrid_%s", tname.c_str()));
    result.hybrid_tot.dt23 = FitGaussianCore(h_dt23_hy, 0.90, Form("fit_dt23_hybrid_%s", tname.c_str()));
    result.hybrid_tot.C    = FitGaussianCore(h_C_hy,    0.90, Form("fit_C_hybrid_%s",    tname.c_str()));

    result.pure_tot.dt12 = FitGaussianCore(h_dt12_pu, 0.90, Form("fit_dt12_pure_%s", tname.c_str()));
    result.pure_tot.dt13 = FitGaussianCore(h_dt13_pu, 0.90, Form("fit_dt13_pure_%s", tname.c_str()));
    result.pure_tot.dt23 = FitGaussianCore(h_dt23_pu, 0.90, Form("fit_dt23_pure_%s", tname.c_str()));
    result.pure_tot.C    = FitGaussianCore(h_C_pu,    0.90, Form("fit_C_pure_%s",    tname.c_str()));
    // Fit del metodo "corrected"
    result.corrected.dt12 = FitGaussianCore(h_dt12_co, 0.90, Form("fit_dt12_corr_%s", tname.c_str()));
    result.corrected.dt13 = FitGaussianCore(h_dt13_co, 0.90, Form("fit_dt13_corr_%s", tname.c_str()));
    result.corrected.dt23 = FitGaussianCore(h_dt23_co, 0.90, Form("fit_dt23_corr_%s", tname.c_str()));
    result.corrected.C    = FitGaussianCore(h_C_co,    0.90, Form("fit_C_corr_%s",    tname.c_str()));
    // ---- Salvataggio TTree + istogrammi ----
    outfile->cd();
    tree->Write();
    h_dt12_nc->Write(); h_dt13_nc->Write(); h_dt23_nc->Write(); h_C_nc->Write();
    h_dt12_hy->Write(); h_dt13_hy->Write(); h_dt23_hy->Write(); h_C_hy->Write();
    h_dt12_pu->Write(); h_dt13_pu->Write(); h_dt23_pu->Write(); h_C_pu->Write();
    h_dt12_co->Write(); h_dt13_co->Write(); h_dt23_co->Write(); h_C_co->Write(); // salvataggio istogrammi metodo corrected

    // ====================================================================
    //  CANVAS COMPARATIVO PER QUESTO DATASET — raw vs corrected (P&R)
    // ====================================================================
    //  Per ogni posizione x_k, mostriamo in 4 pad gli istogrammi di
    //  Dt12/Dt13/Dt23/C calcolati col metodo hybrid_tot (rosso) e con la
    //  correzione di Pietro&Rick (blu). Sovrapponendoli si vede:
    //   - lo spostamento del centro (sistema di rimozione del bias)
    //   - l'eventuale riduzione (o aumento) della larghezza
    //
    //  Per un confronto onesto le due distribuzioni sono mostrate con la
    //  STESSA scala dell'asse X (range adattivo che ingloba entrambe).
    //  La legenda riporta media e sigma del fit gaussiano sul core 90%.
    {
        TCanvas* c_cmp = new TCanvas(
            Form("c_cmp_raw_corr_%s", tname.c_str()),
            Form("Confronto raw vs corrected (P&R) @ %s", tname.c_str()),
            1400, 1000);
        c_cmp->Divide(2, 2);

        // Lambda per disegnare una coppia (raw, corr) sovrapposti su un pad,
        // adattando il range X all'unione delle due distribuzioni.
        auto drawPair = [&](int pad, TH1D* h_raw, TH1D* h_cor,
                            const FitResult &fr_raw, const FitResult &fr_cor,
                            const char* var_label) {
            c_cmp->cd(pad);
            gPad->SetGrid(1, 1);

            // Stili dei due istogrammi
            h_raw->SetLineColor(kRed + 1);
            h_raw->SetLineWidth(2);
            h_raw->SetMarkerStyle(0);
            h_raw->SetStats(0);
            h_cor->SetLineColor(kBlue + 1);
            h_cor->SetLineWidth(2);
            h_cor->SetMarkerStyle(0);
            h_cor->SetStats(0);

            // Range X comune: unione dei due range, allargato per mostrare le code
            double xlo = std::min(h_raw->GetXaxis()->GetXmin(),
                                  h_cor->GetXaxis()->GetXmin());
            double xhi = std::max(h_raw->GetXaxis()->GetXmax(),
                                  h_cor->GetXaxis()->GetXmax());
            // Clonazione per evitare di toccare gli istogrammi originali
            TH1D* h_raw_d = (TH1D*)h_raw->Clone(Form("%s_d", h_raw->GetName()));
            TH1D* h_cor_d = (TH1D*)h_cor->Clone(Form("%s_d", h_cor->GetName()));
            h_raw_d->SetTitle(Form("%s @ %s;%s [ns];Conteggi",
                                   var_label, tname.c_str(), var_label));
            h_raw_d->GetXaxis()->SetRangeUser(xlo, xhi);

            double y_max = std::max(h_raw_d->GetMaximum(), h_cor_d->GetMaximum());
            h_raw_d->GetYaxis()->SetRangeUser(0, 1.15 * y_max);

            h_raw_d->Draw("HIST");
            h_cor_d->Draw("HIST SAME");

            // Legenda con risultati dei fit
            TLegend* leg = new TLegend(0.62, 0.65, 0.97, 0.92);
            leg->SetFillColorAlpha(kWhite, 0.85);
            leg->SetTextSize(0.030);
            leg->AddEntry(h_raw_d,
                Form("raw (hybrid): #mu=%.3f, #sigma=%.3f ns",
                     fr_raw.center, fabs(fr_raw.width)), "l");
            leg->AddEntry(h_cor_d,
                Form("corr (P&R):  #mu=%.3f, #sigma=%.3f ns",
                     fr_cor.center, fabs(fr_cor.width)), "l");
            // Shift della media dovuto alla correzione
            leg->AddEntry((TObject*)nullptr,
                Form("#Delta#mu = %+.3f ns",
                     fr_cor.center - fr_raw.center), "");
            leg->Draw();
        };

        drawPair(1, h_dt12_hy, h_dt12_co,
                 result.hybrid_tot.dt12, result.corrected.dt12,
                 "#Delta t_{12}");
        drawPair(2, h_dt13_hy, h_dt13_co,
                 result.hybrid_tot.dt13, result.corrected.dt13,
                 "#Delta t_{13}");
        drawPair(3, h_dt23_hy, h_dt23_co,
                 result.hybrid_tot.dt23, result.corrected.dt23,
                 "#Delta t_{23}");
        drawPair(4, h_C_hy,    h_C_co,
                 result.hybrid_tot.C,    result.corrected.C,
                 "C");

        c_cmp->Update();
        outfile->cd();
        c_cmp->Write(Form("Compare_raw_corr_%s", tname.c_str()));
    }
// ---- Istogrammi rise time per PMT (range adattivo + fit gaussiano) ----
    // Per ogni PMT creiamo un istogramma del rise time a questa posizione,
    // con range adattivo e fit gaussiano sul 90% centrale sovrapposto.
    const char* pmt_rt_names[3] = {"PMT1", "PMT2", "PMT3"};
    const int   pmt_rt_cols[3]  = {kBlue + 1, kRed + 1, kGreen + 2};
    for (int k = 0; k < 3; k++) {
        double rt_lo, rt_hi;
        computeRange(v_rt_pmt[k], DT_RANGE_NSIGMA, 2.0, rt_lo, rt_hi);
        TH1D* h_rt = new TH1D(
            Form("h_risetime_%s_%s", pmt_rt_names[k], tname.c_str()),
            Form("Rise Time T_{90}-T_{10} %s @ %s;"
                 "Rise Time [ns];Conteggi", pmt_rt_names[k], tname.c_str()),
            NBINS_DT, rt_lo, rt_hi);
        for (size_t i = 0; i < v_rt_pmt[k].size(); i++) h_rt->Fill(v_rt_pmt[k][i]);

// Fit del rise time con EMG (gaussiana + coda esponenziale destra),
        // vedi Punto 2: la coda a destra non e' descritta da una gaussiana.
        // La TF1 resta agganciata all'istogramma per essere disegnata; il
        // colore segue il PMT. Con pochi eventi/non convergenza c'e' fallback.
        if (h_rt->GetEntries() > 20) {
            FitRiseTimeEMG(
                h_rt,
                Form("f_rt_emg_%s_%s", pmt_rt_names[k], tname.c_str()),
                pmt_rt_cols[k],   // colore della curva = colore del PMT
                true);            // attach = true: la curva viene disegnata
        }

        outfile->cd();
        h_rt->Write();
    }
    return result;
}


// ==========================================================================
//  SEZIONE 10: FUNZIONE PRINCIPALE DI CALIBRAZIONE
// ==========================================================================

/// TOF_Calibration(): orchestratore della calibrazione.
///
///   1. Apre il file ROOT di output.
///   2. Calibra il modello A(TOT) per PMT1/2/3 (CalibrateTOT).
///   3. Crea gli istogrammi diagnostici dA aggregati.
///   4. Cicla sulle posizioni: per ognuna chiama ProcessDataset, raccoglie i
///      risultati dei 3 metodi e riempie il TTree "summary" (branch espliciti).
///   5. Costruisce la retta di calibrazione Dt12_hybrid(x) -> m, q, v_eff e il
///      fit lineare C_hybrid(x) -> C_lin_p0, C_lin_p1.
///   6. Disegna i canvas di output.
///   7. Scrive il TTree "fit_params" (CONTRATTO con TOF_Analysis_v9.cpp).
///   8. Chiude il file.
///
/// ARGOMENTI:
///   folder        — cartella contenente i file XML (es. "/home/lux_n/TOF_DRS/")
///   outname       — nome del file ROOT di output
///   parallax_file — file ROOT del MC di parallasse (nullptr = usa PARALLAX_FILE)
///
/// UTILIZZO:
///   root -l 'TOF_Calibration_v13.cpp("/percorso/dati/")'
void TOF_Calibration(const char* folder,
                     const char* outname = "TOF_Calibration_v13_output.root",
                     const char* parallax_file = nullptr) {

gROOT->SetBatch(kTRUE);   // disattiva le finestre grafiche interattive

    // Impedisce la registrazione delle TF1 nella lista globale di ROOT.
    // Senza questa riga, allo shutdown (.q) ROOT tenta di distruggere le TF1
    // orfane dalla lista globale: i loro puntatori interni riferiscono
    // istogrammi/canvas gia' cancellati da fout->Close(), causando un
    // segfault in RecursiveRemove. Con AddToGlobalList(false) le TF1 restano
    // associate agli istogrammi (per il disegno) ma ROOT non le tiene in
    // memoria dopo la chiusura del file -> shutdown pulito.
    TF1::DefaultAddToGlobalList(kFALSE);

    std::cout << "=============================================" << std::endl;
    std::cout << "  TOF Calibration v13 (modello A(TOT) cubico) " << std::endl;
    std::cout << "  Folder: " << folder  << std::endl;
    std::cout << "  Output: " << outname << std::endl;
    std::cout << "=============================================" << std::endl;

    // Stile grafico
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetTitleSize(0.045, "xy");
    gStyle->SetTitleOffset(1.0, "y");

    // Risoluzione del file di parallasse: argomento esplicito > costante > vuoto.
    const char* effective_parallax_file = parallax_file ? parallax_file : PARALLAX_FILE;
    std::cout << "  Parallax file: "
              << (effective_parallax_file && strlen(effective_parallax_file) > 0
                  ? effective_parallax_file
                  : "(non specificato — uso fallback costante)")
              << std::endl;
    std::cout << "=============================================" << std::endl;

    // ---- Apri il file di output ROOT ----
    TFile* fout = new TFile(outname, "RECREATE");
    if (!fout->IsOpen()) {
        std::cerr << "[ERRORE] Impossibile creare " << outname << std::endl;
        return;
    }

    // ====================================================================
    //  PASSO 1: CALIBRAZIONE DEL MODELLO A(TOT)
    // ====================================================================
    // Eseguita PRIMA del loop principale: serve a ricostruire i clippati
    // (hybrid_tot) e a popolare la modalita' PURE_TOT durante il processing.
    CalibrateTOT(folder, DEFAULT_POSITIONS, fout);

    // ====================================================================
    //  PASSO 2: ISTOGRAMMI DIAGNOSTICI dA = A_rec(TOT) - A_vera (aggregati)
    // ====================================================================
    // Creati qui (prima del loop), riempiti dentro ProcessDataset, scritti
    // alla fine. Coprono PMT1/2/3.
    const char* pmt_names[3] = {"PMT1", "PMT2", "PMT3"};
    fout->cd();
    for (int k = 0; k < 3; k++) {
        g_h_dA[k] = new TH1D(
            Form("h_dA_%s", pmt_names[k]),
            Form("#DeltaA = A_{rec}(TOT) - A_{vera} %s (non clippati);"
                 "#DeltaA [mV];Conteggi", pmt_names[k]),
            NBINS_DA, DA_LO, DA_HI);

        g_h_dA_vs_A[k] = new TH2F(
            Form("h_dA_vs_A_%s", pmt_names[k]),
            Form("#DeltaA vs A_{vera} %s;A_{vera} [mV];#DeltaA [mV]", pmt_names[k]),
            NBINS_DA2_A,  DA2_A_LO,  DA2_A_HI,
            NBINS_DA2_DA, DA2_DA_LO, DA2_DA_HI);

        g_h_dA_vs_TOT[k] = new TH2F(
            Form("h_dA_vs_TOT_%s", pmt_names[k]),
            Form("#DeltaA vs TOT %s;TOT [ns];#DeltaA [mV]", pmt_names[k]),
            NBINS_DA2_TOT, DA2_TOT_LO, DA2_TOT_HI,
            NBINS_DA2_DA,  DA2_DA_LO,  DA2_DA_HI);
    }

    // ====================================================================
    //  PASSO 3: TTree "summary" — risultati dei fit per posizione
    // ====================================================================
    // Convenzione di naming: SOLO branch espliciti per metodo (nessun nome
    // "legacy" senza suffisso). TOF_Analysis_v9.cpp legge da qui 'x' e
    // 'C_hybrid_mu' per costruire la funzione C(x).
    fout->cd();
    TTree* summary = new TTree("summary",
                               "Risultati calibrazione per posizione (3 metodi)");
    Float_t s_x, s_dx;
    Char_t  s_name[32];
    // hybrid_tot (metodo fisico ufficiale)
    Float_t s_dt12_hy_mu, s_dt12_hy_mu_err, s_dt12_hy_sig, s_dt12_hy_sig_err;
    Float_t s_dt13_hy_mu, s_dt13_hy_mu_err, s_dt13_hy_sig, s_dt13_hy_sig_err;
    Float_t s_dt23_hy_mu, s_dt23_hy_mu_err, s_dt23_hy_sig, s_dt23_hy_sig_err;
    Float_t s_C_hy_mu,    s_C_hy_mu_err,    s_C_hy_sig,    s_C_hy_sig_err;
    Int_t   s_n_good_hy;
    // noclip (riferimento)
    Float_t s_dt12_nc_mu, s_dt12_nc_mu_err, s_dt12_nc_sig, s_dt12_nc_sig_err;
    Float_t s_C_nc_mu,    s_C_nc_mu_err,    s_C_nc_sig,    s_C_nc_sig_err;
    Int_t   s_n_good_nc;
    // pure_tot (diagnostica)
    Float_t s_dt12_pu_mu, s_dt12_pu_mu_err, s_dt12_pu_sig, s_dt12_pu_sig_err;
    Float_t s_C_pu_mu,    s_C_pu_mu_err,    s_C_pu_sig,    s_C_pu_sig_err;
    Int_t   s_n_good_pu;

    summary->Branch("x",    &s_x,    "x/F");
    summary->Branch("dx",   &s_dx,   "dx/F");
    summary->Branch("name", s_name,  "name[32]/C");
    // --- hybrid_tot ---
    summary->Branch("dt12_hybrid_mu",        &s_dt12_hy_mu,      "dt12_hybrid_mu/F");
    summary->Branch("dt12_hybrid_mu_err",    &s_dt12_hy_mu_err,  "dt12_hybrid_mu_err/F");
    summary->Branch("dt12_hybrid_sigma",     &s_dt12_hy_sig,     "dt12_hybrid_sigma/F");
    summary->Branch("dt12_hybrid_sigma_err", &s_dt12_hy_sig_err, "dt12_hybrid_sigma_err/F");
    summary->Branch("dt13_hybrid_mu",        &s_dt13_hy_mu,      "dt13_hybrid_mu/F");
    summary->Branch("dt13_hybrid_mu_err",    &s_dt13_hy_mu_err,  "dt13_hybrid_mu_err/F");
    summary->Branch("dt13_hybrid_sigma",     &s_dt13_hy_sig,     "dt13_hybrid_sigma/F");
    summary->Branch("dt13_hybrid_sigma_err", &s_dt13_hy_sig_err, "dt13_hybrid_sigma_err/F");
    summary->Branch("dt23_hybrid_mu",        &s_dt23_hy_mu,      "dt23_hybrid_mu/F");
    summary->Branch("dt23_hybrid_mu_err",    &s_dt23_hy_mu_err,  "dt23_hybrid_mu_err/F");
    summary->Branch("dt23_hybrid_sigma",     &s_dt23_hy_sig,     "dt23_hybrid_sigma/F");
    summary->Branch("dt23_hybrid_sigma_err", &s_dt23_hy_sig_err, "dt23_hybrid_sigma_err/F");
    summary->Branch("C_hybrid_mu",           &s_C_hy_mu,         "C_hybrid_mu/F");
    summary->Branch("C_hybrid_mu_err",       &s_C_hy_mu_err,     "C_hybrid_mu_err/F");
    summary->Branch("C_hybrid_sigma",        &s_C_hy_sig,        "C_hybrid_sigma/F");
    summary->Branch("C_hybrid_sigma_err",    &s_C_hy_sig_err,    "C_hybrid_sigma_err/F");
    summary->Branch("n_good_hybrid",         &s_n_good_hy,       "n_good_hybrid/I");
    // --- noclip ---
    summary->Branch("dt12_noclip_mu",        &s_dt12_nc_mu,      "dt12_noclip_mu/F");
    summary->Branch("dt12_noclip_mu_err",    &s_dt12_nc_mu_err,  "dt12_noclip_mu_err/F");
    summary->Branch("dt12_noclip_sigma",     &s_dt12_nc_sig,     "dt12_noclip_sigma/F");
    summary->Branch("dt12_noclip_sigma_err", &s_dt12_nc_sig_err, "dt12_noclip_sigma_err/F");
    summary->Branch("C_noclip_mu",           &s_C_nc_mu,         "C_noclip_mu/F");
    summary->Branch("C_noclip_mu_err",       &s_C_nc_mu_err,     "C_noclip_mu_err/F");
    summary->Branch("C_noclip_sigma",        &s_C_nc_sig,        "C_noclip_sigma/F");
    summary->Branch("C_noclip_sigma_err",    &s_C_nc_sig_err,    "C_noclip_sigma_err/F");
    summary->Branch("n_good_noclip",         &s_n_good_nc,       "n_good_noclip/I");
    // --- pure_tot ---
    summary->Branch("dt12_pure_mu",          &s_dt12_pu_mu,      "dt12_pure_mu/F");
    summary->Branch("dt12_pure_mu_err",      &s_dt12_pu_mu_err,  "dt12_pure_mu_err/F");
    summary->Branch("dt12_pure_sigma",       &s_dt12_pu_sig,     "dt12_pure_sigma/F");
    summary->Branch("dt12_pure_sigma_err",   &s_dt12_pu_sig_err, "dt12_pure_sigma_err/F");
    summary->Branch("C_pure_mu",             &s_C_pu_mu,         "C_pure_mu/F");
    summary->Branch("C_pure_mu_err",         &s_C_pu_mu_err,     "C_pure_mu_err/F");
    summary->Branch("C_pure_sigma",          &s_C_pu_sig,        "C_pure_sigma/F");
    summary->Branch("C_pure_sigma_err",      &s_C_pu_sig_err,    "C_pure_sigma_err/F");
    summary->Branch("n_good_pure",           &s_n_good_pu,       "n_good_pure/I");

    // ---- Vettori per i grafici (riempiti posizione per posizione) ----
    std::vector<double> x_pos, dx_pos;
    std::vector<std::string> labels;
    // hybrid_tot
    std::vector<double> dt12_hy, dt12_hy_e, dt12_hy_w, dt12_hy_we, C_hy, C_hy_e;
    // noclip
    std::vector<double> dt12_nc, dt12_nc_e, C_nc, C_nc_e;
    // pure_tot
    std::vector<double> dt12_pu, dt12_pu_e, C_pu, C_pu_e;
    // Vettori per i fit globali del metodo "corrected" (correzione P&R)
    std::vector<double> dt12_co, dt12_co_e, dt12_co_w, dt12_co_we, C_co, C_co_e;
// Diagnostica Dt13(x) e Dt23(x) — equalizzazione PMT1/PMT2
    // Fisicamente: Dt13 = (L/2 + x)/v + delta_1 - delta_3 -> pendenza +1/v = +m/2
    //              Dt23 = (L/2 - x)/v + delta_2 - delta_3 -> pendenza -1/v = -m/2
    // Se le pendenze estratte dai fit lineari sono esattamente +m/2 e -m/2
    // l'equalizzazione temporale dei PMT e' corretta e il modello di propagazione
    // della luce nella barra e' consistente. Deviazioni indicano:
    //   - asimmetria cavi/elettronica non corretta;
    //   - time-walk differenziale tra i 2 PMT (sbilanciamento guadagni);
    //   - non-linearita' nella propagazione (attenuazione, dispersione).
    std::vector<double> dt13_hy, dt13_hy_e, dt23_hy, dt23_hy_e;
    std::vector<double> dt13_hy_w, dt13_hy_we, dt23_hy_w, dt23_hy_we;
    // ---- Vettori rise time medio per PMT in funzione di x ----
    // Raccolgono mediana e MAD del rise time ad ogni posizione, per i 3 PMT.
    std::vector<double> rt_mu[3], rt_mu_e[3];   // media troncata e suo errore
    // ====================================================================
    //  PASSO 4: LOOP SULLE POSIZIONI
    // ====================================================================
    for (size_t ip = 0; ip < DEFAULT_POSITIONS.size(); ip++) {

        const PositionInfo &pos = DEFAULT_POSITIONS[ip];

        // Risoluzione del path al file XML (con o senza slash finale).
        std::string xml_path = std::string(folder) + "/" + pos.name + ".xml";
        std::ifstream test_file(xml_path.c_str());
        if (!test_file.good()) {
            xml_path = std::string(folder) + pos.name + ".xml";
            test_file.open(xml_path.c_str());
            if (!test_file.good()) {
                std::cerr << "[WARNING] File non trovato: " << pos.name
                          << ".xml — salto questa posizione." << std::endl;
                continue;
            }
        }
        test_file.close();

        // Processa il dataset (3 metodi).
        DatasetResults dres = ProcessDataset(xml_path.c_str(), pos, fout);

        // Calcola l'incertezza sulla posizione e la salva nel risultato.
        double sx_k = ComputeSigmaX(pos.x_cm, effective_parallax_file);
        dres.dx = sx_k;

        // Scarta la posizione se il metodo fisico (hybrid_tot) ha troppo poco.
        if (!dres.hybrid_tot.dt12.fit_ok && dres.hybrid_tot.n_good < 30) {
            std::cerr << "[WARNING] Fit Dt12 hybrid_tot fallito per " << pos.name
                      << " — escludo la posizione dalla retta." << std::endl;
            continue;
        }

        // ---- Raccolta nei vettori per i grafici ----
        x_pos.push_back(pos.x_cm);
        dx_pos.push_back(sx_k);
        labels.push_back(pos.name);

        dt12_hy.push_back  (dres.hybrid_tot.dt12.center);
        dt12_hy_e.push_back(dres.hybrid_tot.dt12.center_err);
        dt12_hy_w.push_back(dres.hybrid_tot.dt12.width);
        dt12_hy_we.push_back(dres.hybrid_tot.dt12.width_err);
        C_hy.push_back     (dres.hybrid_tot.C.center);
        C_hy_e.push_back   (dres.hybrid_tot.C.center_err);

        dt12_nc.push_back  (dres.noclip.dt12.center);
        dt12_nc_e.push_back(dres.noclip.dt12.center_err);
        C_nc.push_back     (dres.noclip.C.center);
        C_nc_e.push_back   (dres.noclip.C.center_err);

        dt12_pu.push_back  (dres.pure_tot.dt12.center);
        dt12_pu_e.push_back(dres.pure_tot.dt12.center_err);
        C_pu.push_back     (dres.pure_tot.C.center);
        C_pu_e.push_back   (dres.pure_tot.C.center_err);
// Dt13 e Dt23 (hybrid_tot) per diagnostica equalizzazione
        dt13_hy.push_back   (dres.hybrid_tot.dt13.center);
        dt13_hy_e.push_back (dres.hybrid_tot.dt13.center_err);
        dt13_hy_w.push_back (dres.hybrid_tot.dt13.width);
        dt13_hy_we.push_back(dres.hybrid_tot.dt13.width_err);
        dt23_hy.push_back   (dres.hybrid_tot.dt23.center);
        dt23_hy_e.push_back (dres.hybrid_tot.dt23.center_err);
        dt23_hy_w.push_back (dres.hybrid_tot.dt23.width);
        dt23_hy_we.push_back(dres.hybrid_tot.dt23.width_err);
        // corrected (correzione P&R, valida solo dove non ci sono clippati)
        // Nota: per posizioni estreme dove la statistica corrected e' bassa
        // il fit puo' fallire; in tal caso il dato viene escluso dal grafico
        // globale ma resta presente per posizione nel file ROOT.
        if (dres.corrected.dt12.fit_ok && dres.corrected.C.fit_ok
            && dres.corrected.n_good >= 30) {
            // Allineiamo la raccolta corrected con le altre: prepariamo i
            // punti in vettori PARALLELI a x_pos (gia' aggiunto sopra).
            dt12_co.push_back   (dres.corrected.dt12.center);
            dt12_co_e.push_back (dres.corrected.dt12.center_err);
            dt12_co_w.push_back (dres.corrected.dt12.width);
            dt12_co_we.push_back(dres.corrected.dt12.width_err);
            C_co.push_back      (dres.corrected.C.center);
            C_co_e.push_back    (dres.corrected.C.center_err);
        } else {
            // Sentinel: questa posizione non contribuisce ai fit corrected.
            // Inseriamo NaN cosi' il TGraph li scarta automaticamente.
            dt12_co.push_back   (std::nan(""));
            dt12_co_e.push_back (0.0);
            dt12_co_w.push_back (std::nan(""));
            dt12_co_we.push_back(0.0);
            C_co.push_back      (std::nan(""));
            C_co_e.push_back    (0.0);
            std::cerr << "[INFO] " << pos.name
                      << ": metodo corrected escluso dai grafici globali "
                      << "(n_good = " << dres.corrected.n_good << ")"
                      << std::endl;
        }
        // ---- Raccolta rise time medio per PMT a questa posizione ----
        // Per ogni PMT, calcoliamo la mediana e l'errore robusto (MAD) del
        // rise time a questa posizione x_k, sui segnali non clippati e puliti.
        // Usiamo i dati dal TTree per-posizione appena scritto.
        // Ripetiamo il parsing perche' i vettori v_rt_pmt sono locali a
        // ProcessDataset e non ritornano. L'approccio piu' pulito e' leggere
        // gli istogrammi di rise time appena scritti nel file ROOT.
        for (int kp = 0; kp < 3; kp++) {
            const char* pmt_rt_tag[3] = {"PMT1", "PMT2", "PMT3"};
            TH1D* h_rt_k = (TH1D*)fout->Get(
                Form("h_risetime_%s_%s", pmt_rt_tag[kp], pos.name.c_str()));
            if (h_rt_k && h_rt_k->GetEntries() > 10) {
// Fit EMG (gaussiana + coda destra) per estrarre il valore
                // centrale (= moda del picco) e il suo errore, vedi Punto 2.
                // attach=false: serve solo il numero, non si ridisegna la curva.
                FitResult fr_rt = FitRiseTimeEMG(h_rt_k,
                    Form("fit_rt_%s_%s", pmt_rt_tag[kp], pos.name.c_str()),
                    -1,        // nessun colore
                    false);    // attach = false
                if (fr_rt.fit_ok) {
                    rt_mu[kp].push_back(fr_rt.center);
                    rt_mu_e[kp].push_back(fr_rt.center_err);
                } else {
                    // Fallback: media e errore dalla media dell'istogramma
                    double mean = h_rt_k->GetMean();
                    double rms  = h_rt_k->GetRMS();
                    double n    = h_rt_k->GetEntries();
                    rt_mu[kp].push_back(mean);
                    rt_mu_e[kp].push_back(rms / sqrt(n));
                }
            } else {
                // Nessun dato sufficiente: valore sentinel
                rt_mu[kp].push_back(-999.0);
                rt_mu_e[kp].push_back(0.0);
            }
        }
        // ---- Riempimento del TTree summary ----
        s_x  = (Float_t)pos.x_cm;
        s_dx = (Float_t)sx_k;
        strncpy(s_name, pos.name.c_str(), 31);
        s_name[31] = '\0';
        // hybrid_tot
        s_dt12_hy_mu      = (Float_t)dres.hybrid_tot.dt12.center;
        s_dt12_hy_mu_err  = (Float_t)dres.hybrid_tot.dt12.center_err;
        s_dt12_hy_sig     = (Float_t)dres.hybrid_tot.dt12.width;
        s_dt12_hy_sig_err = (Float_t)dres.hybrid_tot.dt12.width_err;
        s_dt13_hy_mu      = (Float_t)dres.hybrid_tot.dt13.center;
        s_dt13_hy_mu_err  = (Float_t)dres.hybrid_tot.dt13.center_err;
        s_dt13_hy_sig     = (Float_t)dres.hybrid_tot.dt13.width;
        s_dt13_hy_sig_err = (Float_t)dres.hybrid_tot.dt13.width_err;
        s_dt23_hy_mu      = (Float_t)dres.hybrid_tot.dt23.center;
        s_dt23_hy_mu_err  = (Float_t)dres.hybrid_tot.dt23.center_err;
        s_dt23_hy_sig     = (Float_t)dres.hybrid_tot.dt23.width;
        s_dt23_hy_sig_err = (Float_t)dres.hybrid_tot.dt23.width_err;
        s_C_hy_mu         = (Float_t)dres.hybrid_tot.C.center;
        s_C_hy_mu_err     = (Float_t)dres.hybrid_tot.C.center_err;
        s_C_hy_sig        = (Float_t)dres.hybrid_tot.C.width;
        s_C_hy_sig_err    = (Float_t)dres.hybrid_tot.C.width_err;
        s_n_good_hy       = dres.hybrid_tot.n_good;
        // noclip
        s_dt12_nc_mu      = (Float_t)dres.noclip.dt12.center;
        s_dt12_nc_mu_err  = (Float_t)dres.noclip.dt12.center_err;
        s_dt12_nc_sig     = (Float_t)dres.noclip.dt12.width;
        s_dt12_nc_sig_err = (Float_t)dres.noclip.dt12.width_err;
        s_C_nc_mu         = (Float_t)dres.noclip.C.center;
        s_C_nc_mu_err     = (Float_t)dres.noclip.C.center_err;
        s_C_nc_sig        = (Float_t)dres.noclip.C.width;
        s_C_nc_sig_err    = (Float_t)dres.noclip.C.width_err;
        s_n_good_nc       = dres.noclip.n_good;
        // pure_tot
        s_dt12_pu_mu      = (Float_t)dres.pure_tot.dt12.center;
        s_dt12_pu_mu_err  = (Float_t)dres.pure_tot.dt12.center_err;
        s_dt12_pu_sig     = (Float_t)dres.pure_tot.dt12.width;
        s_dt12_pu_sig_err = (Float_t)dres.pure_tot.dt12.width_err;
        s_C_pu_mu         = (Float_t)dres.pure_tot.C.center;
        s_C_pu_mu_err     = (Float_t)dres.pure_tot.C.center_err;
        s_C_pu_sig        = (Float_t)dres.pure_tot.C.width;
        s_C_pu_sig_err    = (Float_t)dres.pure_tot.C.width_err;
        s_n_good_pu       = dres.pure_tot.n_good;
        summary->Fill();

        std::cout << "[RESULT] " << pos.name << ": x = " << pos.x_cm << " cm" << std::endl;
        std::cout << "   hybrid_tot: Dt12 = "
                  << Form("%.3f +/- %.3f", dres.hybrid_tot.dt12.center,
                          dres.hybrid_tot.dt12.center_err)
                  << " ns, sigma = " << Form("%.3f", dres.hybrid_tot.dt12.width)
                  << " ns, C = "
                  << Form("%.3f +/- %.3f", dres.hybrid_tot.C.center,
                          dres.hybrid_tot.C.center_err)
                  << " ns  (N = " << dres.hybrid_tot.n_good << ")" << std::endl;
    }

    // Salva il TTree summary.
    fout->cd();
    summary->Write();

    // ---- Verifica numero minimo di punti ----
    int npts = (int)x_pos.size();
    if (npts < 3) {
        std::cerr << "[ERRORE] Solo " << npts << " posizioni valide — "
                  << "impossibile fare il fit lineare." << std::endl;
        fout->Close();
        return;
    }
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  Calibrazione con " << npts << " posizioni" << std::endl;
    std::cout << "=============================================" << std::endl;

    // Variabili dei fit lineari, dichiarate qui perche' servono sia ai canvas
    // sia al TTree fit_params.
    double m_fit = 0.0, m_err = 0.0, q_fit = 0.0, q_err = 0.0, chi2ndf = -1.0;
    double v_eff = 0.0, v_eff_err = 0.0;
    double C_lin_p0 = 0.0, C_lin_p0_err = 0.0;
    double C_lin_p1 = 0.0, C_lin_p1_err = 0.0, C_lin_chi2ndf = -1.0;

    // Modello quadratico C(x) = C_poly_p0 + C_poly_p1*x + C_poly_p2*x^2
    // Cattura la curvatura parabolica dovuta alla dipendenza ampiezza-distanza
    // di PMT3 e al time-walk residuo. E' il modello C(x) PRINCIPALE passato
    // all'analisi TOF, perche' riduce il residuo a ~50 ps RMS rispetto ai
    // ~110 ps del modello costante e ai ~90 ps del lineare.
    double C_poly_p0 = 0.0, C_poly_p0_err = 0.0;
    double C_poly_p1 = 0.0, C_poly_p1_err = 0.0;
    double C_poly_p2 = 0.0, C_poly_p2_err = 0.0;
    double C_poly_chi2ndf = -1.0;

    // Modello costante C = <C>: media pesata dei C(x_k) misurati.
    // Usato come alternativa semplice per stimare la sistematica introdotta
    // dalla scelta del modello C(x). La differenza poly vs const nel
    // risultato su beta quantifica l'errore sistematico di calibrazione.
    double C_const_val = 0.0, C_const_err = 0.0;
    double C_const_chi2ndf = -1.0;
// Parametri del fit lineare Dt12_corr = m_co * x + q_co (metodo corrected)
double m_co_fit = 0.0, m_co_err = 0.0, q_co_fit = 0.0, q_co_err = 0.0;
    double v_eff_co = 0.0, v_eff_co_err = 0.0, chi2ndf_co = -1.0;

    // ----------------------------------------------------------------------
    //  PARAMETRI DEL RANGE LINEARE RISTRETTO (corrected) — Punto Lin_Range
    // ----------------------------------------------------------------------
    //  Dichiarate qui (scope di funzione) affinche' siano visibili sia dal
    //  blocco che le calcola (canvas linrange) sia dal blocco che scrive
    //  fit_params (CONTRATTO con l'analisi). Inizializzate a sentinella: se il
    //  range lineare ha meno di 3 punti validi, restano a questi valori e
    //  l'analisi capira' che la calibrazione ristretta non e' disponibile.
    //  Retta Dt12_corr ristretta (per pipeline 2a, gia' pronta):
    double m_lin_co = 0.0, m_lin_co_err = 0.0;     // pendenza ristretta [ns/cm]
    double q_lin_co = 0.0, q_lin_co_err = 0.0;     // intercetta ristretta [ns]
    double chi2ndf_lin_co = -1.0;                  // chi2/ndf retta ristretta
    double v_eff_lin_co = 0.0, v_eff_lin_co_err = 0.0;
    //  C(x) ristretto, modello LINEARE C=C0+C1*x (per pipeline 2a, gia' pronto):
    double C_lin_co_p0 = 0.0, C_lin_co_p0_err = 0.0;
    double C_lin_co_p1 = 0.0, C_lin_co_p1_err = 0.0;
    double C_lin_co_chi2ndf = -1.0;
    //  C(x) ristretto, modello COSTANTE C=<C> (per pipeline 2b, quella scelta):
    double C_const_lin_co = 0.0, C_const_lin_co_err = 0.0;
    double C_const_lin_co_chi2ndf = -1.0;
    bool   linrange_co_valid = false;              // true se i fit ristretti ci sono

    // Parametri del fit polinomiale C_corr(x) = p0_co + p1_co*x + p2_co*x^2
    double C_co_poly_p0 = 0.0, C_co_poly_p0_err = 0.0;
    double C_co_poly_p1 = 0.0, C_co_poly_p1_err = 0.0;
    double C_co_poly_p2 = 0.0, C_co_poly_p2_err = 0.0;
    double C_co_poly_chi2ndf = -1.0;
    // ====================================================================
    //  CANVAS 1: RETTA DI CALIBRAZIONE Dt12_hybrid vs x (con residui)
    // ====================================================================
    // Fit lineare Dt12 = m*x + q. La pendenza m = 2/v_eff: v_eff = 2/|m|.
    // Pannello superiore: dati + retta di fit. Pannello inferiore: residui.
    {
        TCanvas *c_cal = new TCanvas("c_calibration",
                                     "Calibrazione #Delta t_{12} vs x (hybrid_tot)",
                                     900, 700);
        c_cal->Divide(1, 2);

        // --- Pad superiore: dati + fit ---
        c_cal->cd(1);
        gPad->SetPad(0.0, 0.35, 1.0, 1.0);
        gPad->SetBottomMargin(0.02);
        gPad->SetTopMargin(0.08);
        gPad->SetGrid(1, 1);

        TGraphErrors *g_cal = new TGraphErrors(npts, x_pos.data(), dt12_hy.data(),
                                               dx_pos.data(), dt12_hy_e.data());
        g_cal->SetMarkerStyle(20);
        g_cal->SetMarkerSize(1.0);
        g_cal->SetMarkerColor(kBlue + 1);
        g_cal->SetLineColor(kBlue + 1);
        g_cal->SetTitle(";Posizione x [cm];#Delta t_{12} = t_{1} - t_{2} [ns]");
        g_cal->GetXaxis()->SetLabelSize(0);   // etichette X mostrate dal pad residui

        TF1 *f_lin = new TF1("f_lin", "[0] + [1]*x", -150, 150);
        f_lin->SetParNames("q (offset)", "m (2/v)");
        g_cal->Fit(f_lin, "S R Q");
        g_cal->Draw("AP");

        m_fit   = f_lin->GetParameter(1);
        m_err   = f_lin->GetParError(1);
        q_fit   = f_lin->GetParameter(0);
        q_err   = f_lin->GetParError(0);
        chi2ndf = (f_lin->GetNDF() > 0) ? f_lin->GetChisquare() / f_lin->GetNDF() : -1.0;

        // v_eff = 2/|m| ; propagazione: sigma_v = 2*sigma_m / m^2
        v_eff     = 2.0 / fabs(m_fit);
        v_eff_err = 2.0 * m_err / (m_fit * m_fit);

        TPaveText *pt_cal = new TPaveText(0.15, 0.62, 0.55, 0.92, "NDC");
        pt_cal->SetFillColorAlpha(kWhite, 0.85);
        pt_cal->SetLineColor(kBlack);
        pt_cal->SetTextFont(42);
        pt_cal->SetTextSize(0.045);
        pt_cal->SetTextAlign(12);
        pt_cal->AddText("Fit: #Delta t_{12} = m #upoint x + q  (hybrid_tot)");
        pt_cal->AddText(Form("m = %.5f #pm %.5f ns/cm", m_fit, m_err));
        pt_cal->AddText(Form("q = %.3f #pm %.3f ns", q_fit, q_err));
        pt_cal->AddText(Form("#chi^{2}/ndf = %.2f", chi2ndf));
        pt_cal->AddText(Form("v_{eff} = %.2f #pm %.2f cm/ns", v_eff, v_eff_err));
        pt_cal->Draw();

        // --- Pad inferiore: residui ---
        c_cal->cd(2);
        gPad->SetPad(0.0, 0.0, 1.0, 0.35);
        gPad->SetTopMargin(0.02);
        gPad->SetBottomMargin(0.25);
        gPad->SetGrid(1, 1);

        std::vector<double> residui(npts);
        for (int i = 0; i < npts; i++)
            residui[i] = dt12_hy[i] - f_lin->Eval(x_pos[i]);

        TGraphErrors *g_res = new TGraphErrors(npts, x_pos.data(), residui.data(),
                                               dx_pos.data(), dt12_hy_e.data());
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

        double xmin_r = *std::min_element(x_pos.begin(), x_pos.end()) - 10;
        double xmax_r = *std::max_element(x_pos.begin(), x_pos.end()) + 10;
        TLine *l_zero = new TLine(xmin_r, 0, xmax_r, 0);
        l_zero->SetLineColor(kRed);
        l_zero->SetLineStyle(2);
        l_zero->Draw("same");

        c_cal->Update();
        fout->cd();
        c_cal->Write("Calibration_Dt12_vs_x");
    }
// ====================================================================
    //  CANVAS 1-bis: RETTA DI CALIBRAZIONE Dt12_corr vs x (correzione P&R)
    // ====================================================================
    // Stessa logica del CANVAS 1, ma sui Dt12 calcolati con la correzione
    // del rise time variabile alla Pietro&Rick. Confronto critico:
    //   - la pendenza m_co (-> v_eff_co) DEVE essere consistente con la
    //     pendenza m del fit hybrid_tot entro le incertezze: la
    //     correzione non deve modificare la velocita' di propagazione,
    //     che e' una proprieta' fisica intrinseca della barra;
    //   - l'intercetta q_co PUO' differire da q a causa dello shift
    //     globale introdotto dalla correzione (dT non si annulla in media).
    {
        // Filtra i punti validi (NaN scartati: vedi sentinel nel BLOCCO L)
        std::vector<double> xv, yv, xev, yev;
        for (int i = 0; i < npts; i++) {
            if (std::isnan(dt12_co[i])) continue;
            xv.push_back(x_pos[i]);
            yv.push_back(dt12_co[i]);
            xev.push_back(dx_pos[i]);
            yev.push_back(dt12_co_e[i]);
        }
        int npts_co = (int)xv.size();

        if (npts_co < 3) {
            std::cerr << "[WARNING] Solo " << npts_co
                      << " punti corrected validi: salto il fit Dt12_corr vs x"
                      << std::endl;
        } else {
            TCanvas *c_cal_co = new TCanvas("c_calibration_corr",
                "Calibrazione #Delta t_{12} vs x (corrected, P&R)",
                900, 700);
            c_cal_co->Divide(1, 2);

            // --- Pad superiore: dati + fit ---
            c_cal_co->cd(1);
            gPad->SetPad(0.0, 0.35, 1.0, 1.0);
            gPad->SetBottomMargin(0.02);
            gPad->SetTopMargin(0.08);
            gPad->SetGrid(1, 1);

            TGraphErrors *g_cal_co = new TGraphErrors(
                npts_co, xv.data(), yv.data(), xev.data(), yev.data());
            g_cal_co->SetMarkerStyle(20);
            g_cal_co->SetMarkerSize(1.0);
            g_cal_co->SetMarkerColor(kBlue + 1);
            g_cal_co->SetLineColor(kBlue + 1);
            g_cal_co->SetTitle(";Posizione x [cm];#Delta t_{12}^{corr} [ns]");
            g_cal_co->GetXaxis()->SetLabelSize(0);

            TF1 *f_lin_co = new TF1("f_lin_co", "[0] + [1]*x", -150, 150);
            f_lin_co->SetParNames("q_{corr}", "m_{corr}");
            g_cal_co->Fit(f_lin_co, "S R Q");
            g_cal_co->Draw("AP");

            m_co_fit = f_lin_co->GetParameter(1);
            m_co_err = f_lin_co->GetParError(1);
            q_co_fit = f_lin_co->GetParameter(0);
            q_co_err = f_lin_co->GetParError(0);
            chi2ndf_co = (f_lin_co->GetNDF() > 0)
                       ? f_lin_co->GetChisquare() / f_lin_co->GetNDF() : -1.0;

            v_eff_co     = 2.0 / fabs(m_co_fit);
            v_eff_co_err = 2.0 * m_co_err / (m_co_fit * m_co_fit);

            TPaveText *pt_co = new TPaveText(0.15, 0.55, 0.55, 0.92, "NDC");
            pt_co->SetFillColorAlpha(kWhite, 0.85);
            pt_co->SetLineColor(kBlack);
            pt_co->SetTextFont(42);
            pt_co->SetTextSize(0.045);
            pt_co->SetTextAlign(12);
            pt_co->AddText("Fit: #Delta t_{12}^{corr} = m_{c}#upointx + q_{c}  (corrected, P&R)");
            pt_co->AddText(Form("m_{c} = %.5f #pm %.5f ns/cm", m_co_fit, m_co_err));
            pt_co->AddText(Form("q_{c} = %.3f #pm %.3f ns",    q_co_fit, q_co_err));
            pt_co->AddText(Form("#chi^{2}/ndf = %.2f", chi2ndf_co));
            pt_co->AddText(Form("v_{eff,c} = %.2f #pm %.2f cm/ns",
                                v_eff_co, v_eff_co_err));
            // Confronto numerico con il fit raw (hybrid_tot)
            pt_co->AddText(Form("#Deltam = m_{c} - m = %+.5f ns/cm",
                                m_co_fit - m_fit));
            pt_co->Draw();

            // --- Pad inferiore: residui ---
            c_cal_co->cd(2);
            gPad->SetPad(0.0, 0.0, 1.0, 0.35);
            gPad->SetTopMargin(0.02);
            gPad->SetBottomMargin(0.25);
            gPad->SetGrid(1, 1);

            std::vector<double> residui_co(npts_co);
            for (int i = 0; i < npts_co; i++)
                residui_co[i] = yv[i] - f_lin_co->Eval(xv[i]);

            TGraphErrors *g_res_co = new TGraphErrors(
                npts_co, xv.data(), residui_co.data(), xev.data(), yev.data());
            g_res_co->SetMarkerStyle(20);
            g_res_co->SetMarkerSize(1.0);
            g_res_co->SetMarkerColor(kBlue + 1);
            g_res_co->SetLineColor(kBlue + 1);
            g_res_co->SetTitle(";Posizione x [cm];Residui [ns]");
            g_res_co->GetXaxis()->SetTitleSize(0.08);
            g_res_co->GetXaxis()->SetLabelSize(0.07);
            g_res_co->GetYaxis()->SetTitleSize(0.08);
            g_res_co->GetYaxis()->SetLabelSize(0.07);
            g_res_co->GetYaxis()->SetTitleOffset(0.5);
            g_res_co->Draw("AP");

            double xmin_r = *std::min_element(xv.begin(), xv.end()) - 10;
            double xmax_r = *std::max_element(xv.begin(), xv.end()) + 10;
            TLine *l_zero_co = new TLine(xmin_r, 0, xmax_r, 0);
            l_zero_co->SetLineColor(kRed);
            l_zero_co->SetLineStyle(2);
            l_zero_co->Draw("same");

            c_cal_co->Update();
            fout->cd();
            c_cal_co->Write("Calibration_Dt12_corr_vs_x");
        }
    }
    // ====================================================================
    //  CANVAS 2: RISOLUZIONE TEMPORALE sigma(Dt12) vs x (hybrid_tot)
    // ====================================================================
    {
        TCanvas *c_res = new TCanvas("c_resolution",
                                     "Risoluzione #sigma(#Delta t_{12}) vs x (hybrid_tot)",
                                     800, 500);
        c_res->SetGrid(1, 1);

        TGraphErrors *g_width = new TGraphErrors(npts, x_pos.data(), dt12_hy_w.data(),
                                                 dx_pos.data(), dt12_hy_we.data());
        g_width->SetMarkerStyle(21);
        g_width->SetMarkerSize(1.0);
        g_width->SetMarkerColor(kRed + 1);
        g_width->SetLineColor(kRed + 1);
        g_width->SetTitle("Risoluzione temporale vs posizione (hybrid_tot);"
                          "Posizione x [cm];#sigma(#Delta t_{12}) [ns]");
        g_width->Draw("AP");
        c_res->Update();
        fout->cd();
        c_res->Write("Resolution_Dt12_vs_x");
    }

    // ====================================================================
    //  CANVAS 3: FIT LINEARE C_hybrid(x) = C0 + C1*x (con residui)
    // ====================================================================
    // C(x) dovrebbe essere ~ costante; un termine lineare C1 != 0 segnala una
    // sistematica residua. I parametri C_lin_p0, C_lin_p1 vanno nel TTree
    // fit_params e sono il modello C(x) usato da TOF_Analysis_v9.cpp.
    {
        TCanvas *c_C = new TCanvas("c_C_lin",
                                   "Fit lineare C(x) (hybrid_tot)", 900, 700);
        c_C->Divide(1, 2);

        // --- Pad superiore: dati + fit ---
        c_C->cd(1);
        gPad->SetPad(0.0, 0.35, 1.0, 1.0);
        gPad->SetBottomMargin(0.02);
        gPad->SetGrid(1, 1);

        TGraphErrors *g_C = new TGraphErrors(npts, x_pos.data(), C_hy.data(),
                                             dx_pos.data(), C_hy_e.data());
        g_C->SetMarkerStyle(22);
        g_C->SetMarkerSize(1.3);
        g_C->SetMarkerColor(kGreen + 2);
        g_C->SetLineColor(kGreen + 2);
        g_C->SetTitle(";Posizione x [cm];C = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]");
        g_C->GetXaxis()->SetLabelSize(0);

TF1 *f_C_lin = new TF1("f_C_lin", "[0] + [1]*x", -150, 150);
        f_C_lin->SetParNames("C_{0}", "C_{1}");
        // Fit con opzione "S R Q" (NO opzione "E"/MINOS: con chi2/ndf elevato,
        // MINOS gonfia paradossalmente gli errori dei parametri perche' la
        // chi2 e' piatta nel parametro -> servirebbe variarlo molto per
        // alzare chi2 di 1. Qui usiamo gli errori Hessiani standard, poi li
        // re-scaliamo per sqrt(chi2/ndf) se >1: e' la pratica "consistency
        // rescaling", che riconosce esplicitamente che il modello lineare non
        // descrive perfettamente i dati e gonfia gli errori in modo realistico.
        g_C->Fit(f_C_lin, "S R Q");

        C_lin_p0      = f_C_lin->GetParameter(0);
        C_lin_p0_err  = f_C_lin->GetParError(0);
        C_lin_p1      = f_C_lin->GetParameter(1);
        C_lin_p1_err  = f_C_lin->GetParError(1);
        C_lin_chi2ndf = (f_C_lin->GetNDF() > 0)
                      ? f_C_lin->GetChisquare() / f_C_lin->GetNDF() : -1.0;

        // --- Consistency rescaling degli errori ---
        // Se chi2/ndf > 1, il modello lineare e' sotto-descritto dai dati e
        // gli errori statistici sui punti sono ottimistici. Si rinormalizzano
        // gli errori dei parametri per sqrt(chi2/ndf) in modo da avere errori
        // tali che il chi2 ridotto del fit "scalato" sia consistente con 1.
        double rescale_factor = 1.0;
        if (C_lin_chi2ndf > 1.0) {
            rescale_factor = sqrt(C_lin_chi2ndf);
            C_lin_p0_err *= rescale_factor;
            C_lin_p1_err *= rescale_factor;
            std::cout << "[INFO] C(x) fit: chi2/ndf = " << C_lin_chi2ndf
                      << " > 1 -> errori parametri rinormalizzati per sqrt("
                      << C_lin_chi2ndf << ") = " << rescale_factor << std::endl;
        }

        // Disegno manuale di dati e retta (il "draw automatico" del Fit non
        // viene usato per controllare colore/stile e per non sovrapporre
        // linee fantasma del fit precedente in canvas riusati).
        g_C->Draw("AP");
        f_C_lin->SetLineColor(kRed + 1);
        f_C_lin->SetLineWidth(2);
        f_C_lin->Draw("same");

        // Banda +/-1 sigma del fit (propagazione: sigma_C^2 = sigma_p0^2 + x^2*sigma_p1^2)
        const int nband = 200;
        std::vector<double> xb(nband), yb_up(nband), yb_dn(nband);
        for (int i = 0; i < nband; i++) {
            xb[i] = -150.0 + 300.0 * i / (nband - 1);
            double c_val = C_lin_p0 + C_lin_p1 * xb[i];
            double c_err = sqrt(C_lin_p0_err * C_lin_p0_err
                              + xb[i] * xb[i] * C_lin_p1_err * C_lin_p1_err);
            yb_up[i] = c_val + c_err;
            yb_dn[i] = c_val - c_err;
        }
        TGraph *g_band_up = new TGraph(nband, xb.data(), yb_up.data());
        TGraph *g_band_dn = new TGraph(nband, xb.data(), yb_dn.data());
        g_band_up->SetLineColor(kGreen - 7); g_band_up->SetLineStyle(7);
        g_band_dn->SetLineColor(kGreen - 7); g_band_dn->SetLineStyle(7);
        g_band_up->Draw("L same");
        g_band_dn->Draw("L same");
        g_C->Draw("P same");

        TPaveText *pt_C = new TPaveText(0.15, 0.15, 0.65, 0.42, "NDC");
        pt_C->SetFillColorAlpha(kWhite, 0.85);
        pt_C->SetTextFont(42);
        pt_C->SetTextSize(0.045);
        pt_C->SetTextAlign(12);
        pt_C->AddText("C(x) = C_{0} + C_{1} #upoint x  (hybrid_tot)");
        pt_C->AddText(Form("C_{0} = %.4f #pm %.4f ns", C_lin_p0, C_lin_p0_err));
        pt_C->AddText(Form("C_{1} = %.6f #pm %.6f ns/cm", C_lin_p1, C_lin_p1_err));
        pt_C->AddText(Form("#chi^{2}/ndf = %.2f", C_lin_chi2ndf));
        pt_C->Draw();

        // --- Pad inferiore: residui (pull) ---
        c_C->cd(2);
        gPad->SetPad(0.0, 0.0, 1.0, 0.35);
        gPad->SetTopMargin(0.02);
        gPad->SetBottomMargin(0.25);
        gPad->SetGrid(1, 1);

// Pull calcolati con gli errori ORIGINALI (non rescalati): cosi' il
        // pull plot mostra effettivamente quanto i punti deviano in unita' di
        // sigma_statistico, ed e' la diagnosi visiva della bonta' del modello.
        std::vector<double> pull_x(npts), pull_y(npts), pull_ex(npts), pull_ey(npts);
        double pull_max_abs = 0.0;
        for (int i = 0; i < npts; i++) {
            double C_fit_i = C_lin_p0 + C_lin_p1 * x_pos[i];
            double res_i   = C_hy[i] - C_fit_i;
            // Errore originale (pre-rescaling) per il pull statistico.
            double sigma_orig = (rescale_factor > 0.0) ? C_hy_e[i] : 1.0;
            pull_x[i]  = x_pos[i];
            pull_ex[i] = dx_pos[i];
            pull_y[i]  = (sigma_orig > 0) ? res_i / sigma_orig : 0.0;
            pull_ey[i] = (sigma_orig > 0) ? 1.0 : 0.0;
            if (fabs(pull_y[i]) > pull_max_abs) pull_max_abs = fabs(pull_y[i]);
        }
        TGraphErrors *g_pull = new TGraphErrors(npts, pull_x.data(), pull_y.data(),
                                                pull_ex.data(), pull_ey.data());
        g_pull->SetMarkerStyle(22);
        g_pull->SetMarkerSize(1.1);
        g_pull->SetMarkerColor(kGreen + 2);
        g_pull->SetLineColor(kGreen + 2);
        g_pull->SetTitle(";Posizione x [cm];Pull (C_{mis} - C_{fit}) / #sigma");
        g_pull->GetXaxis()->SetTitleSize(0.09);
        g_pull->GetXaxis()->SetLabelSize(0.08);
        g_pull->GetYaxis()->SetTitleSize(0.09);
        g_pull->GetYaxis()->SetTitleOffset(0.5);
        g_pull->GetYaxis()->SetLabelSize(0.08);
        // Range Y dinamico: di default +/-3.5 (caso "buono"), ma allargato
        // a +/-1.2 * max(|pull|) se i pull sono fuori scala (modello sotto-fit).
        double y_range = std::max(3.5, 1.2 * pull_max_abs);
        g_pull->GetYaxis()->SetRangeUser(-y_range, y_range);
        g_pull->GetXaxis()->SetLimits(-150, 150);
        g_pull->Draw("AP");

        TLine *l0 = new TLine(-150, 0, 150, 0);
        l0->SetLineColor(kBlack); l0->SetLineWidth(2); l0->SetLineStyle(2);
        l0->Draw("same");
        TLine *lp1 = new TLine(-150, 1, 150, 1);
        lp1->SetLineColor(kGray + 1); lp1->SetLineStyle(3); lp1->Draw("same");
        TLine *lm1 = new TLine(-150, -1, 150, -1);
        lm1->SetLineColor(kGray + 1); lm1->SetLineStyle(3); lm1->Draw("same");

c_C->Update();
        fout->cd();
        c_C->Write("Calibration_C_vs_x");
    }

   // ====================================================================
    //  CANVAS 3-ter / 3-quater: ANALISI NEL RANGE LINEARE [-84, +70] cm
    //  (metodo CORRECTED, correzione P&R del rise time variabile)
    // ====================================================================
    //  Vedi Punto 3, versione aggiornata: i punti usati per i fit ristretti
    //  sono quelli del metodo "corrected" (Dt12_corr, C_corr), NON hybrid_tot.
    //  Motivazione fisica: la correzione di Pietro&Rick rimuove il time-walk
    //  dovuto al rise time variabile con x. Proprio nella finestra [-84,+70] cm
    //  il rise time e' ~lineare in x, quindi la correzione e' ben definita:
    //  combinare "zona lineare" + "correzione rise-time" da' la stima di v_eff
    //  piu' libera sia dalle non-linearita' di bordo sia dal time-walk residuo.
    //  I vettori dt12_co/C_co sono PARALLELI a x_pos ma contengono NaN dove il
    //  fit corrected per posizione e' stato escluso: tali punti vanno saltati.
    //  Si producono due canvas diagnostici:
    //    - Calibration_Dt12_corr_vs_x_linrange : fit lineare -> v_eff ristretta
    //    - Calibration_C_corr_vs_x_linrange    : andamento di C^corr(x)
    //  Nessun parametro viene salvato in fit_params: sono diagnostiche.
    {
        // ---- Filtro: dentro la finestra lineare E con valore corrected valido ----
        std::vector<double> xl, dxl, dt12l, dt12l_e, Cl, Cl_e;
        for (int i = 0; i < npts; i++) {
            if (x_pos[i] < LINRANGE_X_LO || x_pos[i] > LINRANGE_X_HI) continue;
            // Scarta i punti dove il metodo corrected non e' disponibile (NaN).
            if (std::isnan(dt12_co[i]) || std::isnan(C_co[i]))        continue;
            xl.push_back(x_pos[i]);
            dxl.push_back(dx_pos[i]);
            dt12l.push_back(dt12_co[i]);
            dt12l_e.push_back(dt12_co_e[i]);
            Cl.push_back(C_co[i]);
            Cl_e.push_back(C_co_e[i]);
        }
        int nl = (int)xl.size();

        if (nl < 3) {
            std::cerr << "[WARNING] Range lineare [" << LINRANGE_X_LO << ","
                      << LINRANGE_X_HI << "] cm (corrected): solo " << nl
                      << " punti validi — canvas ristretti non generati." << std::endl;
        } else {
            std::cout << "\n[INFO] Analisi nel range lineare ["
                      << LINRANGE_X_LO << ", " << LINRANGE_X_HI
                      << "] cm (metodo corrected): " << nl << " posizioni." << std::endl;

            // Estremi grafici comodi per assi/fit.
            double xg_lo = LINRANGE_X_LO - 8.0;
            double xg_hi = LINRANGE_X_HI + 8.0;

            // ============================================================
            //  CANVAS 3-ter: RETTA Dt12_corr(x) NEL RANGE RISTRETTO + RESIDUI
            // ============================================================
            double m_lin = 0.0, m_lin_err = 0.0;
            double q_lin = 0.0, q_lin_err = 0.0;
            double chi2ndf_lin = -1.0, v_eff_lin = 0.0, v_eff_lin_err = 0.0;
            {
                TCanvas *c_lr = new TCanvas("c_dt12_corr_linrange",
                    "Calibrazione #Delta t_{12}^{corr} vs x (range lineare)", 900, 700);
                c_lr->Divide(1, 2);

                // --- Pad superiore: dati + fit ---
                c_lr->cd(1);
                gPad->SetPad(0.0, 0.35, 1.0, 1.0);
                gPad->SetBottomMargin(0.02);
                gPad->SetTopMargin(0.08);
                gPad->SetGrid(1, 1);

                TGraphErrors *g_lr = new TGraphErrors(nl, xl.data(), dt12l.data(),
                                                      dxl.data(), dt12l_e.data());
                g_lr->SetName("g_dt12_corr_linrange");
                g_lr->SetMarkerStyle(20);
                g_lr->SetMarkerSize(1.0);
                g_lr->SetMarkerColor(kMagenta + 2);
                g_lr->SetLineColor(kMagenta + 2);
                g_lr->SetTitle(";Posizione x [cm];"
                               "#Delta t_{12}^{corr} = t_{1} - t_{2} [ns]");
                g_lr->GetXaxis()->SetLabelSize(0);
                g_lr->GetXaxis()->SetLimits(xg_lo, xg_hi);

                TF1 *f_lr = new TF1("f_lin_corr_linrange", "[0] + [1]*x", xg_lo, xg_hi);
                f_lr->SetParNames("q (offset)", "m (2/v)");
                g_lr->Fit(f_lr, "S R Q");
                g_lr->Draw("AP");

                m_lin     = f_lr->GetParameter(1);
                m_lin_err = f_lr->GetParError(1);
                q_lin     = f_lr->GetParameter(0);
                q_lin_err = f_lr->GetParError(0);
                chi2ndf_lin = (f_lr->GetNDF() > 0)
                            ? f_lr->GetChisquare() / f_lr->GetNDF() : -1.0;
                v_eff_lin     = 2.0 / fabs(m_lin);
                v_eff_lin_err = 2.0 * m_lin_err / (m_lin * m_lin);

                // Copia nelle variabili di scope funzione -> fit_params (CONTRATTO).
                // Servono alla pipeline 2a (retta ristretta), gia' pronta per
                // l'analisi anche se la pipeline scelta ora e' la 2b.
                m_lin_co         = m_lin;
                m_lin_co_err     = m_lin_err;
                q_lin_co         = q_lin;
                q_lin_co_err     = q_lin_err;
                chi2ndf_lin_co   = chi2ndf_lin;
                v_eff_lin_co     = v_eff_lin;
                v_eff_lin_co_err = v_eff_lin_err;
                linrange_co_valid = true;   // almeno la retta ristretta c'e'   

                // Box di confronto: TRE valori di v_eff per disaccoppiare
                //   - effetto bordi   : ristretto (corr) vs globale (corr)
                //   - effetto timewalk: corrected      vs hybrid
                TPaveText *pt_lr = new TPaveText(0.13, 0.55, 0.62, 0.92, "NDC");
                pt_lr->SetFillColorAlpha(kWhite, 0.85);
                pt_lr->SetLineColor(kBlack);
                pt_lr->SetTextFont(42);
                pt_lr->SetTextSize(0.038);
                pt_lr->SetTextAlign(12);
                pt_lr->AddText(Form("Fit #Delta t_{12}^{corr} = m x + q  (x #in [%.0f, %.0f] cm)",
                                    LINRANGE_X_LO, LINRANGE_X_HI));
                pt_lr->AddText(Form("m = %.5f #pm %.5f ns/cm", m_lin, m_lin_err));
                pt_lr->AddText(Form("q = %.3f #pm %.3f ns", q_lin, q_lin_err));
                pt_lr->AddText(Form("#chi^{2}/ndf = %.2f", chi2ndf_lin));
                pt_lr->AddText(Form("v_{eff} = %.2f #pm %.2f cm/ns  (corr, ristretto)",
                                    v_eff_lin, v_eff_lin_err));
                // Riferimento: corrected globale (puo' essere n/d se il fit globale fallisce).
                if (v_eff_co > 0.0)
                    pt_lr->AddText(Form("globale corr:   v_{eff} = %.2f #pm %.2f cm/ns",
                                        v_eff_co, v_eff_co_err));
                else
                    pt_lr->AddText("globale corr:   v_{eff} = n/d (fit corrected non valido)");
                // Riferimento: hybrid globale (metodo fisico ufficiale).
                pt_lr->AddText(Form("globale hybrid: v_{eff} = %.2f #pm %.2f cm/ns",
                                    v_eff, v_eff_err));
                pt_lr->Draw();

                // --- Pad inferiore: residui ---
                c_lr->cd(2);
                gPad->SetPad(0.0, 0.0, 1.0, 0.35);
                gPad->SetTopMargin(0.02);
                gPad->SetBottomMargin(0.25);
                gPad->SetGrid(1, 1);

                std::vector<double> res_lr(nl);
                for (int i = 0; i < nl; i++)
                    res_lr[i] = dt12l[i] - f_lr->Eval(xl[i]);

                TGraphErrors *g_res_lr = new TGraphErrors(nl, xl.data(), res_lr.data(),
                                                          dxl.data(), dt12l_e.data());
                g_res_lr->SetMarkerStyle(20);
                g_res_lr->SetMarkerSize(1.0);
                g_res_lr->SetMarkerColor(kMagenta + 2);
                g_res_lr->SetLineColor(kMagenta + 2);
                g_res_lr->SetTitle(";Posizione x [cm];Residui [ns]");
                g_res_lr->GetXaxis()->SetTitleSize(0.08);
                g_res_lr->GetXaxis()->SetLabelSize(0.07);
                g_res_lr->GetYaxis()->SetTitleSize(0.08);
                g_res_lr->GetYaxis()->SetLabelSize(0.07);
                g_res_lr->GetYaxis()->SetTitleOffset(0.5);
                g_res_lr->GetXaxis()->SetLimits(xg_lo, xg_hi);
                g_res_lr->Draw("AP");

                TLine *l0_lr = new TLine(xg_lo, 0, xg_hi, 0);
                l0_lr->SetLineColor(kRed);
                l0_lr->SetLineStyle(2);
                l0_lr->Draw("same");

                c_lr->Update();
                fout->cd();
                c_lr->Write("Calibration_Dt12_corr_vs_x_linrange");
            }

            // ============================================================
            //  CANVAS 3-quater: C_corr(x) NEL RANGE RISTRETTO + RESIDUI (lineare)
            // ============================================================
            //  Pad superiore: punti C_corr(x_k) con (a) fit lineare C=C0+C1*x e
            //  (b) media costante <C^corr> sovrapposte. Pad inferiore: residui
            //  rispetto al FIT LINEARE, per vederne l'andamento. Il fit lineare
            //  ristretto serve alla pipeline 2a (gia' pronto); la media e' la
            //  base del fit costante 2b realizzato nel canvas dedicato sotto.
            {
                TCanvas *c_Clr = new TCanvas("c_C_corr_linrange",
                    "C^{corr}(x) nel range lineare (lineare + residui)", 900, 700);
                c_Clr->Divide(1, 2);

                // --- Pad superiore: dati + fit lineare + media ---
                c_Clr->cd(1);
                gPad->SetPad(0.0, 0.35, 1.0, 1.0);
                gPad->SetBottomMargin(0.02);
                gPad->SetTopMargin(0.08);
                gPad->SetGrid(1, 1);

                TGraphErrors *g_Clr = new TGraphErrors(nl, xl.data(), Cl.data(),
                                                       dxl.data(), Cl_e.data());
                g_Clr->SetName("g_C_corr_linrange");
                g_Clr->SetMarkerStyle(22);
                g_Clr->SetMarkerSize(1.3);
                g_Clr->SetMarkerColor(kGreen + 2);
                g_Clr->SetLineColor(kGreen + 2);
                g_Clr->SetTitle(Form("C^{corr}(x) nel range [%.0f, %.0f] cm;"
                    ";C^{corr} = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]",
                    LINRANGE_X_LO, LINRANGE_X_HI));
                g_Clr->GetXaxis()->SetLabelSize(0);
                g_Clr->GetXaxis()->SetLimits(xg_lo, xg_hi);
                g_Clr->Draw("AP");

                // (a) Fit lineare C^corr(x) = C0 + C1*x.
                TF1 *f_Clr = new TF1("f_C_lin_corr_linrange", "[0] + [1]*x", xg_lo, xg_hi);
                f_Clr->SetParNames("C_{0}", "C_{1}");
                {
                    double sum_Clr = 0.0;
                    for (int i = 0; i < nl; i++) sum_Clr += Cl[i];
                    f_Clr->SetParameters(sum_Clr / nl, 0.0);
                }
                g_Clr->Fit(f_Clr, "S R Q");
                double C0_lr     = f_Clr->GetParameter(0);
                double C0_lr_err = f_Clr->GetParError(0);
                double C1_lr     = f_Clr->GetParameter(1);
                double C1_lr_err = f_Clr->GetParError(1);
                double C_lr_chi2 = (f_Clr->GetNDF() > 0)
                                 ? f_Clr->GetChisquare() / f_Clr->GetNDF() : -1.0;
                f_Clr->SetLineColor(kRed + 1);
                f_Clr->SetLineWidth(2);
                f_Clr->Draw("same");

                // (b) Media costante <C^corr> (media pesata dei punti nel range).
                double sw = 0.0, swc = 0.0;
                for (int i = 0; i < nl; i++) {
                    double w = (Cl_e[i] > 0.0) ? 1.0 / (Cl_e[i] * Cl_e[i]) : 0.0;
                    sw  += w;
                    swc += w * Cl[i];
                }
                double C_mean     = (sw > 0.0) ? swc / sw : 0.0;
                double C_mean_err = (sw > 0.0) ? sqrt(1.0 / sw) : 0.0;
                TLine *l_mean = new TLine(xg_lo, C_mean, xg_hi, C_mean);
                l_mean->SetLineColor(kBlue + 1);
                l_mean->SetLineStyle(2);
                l_mean->SetLineWidth(2);
                l_mean->Draw("same");

                g_Clr->Draw("P same");

                TLegend *leg_Clr = new TLegend(0.55, 0.72, 0.89, 0.89);
                leg_Clr->SetTextFont(42);
                leg_Clr->SetTextSize(0.032);
                leg_Clr->AddEntry(g_Clr, "C^{corr}(x_{k}) misurati", "pe");
                leg_Clr->AddEntry(f_Clr,
                    Form("Fit lin.: C_{1} = %.5f #pm %.5f ns/cm", C1_lr, C1_lr_err), "l");
                leg_Clr->AddEntry(l_mean,
                    Form("<C^{corr}> = %.4f #pm %.4f ns", C_mean, C_mean_err), "l");
                leg_Clr->Draw();

                TPaveText *pt_Clr = new TPaveText(0.13, 0.13, 0.55, 0.30, "NDC");
                pt_Clr->SetFillColorAlpha(kWhite, 0.85);
                pt_Clr->SetTextFont(42);
                pt_Clr->SetTextSize(0.032);
                pt_Clr->SetTextAlign(12);
                pt_Clr->AddText(Form("C_{0} = %.4f #pm %.4f ns", C0_lr, C0_lr_err));
                pt_Clr->AddText(Form("#chi^{2}/ndf (lin.) = %.2f", C_lr_chi2));
                pt_Clr->Draw();

                // --- Pad inferiore: residui rispetto al FIT LINEARE ---
                c_Clr->cd(2);
                gPad->SetPad(0.0, 0.0, 1.0, 0.35);
                gPad->SetTopMargin(0.02);
                gPad->SetBottomMargin(0.25);
                gPad->SetGrid(1, 1);

                std::vector<double> res_Clin(nl);
                for (int i = 0; i < nl; i++)
                    res_Clin[i] = Cl[i] - f_Clr->Eval(xl[i]);

                TGraphErrors *g_res_Clin = new TGraphErrors(nl, xl.data(),
                                              res_Clin.data(), dxl.data(), Cl_e.data());
                g_res_Clin->SetMarkerStyle(22);
                g_res_Clin->SetMarkerSize(1.1);
                g_res_Clin->SetMarkerColor(kRed + 1);
                g_res_Clin->SetLineColor(kRed + 1);
                g_res_Clin->SetTitle(";Posizione x [cm];Residui lin. [ns]");
                g_res_Clin->GetXaxis()->SetTitleSize(0.08);
                g_res_Clin->GetXaxis()->SetLabelSize(0.07);
                g_res_Clin->GetYaxis()->SetTitleSize(0.08);
                g_res_Clin->GetYaxis()->SetLabelSize(0.07);
                g_res_Clin->GetYaxis()->SetTitleOffset(0.5);
                g_res_Clin->GetXaxis()->SetLimits(xg_lo, xg_hi);
                g_res_Clin->Draw("AP");
                TLine *l0_Clin = new TLine(xg_lo, 0, xg_hi, 0);
                l0_Clin->SetLineColor(kBlack);
                l0_Clin->SetLineStyle(2);
                l0_Clin->Draw("same");

                c_Clr->Update();
                fout->cd();
                c_Clr->Write("Calibration_C_corr_vs_x_linrange");

                // Copia nelle variabili di scope funzione: C lineare ristretto
                // (pipeline 2a, gia' pronto per l'analisi).
                C_lin_co_p0      = C0_lr;
                C_lin_co_p0_err  = C0_lr_err;
                C_lin_co_p1      = C1_lr;
                C_lin_co_p1_err  = C1_lr_err;
                C_lin_co_chi2ndf = C_lr_chi2;

                std::cout << "[INFO] Range lineare (corrected): v_eff = "
                          << Form("%.2f +/- %.2f", v_eff_lin, v_eff_lin_err)
                          << " cm/ns" << std::endl;
                if (v_eff_co > 0.0)
                    std::cout << "       confronto: globale corr = "
                              << Form("%.2f +/- %.2f", v_eff_co, v_eff_co_err)
                              << " cm/ns, globale hybrid = "
                              << Form("%.2f +/- %.2f", v_eff, v_eff_err) << " cm/ns"
                              << std::endl;
                else
                    std::cout << "       confronto: globale corr = n/d, "
                              << "globale hybrid = "
                              << Form("%.2f +/- %.2f", v_eff, v_eff_err) << " cm/ns"
                              << std::endl;
                std::cout << "       C^corr(x) range lineare: C1 = "
                          << Form("%.5f +/- %.5f", C1_lr, C1_lr_err)
                          << " ns/cm, <C^corr> = "
                          << Form("%.4f +/- %.4f", C_mean, C_mean_err) << " ns"
                          << std::endl;

                // ========================================================
                //  CANVAS 3-quinquies: FIT COSTANTE di C_corr(x) + RESIDUI (2b)
                // ========================================================
                //  Pipeline SCELTA (2b): C(x) nel range ristretto modellato da
                //  una COSTANTE. Qui realizziamo un vero fit chi^2 con TF1 "[0]"
                //  inizializzata alla media semplice dei C_corr nel range. Il
                //  pad inferiore mostra i residui (C_k - costante): se piatti e
                //  centrati su zero, il modello costante e' adeguato.
                {
                    TCanvas *c_Ccst = new TCanvas("c_C_corr_linrange_const",
                        "C^{corr}(x) range lineare: fit costante + residui", 900, 700);
                    c_Ccst->Divide(1, 2);

                    // --- Pad superiore: dati + fit costante ---
                    c_Ccst->cd(1);
                    gPad->SetPad(0.0, 0.35, 1.0, 1.0);
                    gPad->SetBottomMargin(0.02);
                    gPad->SetTopMargin(0.08);
                    gPad->SetGrid(1, 1);

                    TGraphErrors *g_Ccst = new TGraphErrors(nl, xl.data(), Cl.data(),
                                                           dxl.data(), Cl_e.data());
                    g_Ccst->SetName("g_C_corr_linrange_const");
                    g_Ccst->SetMarkerStyle(22);
                    g_Ccst->SetMarkerSize(1.3);
                    g_Ccst->SetMarkerColor(kGreen + 2);
                    g_Ccst->SetLineColor(kGreen + 2);
                    g_Ccst->SetTitle(Form("C^{corr}(x) range [%.0f, %.0f] cm: fit costante;"
                        ";C^{corr} = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]",
                        LINRANGE_X_LO, LINRANGE_X_HI));
                    g_Ccst->GetXaxis()->SetLabelSize(0);
                    g_Ccst->GetXaxis()->SetLimits(xg_lo, xg_hi);
                    g_Ccst->Draw("AP");

                    // Fit costante chi^2: TF1 "[0]", init = media semplice.
                    TF1 *f_const = new TF1("f_C_const_corr_linrange", "[0]", xg_lo, xg_hi);
                    f_const->SetParName(0, "C_{const}");
                    {
                        double sum_simple = 0.0;
                        for (int i = 0; i < nl; i++) sum_simple += Cl[i];
                        f_const->SetParameter(0, (nl > 0) ? sum_simple / nl : 0.0);
                    }
                    g_Ccst->Fit(f_const, "S R Q");
                    double Cc      = f_const->GetParameter(0);
                    double Cc_err  = f_const->GetParError(0);
                    double Cc_chi2 = (f_const->GetNDF() > 0)
                                   ? f_const->GetChisquare() / f_const->GetNDF() : -1.0;
                    f_const->SetLineColor(kBlue + 1);
                    f_const->SetLineWidth(2);
                    f_const->Draw("same");
                    g_Ccst->Draw("P same");

                    TPaveText *pt_Ccst = new TPaveText(0.50, 0.70, 0.89, 0.90, "NDC");
                    pt_Ccst->SetFillColorAlpha(kWhite, 0.85);
                    pt_Ccst->SetLineColor(kBlack);
                    pt_Ccst->SetTextFont(42);
                    pt_Ccst->SetTextSize(0.034);
                    pt_Ccst->SetTextAlign(12);
                    pt_Ccst->AddText(Form("Fit costante (x #in [%.0f, %.0f] cm)",
                                          LINRANGE_X_LO, LINRANGE_X_HI));
                    pt_Ccst->AddText(Form("C_{const} = %.4f #pm %.4f ns", Cc, Cc_err));
                    pt_Ccst->AddText(Form("#chi^{2}/ndf = %.2f", Cc_chi2));
                    pt_Ccst->Draw();

                    // --- Pad inferiore: residui rispetto alla costante ---
                    c_Ccst->cd(2);
                    gPad->SetPad(0.0, 0.0, 1.0, 0.35);
                    gPad->SetTopMargin(0.02);
                    gPad->SetBottomMargin(0.25);
                    gPad->SetGrid(1, 1);

                    std::vector<double> res_Cc(nl);
                    for (int i = 0; i < nl; i++) res_Cc[i] = Cl[i] - Cc;

                    TGraphErrors *g_res_Cc = new TGraphErrors(nl, xl.data(),
                                                res_Cc.data(), dxl.data(), Cl_e.data());
                    g_res_Cc->SetMarkerStyle(22);
                    g_res_Cc->SetMarkerSize(1.1);
                    g_res_Cc->SetMarkerColor(kBlue + 1);
                    g_res_Cc->SetLineColor(kBlue + 1);
                    g_res_Cc->SetTitle(";Posizione x [cm];Residui cost. [ns]");
                    g_res_Cc->GetXaxis()->SetTitleSize(0.08);
                    g_res_Cc->GetXaxis()->SetLabelSize(0.07);
                    g_res_Cc->GetYaxis()->SetTitleSize(0.08);
                    g_res_Cc->GetYaxis()->SetLabelSize(0.07);
                    g_res_Cc->GetYaxis()->SetTitleOffset(0.5);
                    g_res_Cc->GetXaxis()->SetLimits(xg_lo, xg_hi);
                    g_res_Cc->Draw("AP");
                    TLine *l0_Cc = new TLine(xg_lo, 0, xg_hi, 0);
                    l0_Cc->SetLineColor(kBlack);
                    l0_Cc->SetLineStyle(2);
                    l0_Cc->Draw("same");

                    c_Ccst->Update();
                    fout->cd();
                    c_Ccst->Write("Calibration_C_corr_vs_x_linrange_const");

                    // Copia nelle variabili di scope funzione: C COSTANTE
                    // ristretto (pipeline 2b, quella scelta).
                    C_const_lin_co        = Cc;
                    C_const_lin_co_err    = Cc_err;
                    C_const_lin_co_chi2ndf = Cc_chi2;

                    std::cout << "       [2b] C costante ristretto: C_const = "
                              << Form("%.4f +/- %.4f", Cc, Cc_err)
                              << " ns  (chi2/ndf = " << Form("%.2f", Cc_chi2) << ")"
                              << std::endl;
                }
            }
        }
    }
// ====================================================================
    //  CANVAS 3-bis: FIT POLINOMIALE C_corr(x) = p0+p1*x+p2*x^2 (corrected)
    // ====================================================================
    // Test critico della correzione di Pietro&Rick: se il rise time
    // variabile e' la causa dominante della curvatura parabolica di C(x),
    // dopo aver sottratto il bias dT_n event-by-event il coefficiente
    // quadratico p2 deve scendere significativamente verso zero.
    //
    // Il valore di p2 nel raw e' tipicamente ~ -1.5e-5 ns/cm^2 (parabola
    // capovolta, vedi grafico c_xK.png con curvatura di ~0.4 ns su 2.8 m).
    // Aspettativa P&R: |p2_corr| < |p2_raw| / 2 e/o sigma_residuo
    // notevolmente piu' piccolo. Se non avviene -> la correzione P&R con
    // coefficiente 1 e' inadeguata e serve un modello piu' raffinato
    // (es. coefficiente determinato dai dati, BLOCCO discusso in chat).
    {
        // Filtra i punti validi (NaN scartati)
        std::vector<double> xv, yv, xev, yev;
        for (int i = 0; i < npts; i++) {
            if (std::isnan(C_co[i])) continue;
            xv.push_back(x_pos[i]);
            yv.push_back(C_co[i]);
            xev.push_back(dx_pos[i]);
            yev.push_back(C_co_e[i]);
        }
        int npts_co = (int)xv.size();

        if (npts_co < 4) {
            std::cerr << "[WARNING] Solo " << npts_co
                      << " punti corrected validi: salto il fit C_corr(x)"
                      << std::endl;
        } else {
            TCanvas *c_Cpoly_co = new TCanvas("c_C_poly_corr",
                "Fit polinomiale C(x) corrected (P&R)",
                1200, 800);
            TPad *pad_fit_co = new TPad("pad_Cpoly_co_fit", "", 0.0, 0.32, 1.0, 1.0);
            TPad *pad_res_co = new TPad("pad_Cpoly_co_res", "", 0.0, 0.00, 1.0, 0.32);
            pad_fit_co->SetBottomMargin(0.025);
            pad_fit_co->SetTopMargin(0.070);
            pad_fit_co->SetLeftMargin(0.095);
            pad_fit_co->SetRightMargin(0.035);
            pad_fit_co->SetGrid(1, 1);
            pad_res_co->SetTopMargin(0.030);
            pad_res_co->SetBottomMargin(0.300);
            pad_res_co->SetLeftMargin(0.095);
            pad_res_co->SetRightMargin(0.035);
            pad_res_co->SetGrid(1, 1);
            c_Cpoly_co->cd();
            pad_fit_co->Draw();
            pad_res_co->Draw();

            // --- Fit polinomiale ---
            pad_fit_co->cd();
            TGraphErrors *g_Cco = new TGraphErrors(
                npts_co, xv.data(), yv.data(), xev.data(), yev.data());
            g_Cco->SetName("g_C_poly_corr_points");
            g_Cco->SetMarkerStyle(22);
            g_Cco->SetMarkerSize(1.3);
            g_Cco->SetMarkerColor(kBlue + 1);
            g_Cco->SetLineColor(kBlue + 1);
            g_Cco->SetTitle(";Posizione x [cm];"
                "C^{corr} = t_{3}^{c} - #frac{t_{1}^{c}+t_{2}^{c}}{2} [ns]");
            g_Cco->GetXaxis()->SetLabelSize(0);

            TF1 *f_C_poly_co = new TF1("f_C_poly_co",
                "[0] + [1]*x + [2]*x*x", -150, 150);
            f_C_poly_co->SetParNames("p_{0}^{c}", "p_{1}^{c}", "p_{2}^{c}");
            // Parametri iniziali: p0 = media dei C_corr misurati (il termine
            // dominante), p1 = p2 = 0 (C(x) e' quasi costante). Senza questa
            // inizializzazione ROOT parte da valori di default lontanissimi
            // dal minimo e il fit non converge (residui e parametri assurdi).
            {
                double sum_co = 0.0;
                for (int i = 0; i < npts_co; i++) sum_co += yv[i];
                double mean_co = sum_co / npts_co;
                f_C_poly_co->SetParameters(mean_co, 0.0, 0.0);
            }
            g_Cco->Fit(f_C_poly_co, "S R Q");
            g_Cco->Draw("AP");

            C_co_poly_p0     = f_C_poly_co->GetParameter(0);
            C_co_poly_p0_err = f_C_poly_co->GetParError(0);
            C_co_poly_p1     = f_C_poly_co->GetParameter(1);
            C_co_poly_p1_err = f_C_poly_co->GetParError(1);
            C_co_poly_p2     = f_C_poly_co->GetParameter(2);
            C_co_poly_p2_err = f_C_poly_co->GetParError(2);
            C_co_poly_chi2ndf = (f_C_poly_co->GetNDF() > 0)
                              ? f_C_poly_co->GetChisquare() / f_C_poly_co->GetNDF()
                              : -1.0;

            // Rescaling errori se chi2/ndf > 1 (consistenza con CANVAS 3b)
            double rescale_co = 1.0;
            if (C_co_poly_chi2ndf > 1.0) {
                rescale_co = sqrt(C_co_poly_chi2ndf);
                C_co_poly_p0_err *= rescale_co;
                C_co_poly_p1_err *= rescale_co;
                C_co_poly_p2_err *= rescale_co;
            }

            f_C_poly_co->SetLineColor(kRed + 1);
            f_C_poly_co->SetLineWidth(2);
            f_C_poly_co->Draw("same");

            TPaveText *pt_Cco = new TPaveText(0.55, 0.62, 0.97, 0.92, "NDC");
            pt_Cco->SetFillColorAlpha(kWhite, 0.85);
            pt_Cco->SetTextFont(42);
            pt_Cco->SetTextSize(0.038);
            pt_Cco->SetTextAlign(12);
            pt_Cco->AddText("C^{corr}(x) = p_{0} + p_{1} x + p_{2} x^{2}  (P&R)");
            pt_Cco->AddText(Form("p_{0} = %.4f #pm %.4f ns",
                                 C_co_poly_p0, C_co_poly_p0_err));
            pt_Cco->AddText(Form("p_{1} = %.6f #pm %.6f ns/cm",
                                 C_co_poly_p1, C_co_poly_p1_err));
            pt_Cco->AddText(Form("p_{2} = %.3e #pm %.3e ns/cm^{2}",
                                 C_co_poly_p2, C_co_poly_p2_err));
            pt_Cco->AddText(Form("#chi^{2}/ndf = %.2f", C_co_poly_chi2ndf));
            // Confronto |p2| con il raw: indicatore quantitativo dell'efficacia
            pt_Cco->AddText(Form("|p_{2}^{corr}/p_{2}^{raw}| = %.2f",
                                 (C_poly_p2 != 0.0)
                                 ? fabs(C_co_poly_p2 / C_poly_p2)
                                 : 0.0));
            pt_Cco->Draw();

            // --- Pad residui ---
            pad_res_co->cd();
            std::vector<double> res_co(npts_co), pull_x(npts_co);
            for (int i = 0; i < npts_co; i++) {
                res_co[i] = yv[i] - f_C_poly_co->Eval(xv[i]);
                pull_x[i] = xv[i];
            }
            TGraphErrors *g_res_Cco = new TGraphErrors(
                npts_co, pull_x.data(), res_co.data(), xev.data(), yev.data());
            g_res_Cco->SetMarkerStyle(22);
            g_res_Cco->SetMarkerSize(1.1);
            g_res_Cco->SetMarkerColor(kBlue + 1);
            g_res_Cco->SetLineColor(kBlue + 1);
            g_res_Cco->SetTitle(";Posizione x [cm];"
                "C^{corr}_{mis} - C^{corr}_{fit} [ns]");
            g_res_Cco->GetXaxis()->SetTitleSize(0.09);
            g_res_Cco->GetXaxis()->SetLabelSize(0.08);
            g_res_Cco->GetYaxis()->SetTitleSize(0.09);
            g_res_Cco->GetYaxis()->SetTitleOffset(0.5);
            g_res_Cco->GetYaxis()->SetLabelSize(0.08);
            g_res_Cco->GetXaxis()->SetLimits(-150, 150);
            g_res_Cco->Draw("AP");
            TLine *l0_co = new TLine(-150, 0, 150, 0);
            l0_co->SetLineColor(kBlack);
            l0_co->SetLineWidth(2);
            l0_co->SetLineStyle(2);
            l0_co->Draw("same");

            // Statistica residui
            double res_mean = 0.0, res_rms = 0.0;
            for (double r : res_co) { res_mean += r; res_rms += r * r; }
            res_mean /= npts_co;
            res_rms = sqrt(res_rms / npts_co - res_mean * res_mean);
            TPaveText *pt_res = new TPaveText(0.10, 0.78, 0.45, 0.95, "NDC");
            pt_res->SetFillColorAlpha(kWhite, 0.85);
            pt_res->SetTextFont(42);
            pt_res->SetTextSize(0.07);
            pt_res->SetTextAlign(12);
            pt_res->AddText(Form("<res> = %+.4f ns", res_mean));
            pt_res->AddText(Form("RMS(res) = %.4f ns", res_rms));
            pt_res->Draw();

            c_Cpoly_co->Update();
            fout->cd();
            c_Cpoly_co->Write("Calibration_C_poly_corr_vs_x");
        }
    }
    // ====================================================================
//  CANVAS 3b: FIT POLINOMIALE C(x) = p0 + p1*x + p2*x^2  (hybrid_tot)
// ====================================================================
//  Questo canvas contiene SOLO il modello polinomiale usato per C(x)
//  e, nel pad inferiore, i residui punto-per-punto:
//
//      residuo_i = C_i - C_poly(x_i)
//
//  La vecchia divisione laterale del canvas e il confronto grafico con
//  i modelli costante/lineare sono stati rimossi. La media pesata
//  C_const viene comunque calcolata e salvata nelle variabili gia'
//  previste, cosi' da non rompere la compatibilita' con il TTree
//  fit_params o con eventuali controlli sistematici successivi.
{
    TCanvas *c_Cpoly = new TCanvas("c_C_poly",
                                   "Fit polinomiale C(x) con residui (hybrid_tot)",
                                   1200, 800);

    // Due pad verticali a tutta larghezza: fit sopra, residui sotto.
    // Non si usa c_Cpoly->Divide(), cosi' si evita il problema di cd(3)
    // e non resta nessun pannello laterale vuoto nel file ROOT.
    TPad *pad_Cpoly_fit = new TPad("pad_Cpoly_fit", "", 0.0, 0.32, 1.0, 1.0);
    TPad *pad_Cpoly_res = new TPad("pad_Cpoly_res", "", 0.0, 0.00, 1.0, 0.32);

    pad_Cpoly_fit->SetBottomMargin(0.025);
    pad_Cpoly_fit->SetTopMargin(0.070);
    pad_Cpoly_fit->SetLeftMargin(0.095);
    pad_Cpoly_fit->SetRightMargin(0.035);
    pad_Cpoly_fit->SetGrid(1, 1);

    pad_Cpoly_res->SetTopMargin(0.030);
    pad_Cpoly_res->SetBottomMargin(0.300);
    pad_Cpoly_res->SetLeftMargin(0.095);
    pad_Cpoly_res->SetRightMargin(0.035);
    pad_Cpoly_res->SetGrid(1, 1);

    c_Cpoly->cd();
    pad_Cpoly_fit->Draw();
    pad_Cpoly_res->Draw();

    // ------------------------------------------------------------------
    //  (a) FIT POLINOMIALE C(x) = [p0] + [p1]*x + [p2]*x^2
    // ------------------------------------------------------------------
    pad_Cpoly_fit->cd();

    TGraphErrors *g_Cpoly = new TGraphErrors(npts, x_pos.data(), C_hy.data(),
                                              dx_pos.data(), C_hy_e.data());
    g_Cpoly->SetName("g_C_poly_points");
    g_Cpoly->SetMarkerStyle(22);
    g_Cpoly->SetMarkerSize(1.3);
    g_Cpoly->SetMarkerColor(kGreen + 2);
    g_Cpoly->SetLineColor(kGreen + 2);
    g_Cpoly->SetTitle(";Posizione x [cm];C = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]");
    g_Cpoly->GetXaxis()->SetLimits(-150, 150);
    g_Cpoly->GetXaxis()->SetLabelSize(0.0);
    g_Cpoly->GetXaxis()->SetTitleSize(0.0);
    g_Cpoly->GetYaxis()->SetTitleSize(0.055);
    g_Cpoly->GetYaxis()->SetLabelSize(0.047);
    g_Cpoly->GetYaxis()->SetTitleOffset(0.78);

    TF1 *f_C_poly = new TF1("f_C_poly", "[0] + [1]*x + [2]*x*x", -150, 150);
    f_C_poly->SetParNames("p_{0}", "p_{1}", "p_{2}");
    // Inizializzazione: p0 ~ media C, p1 ~ pendenza grezza, p2 ~ curvatura piccola.
    f_C_poly->SetParameters(-12.2, 1e-3, -1e-5);
    g_Cpoly->Fit(f_C_poly, "S R Q");

    C_poly_p0      = f_C_poly->GetParameter(0);
    C_poly_p0_err  = f_C_poly->GetParError(0);
    C_poly_p1      = f_C_poly->GetParameter(1);
    C_poly_p1_err  = f_C_poly->GetParError(1);
    C_poly_p2      = f_C_poly->GetParameter(2);
    C_poly_p2_err  = f_C_poly->GetParError(2);
    C_poly_chi2ndf = (f_C_poly->GetNDF() > 0)
                   ? f_C_poly->GetChisquare() / f_C_poly->GetNDF() : -1.0;

    // Consistency rescaling se chi2/ndf > 1.
    // Nota: i parametri centrali restano quelli del fit; vengono gonfiati
    // solo gli errori, come nel codice precedente.
    double rescale_poly = 1.0;
    if (C_poly_chi2ndf > 1.0) {
        rescale_poly = sqrt(C_poly_chi2ndf);
        C_poly_p0_err *= rescale_poly;
        C_poly_p1_err *= rescale_poly;
        C_poly_p2_err *= rescale_poly;
        std::cout << "[INFO] C(x) poly fit: chi2/ndf = " << C_poly_chi2ndf
                  << " > 1 -> errori parametri rinormalizzati per sqrt("
                  << C_poly_chi2ndf << ") = " << rescale_poly << std::endl;
    }

    g_Cpoly->Draw("AP");
    f_C_poly->SetLineColor(kRed + 1);
    f_C_poly->SetLineWidth(2);
    f_C_poly->Draw("same");

    // Banda +/-1 sigma del fit polinomiale.
    // Si mantiene la stessa approssimazione diagonale del codice originale,
    // usando gli errori eventualmente gia' riscalati.
    {
        const int nband2 = 200;
        std::vector<double> xb2(nband2), yb2_up(nband2), yb2_dn(nband2);
        for (int i = 0; i < nband2; i++) {
            xb2[i] = -150.0 + 300.0 * i / (nband2 - 1);
            double xx = xb2[i];
            double c_val = C_poly_p0 + C_poly_p1 * xx + C_poly_p2 * xx * xx;
            double c_err = sqrt(C_poly_p0_err * C_poly_p0_err
                              + xx * xx * C_poly_p1_err * C_poly_p1_err
                              + xx * xx * xx * xx * C_poly_p2_err * C_poly_p2_err);
            yb2_up[i] = c_val + c_err;
            yb2_dn[i] = c_val - c_err;
        }
        TGraph *g_band2_up = new TGraph(nband2, xb2.data(), yb2_up.data());
        TGraph *g_band2_dn = new TGraph(nband2, xb2.data(), yb2_dn.data());
        g_band2_up->SetName("g_C_poly_band_up");
        g_band2_dn->SetName("g_C_poly_band_down");
        g_band2_up->SetLineColor(kOrange - 3);
        g_band2_dn->SetLineColor(kOrange - 3);
        g_band2_up->SetLineStyle(7);
        g_band2_dn->SetLineStyle(7);
        g_band2_up->Draw("L same");
        g_band2_dn->Draw("L same");
        g_Cpoly->Draw("P same");
    }

    TPaveText *pt_Cp = new TPaveText(0.135, 0.150, 0.640, 0.455, "NDC");
    pt_Cp->SetFillColorAlpha(kWhite, 0.85);
    pt_Cp->SetTextFont(42);
    pt_Cp->SetTextSize(0.040);
    pt_Cp->SetTextAlign(12);
    pt_Cp->AddText("C(x) = p_{0} + p_{1} #upoint x + p_{2} #upoint x^{2}  (hybrid_tot)");
    pt_Cp->AddText(Form("p_{0} = %.4f #pm %.4f ns", C_poly_p0, C_poly_p0_err));
    pt_Cp->AddText(Form("p_{1} = %.6f #pm %.6f ns/cm", C_poly_p1, C_poly_p1_err));
    pt_Cp->AddText(Form("p_{2} = %.3e #pm %.3e ns/cm^{2}", C_poly_p2, C_poly_p2_err));
    pt_Cp->AddText(Form("#chi^{2}/ndf = %.2f", C_poly_chi2ndf));
    pt_Cp->Draw();

    // ------------------------------------------------------------------
    //  (b) RESIDUI DEL FIT POLINOMIALE
    // ------------------------------------------------------------------
    pad_Cpoly_res->cd();

    std::vector<double> res2_x(npts), res2_y(npts), res2_ex(npts), res2_ey(npts);
    double res2_max = 0.0;
    double res2_sum = 0.0;
    double res2_sum2 = 0.0;
    int    res2_n = 0;

    for (int i = 0; i < npts; i++) {
        double C_fit_i = C_poly_p0 + C_poly_p1 * x_pos[i]
                       + C_poly_p2 * x_pos[i] * x_pos[i];
        double res_i = C_hy[i] - C_fit_i;

        res2_x[i]  = x_pos[i];
        res2_ex[i] = dx_pos[i];
        res2_y[i]  = res_i;
        res2_ey[i] = C_hy_e[i];

        double half_height_i = fabs(res_i) + ((C_hy_e[i] > 0.0) ? C_hy_e[i] : 0.0);
        if (half_height_i > res2_max) res2_max = half_height_i;

        res2_sum  += res_i;
        res2_sum2 += res_i * res_i;
        res2_n++;
    }

    double res2_mean = (res2_n > 0) ? res2_sum / res2_n : 0.0;
    double res2_rms  = 0.0;
    if (res2_n > 0) {
        double var = res2_sum2 / res2_n - res2_mean * res2_mean;
        res2_rms = (var > 0.0) ? sqrt(var) : 0.0;
    }

    TGraphErrors *g_res2 = new TGraphErrors(npts, res2_x.data(), res2_y.data(),
                                             res2_ex.data(), res2_ey.data());
    g_res2->SetName("g_C_poly_residuals");
    g_res2->SetMarkerStyle(22);
    g_res2->SetMarkerSize(1.05);
    g_res2->SetMarkerColor(kGreen + 2);
    g_res2->SetLineColor(kGreen + 2);
    g_res2->SetTitle(";Posizione x [cm];C_{mis} - C_{fit} [ns]");
    g_res2->GetXaxis()->SetLimits(-150, 150);
    g_res2->GetXaxis()->SetTitleSize(0.095);
    g_res2->GetXaxis()->SetLabelSize(0.082);
    g_res2->GetYaxis()->SetTitleSize(0.086);
    g_res2->GetYaxis()->SetTitleOffset(0.48);
    g_res2->GetYaxis()->SetLabelSize(0.075);
    g_res2->GetYaxis()->SetNdivisions(505);

    double y_range2 = std::max(0.10, 1.25 * res2_max);
    g_res2->GetYaxis()->SetRangeUser(-y_range2, y_range2);
    g_res2->Draw("AP");

    TLine *l0p = new TLine(-150, 0, 150, 0);
    l0p->SetLineColor(kBlack);
    l0p->SetLineWidth(2);
    l0p->SetLineStyle(2);
    l0p->Draw("same");

    TPaveText *pt_res = new TPaveText(0.135, 0.740, 0.410, 0.925, "NDC");
    pt_res->SetFillColorAlpha(kWhite, 0.85);
    pt_res->SetTextFont(42);
    pt_res->SetTextSize(0.070);
    pt_res->SetTextAlign(12);
    pt_res->AddText(Form("#LTres#GT = %.4f ns", res2_mean));
    pt_res->AddText(Form("RMS(res) = %.4f ns", res2_rms));
    pt_res->Draw();

    // ------------------------------------------------------------------
    //  (c) MEDIA PESATA C_const = <C>  [calcolata ma non disegnata]
    // ------------------------------------------------------------------
    //  Manteniamo il calcolo per compatibilita' con il salvataggio dei
    //  parametri e con l'analisi sistematica, ma non lo mostriamo piu'
    //  nel canvas richiesto.
    double sum_w = 0.0, sum_wC = 0.0;
    for (int i = 0; i < npts; i++) {
        if (C_hy_e[i] <= 0.0) continue;
        double w_i = 1.0 / (C_hy_e[i] * C_hy_e[i]);
        sum_w  += w_i;
        sum_wC += w_i * C_hy[i];
    }

    if (sum_w > 0.0) {
        C_const_val = sum_wC / sum_w;
        C_const_err = 1.0 / sqrt(sum_w);
    } else {
        C_const_val = 0.0;
        C_const_err = 0.0;
    }

    double C_const_chi2 = 0.0;
    for (int i = 0; i < npts; i++) {
        if (C_hy_e[i] <= 0.0) continue;
        double dC = C_hy[i] - C_const_val;
        C_const_chi2 += (dC * dC) / (C_hy_e[i] * C_hy_e[i]);
    }
    C_const_chi2ndf = (npts > 1) ? C_const_chi2 / (npts - 1) : -1.0;

    if (C_const_chi2ndf > 1.0) {
        C_const_err *= sqrt(C_const_chi2ndf);
        std::cout << "[INFO] C_const: chi2/ndf = " << C_const_chi2ndf
                  << " > 1 -> errore media rinormalizzato per sqrt("
                  << C_const_chi2ndf << ")" << std::endl;
    }

    c_Cpoly->cd();
    c_Cpoly->Update();
    fout->cd();
    c_Cpoly->Write("Calibration_C_poly_and_const");

    std::cout << "\n[INFO] Modello polinomiale C(x):" << std::endl;
    std::cout << "       C(x) = " << Form("%.4f", C_poly_p0)
              << " + (" << Form("%.6f", C_poly_p1) << ")*x + ("
              << Form("%.3e", C_poly_p2) << ")*x^2  ns" << std::endl;
    std::cout << "       chi2/ndf = " << Form("%.2f", C_poly_chi2ndf) << std::endl;
    std::cout << "       <res> = " << Form("%.4f", res2_mean)
              << " ns, RMS(res) = " << Form("%.4f", res2_rms) << " ns" << std::endl;
    std::cout << "[INFO] Modello costante C calcolato ma non disegnato:" << std::endl;
    std::cout << "       <C> = " << Form("%.4f +/- %.4f", C_const_val, C_const_err)
              << " ns   (chi2/ndf = " << Form("%.1f", C_const_chi2ndf)
              << ")" << std::endl;
}

    // ====================================================================
    //  CANVAS 4: CONFRONTO Dt12(x) — noclip vs hybrid_tot vs pure_tot
    // ====================================================================
    {
        TCanvas *c_cmp = new TCanvas("c_dt12_3methods",
                                     "#Delta t_{12} vs x: confronto 3 metodi",
                                     1000, 700);
        c_cmp->SetGrid(1, 1);

        // noclip (blu, cerchi)
        TGraphErrors *g_nc = new TGraphErrors(npts, x_pos.data(), dt12_nc.data(),
                                              dx_pos.data(), dt12_nc_e.data());
        g_nc->SetMarkerStyle(20); g_nc->SetMarkerSize(1.1);
        g_nc->SetMarkerColor(kBlue + 1); g_nc->SetLineColor(kBlue + 1);
        g_nc->SetTitle("Confronto #Delta t_{12}(x) tra metodi;"
                       "Posizione x [cm];#Delta t_{12} = t_{1} - t_{2} [ns]");
        TF1 *f_nc = new TF1("f_dt12_3_nc", "[0]+[1]*x", -150, 150);
        g_nc->Fit(f_nc, "S Q R N");
        g_nc->Draw("AP");
        f_nc->SetLineColor(kBlue + 1); f_nc->SetLineWidth(2); f_nc->Draw("same");

        // hybrid_tot (rosso, quadrati)
        TGraphErrors *g_hy = new TGraphErrors(npts, x_pos.data(), dt12_hy.data(),
                                              dx_pos.data(), dt12_hy_e.data());
        g_hy->SetMarkerStyle(21); g_hy->SetMarkerSize(1.1);
        g_hy->SetMarkerColor(kRed + 1); g_hy->SetLineColor(kRed + 1);
        TF1 *f_hy = new TF1("f_dt12_3_hy", "[0]+[1]*x", -150, 150);
        g_hy->Fit(f_hy, "S Q R N");
        g_hy->Draw("P same");
        f_hy->SetLineColor(kRed + 1); f_hy->SetLineWidth(2);
        f_hy->SetLineStyle(2); f_hy->Draw("same");

        // pure_tot (viola, diamanti)
        TGraphErrors *g_pu = new TGraphErrors(npts, x_pos.data(), dt12_pu.data(),
                                              dx_pos.data(), dt12_pu_e.data());
        g_pu->SetMarkerStyle(33); g_pu->SetMarkerSize(1.5);
        g_pu->SetMarkerColor(kViolet + 1); g_pu->SetLineColor(kViolet + 1);
        TF1 *f_pu = new TF1("f_dt12_3_pu", "[0]+[1]*x", -150, 150);
        g_pu->Fit(f_pu, "S Q R N");
        g_pu->Draw("P same");
        f_pu->SetLineColor(kViolet + 1); f_pu->SetLineWidth(2);
        f_pu->SetLineStyle(9); f_pu->Draw("same");

        // Box riepilogativo con le 3 velocita' derivate.
        auto vEff = [](TF1* f) { return 2.0 / fabs(f->GetParameter(1)); };
        auto vErr = [](TF1* f) {
            double m = f->GetParameter(1);
            return 2.0 * f->GetParError(1) / (m * m);
        };
        TPaveText *pt = new TPaveText(0.13, 0.64, 0.57, 0.90, "NDC");
        pt->SetFillColorAlpha(kWhite, 0.85);
        pt->SetTextFont(42); pt->SetTextSize(0.032);
        pt->SetTextAlign(12);
        pt->AddText("v_{eff} = 2/|m| [cm/ns]:");
        pt->AddText(Form("#color[4]{noclip:     %.2f #pm %.2f}", vEff(f_nc), vErr(f_nc)));
        pt->AddText(Form("#color[2]{hybrid_tot: %.2f #pm %.2f}", vEff(f_hy), vErr(f_hy)));
        pt->AddText(Form("#color[%d]{pure_tot:   %.2f #pm %.2f}",
                         kViolet + 1, vEff(f_pu), vErr(f_pu)));
        pt->Draw();

        TLegend *leg = new TLegend(0.62, 0.74, 0.90, 0.92);
        leg->SetTextFont(42); leg->SetTextSize(0.035);
        leg->AddEntry(g_nc, "noclip",     "p");
        leg->AddEntry(g_hy, "hybrid_tot", "p");
        leg->AddEntry(g_pu, "pure_tot",   "p");
        leg->Draw();

        c_cmp->Update();
        fout->cd();
        c_cmp->Write("Compare_Dt12_3_methods");
    }

    // ====================================================================
    //  CANVAS 5: CONFRONTO C(x) — noclip vs hybrid_tot vs pure_tot
    // ====================================================================
    // Confronto cruciale: se hybrid_tot ha una pendenza C1 != 0 ma pure_tot
    // (ricostruzione simmetrica) la attenua, l'effetto e' un artefatto del
    // trattamento asimmetrico clip/non-clip e non fisica della barra.
    {
        TCanvas *c_C3 = new TCanvas("c_C_3methods",
                                    "C(x): confronto 3 metodi", 1000, 700);
        c_C3->SetGrid(1, 1);

 // Lambda che esegue il fit lineare con consistency rescaling degli
        // errori dei parametri (stessa logica del fit principale del Canvas 3).
        // Ritorna (slope, slope_err) gia' rescalato, e l'intercetta rescalata
        // resta nel TF1 (ci serve per disegnare la retta col color del metodo).
        auto fitC = [](TGraphErrors* g, TF1* f) -> std::pair<double,double> {
            g->Fit(f, "S R Q");
            double chi2ndf = (f->GetNDF() > 0)
                           ? f->GetChisquare() / f->GetNDF() : -1.0;
            double rs = (chi2ndf > 1.0) ? sqrt(chi2ndf) : 1.0;
            return { f->GetParameter(1), f->GetParError(1) * rs };
        };

        TGraphErrors *g_nc = new TGraphErrors(npts, x_pos.data(), C_nc.data(),
                                              dx_pos.data(), C_nc_e.data());
        g_nc->SetMarkerStyle(20); g_nc->SetMarkerSize(1.1);
        g_nc->SetMarkerColor(kBlue + 1); g_nc->SetLineColor(kBlue + 1);
        g_nc->SetTitle("Confronto C(x) tra metodi;"
                       "Posizione x [cm];C = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]");
        TF1 *f_nc = new TF1("f_C_3_nc", "[0]+[1]*x", -150, 150);
        auto res_nc = fitC(g_nc, f_nc);
        g_nc->Draw("AP");
        f_nc->SetLineColor(kBlue + 1); f_nc->SetLineWidth(2); f_nc->Draw("same");

        TGraphErrors *g_hy = new TGraphErrors(npts, x_pos.data(), C_hy.data(),
                                              dx_pos.data(), C_hy_e.data());
        g_hy->SetMarkerStyle(21); g_hy->SetMarkerSize(1.1);
        g_hy->SetMarkerColor(kRed + 1); g_hy->SetLineColor(kRed + 1);
        TF1 *f_hy = new TF1("f_C_3_hy", "[0]+[1]*x", -150, 150);
        auto res_hy = fitC(g_hy, f_hy);
        g_hy->Draw("P same");
        f_hy->SetLineColor(kRed + 1); f_hy->SetLineWidth(2);
        f_hy->SetLineStyle(2); f_hy->Draw("same");

        TGraphErrors *g_pu = new TGraphErrors(npts, x_pos.data(), C_pu.data(),
                                              dx_pos.data(), C_pu_e.data());
        g_pu->SetMarkerStyle(33); g_pu->SetMarkerSize(1.5);
        g_pu->SetMarkerColor(kViolet + 1); g_pu->SetLineColor(kViolet + 1);
        TF1 *f_pu = new TF1("f_C_3_pu", "[0]+[1]*x", -150, 150);
        auto res_pu = fitC(g_pu, f_pu);
        g_pu->Draw("P same");
        f_pu->SetLineColor(kViolet + 1); f_pu->SetLineWidth(2);
        f_pu->SetLineStyle(9); f_pu->Draw("same");

        TPaveText *pt = new TPaveText(0.13, 0.64, 0.57, 0.90, "NDC");
        pt->SetFillColorAlpha(kWhite, 0.85);
        pt->SetTextFont(42); pt->SetTextSize(0.032);
        pt->SetTextAlign(12);
        pt->AddText("Pendenza lineare C_{1} [ns/cm] (errori rescalati):");
        pt->AddText(Form("#color[4]{noclip:     %.6f #pm %.6f}",
                         res_nc.first, res_nc.second));
        pt->AddText(Form("#color[2]{hybrid_tot: %.6f #pm %.6f}",
                         res_hy.first, res_hy.second));
        pt->AddText(Form("#color[%d]{pure_tot:   %.6f #pm %.6f}",
                         kViolet + 1, res_pu.first, res_pu.second));
        pt->Draw();

        TLegend *leg = new TLegend(0.62, 0.74, 0.90, 0.92);
        leg->SetTextFont(42); leg->SetTextSize(0.035);
        leg->AddEntry(g_nc, "noclip",     "p");
        leg->AddEntry(g_hy, "hybrid_tot", "p");
        leg->AddEntry(g_pu, "pure_tot",   "p");
        leg->Draw();

        c_C3->Update();
        fout->cd();
        c_C3->Write("Compare_C_3_methods");
    }

    // ====================================================================
    //  CANVAS 6: DIAGNOSTICA dA = A_rec(TOT) - A_vera  (PMT1, PMT2, PMT3)
    // ====================================================================
    // Riga 1: istogrammi 1D di dA (bias e larghezza del recupero TOT).
    // Riga 2: scatter 2D dA vs A_vera con profilo (dipendenza dall'ampiezza).
    // Un dA centrato su 0 e senza pendenza in A indica un modello A(TOT) ben
    // calibrato; una pendenza segnala bias residuo dipendente dall'ampiezza.
    {
        TCanvas *c_dA = new TCanvas("c_diagnostic_dA",
                                    "Diagnostica #DeltaA = A_{rec}(TOT) - A_{vera}",
                                    1500, 800);
        c_dA->Divide(3, 2);
        int colors_d[3] = {kBlue + 1, kRed + 1, kGreen + 2};

        for (int k = 0; k < 3; k++) {
            // --- Pad superiore: istogramma 1D di dA ---
            c_dA->cd(k + 1);
            gPad->SetGrid(1, 1);
            if (g_h_dA[k] && g_h_dA[k]->GetEntries() > 10) {
                g_h_dA[k]->SetLineColor(colors_d[k]);
                g_h_dA[k]->SetFillColorAlpha(colors_d[k] - 9, 0.3);
                g_h_dA[k]->Draw();

                TPaveText *ptd = new TPaveText(0.60, 0.66, 0.92, 0.92, "NDC");
                ptd->SetFillColorAlpha(kWhite, 0.85);
                ptd->SetTextFont(42); ptd->SetTextSize(0.045);
                ptd->SetTextAlign(12);
                ptd->AddText(Form("%s", pmt_names[k]));
                ptd->AddText(Form("media = %.2f mV", g_h_dA[k]->GetMean()));
                ptd->AddText(Form("RMS = %.2f mV", g_h_dA[k]->GetRMS()));
                ptd->AddText(Form("N = %.0f", g_h_dA[k]->GetEntries()));
                ptd->Draw();
            }

            // --- Pad inferiore: scatter 2D dA vs A_vera + profilo ---
            c_dA->cd(k + 4);
            gPad->SetGrid(1, 1);
            if (g_h_dA_vs_A[k] && g_h_dA_vs_A[k]->GetEntries() > 10) {
                g_h_dA_vs_A[k]->SetStats(0);
                g_h_dA_vs_A[k]->Draw("COLZ");
                TProfile *p_prof = g_h_dA_vs_A[k]->ProfileX(
                    Form("p_dA_vs_A_%s", pmt_names[k]));
                p_prof->SetLineColor(kBlack);
                p_prof->SetLineWidth(2);
                p_prof->SetMarkerStyle(20);
                p_prof->SetMarkerSize(0.8);
                p_prof->Draw("same");

                TLine *l_z = new TLine(DA2_A_LO, 0, DA2_A_HI, 0);
                l_z->SetLineColor(kRed); l_z->SetLineStyle(2);
                l_z->SetLineWidth(2); l_z->Draw("same");
            }
        }

        c_dA->Update();
        fout->cd();
        c_dA->Write("Diagnostic_dA_TOT");

        // Scrivi anche gli istogrammi singoli (utili per analisi successive).
        for (int k = 0; k < 3; k++) {
            if (g_h_dA[k])        g_h_dA[k]->Write();
            if (g_h_dA_vs_A[k])   g_h_dA_vs_A[k]->Write();
            if (g_h_dA_vs_TOT[k]) g_h_dA_vs_TOT[k]->Write();
        }
    }

    // ====================================================================
    //  CANVAS 7-8: GALLERIE ISTOGRAMMI Dt12 e C (metodo hybrid_tot)
    // ====================================================================
    // Una griglia di istogrammi, uno per posizione: controllo visivo della
    // forma delle distribuzioni e della qualita' dei fit.
    {
        int ncols = (npts <= 4) ? npts : (npts <= 9 ? 3 : 4);
        int nrows = (npts + ncols - 1) / ncols;

        // --- Galleria Dt12 (hybrid_tot) ---
        TCanvas *c_gal_dt12 = new TCanvas("c_dt12_gallery",
                                          "Distribuzioni #Delta t_{12} (hybrid_tot)",
                                          350 * ncols, 300 * nrows);
        c_gal_dt12->Divide(ncols, nrows);
        for (int i = 0; i < npts; i++) {
            c_gal_dt12->cd(i + 1);
            gPad->SetGrid(1, 1);
            TH1D *h = (TH1D*)fout->Get(Form("h_dt12_hybrid_%s", labels[i].c_str()));
            if (h) {
                h->SetLineColor(kRed + 1);
                h->SetFillColorAlpha(kRed - 9, 0.3);
                h->Draw();
            }
        }
        c_gal_dt12->Update();
        fout->cd();
        c_gal_dt12->Write("Gallery_Dt12_hybrid_tot");

        // --- Galleria C (hybrid_tot) ---
        TCanvas *c_gal_C = new TCanvas("c_C_gallery",
                                       "Distribuzioni C (hybrid_tot)",
                                       350 * ncols, 300 * nrows);
        c_gal_C->Divide(ncols, nrows);
        for (int i = 0; i < npts; i++) {
            c_gal_C->cd(i + 1);
            gPad->SetGrid(1, 1);
            TH1D *h = (TH1D*)fout->Get(Form("h_C_hybrid_%s", labels[i].c_str()));
            if (h) {
                h->SetLineColor(kGreen + 2);
                h->SetFillColorAlpha(kGreen - 9, 0.3);
                h->Draw();
            }
        }
        c_gal_C->Update();
        fout->cd();
        c_gal_C->Write("Gallery_C_hybrid_tot");
    }
// ====================================================================
    //  CANVAS 9: DIAGNOSTICA Dt13(x) e Dt23(x) — equalizzazione PMT1/PMT2
    // ====================================================================
    //  In un apparato con PMT1 a x = -L/2 e PMT2 a x = +L/2, i tempi di
    //  arrivo della luce nei due PMT dipendono linearmente dalla posizione
    //  di impatto x:
    //
    //     t1 = t0 + (L/2 + x)/v + delta_1
    //     t2 = t0 + (L/2 - x)/v + delta_2
    //     t3 = t0              + delta_3    (Guida A: contatto con la barra)
    //
    //  Di conseguenza:
    //     Dt13 = t1 - t3 = (L/2 + x)/v + (delta_1 - delta_3)
    //                     = x/v + const_1    -> pendenza attesa = +1/v = +|m|/2
    //     Dt23 = t2 - t3 = (L/2 - x)/v + (delta_2 - delta_3)
    //                     = -x/v + const_2   -> pendenza attesa = -1/v = -|m|/2
    //
    //  dove m = 2/v_eff e' la pendenza della retta Dt12(x) gia' fittata.
    //
    //  PANNELLI:
    //   [1] In alto a sinistra: Dt13(x) dati + fit lineare
    //   [2] In alto a destra:   Dt23(x) dati + fit lineare
    //   [3] In basso a sinistra:  residui Dt13 dal fit lineare
    //   [4] In basso a destra:    residui Dt23 dal fit lineare
    //
    //  La legenda riporta: pendenza misurata, pendenza attesa |m|/2,
    //  deviazione in sigma, e il rapporto m13/m23 (atteso = -1).
    //
    //  Un canvas aggiuntivo mostra le risoluzioni sigma(Dt13) e sigma(Dt23)
    //  in funzione di x, sovrapposte: la risoluzione dovrebbe essere simmetrica
    //  attorno al centro barra se i PMT sono equalizzati.
    {
        TCanvas *c_dt13_23 = new TCanvas("c_dt13_dt23",
            "Diagnostica #Delta t_{13}(x) e #Delta t_{23}(x)", 1200, 800);
        c_dt13_23->Divide(2, 2);

        // Pendenza attesa dal fit Dt12: m_fit = 2/v, quindi 1/v = m_fit/2.
        // Dt13: pendenza attesa = +|m_fit|/2
        // Dt23: pendenza attesa = -|m_fit|/2
        double m_expected_13 = +fabs(m_fit) / 2.0;
        double m_expected_23 = -fabs(m_fit) / 2.0;
        // Errore propagato: sigma(m/2) = sigma_m / 2
        double m_expected_err = m_err / 2.0;

        // ---- (1) Dt13(x): fit lineare ----
        c_dt13_23->cd(1);
        gPad->SetPad(0.0, 0.4, 0.5, 1.0);
        gPad->SetBottomMargin(0.02);
        gPad->SetGrid(1, 1);

        TGraphErrors *g_dt13 = new TGraphErrors(npts, x_pos.data(), dt13_hy.data(),
                                                 dx_pos.data(), dt13_hy_e.data());
        g_dt13->SetMarkerStyle(20); g_dt13->SetMarkerSize(1.0);
        g_dt13->SetMarkerColor(kBlue + 1); g_dt13->SetLineColor(kBlue + 1);
        g_dt13->SetTitle(";Posizione x [cm];#Delta t_{13} = t_{1} - t_{3} [ns]");
        g_dt13->GetXaxis()->SetLabelSize(0);

        TF1 *f_dt13 = new TF1("f_dt13_lin", "[0] + [1]*x", -150, 150);
        f_dt13->SetParNames("q_{13}", "m_{13}");
        {
            double q13_init = 0.0;
            for (int i = 0; i < npts; i++)
                q13_init += dt13_hy[i] - m_expected_13 * x_pos[i];
            q13_init /= (npts > 0 ? npts : 1);
            f_dt13->SetParameters(q13_init, m_expected_13);
        }
        g_dt13->Fit(f_dt13, "S R Q");
        g_dt13->Draw("AP");
        f_dt13->SetLineColor(kRed + 1); f_dt13->SetLineWidth(2);
        f_dt13->Draw("same");

        double m13     = f_dt13->GetParameter(1);
        double m13_err = f_dt13->GetParError(1);
        double q13     = f_dt13->GetParameter(0);
        double q13_err = f_dt13->GetParError(0);
        double chi2_13 = (f_dt13->GetNDF() > 0) ?
                          f_dt13->GetChisquare() / f_dt13->GetNDF() : -1.0;
        // Deviazione dalla pendenza attesa in unita' di sigma combinata
        double sigma_comb_13 = sqrt(m13_err * m13_err
                                  + m_expected_err * m_expected_err);
        double dev_sigma_13  = (sigma_comb_13 > 0)
                             ? (m13 - m_expected_13) / sigma_comb_13 : 0.0;

        TPaveText *pt13 = new TPaveText(0.14, 0.14, 0.62, 0.48, "NDC");
        pt13->SetFillColorAlpha(kWhite, 0.85);
        pt13->SetTextFont(42); pt13->SetTextSize(0.038);
        pt13->SetTextAlign(12);
        pt13->AddText(Form("#Delta t_{13} = q_{13} + m_{13} #upoint x"));
        pt13->AddText(Form("m_{13} = %.5f #pm %.5f ns/cm", m13, m13_err));
        pt13->AddText(Form("m_{att} = +|m|/2 = %.5f #pm %.5f", m_expected_13, m_expected_err));
        pt13->AddText(Form("Deviazione: %.2f #sigma", dev_sigma_13));
        pt13->AddText(Form("#chi^{2}/ndf = %.2f", chi2_13));
        pt13->Draw();

        // ---- (2) Dt23(x): fit lineare ----
        c_dt13_23->cd(2);
        gPad->SetPad(0.5, 0.4, 1.0, 1.0);
        gPad->SetBottomMargin(0.02);
        gPad->SetGrid(1, 1);

        TGraphErrors *g_dt23 = new TGraphErrors(npts, x_pos.data(), dt23_hy.data(),
                                                 dx_pos.data(), dt23_hy_e.data());
        g_dt23->SetMarkerStyle(21); g_dt23->SetMarkerSize(1.0);
        g_dt23->SetMarkerColor(kGreen + 2); g_dt23->SetLineColor(kGreen + 2);
        g_dt23->SetTitle(";Posizione x [cm];#Delta t_{23} = t_{2} - t_{3} [ns]");
        g_dt23->GetXaxis()->SetLabelSize(0);

        TF1 *f_dt23 = new TF1("f_dt23_lin", "[0] + [1]*x", -150, 150);
        f_dt23->SetParNames("q_{23}", "m_{23}");
        // Stima dei parametri iniziali dalla fisica nota e dai dati stessi.
        // La pendenza attesa e' -|m_fit|/2; l'intercetta si stima come
        // media di (dt23[i] - m_att * x[i]), cioe' il valore che il fit
        // dovrebbe dare a x=0 se il modello lineare fosse perfetto.
        // Senza questa inizializzazione ROOT parte da (1,1) e il minimizzatore
        // finisce in una valle sbagliata, restituendo parametri assurdi.
        {
            double q23_init = 0.0;
            for (int i = 0; i < npts; i++)
                q23_init += dt23_hy[i] - m_expected_23 * x_pos[i];
            q23_init /= (npts > 0 ? npts : 1);
            f_dt23->SetParameters(q23_init, m_expected_23);
        }
        g_dt23->Fit(f_dt23, "S R Q");
        g_dt23->Draw("AP");
        f_dt23->SetLineColor(kRed + 1); f_dt23->SetLineWidth(2);
        f_dt23->Draw("same");

        double m23     = f_dt23->GetParameter(1);
        double m23_err = f_dt23->GetParError(1);
        double q23     = f_dt23->GetParameter(0);
        double q23_err = f_dt23->GetParError(0);
        double chi2_23 = (f_dt23->GetNDF() > 0) ?
                          f_dt23->GetChisquare() / f_dt23->GetNDF() : -1.0;
        double sigma_comb_23 = sqrt(m23_err * m23_err
                                  + m_expected_err * m_expected_err);
        double dev_sigma_23  = (sigma_comb_23 > 0)
                             ? (m23 - m_expected_23) / sigma_comb_23 : 0.0;

        TPaveText *pt23 = new TPaveText(0.42, 0.14, 0.90, 0.48, "NDC");
        pt23->SetFillColorAlpha(kWhite, 0.85);
        pt23->SetTextFont(42); pt23->SetTextSize(0.038);
        pt23->SetTextAlign(12);
        pt23->AddText(Form("#Delta t_{23} = q_{23} + m_{23} #upoint x"));
        pt23->AddText(Form("m_{23} = %.5f #pm %.5f ns/cm", m23, m23_err));
        pt23->AddText(Form("m_{att} = -|m|/2 = %.5f #pm %.5f", m_expected_23, m_expected_err));
        pt23->AddText(Form("Deviazione: %.2f #sigma", dev_sigma_23));
        pt23->AddText(Form("#chi^{2}/ndf = %.2f", chi2_23));
        pt23->Draw();

        // ---- (3) Residui Dt13 ----
        c_dt13_23->cd(3);
        gPad->SetPad(0.0, 0.0, 0.5, 0.4);
        gPad->SetTopMargin(0.02);
        gPad->SetBottomMargin(0.22);
        gPad->SetGrid(1, 1);

        std::vector<double> res13(npts);
        for (int i = 0; i < npts; i++)
            res13[i] = dt13_hy[i] - f_dt13->Eval(x_pos[i]);

        TGraphErrors *g_res13 = new TGraphErrors(npts, x_pos.data(), res13.data(),
                                                  dx_pos.data(), dt13_hy_e.data());
        g_res13->SetMarkerStyle(20); g_res13->SetMarkerSize(1.0);
        g_res13->SetMarkerColor(kBlue + 1); g_res13->SetLineColor(kBlue + 1);
        g_res13->SetTitle(";Posizione x [cm];Residui #Delta t_{13} [ns]");
        g_res13->GetXaxis()->SetTitleSize(0.07);
        g_res13->GetXaxis()->SetLabelSize(0.06);
        g_res13->GetYaxis()->SetTitleSize(0.07);
        g_res13->GetYaxis()->SetLabelSize(0.06);
        g_res13->GetYaxis()->SetTitleOffset(0.6);
        g_res13->Draw("AP");
        TLine *lz13 = new TLine(-150, 0, 150, 0);
        lz13->SetLineColor(kRed); lz13->SetLineStyle(2); lz13->Draw("same");

        // ---- (4) Residui Dt23 ----
        c_dt13_23->cd(4);
        gPad->SetPad(0.5, 0.0, 1.0, 0.4);
        gPad->SetTopMargin(0.02);
        gPad->SetBottomMargin(0.22);
        gPad->SetGrid(1, 1);

        std::vector<double> res23(npts);
        for (int i = 0; i < npts; i++)
            res23[i] = dt23_hy[i] - f_dt23->Eval(x_pos[i]);

        TGraphErrors *g_res23 = new TGraphErrors(npts, x_pos.data(), res23.data(),
                                                  dx_pos.data(), dt23_hy_e.data());
        g_res23->SetMarkerStyle(21); g_res23->SetMarkerSize(1.0);
        g_res23->SetMarkerColor(kGreen + 2); g_res23->SetLineColor(kGreen + 2);
        g_res23->SetTitle(";Posizione x [cm];Residui #Delta t_{23} [ns]");
        g_res23->GetXaxis()->SetTitleSize(0.07);
        g_res23->GetXaxis()->SetLabelSize(0.06);
        g_res23->GetYaxis()->SetTitleSize(0.07);
        g_res23->GetYaxis()->SetLabelSize(0.06);
        g_res23->GetYaxis()->SetTitleOffset(0.6);
        g_res23->Draw("AP");
        TLine *lz23 = new TLine(-150, 0, 150, 0);
        lz23->SetLineColor(kRed); lz23->SetLineStyle(2); lz23->Draw("same");

        c_dt13_23->Update();
        fout->cd();
        c_dt13_23->Write("Diagnostic_Dt13_Dt23_vs_x");

        // ---- Stampa a schermo del riepilogo equalizzazione ----
        std::cout << "\n[INFO] Diagnostica equalizzazione PMT1/PMT2:" << std::endl;
        std::cout << "       Dt13(x): m13 = " << Form("%.5f +/- %.5f", m13, m13_err)
                  << "  atteso = +" << Form("%.5f", m_expected_13)
                  << "  dev = " << Form("%.1f sigma", dev_sigma_13)
                  << "  chi2/ndf = " << Form("%.2f", chi2_13) << std::endl;
        std::cout << "       Dt23(x): m23 = " << Form("%.5f +/- %.5f", m23, m23_err)
                  << "  atteso = " << Form("%.5f", m_expected_23)
                  << "  dev = " << Form("%.1f sigma", dev_sigma_23)
                  << "  chi2/ndf = " << Form("%.2f", chi2_23) << std::endl;
        // Rapporto m13/m23: atteso = -1. Deviazione = asimmetria nei PMT.
        double ratio_m = (fabs(m23) > 1e-9) ? m13 / m23 : 0.0;
        double ratio_m_err = (fabs(m23) > 1e-9)
                           ? fabs(ratio_m) * sqrt((m13_err/m13)*(m13_err/m13)
                                                + (m23_err/m23)*(m23_err/m23)) : 0.0;
        std::cout << "       m13/m23 = " << Form("%.4f +/- %.4f", ratio_m, ratio_m_err)
                  << "  (atteso: -1.000)" << std::endl;
        // Intercette: la differenza q13 - q23 = (delta_1 - delta_2) + L/v e' legata
        // alla differenza dei ritardi dell'elettronica. q13 + q23 = L/v + delta_1 + delta_2 - 2*delta_3.
        std::cout << "       q13 = " << Form("%.3f +/- %.3f", q13, q13_err)
                  << "  q23 = " << Form("%.3f +/- %.3f", q23, q23_err) << " ns"
                  << std::endl;
std::cout << "       q13 - q23 = " << Form("%.3f", q13 - q23)
                  << " ns  (= delta_1 - delta_2 = q di Dt12: "
                  << Form("%.3f", q_fit) << " ns)" << std::endl;

        // ================================================================
        //  CANVAS 9b: CROSS-CHECK C(x) da Dt13 + Dt23
        // ================================================================
        //  Verifica di consistenza interna della calibrazione.
        //
        //  DERIVAZIONE:
        //    Nel codice:  Dt13 = t1 - t3,   Dt23 = t2 - t3
        //    Per definizione:  C = t3 - (t1 + t2)/2
        //    Sommando:  Dt13 + Dt23 = (t1 - t3) + (t2 - t3) = t1 + t2 - 2*t3 = -2C
        //    Quindi:    C_cross = -(Dt13 + Dt23) / 2
        //
        //  Si confrontano TRE stime di C(x):
        //    1) C_direct(x_k):  media della distribuzione t3-(t1+t2)/2 per
        //                       ogni posizione (vettore C_hy gia' disponibile)
        //    2) C_cross(x_k):   -(dt13_hy[k] + dt23_hy[k])/2 punto per punto
        //                       dai centri delle distribuzioni Dt13, Dt23
        //    3) C_cross_fit(x): -(q13 + q23)/2 - (m13 + m23)/2 * x
        //                       dai parametri dei fit lineari di Dt13(x), Dt23(x)
        //
        //  Se le tre stime coincidono, i fit delle singole distribuzioni
        //  (Dt12, Dt13, Dt23, C) sono tutti mutualmente consistenti e la
        //  pipeline di analisi waveform non introduce bias differenziali.
        //  Deviazioni indicano: bias nella ricostruzione del tempo di un
        //  canale specifico (es. clipping recovery differente su PMT1 vs PMT2),
        //  oppure un errore nella definizione dei Dt (segni, canali invertiti).
        {
            TCanvas *c_Ccross = new TCanvas("c_C_crosscheck",
                "Cross-check C(x) da #Delta t_{13} + #Delta t_{23}", 1200, 700);
            c_Ccross->Divide(2, 1);

            // ---- Calcolo punto per punto di C_cross = -(Dt13 + Dt23)/2 ----
            std::vector<double> C_cross(npts), C_cross_e(npts);
            for (int i = 0; i < npts; i++) {
                C_cross[i] = -(dt13_hy[i] + dt23_hy[i]) / 2.0;
                // Propagazione errori (indipendenti):
                // sigma_C_cross = sqrt(sigma_13^2 + sigma_23^2) / 2
                C_cross_e[i] = sqrt(dt13_hy_e[i] * dt13_hy_e[i]
                                  + dt23_hy_e[i] * dt23_hy_e[i]) / 2.0;
            }

            // ---- Curva analitica dal fit lineare:
            //      C_cross_fit(x) = -(q13 + q23)/2 - (m13 + m23)/2 * x ----
            double Ccf_q     = -(q13 + q23) / 2.0;
            double Ccf_m     = -(m13 + m23) / 2.0;
            // Propagazione errori sui parametri della retta combinata
            double Ccf_q_err = sqrt(q13_err * q13_err + q23_err * q23_err) / 2.0;
            double Ccf_m_err = sqrt(m13_err * m13_err + m23_err * m23_err) / 2.0;

// ---- Pannello sinistro: C_direct vs C_cross vs fit quadratici ----
            c_Ccross->cd(1);
            gPad->SetGrid(1, 1);

            // C diretto (misurato dalla distribuzione (t3-(t1+t2)/2))
            TGraphErrors *g_Cdir = new TGraphErrors(npts, x_pos.data(), C_hy.data(),
                                                     dx_pos.data(), C_hy_e.data());
            g_Cdir->SetMarkerStyle(20); g_Cdir->SetMarkerSize(1.1);
            g_Cdir->SetMarkerColor(kRed + 1); g_Cdir->SetLineColor(kRed + 1);
            g_Cdir->SetTitle("Cross-check C(x);"
                             "Posizione x [cm];"
                             "C = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]");

            // C cross-check punto per punto: -(Dt13 + Dt23)/2
            TGraphErrors *g_Ccrs = new TGraphErrors(npts, x_pos.data(), C_cross.data(),
                                                     dx_pos.data(), C_cross_e.data());
            g_Ccrs->SetMarkerStyle(22); g_Ccrs->SetMarkerSize(1.2);
            g_Ccrs->SetMarkerColor(kBlue + 1); g_Ccrs->SetLineColor(kBlue + 1);

            // ---- Fit polinomiale di grado 2 sui punti C_cross ----
            // Si usa la stessa procedura del Canvas 3b per C_direct:
            // fit pesato con 1/sigma^2, con consistency rescaling se chi2/ndf > 1.
            TF1 *f_Ccrs_poly = new TF1("f_C_cross_poly",
                                        "[0] + [1]*x + [2]*x*x", -150, 150);
            f_Ccrs_poly->SetParNames("p_{0}^{cross}", "p_{1}^{cross}", "p_{2}^{cross}");
            // Inizializzazione dai parametri del fit analitico (Ccf_q, Ccf_m)
            // e termine quadratico vicino a zero come stima iniziale
            f_Ccrs_poly->SetParameters(Ccf_q, Ccf_m, 0.0);
            TFitResultPtr r_cross = g_Ccrs->Fit(f_Ccrs_poly, "S R Q");

            double Ccrs_p0     = f_Ccrs_poly->GetParameter(0);
            double Ccrs_p1     = f_Ccrs_poly->GetParameter(1);
            double Ccrs_p2     = f_Ccrs_poly->GetParameter(2);
            double Ccrs_p0_err = f_Ccrs_poly->GetParError(0);
            double Ccrs_p1_err = f_Ccrs_poly->GetParError(1);
            double Ccrs_p2_err = f_Ccrs_poly->GetParError(2);
            double Ccrs_chi2ndf = (f_Ccrs_poly->GetNDF() > 0)
                                 ? f_Ccrs_poly->GetChisquare() / f_Ccrs_poly->GetNDF()
                                 : -1.0;
            // Consistency rescaling: se chi2/ndf > 1 gli errori dei parametri
            // vengono moltiplicati per sqrt(chi2/ndf) (stessa convenzione del Canvas 3b)
            if (Ccrs_chi2ndf > 1.0) {
                double sf = sqrt(Ccrs_chi2ndf);
                Ccrs_p0_err *= sf; Ccrs_p1_err *= sf; Ccrs_p2_err *= sf;
            }

            // Fit quadratico C_poly gia' eseguito in Canvas 3b (C_direct)
            TF1 *f_Cpoly_show = new TF1("f_Cpoly_show_cc",
                                         "[0] + [1]*x + [2]*x*x", -150, 150);
            f_Cpoly_show->SetParameters(C_poly_p0, C_poly_p1, C_poly_p2);
            f_Cpoly_show->SetLineColor(kRed + 1);
            f_Cpoly_show->SetLineStyle(2);
            f_Cpoly_show->SetLineWidth(2);

            f_Ccrs_poly->SetLineColor(kBlue + 1);
            f_Ccrs_poly->SetLineStyle(2);
            f_Ccrs_poly->SetLineWidth(2);

            g_Cdir->Draw("AP");
            g_Ccrs->Draw("P same");
            f_Cpoly_show->Draw("same");
            f_Ccrs_poly->Draw("same");

            TLegend *lg1 = new TLegend(0.12, 0.11, 0.88, 0.42);
            lg1->SetTextFont(42); lg1->SetTextSize(0.030);
            lg1->SetNColumns(1);
            lg1->AddEntry(g_Cdir,
                "C_{diretto}: #LT t_{3}-(t_{1}+t_{2})/2 #GT", "p");
            lg1->AddEntry(f_Cpoly_show,
                Form("  fit poly: %.4f + (%.2e)x + (%.2e)x^{2}  [#chi^{2}/ndf=%.2f]",
                     C_poly_p0, C_poly_p1, C_poly_p2, C_poly_chi2ndf), "l");
            lg1->AddEntry(g_Ccrs,
                "C_{cross}: -(#Delta t_{13}+#Delta t_{23})/2", "p");
            lg1->AddEntry(f_Ccrs_poly,
                Form("  fit poly: %.4f + (%.2e)x + (%.2e)x^{2}  [#chi^{2}/ndf=%.2f]",
                     Ccrs_p0, Ccrs_p1, Ccrs_p2, Ccrs_chi2ndf), "l");
            lg1->Draw();

            // Stampa a schermo per confronto diretto dei parametri
            std::cout << "\n[INFO] Fit quadratico C_cross:" << std::endl;
            std::cout << "       C_cross(x) = "
                      << Form("%.4f +/- %.4f", Ccrs_p0, Ccrs_p0_err)
                      << " + (" << Form("%.6f +/- %.6f", Ccrs_p1, Ccrs_p1_err)
                      << ")*x + (" << Form("%.3e +/- %.3e", Ccrs_p2, Ccrs_p2_err)
                      << ")*x^2  ns" << std::endl;
            std::cout << "       chi2/ndf = " << Form("%.2f", Ccrs_chi2ndf) << std::endl;
            std::cout << "       Confronto p2: C_direct = "
                      << Form("%.3e", C_poly_p2) << "  C_cross = "
                      << Form("%.3e", Ccrs_p2)
                      << "  (compatibili se delta < 2 sigma)" << std::endl;

            // ---- Pannello destro: residuo C_direct - C_cross ----
            // Se la pipeline e' consistente, questo residuo deve essere
            // compatibile con zero entro le barre d'errore per ogni posizione.
            c_Ccross->cd(2);
            gPad->SetGrid(1, 1);

            std::vector<double> dC(npts), dC_e(npts);
            double dC_chi2 = 0.0;
            for (int i = 0; i < npts; i++) {
                dC[i]   = C_hy[i] - C_cross[i];
                // Errore combinato (C_hy e C_cross sono costruiti dagli stessi
                // eventi ma da distribuzioni diverse, quindi sono correlati;
                // l'errore sulla differenza e' conservativo)
                dC_e[i] = sqrt(C_hy_e[i] * C_hy_e[i]
                             + C_cross_e[i] * C_cross_e[i]);
                if (dC_e[i] > 0)
                    dC_chi2 += (dC[i] / dC_e[i]) * (dC[i] / dC_e[i]);
            }
            double dC_chi2ndf = (npts > 0) ? dC_chi2 / npts : -1.0;

            TGraphErrors *g_dC = new TGraphErrors(npts, x_pos.data(), dC.data(),
                                                   dx_pos.data(), dC_e.data());
            g_dC->SetMarkerStyle(20); g_dC->SetMarkerSize(1.1);
            g_dC->SetMarkerColor(kViolet + 1); g_dC->SetLineColor(kViolet + 1);
            g_dC->SetTitle(Form("Residuo C_{diretto} - C_{cross}  "
                                "(#chi^{2}/n = %.2f);Posizione x [cm];"
                                "#Delta C [ns]", dC_chi2ndf));
            g_dC->Draw("AP");

            TLine *lz = new TLine(-150, 0, 150, 0);
            lz->SetLineColor(kBlack); lz->SetLineStyle(2); lz->SetLineWidth(2);
            lz->Draw("same");

            // Bande +/- 50 ps come riferimento visivo
            TLine *lp = new TLine(-150, 0.05, 150, 0.05);
            lp->SetLineColor(kGray + 1); lp->SetLineStyle(3); lp->Draw("same");
            TLine *lm = new TLine(-150, -0.05, 150, -0.05);
            lm->SetLineColor(kGray + 1); lm->SetLineStyle(3); lm->Draw("same");

            c_Ccross->Update();
            fout->cd();
            c_Ccross->Write("Crosscheck_C_from_Dt13_Dt23");

            std::cout << "\n[INFO] Cross-check C(x) da Dt13+Dt23:" << std::endl;
            std::cout << "       Retta C_cross(x) = "
                      << Form("%.4f +/- %.4f", Ccf_q, Ccf_q_err)
                      << " + (" << Form("%.6f +/- %.6f", Ccf_m, Ccf_m_err)
                      << ")*x  ns" << std::endl;
            std::cout << "       Pendenza residua (m13+m23)/2 = "
                      << Form("%.6f +/- %.6f", -Ccf_m, Ccf_m_err)
                      << " ns/cm  (attesa 0 se v uguale nei 2 sensi)"
                      << std::endl;
            std::cout << "       chi2/n (C_dir - C_cross = 0): "
                      << Form("%.2f", dC_chi2ndf) << std::endl;
        }
    }

    // ====================================================================
    //  CANVAS 10: RISOLUZIONE sigma(Dt13) e sigma(Dt23) vs x
    // ====================================================================
    //  Se i due PMT hanno risoluzioni intrinseche uguali (sigma_1 = sigma_2),
    //  le curve sigma(Dt13) e sigma(Dt23) dovrebbero essere simmetriche
    //  l'una dell'altra rispetto al centro barra: sigma(Dt13) a x = +D e'
    //  uguale a sigma(Dt23) a x = -D, perche' in entrambi i casi il PMT
    //  della barra vede un percorso ottico di uguale lunghezza.
    //  Se una delle due curve e' sistematicamente piu' larga, il PMT
    //  corrispondente ha risoluzione peggiore (ringing, time-walk, guadagno).
    {
        TCanvas *c_res13_23 = new TCanvas("c_sigma_dt13_23",
            "Risoluzione #sigma(#Delta t_{13}) e #sigma(#Delta t_{23}) vs x",
            900, 550);
        c_res13_23->SetGrid(1, 1);

        TGraphErrors *g_s13 = new TGraphErrors(npts, x_pos.data(), dt13_hy_w.data(),
                                                dx_pos.data(), dt13_hy_we.data());
        g_s13->SetMarkerStyle(20); g_s13->SetMarkerSize(1.0);
        g_s13->SetMarkerColor(kBlue + 1); g_s13->SetLineColor(kBlue + 1);
        g_s13->SetTitle("Risoluzione vs posizione;"
                        "Posizione x [cm];#sigma [ns]");

        TGraphErrors *g_s23 = new TGraphErrors(npts, x_pos.data(), dt23_hy_w.data(),
                                                dx_pos.data(), dt23_hy_we.data());
        g_s23->SetMarkerStyle(21); g_s23->SetMarkerSize(1.0);
        g_s23->SetMarkerColor(kGreen + 2); g_s23->SetLineColor(kGreen + 2);

        // Anche sigma(Dt12) come riferimento
        TGraphErrors *g_s12 = new TGraphErrors(npts, x_pos.data(), dt12_hy_w.data(),
                                                dx_pos.data(), dt12_hy_we.data());
        g_s12->SetMarkerStyle(22); g_s12->SetMarkerSize(1.0);
        g_s12->SetMarkerColor(kRed + 1); g_s12->SetLineColor(kRed + 1);

        // Trova il range Y adeguato
        double y_min_s = 1e9, y_max_s = 0;
        for (int i = 0; i < npts; i++) {
            if (dt13_hy_w[i] < y_min_s) y_min_s = dt13_hy_w[i];
            if (dt23_hy_w[i] < y_min_s) y_min_s = dt23_hy_w[i];
            if (dt12_hy_w[i] < y_min_s) y_min_s = dt12_hy_w[i];
            if (dt13_hy_w[i] > y_max_s) y_max_s = dt13_hy_w[i];
            if (dt23_hy_w[i] > y_max_s) y_max_s = dt23_hy_w[i];
            if (dt12_hy_w[i] > y_max_s) y_max_s = dt12_hy_w[i];
        }
        g_s13->GetYaxis()->SetRangeUser(y_min_s * 0.85, y_max_s * 1.15);
        g_s13->Draw("AP");
        g_s23->Draw("P same");
        g_s12->Draw("P same");

        TLegend *leg_s = new TLegend(0.55, 0.72, 0.89, 0.89);
        leg_s->SetTextFont(42); leg_s->SetTextSize(0.035);
        leg_s->AddEntry(g_s13, "#sigma(#Delta t_{13})", "p");
        leg_s->AddEntry(g_s23, "#sigma(#Delta t_{23})", "p");
        leg_s->AddEntry(g_s12, "#sigma(#Delta t_{12})", "p");
        leg_s->Draw();

        c_res13_23->Update();
        fout->cd();
        c_res13_23->Write("Resolution_Dt13_Dt23_vs_x");
    }

    // ====================================================================
    //  CANVAS 11a/b/c: RISE TIME vs POSIZIONE — un canvas per PMT
    // ====================================================================
    //  Per ogni PMT viene generato un canvas dedicato con:
    //    - pannello superiore: dati + fit (lineare per PMT1/2, costante per PMT3)
    //    - pannello inferiore: residui dal fit
    //
    //  PMT1: distanza ottica d1 = L/2 + x -> rise time cresce con x -> pendenza > 0
    //  PMT2: distanza ottica d2 = L/2 - x -> rise time decresce con x -> pendenza < 0
    //  PMT3: percorso ottico fisso -> rise time costante, pendenza ~ 0
    //
    //  CANVAS 11d: RECAP — sovrapposizione dei 3 PMT senza fit e senza residui
    //  per visualizzare in modo pulito le differenze di rise time.
    {
        // Filtra le posizioni con dati validi (rise time != -999)
        std::vector<double> rt_x_valid, rt_dx_valid;
        std::vector<double> rt_y[3], rt_ey[3];

        for (int i = 0; i < npts; i++) {
            bool all_valid = true;
            for (int kp = 0; kp < 3; kp++) {
                if (rt_mu[kp][i] < -900.0) { all_valid = false; break; }
            }
            if (!all_valid) continue;
            rt_x_valid.push_back(x_pos[i]);
            rt_dx_valid.push_back(dx_pos[i]);
            for (int kp = 0; kp < 3; kp++) {
                rt_y[kp].push_back(rt_mu[kp][i]);
                rt_ey[kp].push_back(rt_mu_e[kp][i]);
            }
        }
        int nrt = (int)rt_x_valid.size();

        if (nrt < 3) {
            std::cerr << "[WARNING] Meno di 3 posizioni con rise time valido — "
                      << "canvas Rise Time non generati." << std::endl;
        } else {

            const int   rt_colors[3]  = {kBlue + 1, kRed + 1, kGreen + 2};
            const int   rt_markers[3] = {20, 21, 22};
            const char* rt_pmt_lab[3] = {"PMT1", "PMT2", "PMT3"};
            // Nomi dei canvas e titoli ROOT per i 3 PMT
            const char* rt_canvas_name[3] = {
                "c_risetime_PMT1", "c_risetime_PMT2", "c_risetime_PMT3"};
            const char* rt_canvas_title[3] = {
                "Rise Time T_{90}-T_{10} vs Posizione — PMT1",
                "Rise Time T_{90}-T_{10} vs Posizione — PMT2",
                "Rise Time T_{90}-T_{10} vs Posizione — PMT3"};
            const char* rt_write_name[3] = {
                "RiseTime_vs_Position_PMT1",
                "RiseTime_vs_Position_PMT2",
                "RiseTime_vs_Position_PMT3"};

            // Array per salvare i grafici e i fit (servono per il recap e per i printout)
            TGraphErrors *g_rt_all[3] = {nullptr, nullptr, nullptr};
            TF1          *f_rt_all[3] = {nullptr, nullptr, nullptr};

            // Stima delle pendenze iniziali dai dati: rise time approssimativo
            // al primo e ultimo punto, diviso per l'escursione in x.
            // PMT1: pendenza positiva, PMT2: pendenza negativa, PMT3: zero.
            double rt_slope_init[3] = {0.0, 0.0, 0.0};
            double rt_intercept_init[3] = {0.0, 0.0, 0.0};
            for (int kp = 0; kp < 3; kp++) {
                // Media di tutti i rise time come stima dell'intercetta (valore a x=0)
                double sum_rt = 0.0;
                for (int i = 0; i < nrt; i++) sum_rt += rt_y[kp][i];
                rt_intercept_init[kp] = sum_rt / nrt;
                // Pendenza: differenza primo-ultimo diviso escursione x
                if (nrt > 1) {
                    rt_slope_init[kp] = (rt_y[kp][nrt-1] - rt_y[kp][0])
                                      / (rt_x_valid[nrt-1] - rt_x_valid[0]);
                }
            }

            // ---- Loop sui 3 PMT: un canvas ciascuno ----
            for (int kp = 0; kp < 3; kp++) {

                TCanvas *c_rt_k = new TCanvas(rt_canvas_name[kp],
                    rt_canvas_title[kp], 1000, 750);
                c_rt_k->Divide(1, 2);

                // --- Pad superiore: dati + fit ---
                TPad *pad_top = (TPad*)c_rt_k->cd(1);
                pad_top->SetPad(0.0, 0.35, 1.0, 1.0);
                pad_top->SetBottomMargin(0.02);
                pad_top->SetTopMargin(0.08);
                pad_top->SetGrid(1, 1);

                TGraphErrors *g_rt = new TGraphErrors(nrt, rt_x_valid.data(),
                    rt_y[kp].data(), rt_dx_valid.data(), rt_ey[kp].data());
                g_rt->SetName(Form("g_risetime_%s", rt_pmt_lab[kp]));
                g_rt->SetMarkerStyle(rt_markers[kp]);
                g_rt->SetMarkerSize(1.1);
                g_rt->SetMarkerColor(rt_colors[kp]);
                g_rt->SetLineColor(rt_colors[kp]);
                g_rt->SetTitle(Form("Rise Time vs Posizione %s;"
                    "Posizione x [cm];Rise Time T_{90}-T_{10} [ns]",
                    rt_pmt_lab[kp]));
                // L'asse X viene mostrato nel pad superiore con le sue etichette.
                // Nota: SetLabelSize(0) era corretto solo nel vecchio canvas unificato
                // con residui, dove l'asse X era nascosto per evitare duplicati col pad
                // inferiore. Con canvas separati per PMT va rimosso.
                g_rt->GetXaxis()->SetLabelSize(0.045);
                g_rt->GetXaxis()->SetTitleSize(0.048);

                // --- Definizione del modello di fit ---
                TF1 *f_rt = nullptr;
                if (kp < 2) {
                    // PMT1 e PMT2: fit lineare  tau(x) = q + m*x
                    f_rt = new TF1(Form("f_rt_%s", rt_pmt_lab[kp]),
                                   "[0] + [1]*x", -150, 150);
                    f_rt->SetParNames("q (intercetta)", "m (pendenza)");
                    f_rt->SetParameters(rt_intercept_init[kp], rt_slope_init[kp]);
                } else {
                    // PMT3: fit costante  tau = C
                    f_rt = new TF1(Form("f_rt_%s", rt_pmt_lab[kp]),
                                   "[0]", -150, 150);
                    f_rt->SetParNames("#LT#tau#GT");
                    f_rt->SetParameter(0, rt_intercept_init[kp]);
                }
                f_rt->SetLineColor(rt_colors[kp]);
                f_rt->SetLineWidth(2);
                f_rt->SetLineStyle(2);

                // Fit con errori pesati
                g_rt->Fit(f_rt, "S R Q");

                // Consistency rescaling: se chi2/ndf > 1, gli errori dei
                // parametri vengono inflazionati per sqrt(chi2/ndf)
                double rt_chi2ndf = (f_rt->GetNDF() > 0)
                    ? f_rt->GetChisquare() / f_rt->GetNDF() : -1.0;
                double rt_rescale = (rt_chi2ndf > 1.0) ? sqrt(rt_chi2ndf) : 1.0;

                g_rt->Draw("AP");
                f_rt->Draw("same");

                // Box informativo con i risultati del fit
                TPaveText *pt_rt = new TPaveText(0.13, 0.13, 0.60, 0.42, "NDC");
                pt_rt->SetFillColorAlpha(kWhite, 0.85);
                pt_rt->SetTextFont(42);
                pt_rt->SetTextSize(0.042);
                pt_rt->SetTextAlign(12);
                if (kp < 2) {
                    double p0     = f_rt->GetParameter(0);
                    double p0_err = f_rt->GetParError(0) * rt_rescale;
                    double p1     = f_rt->GetParameter(1);
                    double p1_err = f_rt->GetParError(1) * rt_rescale;
                    pt_rt->AddText(Form("#tau(x) = q + m #upoint x   (%s)", rt_pmt_lab[kp]));
                    pt_rt->AddText(Form("q = %.3f #pm %.3f ns", p0, p0_err));
                    pt_rt->AddText(Form("m = %.5f #pm %.5f ns/cm", p1, p1_err));
                    pt_rt->AddText(Form("#chi^{2}/ndf = %.2f", rt_chi2ndf));
                } else {
                    double c0     = f_rt->GetParameter(0);
                    double c0_err = f_rt->GetParError(0) * rt_rescale;
                    pt_rt->AddText(Form("#tau = costante   (%s)", rt_pmt_lab[kp]));
                    pt_rt->AddText(Form("#LT#tau#GT = %.3f #pm %.3f ns", c0, c0_err));
                    pt_rt->AddText(Form("#chi^{2}/ndf = %.2f", rt_chi2ndf));
                }
                pt_rt->Draw();

                // Salva puntatori per recap e printout
                g_rt_all[kp] = g_rt;
                f_rt_all[kp] = f_rt;

                // --- Pad inferiore: residui dal fit ---
                TPad *pad_bot = (TPad*)c_rt_k->cd(2);
                pad_bot->SetPad(0.0, 0.0, 1.0, 0.35);
                pad_bot->SetTopMargin(0.02);
                pad_bot->SetBottomMargin(0.28);
                pad_bot->SetGrid(1, 1);

                std::vector<double> res_y(nrt), res_ey(nrt);
                double res_max = 0.0;
                for (int i = 0; i < nrt; i++) {
                    res_y[i]  = rt_y[kp][i] - f_rt->Eval(rt_x_valid[i]);
                    res_ey[i] = rt_ey[kp][i];
                    double hh = fabs(res_y[i]) + res_ey[i];
                    if (hh > res_max) res_max = hh;
                }
                TGraphErrors *g_res = new TGraphErrors(nrt, rt_x_valid.data(),
                    res_y.data(), rt_dx_valid.data(), res_ey.data());
                g_res->SetMarkerStyle(rt_markers[kp]);
                g_res->SetMarkerSize(1.0);
                g_res->SetMarkerColor(rt_colors[kp]);
                g_res->SetLineColor(rt_colors[kp]);
                g_res->SetTitle(Form(";Posizione x [cm];"
                    "#tau_{mis} - #tau_{fit} [ns]  (%s)", rt_pmt_lab[kp]));
                g_res->GetXaxis()->SetLimits(-150, 150);
                double y_res = std::max(0.20, 1.3 * res_max);
                g_res->GetYaxis()->SetRangeUser(-y_res, y_res);
                g_res->GetXaxis()->SetTitleSize(0.095);
                g_res->GetXaxis()->SetLabelSize(0.082);
                g_res->GetYaxis()->SetTitleSize(0.075);
                g_res->GetYaxis()->SetTitleOffset(0.55);
                g_res->GetYaxis()->SetLabelSize(0.075);
                g_res->GetYaxis()->SetNdivisions(505);
                g_res->Draw("AP");

                TLine *lz = new TLine(-150, 0, 150, 0);
                lz->SetLineColor(kBlack);
                lz->SetLineWidth(2);
                lz->SetLineStyle(2);
                lz->Draw("same");

                c_rt_k->Update();
                fout->cd();
                c_rt_k->Write(rt_write_name[kp]);
                g_rt->Write();
            }

            // ---- Canvas 11d: RECAP — 3 PMT sovrapposti, solo dati ----
            // Nessun fit, nessun residuo: visualizzazione pulita delle
            // differenze di rise time tra i 3 PMT in funzione della posizione.
            {
                TCanvas *c_recap = new TCanvas("c_risetime_recap",
                    "Rise Time vs Posizione — Confronto PMT1 / PMT2 / PMT3",
                    1000, 600);
                c_recap->SetGrid(1, 1);

                // Range Y globale
                double y_lo = 1e9, y_hi = 0;
                for (int kp = 0; kp < 3; kp++) {
                    for (int i = 0; i < nrt; i++) {
                        double lo = rt_y[kp][i] - rt_ey[kp][i];
                        double hi = rt_y[kp][i] + rt_ey[kp][i];
                        if (lo < y_lo) y_lo = lo;
                        if (hi > y_hi) y_hi = hi;
                    }
                }
                double y_margin = 0.08 * (y_hi - y_lo);

                // Crea copie dei grafici per il recap (evita conflitti di ownership)
                TGraphErrors *g_recap[3];
                for (int kp = 0; kp < 3; kp++) {
                    g_recap[kp] = new TGraphErrors(nrt, rt_x_valid.data(),
                        rt_y[kp].data(), rt_dx_valid.data(), rt_ey[kp].data());
                    g_recap[kp]->SetMarkerStyle(rt_markers[kp]);
                    g_recap[kp]->SetMarkerSize(1.1);
                    g_recap[kp]->SetMarkerColor(rt_colors[kp]);
                    g_recap[kp]->SetLineColor(rt_colors[kp]);
                }
                g_recap[0]->SetTitle("Rise Time vs Posizione — Confronto PMT;"
                    "Posizione x [cm];Rise Time T_{90}-T_{10} [ns]");
                g_recap[0]->GetYaxis()->SetRangeUser(y_lo - y_margin, y_hi + y_margin);
                g_recap[0]->Draw("AP");
                g_recap[1]->Draw("P same");
                g_recap[2]->Draw("P same");

                TLegend *leg_rc = new TLegend(0.62, 0.75, 0.89, 0.90);
                leg_rc->SetTextFont(42);
                leg_rc->SetTextSize(0.038);
                leg_rc->AddEntry(g_recap[0], "PMT1", "p");
                leg_rc->AddEntry(g_recap[1], "PMT2", "p");
                leg_rc->AddEntry(g_recap[2], "PMT3", "p");
                leg_rc->Draw();

                c_recap->Update();
                fout->cd();
                c_recap->Write("RiseTime_vs_Position_Recap");
            }

            // ---- Stampa riepilogativa a schermo ----
            std::cout << "\n[INFO] Rise Time T90-T10 vs posizione:" << std::endl;
            for (int kp = 0; kp < 3; kp++) {
                if (!f_rt_all[kp]) continue;
                if (kp < 2) {
                    double rt_chi2 = (f_rt_all[kp]->GetNDF() > 0)
                        ? f_rt_all[kp]->GetChisquare() / f_rt_all[kp]->GetNDF() : -1.0;
                    double rs = (rt_chi2 > 1.0) ? sqrt(rt_chi2) : 1.0;
                    std::cout << "       " << rt_pmt_lab[kp]
                        << ": q = "
                        << Form("%.3f +/- %.3f", f_rt_all[kp]->GetParameter(0),
                                f_rt_all[kp]->GetParError(0) * rs)
                        << " ns,  m = "
                        << Form("%.5f +/- %.5f", f_rt_all[kp]->GetParameter(1),
                                f_rt_all[kp]->GetParError(1) * rs)
                        << " ns/cm  (chi2/ndf = "
                        << Form("%.2f", rt_chi2) << ")" << std::endl;
                } else {
                    double rt_chi2 = (f_rt_all[kp]->GetNDF() > 0)
                        ? f_rt_all[kp]->GetChisquare() / f_rt_all[kp]->GetNDF() : -1.0;
                    double rs = (rt_chi2 > 1.0) ? sqrt(rt_chi2) : 1.0;
                    std::cout << "       " << rt_pmt_lab[kp]
                        << ": <tau> = "
                        << Form("%.3f +/- %.3f", f_rt_all[kp]->GetParameter(0),
                                f_rt_all[kp]->GetParError(0) * rs)
                        << " ns  (costante, chi2/ndf = "
                        << Form("%.2f", rt_chi2) << ")" << std::endl;
                }
            }
            if (f_rt_all[0] && f_rt_all[1]) {
                double s1 = f_rt_all[0]->GetParameter(1);
                double s2 = f_rt_all[1]->GetParameter(1);
                std::cout << "       Rapporto pendenze m1/m2 = "
                    << Form("%.3f", (s2 != 0 ? s1 / s2 : 0))
                    << "  (atteso ~ -1 se PMT simmetrici)" << std::endl;
            }
            
        }
    }
                     
    // ====================================================================
    //  RIEPILOGO FINALE A SCHERMO
    // ====================================================================
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  RISULTATI CALIBRAZIONE (metodo hybrid_tot)  " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Velocita' effettiva: v_eff = "
              << Form("%.2f +/- %.2f", v_eff, v_eff_err) << " cm/ns" << std::endl;
    std::cout << "  Retta Dt12: m = " << Form("%.5f +/- %.5f", m_fit, m_err)
              << " ns/cm,  q = " << Form("%.3f +/- %.3f", q_fit, q_err) << " ns"
              << std::endl;
    std::cout << "  chi2/ndf (retta Dt12): " << Form("%.2f", chi2ndf) << std::endl;
    std::cout << "  C(x) = " << Form("%.4f +/- %.4f", C_lin_p0, C_lin_p0_err)
              << " + (" << Form("%.6f +/- %.6f", C_lin_p1, C_lin_p1_err)
              << ") * x  ns" << std::endl;
    std::cout << "  chi2/ndf (C lineare): " << Form("%.2f", C_lin_chi2ndf) << std::endl;
    std::cout << "  C(x) = " << Form("%.4f", C_poly_p0)
              << " + (" << Form("%.6f", C_poly_p1) << ")*x + ("
              << Form("%.3e", C_poly_p2) << ")*x^2  ns  [PRINCIPALE]" << std::endl;
    std::cout << "  chi2/ndf (C quadratico): " << Form("%.2f", C_poly_chi2ndf) << std::endl;
    std::cout << "  <C> costante = " << Form("%.4f +/- %.4f", C_const_val, C_const_err)
              << " ns   (chi2/ndf = " << Form("%.1f", C_const_chi2ndf) << ")" << std::endl;
    std::cout << "  Output salvato in: " << outname << std::endl;
    std::cout << "=============================================" << std::endl;

    // ====================================================================
    //  PASSO 5: TTree "fit_params" — CONTRATTO con TOF_Analysis_v9.cpp
    // ====================================================================
    // Singolo entry con tutti i parametri numerici che TOF_Analysis_v9.cpp
    // deve caricare automaticamente. NON contiene piu' parametri slew o
    // esponenziali: solo retta Dt12, modello lineare C(x) e i coefficienti
    // del polinomio cubico A(TOT) per PMT1/2/3 (tutti dal metodo hybrid_tot).
    {
        fout->cd();
        TTree *fit_params = new TTree("fit_params",
                                      "Parametri di calibrazione (retta + C(x) + A(TOT))");

        // --- Retta di calibrazione Dt12 = m*x + q ---
        Float_t fp_m, fp_m_err, fp_q, fp_q_err;
        Float_t fp_v_eff, fp_v_eff_err, fp_chi2ndf;
        fit_params->Branch("m",         &fp_m,         "m/F");
        fit_params->Branch("m_err",     &fp_m_err,     "m_err/F");
        fit_params->Branch("q",         &fp_q,         "q/F");
        fit_params->Branch("q_err",     &fp_q_err,     "q_err/F");
        fit_params->Branch("v_eff",     &fp_v_eff,     "v_eff/F");
        fit_params->Branch("v_eff_err", &fp_v_eff_err, "v_eff_err/F");
        fit_params->Branch("chi2ndf",   &fp_chi2ndf,   "chi2ndf/F");

        // --- Modello lineare C(x) = C_lin_p0 + C_lin_p1*x ---
        Float_t fp_C_lin_p0, fp_C_lin_p1, fp_C_lin_p0_err, fp_C_lin_p1_err;
        Float_t fp_C_lin_chi2ndf;
        fit_params->Branch("C_lin_p0",      &fp_C_lin_p0,      "C_lin_p0/F");
        fit_params->Branch("C_lin_p1",      &fp_C_lin_p1,      "C_lin_p1/F");
        fit_params->Branch("C_lin_p0_err",  &fp_C_lin_p0_err,  "C_lin_p0_err/F");
        fit_params->Branch("C_lin_p1_err",  &fp_C_lin_p1_err,  "C_lin_p1_err/F");
        fit_params->Branch("C_lin_chi2ndf", &fp_C_lin_chi2ndf, "C_lin_chi2ndf/F");

        // --- Modello quadratico C(x) = C_poly_p0 + C_poly_p1*x + C_poly_p2*x^2 ---
        // MODELLO PRINCIPALE per C(x): cattura la curvatura parabolica.
        Float_t fp_C_poly_p0, fp_C_poly_p1, fp_C_poly_p2;
        Float_t fp_C_poly_p0_err, fp_C_poly_p1_err, fp_C_poly_p2_err;
        Float_t fp_C_poly_chi2ndf;
        fit_params->Branch("C_poly_p0",      &fp_C_poly_p0,      "C_poly_p0/F");
        fit_params->Branch("C_poly_p1",      &fp_C_poly_p1,      "C_poly_p1/F");
        fit_params->Branch("C_poly_p2",      &fp_C_poly_p2,      "C_poly_p2/F");
        fit_params->Branch("C_poly_p0_err",  &fp_C_poly_p0_err,  "C_poly_p0_err/F");
        fit_params->Branch("C_poly_p1_err",  &fp_C_poly_p1_err,  "C_poly_p1_err/F");
        fit_params->Branch("C_poly_p2_err",  &fp_C_poly_p2_err,  "C_poly_p2_err/F");
        fit_params->Branch("C_poly_chi2ndf", &fp_C_poly_chi2ndf, "C_poly_chi2ndf/F");
// --- Retta di calibrazione Dt12_corr = m_co * x + q_co (P&R) ---
        Float_t fp_m_co, fp_m_co_err, fp_q_co, fp_q_co_err;
        Float_t fp_v_eff_co, fp_v_eff_co_err, fp_chi2ndf_co;
        fit_params->Branch("m_corr",         &fp_m_co,         "m_corr/F");
        fit_params->Branch("m_corr_err",     &fp_m_co_err,     "m_corr_err/F");
        fit_params->Branch("q_corr",         &fp_q_co,         "q_corr/F");
        fit_params->Branch("q_corr_err",     &fp_q_co_err,     "q_corr_err/F");
        fit_params->Branch("v_eff_corr",     &fp_v_eff_co,     "v_eff_corr/F");
        fit_params->Branch("v_eff_corr_err", &fp_v_eff_co_err, "v_eff_corr_err/F");
        fit_params->Branch("chi2ndf_corr",   &fp_chi2ndf_co,   "chi2ndf_corr/F");

        // --- Modello quadratico C_corr(x) = p0 + p1*x + p2*x^2 (P&R) ---
        Float_t fp_C_co_poly_p0, fp_C_co_poly_p1, fp_C_co_poly_p2;
        Float_t fp_C_co_poly_p0_err, fp_C_co_poly_p1_err, fp_C_co_poly_p2_err;
        Float_t fp_C_co_poly_chi2ndf;
        fit_params->Branch("C_poly_corr_p0",      &fp_C_co_poly_p0,      "C_poly_corr_p0/F");
        fit_params->Branch("C_poly_corr_p1",      &fp_C_co_poly_p1,      "C_poly_corr_p1/F");
        fit_params->Branch("C_poly_corr_p2",      &fp_C_co_poly_p2,      "C_poly_corr_p2/F");
        fit_params->Branch("C_poly_corr_p0_err",  &fp_C_co_poly_p0_err,  "C_poly_corr_p0_err/F");
        fit_params->Branch("C_poly_corr_p1_err",  &fp_C_co_poly_p1_err,  "C_poly_corr_p1_err/F");
        fit_params->Branch("C_poly_corr_p2_err",  &fp_C_co_poly_p2_err,  "C_poly_corr_p2_err/F");
        fit_params->Branch("C_poly_corr_chi2ndf", &fp_C_co_poly_chi2ndf, "C_poly_corr_chi2ndf/F");
        // --- Modello costante C(x) = <C> (media pesata) ---
Float_t fp_C_const_val, fp_C_const_err, fp_C_const_chi2ndf;
        fit_params->Branch("C_const_val",      &fp_C_const_val,      "C_const_val/F");
        fit_params->Branch("C_const_err",      &fp_C_const_err,      "C_const_err/F");
        fit_params->Branch("C_const_chi2ndf",  &fp_C_const_chi2ndf,  "C_const_chi2ndf/F");

        // ---- Parametri del RANGE LINEARE RISTRETTO (corrected) — Lin_Range ----
        //  Aggiunti come CONTRATTO con TOF_Analysis per la pipeline Lin_Range.
        //  - Retta Dt12_corr ristretta (m,q,v_eff): per la pipeline 2a (pronta).
        //  - C(x) ristretto LINEARE  (C_lin_co_p0/p1): per la pipeline 2a (pronta).
        //  - C(x) ristretto COSTANTE (C_const_lin):    per la pipeline 2b (scelta).
        //  - linrange_valid: 1 se i fit ristretti sono disponibili, 0 altrimenti.
        Float_t fp_m_lin_co, fp_m_lin_co_err, fp_q_lin_co, fp_q_lin_co_err;
        Float_t fp_v_eff_lin_co, fp_v_eff_lin_co_err, fp_chi2ndf_lin_co;
        fit_params->Branch("m_corr_lin",        &fp_m_lin_co,        "m_corr_lin/F");
        fit_params->Branch("m_corr_lin_err",    &fp_m_lin_co_err,    "m_corr_lin_err/F");
        fit_params->Branch("q_corr_lin",        &fp_q_lin_co,        "q_corr_lin/F");
        fit_params->Branch("q_corr_lin_err",    &fp_q_lin_co_err,    "q_corr_lin_err/F");
        fit_params->Branch("v_eff_corr_lin",    &fp_v_eff_lin_co,    "v_eff_corr_lin/F");
        fit_params->Branch("v_eff_corr_lin_err",&fp_v_eff_lin_co_err,"v_eff_corr_lin_err/F");
        fit_params->Branch("chi2ndf_corr_lin",  &fp_chi2ndf_lin_co,  "chi2ndf_corr_lin/F");

        Float_t fp_C_lin_co_p0, fp_C_lin_co_p0_err;
        Float_t fp_C_lin_co_p1, fp_C_lin_co_p1_err, fp_C_lin_co_chi2ndf;
        fit_params->Branch("C_corr_lin_p0",      &fp_C_lin_co_p0,      "C_corr_lin_p0/F");
        fit_params->Branch("C_corr_lin_p0_err",  &fp_C_lin_co_p0_err,  "C_corr_lin_p0_err/F");
        fit_params->Branch("C_corr_lin_p1",      &fp_C_lin_co_p1,      "C_corr_lin_p1/F");
        fit_params->Branch("C_corr_lin_p1_err",  &fp_C_lin_co_p1_err,  "C_corr_lin_p1_err/F");
        fit_params->Branch("C_corr_lin_chi2ndf", &fp_C_lin_co_chi2ndf, "C_corr_lin_chi2ndf/F");

        Float_t fp_C_const_lin_co, fp_C_const_lin_co_err, fp_C_const_lin_co_chi2ndf;
        fit_params->Branch("C_const_corr_lin",        &fp_C_const_lin_co,
                           "C_const_corr_lin/F");
        fit_params->Branch("C_const_corr_lin_err",    &fp_C_const_lin_co_err,
                           "C_const_corr_lin_err/F");
        fit_params->Branch("C_const_corr_lin_chi2ndf",&fp_C_const_lin_co_chi2ndf,
                           "C_const_corr_lin_chi2ndf/F");

        Int_t fp_linrange_valid = linrange_co_valid ? 1 : 0;
        fit_params->Branch("linrange_valid", &fp_linrange_valid, "linrange_valid/I");

        // --- Metadati del modello A(TOT) ---
        Float_t fp_TOT_thr_ref  = (Float_t)TOT_THR_MV[TOT_CALIB_THR_IDX];
        Int_t   fp_TOT_calib_idx = TOT_CALIB_THR_IDX;
        Int_t   fp_poly_degree   = TOT_POLY_DEGREE;
        Int_t   fp_poly_npar     = TOT_POLY_NPAR;
        Int_t   fp_fix_a0        = FIX_A0_TO_QTHR ? 1 : 0;
        fit_params->Branch("TOT_thr_ref",    &fp_TOT_thr_ref,   "TOT_thr_ref/F");
        fit_params->Branch("TOT_calib_idx",  &fp_TOT_calib_idx, "TOT_calib_idx/I");
        fit_params->Branch("poly_degree",    &fp_poly_degree,   "poly_degree/I");
        fit_params->Branch("poly_npar",      &fp_poly_npar,     "poly_npar/I");
        fit_params->Branch("fix_a0_to_qthr", &fp_fix_a0,        "fix_a0_to_qthr/I");

        // --- Cap di estrapolazione per i 3 PMT ---
        Float_t fp_Amax_TOT_PMT1, fp_Amax_TOT_PMT2, fp_Amax_TOT_PMT3;
        fit_params->Branch("Amax_TOT_PMT1", &fp_Amax_TOT_PMT1, "Amax_TOT_PMT1/F");
        fit_params->Branch("Amax_TOT_PMT2", &fp_Amax_TOT_PMT2, "Amax_TOT_PMT2/F");
        fit_params->Branch("Amax_TOT_PMT3", &fp_Amax_TOT_PMT3, "Amax_TOT_PMT3/F");

        // --- Coefficienti del polinomio A(TOT) per i 3 PMT ---
        // Array di dimensione fissa TOT_POLY_MAX_DEG+1 = 6. Gli elementi con
        // indice > poly_degree restano a zero. Indici: [0]=a0, ..., [3]=a3.
        Float_t fp_poly_TOT_PMT1     [TOT_POLY_MAX_DEG + 1] = {0.0f};
        Float_t fp_poly_TOT_PMT2     [TOT_POLY_MAX_DEG + 1] = {0.0f};
        Float_t fp_poly_TOT_PMT3     [TOT_POLY_MAX_DEG + 1] = {0.0f};
        Float_t fp_poly_TOT_err_PMT1 [TOT_POLY_MAX_DEG + 1] = {0.0f};
        Float_t fp_poly_TOT_err_PMT2 [TOT_POLY_MAX_DEG + 1] = {0.0f};
        Float_t fp_poly_TOT_err_PMT3 [TOT_POLY_MAX_DEG + 1] = {0.0f};
        Float_t fp_poly_TOT_chi2ndf_PMT1 = 0.0f;
        Float_t fp_poly_TOT_chi2ndf_PMT2 = 0.0f;
        Float_t fp_poly_TOT_chi2ndf_PMT3 = 0.0f;

        fit_params->Branch("poly_TOT_PMT1", fp_poly_TOT_PMT1,
                           Form("poly_TOT_PMT1[%d]/F", TOT_POLY_MAX_DEG + 1));
        fit_params->Branch("poly_TOT_PMT2", fp_poly_TOT_PMT2,
                           Form("poly_TOT_PMT2[%d]/F", TOT_POLY_MAX_DEG + 1));
        fit_params->Branch("poly_TOT_PMT3", fp_poly_TOT_PMT3,
                           Form("poly_TOT_PMT3[%d]/F", TOT_POLY_MAX_DEG + 1));
        fit_params->Branch("poly_TOT_err_PMT1", fp_poly_TOT_err_PMT1,
                           Form("poly_TOT_err_PMT1[%d]/F", TOT_POLY_MAX_DEG + 1));
        fit_params->Branch("poly_TOT_err_PMT2", fp_poly_TOT_err_PMT2,
                           Form("poly_TOT_err_PMT2[%d]/F", TOT_POLY_MAX_DEG + 1));
        fit_params->Branch("poly_TOT_err_PMT3", fp_poly_TOT_err_PMT3,
                           Form("poly_TOT_err_PMT3[%d]/F", TOT_POLY_MAX_DEG + 1));
        fit_params->Branch("poly_TOT_chi2ndf_PMT1",
                           &fp_poly_TOT_chi2ndf_PMT1, "poly_TOT_chi2ndf_PMT1/F");
        fit_params->Branch("poly_TOT_chi2ndf_PMT2",
                           &fp_poly_TOT_chi2ndf_PMT2, "poly_TOT_chi2ndf_PMT2/F");
        fit_params->Branch("poly_TOT_chi2ndf_PMT3",
                           &fp_poly_TOT_chi2ndf_PMT3, "poly_TOT_chi2ndf_PMT3/F");

        // --- Popolamento ---
        fp_m         = (Float_t)m_fit;
        fp_m_err     = (Float_t)m_err;
        fp_q         = (Float_t)q_fit;
        fp_q_err     = (Float_t)q_err;
        fp_v_eff     = (Float_t)v_eff;
        fp_v_eff_err = (Float_t)v_eff_err;
        fp_chi2ndf   = (Float_t)chi2ndf;

        fp_C_lin_p0      = (Float_t)C_lin_p0;
        fp_C_lin_p1      = (Float_t)C_lin_p1;
        fp_C_lin_p0_err  = (Float_t)C_lin_p0_err;
        fp_C_lin_p1_err  = (Float_t)C_lin_p1_err;
        fp_C_lin_chi2ndf = (Float_t)C_lin_chi2ndf;

        fp_C_poly_p0      = (Float_t)C_poly_p0;
        fp_C_poly_p1      = (Float_t)C_poly_p1;
        fp_C_poly_p2      = (Float_t)C_poly_p2;
        fp_C_poly_p0_err  = (Float_t)C_poly_p0_err;
        fp_C_poly_p1_err  = (Float_t)C_poly_p1_err;
        fp_C_poly_p2_err  = (Float_t)C_poly_p2_err;
        fp_C_poly_chi2ndf = (Float_t)C_poly_chi2ndf;
        
// Popolamento parametri corrected (P&R)
        fp_m_co         = (Float_t)m_co_fit;
        fp_m_co_err     = (Float_t)m_co_err;
        fp_q_co         = (Float_t)q_co_fit;
        fp_q_co_err     = (Float_t)q_co_err;
        fp_v_eff_co     = (Float_t)v_eff_co;
        fp_v_eff_co_err = (Float_t)v_eff_co_err;
        fp_chi2ndf_co   = (Float_t)chi2ndf_co;

        fp_C_co_poly_p0      = (Float_t)C_co_poly_p0;
        fp_C_co_poly_p1      = (Float_t)C_co_poly_p1;
        fp_C_co_poly_p2      = (Float_t)C_co_poly_p2;
        fp_C_co_poly_p0_err  = (Float_t)C_co_poly_p0_err;
        fp_C_co_poly_p1_err  = (Float_t)C_co_poly_p1_err;
        fp_C_co_poly_p2_err  = (Float_t)C_co_poly_p2_err;
        fp_C_co_poly_chi2ndf = (Float_t)C_co_poly_chi2ndf;
        fp_C_const_val      = (Float_t)C_const_val;
        fp_C_const_err      = (Float_t)C_const_err;
        fp_C_const_chi2ndf  = (Float_t)C_const_chi2ndf;

        // --- Parametri del range lineare ristretto (Lin_Range) ---
        fp_m_lin_co         = (Float_t)m_lin_co;
        fp_m_lin_co_err     = (Float_t)m_lin_co_err;
        fp_q_lin_co         = (Float_t)q_lin_co;
        fp_q_lin_co_err     = (Float_t)q_lin_co_err;
        fp_v_eff_lin_co     = (Float_t)v_eff_lin_co;
        fp_v_eff_lin_co_err = (Float_t)v_eff_lin_co_err;
        fp_chi2ndf_lin_co   = (Float_t)chi2ndf_lin_co;

        fp_C_lin_co_p0      = (Float_t)C_lin_co_p0;
        fp_C_lin_co_p0_err  = (Float_t)C_lin_co_p0_err;
        fp_C_lin_co_p1      = (Float_t)C_lin_co_p1;
        fp_C_lin_co_p1_err  = (Float_t)C_lin_co_p1_err;
        fp_C_lin_co_chi2ndf = (Float_t)C_lin_co_chi2ndf;

        fp_C_const_lin_co        = (Float_t)C_const_lin_co;
        fp_C_const_lin_co_err    = (Float_t)C_const_lin_co_err;
        fp_C_const_lin_co_chi2ndf = (Float_t)C_const_lin_co_chi2ndf;

        fp_Amax_TOT_PMT1 = (Float_t)gTOT_Amax_PMT[0];
        fp_Amax_TOT_PMT2 = (Float_t)gTOT_Amax_PMT[1];
        fp_Amax_TOT_PMT3 = (Float_t)gTOT_Amax_PMT[2];

        for (int k = 0; k <= TOT_POLY_MAX_DEG; k++) {
            fp_poly_TOT_PMT1[k]     = (Float_t)gTOT_poly_PMT[0][k];
            fp_poly_TOT_PMT2[k]     = (Float_t)gTOT_poly_PMT[1][k];
            fp_poly_TOT_PMT3[k]     = (Float_t)gTOT_poly_PMT[2][k];
            fp_poly_TOT_err_PMT1[k] = (Float_t)gTOT_poly_err_PMT[0][k];
            fp_poly_TOT_err_PMT2[k] = (Float_t)gTOT_poly_err_PMT[1][k];
            fp_poly_TOT_err_PMT3[k] = (Float_t)gTOT_poly_err_PMT[2][k];
        }
        fp_poly_TOT_chi2ndf_PMT1 = (Float_t)gTOT_poly_chi2ndf_PMT[0];
        fp_poly_TOT_chi2ndf_PMT2 = (Float_t)gTOT_poly_chi2ndf_PMT[1];
        fp_poly_TOT_chi2ndf_PMT3 = (Float_t)gTOT_poly_chi2ndf_PMT[2];

        fit_params->Fill();
        fit_params->Write();

        std::cout << "\n[INFO] TTree 'fit_params' salvato (contratto con l'analisi):"
                  << std::endl;
        std::cout << "       m, q, v_eff (retta Dt12 hybrid_tot)" << std::endl;
        std::cout << "       C_lin_p0, C_lin_p1 (modello C(x) lineare)" << std::endl;
        std::cout << "       C_poly_p0/p1/p2 (modello C(x) QUADRATICO — principale)"
                  << std::endl;
        std::cout << "       C_const_val (modello C costante — confronto)" << std::endl;
        if (linrange_co_valid) {
            std::cout << "       [Lin_Range] m_corr_lin, q_corr_lin (retta ristretta)"
                      << std::endl;
            std::cout << "       [Lin_Range] C_corr_lin_p0/p1 (C lineare ristretto, 2a)"
                      << std::endl;
            std::cout << "       [Lin_Range] C_const_corr_lin (C costante ristretto, 2b) = "
                      << Form("%.4f +/- %.4f", C_const_lin_co, C_const_lin_co_err)
                      << " ns" << std::endl;
        } else {
            std::cout << "       [Lin_Range] parametri ristretti NON disponibili "
                      << "(range lineare con <3 punti validi)" << std::endl;
        }
        std::cout << "       poly_TOT_PMT1/2/3 (polinomio cubico A(TOT), a0 fissato a "
                  << TOT_THR_MV[TOT_CALIB_THR_IDX] << " mV)" << std::endl;
        for (int k = 0; k < 3; k++) {
            std::cout << "       " << pmt_names[k] << ": calibrazione A(TOT) "
                      << (gTOT_calibrated[k] ? "OK" : "FALLITA")
                      << "  (chi2/ndf = " << Form("%.2f", gTOT_poly_chi2ndf_PMT[k])
                      << ")" << std::endl;
        }
        if (!gTOT_calibrated[0] || !gTOT_calibrated[1] || !gTOT_calibrated[2]) {
            std::cerr << "       [WARNING] Almeno un PMT non e' stato calibrato: "
                      << "il recupero clippati per quel canale sara' disabilitato."
                      << std::endl;
        }
    }

fout->Close();
    std::cout << "\n[DONE] Calibrazione v13 completata." << std::endl;
}


// ==========================================================================
//  SEZIONE 11: STUDIO DELLA RISOLUZIONE vs FRAZIONE CFD — vedi Punto 4
// ==========================================================================
//
//  Misura come la risoluzione temporale sigma(Dt12) e sigma(C) dipende dalla
//  frazione CFD usata per definire l'istante di arrivo, in funzione della
//  posizione x. Frazioni testate: {10, 15, 20, 30, 40, 50} %.
//
//  Frazioni basse (vicino alla baseline): fronte ripido ma SNR peggiore.
//  Frazioni alte (vicino al picco): SNR migliore ma fronte che curva.
//  Esiste una frazione che minimizza sigma; questo studio la individua e
//  quantifica quanto e' robusta la scelta CFD_FRACTION = 0.15 usata altrove.
//
//  Sfrutta la cache del Punto 1: ogni posizione viene caricata una sola volta
//  e tutte le frazioni vengono calcolate dagli stessi eventi in memoria.

/// ComputeCFDTimeAtFraction(): ricalcola il tempo CFD di un canale a una
/// frazione ARBITRARIA, riusando baseline e ampiezza gia' calcolate da
/// AnalyzeChannel. La logica di crossing e interpolazione e' IDENTICA a quella
/// di AnalyzeChannel (PASSO 6), parametrizzata sulla frazione.
///
/// INPUT:
///   cd   : canale gia' analizzato (servono baseline, v_min, time[], voltage[])
///   frac : frazione CFD (es. 0.15)
///   ok   : (output) true se il crossing e' stato trovato
/// RITORNA: tempo CFD interpolato [ns], oppure -999.0 se non valido.
double ComputeCFDTimeAtFraction(const ChannelData& cd, double frac, bool& ok) {
    ok = false;
    if (!cd.has_pulse) return -999.0;
    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return -999.0;

    const float* t = cd.time;
    const float* v = cd.voltage;
    double bl   = cd.baseline;
    double vmin = cd.v_min;
    double amp  = bl - vmin;            // = cd.amplitude (impulso negativo)
    if (amp < NOISE_THRESH) return -999.0;

    // Indice del minimo (non e' memorizzato in ChannelData: lo ricalcolo).
    int imin = 0; double vm = v[0];
    for (int i = 1; i < ns; i++) if (v[i] < vm) { vm = v[i]; imin = i; }

    // Soglia CFD: una frazione 'frac' dell'ampiezza sotto la baseline.
    //   v_thr = bl - frac*amp = bl + frac*(vmin - bl)
    double v_thr = bl - frac * amp;

    // Ricerca IN AVANTI dal fondo della baseline al minimo (identica ad
    // AnalyzeChannel): primo intervallo [i,i+1] dove la tensione, scendendo,
    // attraversa la soglia, con interpolazione lineare sub-campione.
    for (int i = NBL_SAMPLES; i < imin; i++) {
        if (v[i] > v_thr && v[i + 1] <= v_thr) {
            double dv = (double)v[i + 1] - v[i];
            if (fabs(dv) > 1e-6) {
                ok = true;
                return t[i] + (v_thr - v[i]) / dv * (t[i + 1] - t[i]);
            }
            break;
        }
    }
    return -999.0;
}

/// CFD_frac_sdy(): studio della risoluzione temporale al variare della frazione
/// CFD. Funzione AUTONOMA, da chiamare separatamente da TOF_Calibration.
///
/// UTILIZZO (ROOT):
///   root -l 'TOF_Calibration_v14.cpp'
///   root [1] CFD_frac_sdy("/percorso/cartella/dati/")
///   root [2] CFD_frac_sdy("/percorso/dati/", "studio_cfd.root")
///
/// ARGOMENTI:
///   folder        — cartella con i file XML di calibrazione (come TOF_Calibration)
///   outname       — file ROOT di output dello studio
///   parallax_file — file MC di parallasse per gli errori su x (nullptr = fallback)
///
/// OUTPUT nel file ROOT:
///   - una sottocartella per frazione (cfd_10_pct, ...) con:
///       * istogrammi h_dt12_cfdXX_<pos> e h_C_cfdXX_<pos> per posizione;
///       * g_sigma_dt12_cfdXX e g_sigma_C_cfdXX (sigma vs x per quella frazione);
///   - due canvas riassuntivi al livello principale:
///       * Resolution_Dt12_vs_CFDfrac : sigma(Dt12) vs x, tutte le frazioni;
///       * Resolution_C_vs_CFDfrac    : sigma(C)    vs x, tutte le frazioni.
void CFD_frac_sdy(const char* folder,
                  const char* outname = "TOF_CFD_frac_study.root",
                  const char* parallax_file = nullptr) {

gROOT->SetBatch(kTRUE);
    TF1::DefaultAddToGlobalList(kFALSE);   // previene segfault allo shutdown (.q)
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // ---- Frazioni CFD da studiare, con colore e marker dedicati ----
    const int    NF = 6;
    const double fracs[NF] = { 0.10, 0.15, 0.20, 0.30, 0.40, 0.50 };
    const int    pct  [NF] = {   10,   15,   20,   30,   40,   50 };
    const int    cols [NF] = { kRed + 1, kOrange + 1, kGreen + 2,
                               kCyan + 2, kBlue + 1,  kViolet + 1 };
    const int    mrks [NF] = { 20, 21, 22, 23, 29, 33 };

    std::cout << "=============================================" << std::endl;
    std::cout << "  STUDIO RISOLUZIONE vs FRAZIONE CFD (Punto 4)" << std::endl;
    std::cout << "  Folder: " << folder  << std::endl;
    std::cout << "  Output: " << outname << std::endl;
    std::cout << "  Frazioni: 10 15 20 30 40 50 %" << std::endl;
    std::cout << "=============================================" << std::endl;

    const char* eff_parallax = parallax_file ? parallax_file : PARALLAX_FILE;

    TFile* fout = new TFile(outname, "RECREATE");
    if (!fout || !fout->IsOpen()) {
        std::cerr << "[ERRORE] Impossibile creare " << outname << std::endl;
        return;
    }

    // ---- Sottocartelle per frazione ----
    TDirectory* dir_frac[NF];
    for (int f = 0; f < NF; f++)
        dir_frac[f] = fout->mkdir(Form("cfd_%02d_pct", pct[f]));

    // ---- Contenitori dei risultati: vettori paralleli, uno per frazione ----
    // res_x[f] = posizioni; res_sdt12[f]/res_sdt12_e[f] = sigma(Dt12) e errore;
    // res_sC[f]/res_sC_e[f] = sigma(C) e errore.
    std::vector<double> res_x[NF], res_dx[NF];
    std::vector<double> res_sdt12[NF], res_sdt12_e[NF];
    std::vector<double> res_sC[NF],    res_sC_e[NF];

    // ---- Lambda: range adattivo da mediana + MAD (come in ProcessDataset) ----
    auto computeRange = [](std::vector<float>& vals, double nsigma,
                           double range_min, double& lo, double& hi) {
        if (vals.empty()) { lo = -10.0; hi = 10.0; return; }
        std::vector<float> v = vals;          // copia: non altero l'ordine originale
        std::sort(v.begin(), v.end());
        int n = (int)v.size();
        double median = (n % 2 == 0)
                      ? 0.5 * ((double)v[n/2 - 1] + (double)v[n/2])
                      : (double)v[n/2];
        std::vector<double> devs(n);
        for (int i = 0; i < n; i++) devs[i] = fabs((double)v[i] - median);
        std::sort(devs.begin(), devs.end());
        double mad = (n % 2 == 0)
                   ? 0.5 * (devs[n/2 - 1] + devs[n/2]) : devs[n/2];
        double sigma_est  = 1.4826 * mad;
        double half_range = std::max(nsigma * sigma_est, range_min);
        lo = median - half_range;
        hi = median + half_range;
    };

    // ====================================================================
    //  LOOP SULLE POSIZIONI (caricate UNA volta) x LOOP SULLE FRAZIONI
    // ====================================================================
    for (size_t ip = 0; ip < DEFAULT_POSITIONS.size(); ip++) {

        const PositionInfo& pos = DEFAULT_POSITIONS[ip];

        // Risoluzione del path XML (con o senza slash, come altrove).
        std::string xml_path = std::string(folder) + "/" + pos.name + ".xml";
        std::ifstream test(xml_path.c_str());
        if (!test.good()) {
            xml_path = std::string(folder) + pos.name + ".xml";
            test.open(xml_path.c_str());
            if (!test.good()) {
                std::cerr << "[WARNING] CFD_frac_sdy: file non trovato per "
                          << pos.name << " — salto." << std::endl;
                continue;
            }
        }
        test.close();

        // Caricamento eventi (cache del Punto 1: una sola lettura per posizione).
        std::vector<EventData> events;
        int n_parsed = LoadOrParseDataset(xml_path.c_str(), pos.name, events);
        if (n_parsed <= 0) continue;

        double dx_k = ComputeSigmaX(pos.x_cm, eff_parallax);

        // Per ogni frazione: vettori dei Dt12 e C ricalcolati a quella frazione.
        std::vector<float> v_dt12[NF], v_C[NF];
        for (int f = 0; f < NF; f++) {
            v_dt12[f].reserve(events.size());
            v_C[f].reserve(events.size());
        }

        // ---- Un solo passaggio sugli eventi: per ognuno, tutte le frazioni ----
        for (const auto& e : events) {
            if (e.nchannels < 3) continue;
            const ChannelData& c1 = e.ch[0];
            const ChannelData& c2 = e.ch[1];
            const ChannelData& c3 = e.ch[2];

            // Selezione PULITA (indipendente dalla frazione): tutti e 3 i PMT
            // con impulso, non clippati, non oscillanti. Garantisce che il
            // campione sia lo stesso per tutte le frazioni -> confronto equo.
            if (!c1.has_pulse || !c2.has_pulse || !c3.has_pulse) continue;
            if (c1.is_clipped || c2.is_clipped || c3.is_clipped) continue;
            if (ENABLE_OSC_CUT &&
                (c1.is_oscillating || c2.is_oscillating || c3.is_oscillating))
                continue;

            for (int f = 0; f < NF; f++) {
                bool o1, o2, o3;
                double t1 = ComputeCFDTimeAtFraction(c1, fracs[f], o1);
                double t2 = ComputeCFDTimeAtFraction(c2, fracs[f], o2);
                double t3 = ComputeCFDTimeAtFraction(c3, fracs[f], o3);
                if (!o1 || !o2 || !o3) continue;   // crossing valido su tutti e 3
                v_dt12[f].push_back((float)(t1 - t2));
                v_C[f]   .push_back((float)(t3 - 0.5 * (t1 + t2)));
            }
        }

        // ---- Per ogni frazione: istogrammi + fit gaussiano di core ----
        for (int f = 0; f < NF; f++) {
            if ((int)v_dt12[f].size() < 30) {
                std::cerr << "[INFO] " << pos.name << " @ CFD " << pct[f]
                          << "%: solo " << v_dt12[f].size()
                          << " eventi — punto escluso." << std::endl;
                continue;
            }

            // --- Istogramma e fit di Dt12 ---
            double lo, hi;
            computeRange(v_dt12[f], DT_RANGE_NSIGMA, DT_RANGE_MIN, lo, hi);
            TH1D* h_dt12 = new TH1D(
                Form("h_dt12_cfd%02d_%s", pct[f], pos.name.c_str()),
                Form("#Delta t_{12} (CFD %d%%) @ %s;#Delta t_{12} [ns];Conteggi",
                     pct[f], pos.name.c_str()),
                NBINS_DT, lo, hi);
            for (float val : v_dt12[f]) h_dt12->Fill(val);
            FitResult fr_dt12 = FitGaussianCore(h_dt12, 0.90,
                Form("fit_dt12_cfd%02d_%s", pct[f], pos.name.c_str()));

            // --- Istogramma e fit di C ---
            computeRange(v_C[f], DT_RANGE_NSIGMA, DT_RANGE_MIN, lo, hi);
            TH1D* h_C = new TH1D(
                Form("h_C_cfd%02d_%s", pct[f], pos.name.c_str()),
                Form("C (CFD %d%%) @ %s;C = t_{3}-#frac{t_{1}+t_{2}}{2} [ns];Conteggi",
                     pct[f], pos.name.c_str()),
                NBINS_DT, lo, hi);
            for (float val : v_C[f]) h_C->Fill(val);
            FitResult fr_C = FitGaussianCore(h_C, 0.90,
                Form("fit_C_cfd%02d_%s", pct[f], pos.name.c_str()));

            // Salva gli istogrammi nella sottocartella della frazione.
            dir_frac[f]->cd();
            h_dt12->Write();
            h_C->Write();

            // Accumula sigma (= risoluzione) e suo errore.
            res_x[f]      .push_back(pos.x_cm);
            res_dx[f]     .push_back(dx_k);
            res_sdt12[f]  .push_back(fr_dt12.width);
            res_sdt12_e[f].push_back(fr_dt12.width_err);
            res_sC[f]     .push_back(fr_C.width);
            res_sC_e[f]   .push_back(fr_C.width_err);
        }

        std::cout << "[INFO] " << pos.name << " elaborata ("
                  << events.size() << " eventi)." << std::endl;
    }

    // ====================================================================
    //  COSTRUZIONE DEI GRAFICI sigma vs x (uno per frazione) E SALVATAGGIO
    // ====================================================================
    // Per ogni frazione si creano due TGraphErrors (Dt12 e C). Vengono salvati
    // nella sottocartella della frazione e raccolti in due TMultiGraph per i
    // canvas riassuntivi.
    TMultiGraph* mg_dt12 = new TMultiGraph();
    TMultiGraph* mg_C    = new TMultiGraph();
    mg_dt12->SetTitle("Risoluzione #sigma(#Delta t_{12}) vs posizione — "
                      "confronto frazione CFD;Posizione x [cm];#sigma(#Delta t_{12}) [ns]");
    mg_C->SetTitle("Risoluzione #sigma(C) vs posizione — "
                   "confronto frazione CFD;Posizione x [cm];#sigma(C) [ns]");

    TLegend* leg_dt12 = new TLegend(0.70, 0.68, 0.89, 0.89);
    TLegend* leg_C    = new TLegend(0.70, 0.68, 0.89, 0.89);
    leg_dt12->SetTextFont(42); leg_dt12->SetTextSize(0.030);
    leg_C->SetTextFont(42);    leg_C->SetTextSize(0.030);

    for (int f = 0; f < NF; f++) {
        int n = (int)res_x[f].size();
        if (n < 1) continue;

        TGraphErrors* g_dt12 = new TGraphErrors(n, res_x[f].data(),
            res_sdt12[f].data(), res_dx[f].data(), res_sdt12_e[f].data());
        g_dt12->SetName(Form("g_sigma_dt12_cfd%02d", pct[f]));
        g_dt12->SetTitle(Form("CFD %d%%", pct[f]));
        g_dt12->SetMarkerStyle(mrks[f]);
        g_dt12->SetMarkerSize(1.2);
        g_dt12->SetMarkerColor(cols[f]);
        g_dt12->SetLineColor(cols[f]);

        TGraphErrors* g_C = new TGraphErrors(n, res_x[f].data(),
            res_sC[f].data(), res_dx[f].data(), res_sC_e[f].data());
        g_C->SetName(Form("g_sigma_C_cfd%02d", pct[f]));
        g_C->SetTitle(Form("CFD %d%%", pct[f]));
        g_C->SetMarkerStyle(mrks[f]);
        g_C->SetMarkerSize(1.2);
        g_C->SetMarkerColor(cols[f]);
        g_C->SetLineColor(cols[f]);

        // Salva i singoli grafici nella sottocartella della frazione.
        dir_frac[f]->cd();
        g_dt12->Write();
        g_C->Write();

        // Aggiungi ai grafici riassuntivi e alle legende.
        mg_dt12->Add(g_dt12, "P");
        mg_C->Add(g_C, "P");
        leg_dt12->AddEntry(g_dt12, Form("CFD %d%%", pct[f]), "pe");
        leg_C->AddEntry(g_C,       Form("CFD %d%%", pct[f]), "pe");
    }

    // ---- Canvas riassuntivo sigma(Dt12) ----
    {
        TCanvas* c1 = new TCanvas("c_res_dt12_cfd",
            "Risoluzione Dt12 vs frazione CFD", 1000, 650);
        c1->SetGrid(1, 1);
        mg_dt12->Draw("A");
        leg_dt12->Draw();
        c1->Update();
        fout->cd();
        c1->Write("Resolution_Dt12_vs_CFDfrac");
    }

    // ---- Canvas riassuntivo sigma(C) ----
    {
        TCanvas* c2 = new TCanvas("c_res_C_cfd",
            "Risoluzione C vs frazione CFD", 1000, 650);
        c2->SetGrid(1, 1);
        mg_C->Draw("A");
        leg_C->Draw();
        c2->Update();
        fout->cd();
        c2->Write("Resolution_C_vs_CFDfrac");
    }

    // ---- Salvataggio dei TMultiGraph (oltre ai canvas) ----
    fout->cd();
    mg_dt12->Write("mg_sigma_dt12_vs_x");
    mg_C->Write("mg_sigma_C_vs_x");

    fout->Close();
    std::cout << "\n[DONE] Studio CFD completato. Output: " << outname << std::endl;
    std::cout << "       Canvas: Resolution_Dt12_vs_CFDfrac, "
              << "Resolution_C_vs_CFDfrac" << std::endl;
}