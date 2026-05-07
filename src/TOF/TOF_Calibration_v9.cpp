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
//  AUTORI: Luca (sviluppo) + Claude (assistenza)
//  DATA:   Aprile 2026
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
const double CLIP_V_LO       = -499.0;
// ============================================================
//  RILEVAMENTO OSCILLAZIONE: discriminatore a isteresi
// ============================================================
//
// MOTIVAZIONE FISICA:
//   La precedente versione del check (vmax > OSC_VMAX_THRESH) marcava
//   come "oscillante" ogni evento con una singola eccursione positiva
//   sopra +100 mV. Questo include due fenomeni completamente diversi:
//
//   (1) OVERSHOOT POST-CLIP: dopo un impulso negativo molto grande
//       (saturato a -500 mV), l'accoppiamento AC della catena di
//       lettura (capacità di disaccoppiamento + impedenza dei cavi)
//       produce UN SINGOLO transiente positivo nella coda, di
//       ampiezza tipica +100 a +300 mV. È un artefatto ripetibile
//       e DETERMINISTICO della catena, NON un segnale patologico.
//       Il fronte di salita iniziale è perfetto e usabile per il
//       timing → l'evento NON deve essere scartato.
//
//   (2) VERA OSCILLAZIONE PATOLOGICA: pickup EMI sui cavi, scariche
//       nel PMT, pile-up massivo che mette in oscillazione l'elettronica.
//       Il segnale presenta una componente sinusoidale persistente
//       che attraversa la soglia +100 mV ripetutamente, salendo e
//       scendendo molte volte → l'evento DEVE essere scartato.
//
// DISCRIMINANTE:
//   Contiamo il numero di volte che il segnale ATTRAVERSA IN SALITA
//   la soglia OSC_V_HIGH dopo il picco negativo. Un solo overshoot
//   conta UNA volta; un'oscillazione con N cicli conta N volte.
//
// ISTERESI (per evitare falsi crossing dovuti a ringing fine):
//   Un nuovo crossing in salita si conta SOLO DOPO che il segnale è
//   prima sceso al di sotto della soglia bassa OSC_V_LOW. Questo
//   implementa un classico "Schmitt trigger" software: il rumore
//   fine intorno alla soglia non genera conteggi spuri.

// Soglia ALTA per definire un crossing di "oscillazione".
// Stesso valore della v5 per mantenere la stessa scala fisica.
const double OSC_V_HIGH      = 50.0;   // [mV] stessa scala di prima

// Soglia BASSA per l'isteresi. Tra OSC_V_LOW e OSC_V_HIGH il
// contatore non cambia stato. Solo quando il segnale scende sotto
// OSC_V_LOW e poi risale sopra OSC_V_HIGH si conta UN nuovo crossing.
// Valore scelto: metà della soglia alta, garantisce robustezza al
// ringing di ampiezza fino a ~25 mV picco-picco.
const double OSC_V_LOW       = 25.0;   // [mV] soglia bassa di isteresi

// Numero di crossing in salita oltre il quale l'evento è dichiarato
// oscillante (e quindi scartato).
//   1 crossing = singolo overshoot post-clip → ammesso (fisiologico)
//   2 crossing = doppio rimbalzo → ammesso (margine di tolleranza)
//   3+ crossing = oscillazione vera, richiede componente periodica → scartato
const int    OSC_NCROSS_MAX  = 3;

// Per RETROCOMPATIBILITÀ: manteniamo la vecchia costante per chi
// volesse riabilitare il vecchio criterio (semplice ma rozzo).
// Non viene più usata di default.
const double OSC_VMAX_THRESH = 100.0;  // [mV] DEPRECATA, usata solo se
                                       // OSC_USE_LEGACY = true (vedi sotto)

// Flag globale per scegliere tra:
//   false → nuovo criterio basato sul conteggio dei crossing (DEFAULT)
//   true  → vecchio criterio basato su V_max > OSC_VMAX_THRESH (legacy)
// Utile per studi di sistematica o per confrontare i risultati.
const bool   OSC_USE_LEGACY  = false;

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
const int    SLEW_NMIN     = 4;

// Soglia di qualità sul χ²/ndf del fit lineare. Eventi con χ²/ndf > soglia
// vengono scartati dalla calibrazione e dalla ricostruzione. Il χ² è
// calcolato usando come errore tipico del campione max(bl_rms, SLEW_SIGMA_FLOOR).
const double SLEW_CHI2_MAX = 1e7;
//PENSO SIA DA ELIMINARE PERCHé controllo troppo rigido che elimina tutti gli eventi clippatirecuperabili
//rumore ed errori troppo piccoli per avere un buon fit lineare
//ma importante avere questi eventi sennò non c'è statistica






// Floor per σ_y nel calcolo del χ²: evita che eventi con baseline RMS
// anomalamente bassa (sotto 2 mV) abbiano χ² gonfiati artificialmente
// dalla pura curvatura fisica del fronte di salita.
const double SLEW_SIGMA_FLOOR = 2.0; // [mV]

// Margine sopra |CLIP_V_LO| sotto cui la soglia CFD ricostruita è
// considerata "troppo vicina al clipping" e l'evento viene scartato.
// Es: se |CLIP_V_LO| = 450 mV e margine = 20 mV, accettiamo solo
// soglie CFD nella regione v_thr > -430 mV (cioè u_thr < 430 mV).
const double SLEW_CLIP_MARGIN = 5.0;  // [mV]

// ============================================================
//  PARAMETRI PER IL RECUPERO TRAMITE TIME-OVER-THRESHOLD (TOT)
// ============================================================
//
// Strategia parallela al recupero via slew rate. Per impulsi a forma
// fissa A·g(t), la durata del segnale sopra una soglia q (TOT) è una
// funzione monotona crescente di A. Calibrando A vs TOT su eventi
// NON clippati, possiamo poi estrapolare ai clippati per stimare A
// e ricostruire il tempo CFD su soglia f·A_TOT.
//
// MODELLO FUNZIONALE:
//      A = p0 · exp(p1 · TOT)
// equivalente a:  ln A = ln p0 + p1·TOT
// È giustificato dal limite asintotico TOT ~ τ·ln(A/q) per impulsi
// con coda esponenziale (τ_eff = 1/p1).

// Insieme di soglie [mV] a cui calcoliamo TOT per ogni evento.
// Le scelgo con un set ragionevolmente ampio per permettere studi
// di sistematica (vedi quale soglia dà miglior risoluzione).
// Vincoli:
//   - tutte ≥ ~5σ_baseline per evitare false detect dovute al rumore
//   - tutte ≤ |CLIP_V_LO| − TOT_CLIP_MARGIN, altrimenti su clippati
//     anche il "fronte sotto soglia" cade nel plateau e il TOT è
//     distorto.
const int    NTOT_THR        = 5;
const double TOT_THR_MV[NTOT_THR] = { 50.0, 100.0, 150.0, 200.0, 300.0 };

// Indice (0-based) della soglia di riferimento usata per il fit di
// calibrazione A vs TOT e per il recupero ampiezza dei clippati.
// 100 mV (indice 1) è un buon default: ~20σ_bl, ~5× sotto il clip.
const int    TOT_CALIB_THR_IDX = 1;   // → TOT_THR_MV[1] = 100 mV

// Margine sopra |CLIP_V_LO| per dichiarare una soglia "troppo vicina
// al clipping". Se q > |CLIP_V_LO| − TOT_CLIP_MARGIN, il TOT a quella
// soglia non viene calcolato neanche sui non clippati (perché il
// fall crossing potrebbe già essere distorto da overshoot post-clip).
const double TOT_CLIP_MARGIN = 20.0;  // [mV]

// Numero minimo di campioni tra rise e fall per accettare il TOT
// (filtra rumore impulsivo da scariche o EMI che potrebbe simulare
// un crossing rapido con TOT di pochi sample).
const int    TOT_NMIN_SAMPLES = 3;

// Soglia massima A_TOT estrapolato accettabile, espressa come
// multiplo dell'ampiezza massima del campione di calibrazione.
// Oltre, l'estrapolazione esponenziale è considerata troppo aggressiva
// e l'evento NON viene recuperato via TOT.
const double TOT_AREC_MAX_FACTOR = 3.0;


// ============================================================
//  COSTANTI DI CALIBRAZIONE TOT (popolate da CalibrateTOT)
// ============================================================
//   A_TOT(TOT) = gTOT_p0_PMT[k] · exp(gTOT_p1_PMT[k] · TOT)
// Indici 0,1 = PMT1, PMT2. PMT3 non viene calibrato (non clippa).
static double gTOT_p0_PMT[MAX_CHANNELS]      = {0.0, 0.0, 0.0, 0.0};
static double gTOT_p1_PMT[MAX_CHANNELS]      = {0.0, 0.0, 0.0, 0.0};
static double gTOT_Amax_PMT[MAX_CHANNELS]    = {0.0, 0.0, 0.0, 0.0};  // A massima nel calib (per il cap)
static bool   gTOT_calibrated                = false;

// --- Range di posizioni usate per la calibrazione di k ---
// Solo dataset con |x_k| <= CALIB_K_X_ABS_MAX contribuiscono al fit
// di k. Giustificazione: a queste posizioni, entrambi i PMT vedono
// abbondanti eventi NON clippati su un range di ampiezze sufficiente
// per determinare la pendenza A = k·s. Posizioni più estreme (±112,
// ±130) sono dominate da clipping su un PMT e quindi non utili per
// la calibrazione (anche se sono i target dell'applicazione).
const double CALIB_K_X_ABS_MAX = 56.0;  // [cm]


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
// Intercetta del fit A = k·s + q. Tipicamente piccola (compatibile
// con 0 entro pochi σ se il fronte è lineare nella finestra slew),
// ma usata nel recupero per coerenza col modello del fit.
static double gQ_PMT[MAX_CHANNELS]      = {0.0, 0.0, 0.0, 0.0};   // q [mV]
static double gQ_PMT_err[MAX_CHANNELS]  = {0.0, 0.0, 0.0, 0.0};   // errore di q
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
const bool TAPE_MODEL_INCREMENTAL = true;

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

// ---- Diagnostica del rilevamento oscillazione ----
    // Numero di crossing in salita della soglia OSC_V_HIGH, contati
    // DOPO il picco negativo dell'impulso, con isteresi su OSC_V_LOW.
    // Vedi PASSO 5 di AnalyzeChannel per la definizione algoritmica.
    int    n_pos_crossings;

    // ---- Time-Over-Threshold multi-soglia (PASSO 8) ----
    // tot_q[k] = durata del segnale sopra la soglia TOT_THR_MV[k] [ns]
    // tot_rise[k], tot_fall[k] = istanti dei due crossing [ns]
    // tot_ok[k]   = true se entrambi i crossing sono validi e
    //               la differenza è > TOT_NMIN_SAMPLES campioni
    // Calcolato per OGNI evento (clippato o no): i non clippati
    // servono per calibrare la curva A vs TOT, i clippati per il
    // recupero in PASSO 9.
    double tot_q   [NTOT_THR];
    double tot_rise[NTOT_THR];
    double tot_fall[NTOT_THR];
    bool   tot_ok  [NTOT_THR];

    // ---- Recupero CFD via TOT (analogo a slew rate) ----
    // Popolati da RecoverClippedCFD_TOT(): A_rec = p0·exp(p1·TOT)
    double amplitude_tot;     // ampiezza ricostruita via TOT [mV]
    double t_cfd_rec_tot;     // tempo CFD su soglia f·A_tot [ns]
    bool   cfd_recovered_tot; // true se CFD ricostruito è valido
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
    cd.n_pos_crossings = 0;          // contatore oscillazione (Schmitt trigger)

    // Inizializzazione campi TOT (PASSO 8)
    for (int kth = 0; kth < NTOT_THR; kth++) {
        cd.tot_q   [kth] = 0.0;
        cd.tot_rise[kth] = 0.0;
        cd.tot_fall[kth] = 0.0;
        cd.tot_ok  [kth] = false;
    }
    cd.amplitude_tot     = 0.0;
    cd.t_cfd_rec_tot     = -999.0;
    cd.cfd_recovered_tot = false;
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

// ---- PASSO 5: RILEVAMENTO OSCILLAZIONE (conteggio dei crossing con isteresi) ----
    //
    // Un impulso di scintillazione del PMT è monopolare negativo per
    // costruzione. La regione di INTERESSE per il timing è il fronte
    // di salita PRIMA del picco negativo (indice imin). La regione
    // DOPO il picco contiene la coda di discesa, eventuale overshoot
    // post-clip, e/o oscillazioni patologiche.
    //
    // ALGORITMO (Schmitt trigger software):
    //   1. Iteriamo sui campioni DOPO il picco (i = imin+1 ... ns-1).
    //   2. Manteniamo uno stato che ricorda se siamo "alti" o "bassi":
    //        state == false → siamo sotto OSC_V_LOW (riposo "basso")
    //        state == true  → abbiamo superato OSC_V_HIGH (riposo "alto")
    //   3. Una transizione false→true (cioè: il segnale era sceso sotto
    //      OSC_V_LOW e ora è risalito sopra OSC_V_HIGH) conta come UN
    //      crossing in salita.
    //   4. La transizione opposta (true→false: salito sopra OSC_V_HIGH e
    //      ora sceso sotto OSC_V_LOW) NON incrementa il contatore, ma
    //      "riarma" lo Schmitt trigger per il prossimo crossing.
    //   5. La banda morta tra OSC_V_LOW e OSC_V_HIGH (= 25 mV di larghezza)
    //      assorbe il ringing fine e il rumore senza generare conteggi.
    //
    // STATO INIZIALE:
    //   Inizializziamo state in base al primo campione dopo imin:
    //   se v[imin+1] è già sopra OSC_V_HIGH (improbabile, perché il
    //   campione subito dopo il picco è ancora vicino a vmin), state=true;
    //   altrimenti state=false. Questo evita di contare un crossing
    //   spurio se la traversata avviene proprio al primo campione.
    //
    // CASO PARTICOLARE: se imin è troppo vicino alla fine del record
    //   (impulso che arriva tardi nella finestra), non c'è abbastanza
    //   coda per analizzare l'oscillazione → n_pos_crossings = 0.
    {
        int n_cross = 0;

        // Stato iniziale dello Schmitt trigger:
        //   true  = siamo nella regione "alta" (sopra OSC_V_HIGH)
        //   false = siamo nella regione "bassa" (sotto OSC_V_LOW)
        // Se siamo nella zona morta intermedia, manteniamo la convenzione
        // "bassa" (significa: il prossimo evento sopra V_HIGH conterà).
        bool above = false;

        // Partiamo dal campione successivo al picco. Nota: vogliamo
        // analizzare la coda, NON il fronte di salita (che per impulsi
        // negativi va verso valori sempre più negativi e non genera
        // crossing positivi).
        for (int i = imin + 1; i < ns; i++) {

            double vi = (double)v[i];

            if (!above) {
                // Stato basso: cerchiamo una salita oltre OSC_V_HIGH.
                // Quando la troviamo, contiamo UN crossing e passiamo
                // allo stato alto.
                if (vi > OSC_V_HIGH) {
                    n_cross++;
                    above = true;
                }
                // Nella zona morta (OSC_V_LOW < vi <= OSC_V_HIGH) non
                // facciamo nulla: aspettiamo un'eccursione decisa.
            } else {
                // Stato alto: cerchiamo una discesa sotto OSC_V_LOW.
                // Quando la troviamo, riarmiamo lo Schmitt trigger
                // SENZA contare nulla (la discesa non è un crossing
                // di oscillazione, è solo il "ritorno alla base").
                if (vi < OSC_V_LOW) {
                    above = false;
                }
            }
        }

        cd.n_pos_crossings = n_cross;

        // ---- Decisione del flag is_oscillating ----
        // Due modi possibili (selezionabili via OSC_USE_LEGACY):
        //   - NUOVO (default): n_pos_crossings >= OSC_NCROSS_MAX
        //   - LEGACY: vmax > OSC_VMAX_THRESH (vecchio criterio v5)
        // Mantenere il legacy permette confronti di sistematica.
        if (OSC_USE_LEGACY) {
            if (vmax > OSC_VMAX_THRESH) cd.is_oscillating = true;
        } else {
            if (n_cross >= OSC_NCROSS_MAX) cd.is_oscillating = true;
        }
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

    // ---- PASSO 8: TIME-OVER-THRESHOLD MULTI-SOGLIA ----
    //
    // Per ciascuna soglia q in TOT_THR_MV[]:
    //   1. cerca il PRIMO crossing sul fronte di salita (rise)
    //      tra NBL_SAMPLES e imin: u[i] < q ≤ u[i+1]
    //   2. cerca il PRIMO crossing sul fronte di discesa (fall)
    //      tra imin e ns-1: u[i] ≥ q > u[i+1]
    //   3. interpolazione lineare per timing sub-campione
    //   4. TOT(q) = t_fall − t_rise
    //
    // Variabile lavorata: u(t) = baseline − v(t), positiva durante
    // l'impulso (analoga alla a(t) del documento di riferimento TOT).
    //
    // VINCOLI DI VALIDITÀ:
    //   - q ≤ |CLIP_V_LO| − TOT_CLIP_MARGIN (altrimenti su clippati il
    //     fall sarebbe contaminato dal plateau saturato — non calcoliamo
    //     proprio il TOT a soglie troppo alte, neanche sui non-clippati,
    //     per coerenza)
    //   - amp ≥ q (il segnale supera la soglia)
    //   - rise e fall trovati entrambi
    //   - t_fall − t_rise ≥ TOT_NMIN_SAMPLES * sample_period
    //
    // RINGING / OSCILLAZIONI:
    //   Cerchiamo il PRIMO rise dopo NBL_SAMPLES e il PRIMO fall dopo
    //   imin. Se ci sono ulteriori crossing nella coda (afterpulse,
    //   overshoot post-clip), vengono ignorati per il TOT principale.
    //   L'evento è già flaggato is_oscillating dal PASSO 5 se la
    //   componente oscillante è patologica.
    {
        // Soglia massima accettabile in tensione, in valore u (positivo)
        const double q_max_safe = fabs(CLIP_V_LO) - TOT_CLIP_MARGIN;

        // Stima del periodo medio di campionamento (per il check sui
        // sample minimi). Usiamo la differenza tra il primo e l'ultimo
        // sample diviso (ns-1).
        double sample_period = (ns > 1) ?
                               ((double)t[ns - 1] - (double)t[0]) / (ns - 1) :
                               0.2;  // fallback 0.2 ns = 5 GS/s
        double dt_min_tot = TOT_NMIN_SAMPLES * sample_period;

        for (int kth = 0; kth < NTOT_THR; kth++) {
            double q = TOT_THR_MV[kth];

            // Skip se la soglia è troppo vicina al clipping
            if (q > q_max_safe) continue;
            // Skip se il segnale non raggiunge nemmeno la soglia
            if (amp < q)        continue;

            // ---- Crossing sul fronte di salita ----
            // Cerca [i, i+1] tale che u[i] < q ≤ u[i+1], i.e. che la
            // tensione v[i] > v_thr ≥ v[i+1] dove v_thr = bl - q.
            double v_thr = bl - q;
            double t_rise = -1.0;
            bool   rise_found = false;

            for (int i = NBL_SAMPLES; i < imin; i++) {
                if (v[i] > v_thr && v[i + 1] <= v_thr) {
                    double dv = (double)v[i + 1] - v[i];
                    if (fabs(dv) > 1e-6) {
                        // Stessa interpolazione del CFD del PASSO 6
                        t_rise = t[i] + (v_thr - v[i]) / dv *
                                 (t[i + 1] - t[i]);
                        rise_found = true;
                    }
                    break;
                }
            }
            if (!rise_found) continue;

            // ---- Crossing sul fronte di discesa ----
            // Cerca [i, i+1] DOPO imin tale che v[i] ≤ v_thr < v[i+1]
            // (equivalente: u[i] ≥ q > u[i+1])
            double t_fall = -1.0;
            bool   fall_found = false;

            for (int i = imin; i < ns - 1; i++) {
                if (v[i] <= v_thr && v[i + 1] > v_thr) {
                    double dv = (double)v[i + 1] - v[i];
                    if (fabs(dv) > 1e-6) {
                        t_fall = t[i] + (v_thr - v[i]) / dv *
                                 (t[i + 1] - t[i]);
                        fall_found = true;
                    }
                    break;
                }
            }
            if (!fall_found) continue;

            // ---- Calcolo TOT e validazione ----
            double tot_val = t_fall - t_rise;
            if (tot_val < dt_min_tot) continue;  // troppo stretto

            cd.tot_rise[kth] = t_rise;
            cd.tot_fall[kth] = t_fall;
            cd.tot_q   [kth] = tot_val;
            cd.tot_ok  [kth] = true;
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
    // Modello completo del fit: A = k·s + q (q tipicamente piccolo;
    // se q = 0 ricadiamo nel vecchio comportamento senza intercetta).
    double A_rec = gK_PMT[ch_index] * cd.slew_rate + gQ_PMT[ch_index];
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

/// RecoverClippedCFD_TOT():
///   Variante del recupero clippati basata sul Time-Over-Threshold
///   alla soglia di riferimento TOT_THR_MV[TOT_CALIB_THR_IDX].
///
/// PRECONDIZIONI:
///   - cd.is_clipped == true
///   - cd.tot_ok[TOT_CALIB_THR_IDX] == true (TOT misurato e valido)
///   - gTOT_calibrated == true e parametri PMT disponibili
///
/// MODELLO: A_rec = p0 · exp(p1 · TOT_ref)
///
/// CAP DI ESTRAPOLAZIONE:
///   Se A_rec > TOT_AREC_MAX_FACTOR × A_max_calib (per quel PMT),
///   l'estrapolazione è considerata non affidabile e il recupero
///   fallisce. Questo protegge contro outlier in cui TOT è
///   anomalamente lungo (pile-up, deriva di baseline).
///
/// VERIFICA SOGLIA CFD: identica a RecoverClippedCFD — la soglia
///   f·A_rec deve cadere SOPRA il livello di clipping con margine.

void RecoverClippedCFD_TOT(ChannelData &cd, int ch_index) {

    // --- Sanity check delle precondizioni ---
    if (!cd.is_clipped)                            return;
    if (ch_index < 0 || ch_index >= MAX_CHANNELS)  return;
    if (!gTOT_calibrated)                          return;
    if (gTOT_p1_PMT[ch_index] <= 0.0)              return;  // p1 non valido

    // Verifica che il TOT alla soglia di riferimento sia stato misurato
    const int kref = TOT_CALIB_THR_IDX;
    if (kref < 0 || kref >= NTOT_THR)              return;
    if (!cd.tot_ok[kref])                          return;

    double tot_ref = cd.tot_q[kref];

    // --- Stima dell'ampiezza ricostruita ---
    //   A_TOT = p0 · exp(p1 · TOT)
double p0 = gTOT_p0_PMT[ch_index];
    double p1 = gTOT_p1_PMT[ch_index];
    double A_tot = p0 * exp(p1 * tot_ref);
    cd.amplitude_tot = A_tot;

    // --- Sanità dell'ampiezza ricostruita ---
    // Richiediamo solo che A_tot sia fisicamente ragionevole:
    //   - A_tot > 0 (ovvio per il modello esponenziale)
    //   - A_tot non ecceda un cap di estrapolazione
    //
    // NON richiediamo più A_tot > cd.amplitude, perché per eventi
    // clippati cd.amplitude è troncata (~499 mV) e questa condizione
    // rigetta sistematicamente gli eventi borderline con ampiezza
    // vera di poco superiore al clipping (la maggioranza dei clippati).
    // La qualità è garantita dal check sulla soglia CFD (sotto).
    if (A_tot <= 0.0) return;

    // Cap di estrapolazione: protegge contro TOT anomalamente lunghi
    // (afterpulse, pile-up) che farebbero esplodere l'esponenziale.
    // Cap: A_tot < TOT_AREC_MAX_FACTOR × A_max_calibrazione
    double A_cap = TOT_AREC_MAX_FACTOR * gTOT_Amax_PMT[ch_index];
    if (A_cap > 0.0 && A_tot > A_cap) return;

    // --- Calcolo della soglia CFD ricostruita ---
    double v_thr = cd.baseline - CFD_FRACTION * A_tot;

    // Verifica CRITICA: la soglia CFD deve cadere SOPRA il livello di
    // clipping (con margine), altrimenti il crossing non esiste sulla
    // forma d'onda. Questo è il vero guardiano della qualità: se la
    // soglia è nel plateau, il tempo è sbagliato per costruzione.
    if (v_thr <= CLIP_V_LO + SLEW_CLIP_MARGIN) return;

    // --- Ricerca del crossing sul fronte di salita ---
    int   ns = cd.nsamples;
    float *t = cd.time;
    float *v = cd.voltage;

    // Trova l'indice del minimo (primo campione del plateau saturato)
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

// Aggiungi una copia visualizzabile della retta libera al TGraph.
        // Così quando il TGraph viene disegnato (o aperto dal TBrowser),
        // la retta del fit libero è visibile come overlay.
        TF1 *f_show = new TF1(Form("%s_show", name),
                              "[0] + [1]*x", s_min, s_max);
        f_show->SetParameters(q_out, k_out);
        f_show->SetLineColor(kRed + 1);
        f_show->SetLineWidth(3);
        gr->GetListOfFunctions()->Add(f_show);

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

        delete f_orig;
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

// ---- Salva nelle globali (usiamo il fit LIBERO con cleaning) ----
    // Il fit con intercetta è il modello più generale: se q ≈ 0, si
    // ricade nel caso del fit per l'origine; se q ≠ 0, l'intercetta
    // cattura offset sistematici (piede non lineare, bias di baseline).
    gK_PMT[0]      = k1;
    gK_PMT_err[0]  = k1_err;
    gK_PMT[1]      = k2;
    gK_PMT_err[1]  = k2_err;
    gQ_PMT[0]      = q1;
    gQ_PMT_err[0]  = q1_err;
    gQ_PMT[1]      = q2;
    gQ_PMT_err[1]  = q2_err;
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
    pt1->AddText("PMT1: A = k_{1} #upoint s + q_{1}");
    pt1->AddText(Form("k_{1} = %.3f #pm %.3f ns", k1, k1_err));
    pt1->AddText(Form("q_{1} = %.1f #pm %.1f mV", q1, q1_err));
    pt1->AddText(Form("#chi^{2}/ndf = %.2f", chi2_lin1));
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
    pt2->AddText("PMT2: A = k_{2} #upoint s + q_{2}");
    pt2->AddText(Form("k_{2} = %.3f #pm %.3f ns", k2, k2_err));
    pt2->AddText(Form("q_{2} = %.1f #pm %.1f mV", q2, q2_err));
    pt2->AddText(Form("#chi^{2}/ndf = %.2f", chi2_lin2));
    pt2->AddText(Form("N punti = %d", n_pts_PMT2));
    pt2->Draw();

    c_calib_k->Update();
    if (fout) {
        fout->cd();
    c_calib_k->Write("Calibration_k_slew_to_amplitude");
    }

    std::cout << "  [DONE] Calibrazione k completata." << std::endl;
    std::cout << "=============================================\n" << std::endl;
}

// ==========================================================================
//  SEZIONE 4c: CALIBRAZIONE TOT → AMPIEZZA
// ==========================================================================

/// CalibrateTOT():
///   Determina i parametri (p0, p1) della relazione
///         A = p0 · exp(p1 · TOT)
///   per ciascun PMT (PMT1 e PMT2), usando i dataset con
///   |x_k| ≤ CALIB_K_X_ABS_MAX e SOLO eventi non clippati.
///
/// MOTIVAZIONE FISICA: per impulso a forma fissa A·g(t) con coda
///   tipo esponenziale, il TOT a soglia q soddisfa
///         TOT ≈ τ_eff · ln(A/q)
///   da cui  A = q · exp(TOT/τ_eff) = p0 · exp(p1 · TOT)
///   con p0 ≈ q (la soglia stessa) e p1 = 1/τ_eff.
///
/// METODO: parsing degli stessi XML usati da CalibrateSlewK (doppio
///   parsing accettabile per modularità). Per ogni evento e ogni PMT,
///   verifica:
///     - has_pulse, !is_clipped, !is_oscillating
///     - tot_ok[TOT_CALIB_THR_IDX] (TOT alla soglia di riferimento valido)
///   Costruisce TGraph (TOT, A) e fitta con TF1 "[0]*exp([1]*x)".
///   Limita p1 > 0 per garantire monotonia crescente.
///
/// EFFETTI: gTOT_p0_PMT[0..1], gTOT_p1_PMT[0..1], gTOT_Amax_PMT[0..1]
///   e gTOT_calibrated vengono aggiornati. Salva due TGraph nel file
///   ROOT di output per ispezione visiva.

void CalibrateTOT(const char* folder,
                  const std::vector<PositionInfo> &positions,
                  TFile* fout) {

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  CALIBRAZIONE TOT → AMPIEZZA                 " << std::endl;
    std::cout << "  Soglia di riferimento: q = "
              << TOT_THR_MV[TOT_CALIB_THR_IDX] << " mV          " << std::endl;
    std::cout << "  Range posizioni: |x| ≤ "
              << CALIB_K_X_ABS_MAX << " cm                     " << std::endl;
    std::cout << "=============================================" << std::endl;

    const int    kref = TOT_CALIB_THR_IDX;
    const double q_ref = TOT_THR_MV[kref];

    // ---- TGraph (TOT, A) per ciascun PMT ----
    TGraph *g_PMT1 = new TGraph();   g_PMT1->SetName("g_calib_TOT_PMT1");
    TGraph *g_PMT2 = new TGraph();   g_PMT2->SetName("g_calib_TOT_PMT2");
    g_PMT1->SetTitle(Form("Calibrazione PMT1: A vs TOT(q=%.0f mV);"
                          "TOT [ns];Ampiezza A [mV]", q_ref));
    g_PMT2->SetTitle(Form("Calibrazione PMT2: A vs TOT(q=%.0f mV);"
                          "TOT [ns];Ampiezza A [mV]", q_ref));

    int n_pts_PMT1 = 0;
    int n_pts_PMT2 = 0;

    int n_files_used     = 0;
    int n_events_total   = 0;
    int n_events_kept_p1 = 0;
    int n_events_kept_p2 = 0;

    // Tracker dei massimi di A per il cap di estrapolazione
    double A_max_PMT1 = 0.0;
    double A_max_PMT2 = 0.0;

    for (const auto &pos : positions) {

        // Filtro |x| ≤ CALIB_K_X_ABS_MAX (stessa logica di CalibrateSlewK)
        if (fabs(pos.x_cm) > CALIB_K_X_ABS_MAX) continue;

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
        int n_parsed = ParseXML(xml_path.c_str(), events);
        if (n_parsed <= 0) continue;

        n_files_used++;
        n_events_total += (int)events.size();

        // ---- Loop sugli eventi ----
        for (size_t ev = 0; ev < events.size(); ev++) {
            const EventData &e = events[ev];
            if (e.nchannels < 2) continue;

            for (int k = 0; k < 2; k++) {
                const ChannelData &cd = e.ch[k];

                // Stessi criteri di CalibrateSlewK + TOT valido alla soglia ref
                if (!cd.has_pulse)        continue;
                if ( cd.is_clipped)       continue;
                if ( cd.is_oscillating)   continue;
                if (!cd.tot_ok[kref])     continue;

                if (k == 0) {
                    g_PMT1->SetPoint(n_pts_PMT1, cd.tot_q[kref], cd.amplitude);
                    n_pts_PMT1++;
                    n_events_kept_p1++;
                    if (cd.amplitude > A_max_PMT1) A_max_PMT1 = cd.amplitude;
                } else {
                    g_PMT2->SetPoint(n_pts_PMT2, cd.tot_q[kref], cd.amplitude);
                    n_pts_PMT2++;
                    n_events_kept_p2++;
                    if (cd.amplitude > A_max_PMT2) A_max_PMT2 = cd.amplitude;
                }
            }
        }
    }

    std::cout << "[INFO] CalibrateTOT: "
              << n_files_used << " dataset, "
              << n_events_total << " eventi totali" << std::endl;
    std::cout << "       PMT1: " << n_events_kept_p1
              << " punti utili (A_max = " << A_max_PMT1 << " mV)" << std::endl;
    std::cout << "       PMT2: " << n_events_kept_p2
              << " punti utili (A_max = " << A_max_PMT2 << " mV)" << std::endl;

    if (n_pts_PMT1 < 100 || n_pts_PMT2 < 100) {
        std::cerr << "[ERRORE] CalibrateTOT: troppi pochi punti per il fit. "
                  << "Calibrazione FALLITA, recupero TOT DISABILITATO."
                  << std::endl;
        gTOT_calibrated = false;
        return;
    }

    // ---- Fit esponenziale: A = p0·exp(p1·TOT) ----
    // ROOT TF1 "expo": exp([0]+[1]*x) → equivalente a A = exp(p0_root)·exp(p1·x)
    // Preferisco la forma esplicita per intuitività dei parametri.
// ---- Fit esponenziale con cleaning robusto: A = p0·exp(p1·TOT) ----
    //
    // MODELLO FISICO:
    //   Per impulso a forma fissa con coda esponenziale:
    //     TOT(q) ≈ τ_eff · ln(A/q)   →   A = q · exp(TOT/τ_eff)
    //   quindi p0 ≈ q (soglia), p1 = 1/τ_eff.
    //
    // PROBLEMI DI CONVERGENZA (v7):
    //   L'inizializzazione usava ComputeRange() sul TGraph grezzo,
    //   che include outlier a TOT > 50 ns (afterpulse, pile-up, crossing
    //   confuso sulla coda). Questo dava p1_init ~ 0.014 invece di ~0.08,
    //   e MINUIT convergeva nel minimo locale sbagliato.
    //
    // SOLUZIONE:
    //   1. Cleaning robusto (mediana + MAD) per rimuovere gli outlier
    //   2. Inizializzazione calcolata sulla nube pulita
    //   3. Taglio ulteriore in ampiezza (A troppo vicine al clip o alla
    //      soglia sono rumorose)
    auto FitOneTOT = [&](TGraph* gr, const char* name,
                         double &p0_out, double &p0_err_out,
                         double &p1_out, double &p1_err_out,
                         double &chi2ndf_out) {

        const int N = gr->GetN();

        // ---- Step 1: estrazione punti ----
        std::vector<double> tot_v(N), A_v(N);
        for (int i = 0; i < N; i++) {
            double xx, yy;
            gr->GetPoint(i, xx, yy);
            tot_v[i] = xx;
            A_v[i]   = yy;
        }

        // ---- Step 2: cleaning robusto su TOT (mediana + MAD) ----
        // La mediana e la MAD (Median Absolute Deviation) sono insensibili
        // agli outlier. Rimuoviamo punti con |TOT - mediana| > N_MAD·σ_MAD
        // dove σ_MAD = MAD · 1.4826 è l'equivalente gaussiano.
        std::vector<double> tot_sorted = tot_v;
        std::sort(tot_sorted.begin(), tot_sorted.end());
        double tot_median = tot_sorted[N / 2];

        std::vector<double> abs_dev(N);
        for (int i = 0; i < N; i++)
            abs_dev[i] = fabs(tot_v[i] - tot_median);
        std::sort(abs_dev.begin(), abs_dev.end());
        double MAD = abs_dev[N / 2];

        // Finestra di accettazione: 4·σ_MAD, con floor a 3 ns per evitare
        // cleaning eccessivo se la distribuzione è molto stretta
        const double N_MAD = 4.0;
        double sigma_MAD = MAD * 1.4826;
        double tot_window = std::max(N_MAD * sigma_MAD, 3.0);
        double tot_lo_cut = tot_median - tot_window;
        double tot_hi_cut = tot_median + tot_window;

        // Soglia A massima: rifiuto eventi vicini al clip (forma deformata)
        const double A_MAX_CLEAN = fabs(CLIP_V_LO) - 20.0;  // ~480 mV
        // Soglia A minima: rifiuto eventi appena sopra soglia TOT (rumorosi)
        const double A_MIN_CLEAN = TOT_THR_MV[TOT_CALIB_THR_IDX] + 10.0;
        // TOT minimo: almeno 1 ns (qualche campione)
        const double TOT_MIN_CLEAN = 1.0;

        // ---- Step 3: costruzione TGraph filtrato ----
        TGraph *gr_clean = new TGraph();
        gr_clean->SetName(Form("%s_clean", gr->GetName()));
        int n_clean = 0;
        for (int i = 0; i < N; i++) {
            if (tot_v[i] < tot_lo_cut || tot_v[i] > tot_hi_cut) continue;
            if (tot_v[i] < TOT_MIN_CLEAN)                       continue;
            if (A_v[i]   < A_MIN_CLEAN || A_v[i] > A_MAX_CLEAN) continue;
            gr_clean->SetPoint(n_clean, tot_v[i], A_v[i]);
            n_clean++;
        }

        std::cout << "  [" << name << "] cleaning: "
                  << N << " → " << n_clean << " punti "
                  << "(TOT median = " << Form("%.1f", tot_median)
                  << " ns, MAD = " << Form("%.1f", MAD)
                  << " ns, finestra [" << Form("%.1f", tot_lo_cut)
                  << ", " << Form("%.1f", tot_hi_cut) << "] ns)"
                  << std::endl;

        if (n_clean < 50) {
            std::cerr << "  [WARNING] " << name << ": troppi pochi punti "
                      << "dopo cleaning (" << n_clean << ")." << std::endl;
            p0_out = 0; p0_err_out = 0;
            p1_out = 0; p1_err_out = 0;
            chi2ndf_out = -1;
            delete gr_clean;
            return;
        }

        // ---- Step 4: calcolo inizializzazione dalla nube pulita ----
        // Usiamo i quantili 10% e 90% della nube pulita per stimare
        // la pendenza esponenziale. Questo è robusto contro i bordi
        // rumorosi della distribuzione.
        std::vector<double> tot_clean(n_clean), A_clean(n_clean);
        for (int i = 0; i < n_clean; i++) {
            double xx, yy;
            gr_clean->GetPoint(i, xx, yy);
            tot_clean[i] = xx;
            A_clean[i]   = yy;
        }

        // Ordina per TOT per estrarre i quantili
        std::vector<int> idx(n_clean);
        for (int i = 0; i < n_clean; i++) idx[i] = i;
        std::sort(idx.begin(), idx.end(),
                  [&tot_clean](int a, int b) {
                      return tot_clean[a] < tot_clean[b];
                  });

        // Quantile 10% e 90% (per TOT e corrispondente A mediana locale)
        int i10 = n_clean / 10;
        int i90 = n_clean - n_clean / 10 - 1;
        // Media locale di A attorno a ogni quantile (±2% del campione)
        auto localMeanA = [&](int icenter, int halfwin) -> double {
            double sum_a = 0.0;
            int count = 0;
            for (int j = std::max(0, icenter - halfwin);
                     j <= std::min(n_clean - 1, icenter + halfwin); j++) {
                sum_a += A_clean[idx[j]];
                count++;
            }
            return (count > 0) ? sum_a / count : A_clean[idx[icenter]];
        };
        int hw = std::max(n_clean / 50, 5);  // mezza-finestra: ~2% del campione
        double TOT_lo = tot_clean[idx[i10]];
        double TOT_hi = tot_clean[idx[i90]];
        double A_lo   = localMeanA(i10, hw);
        double A_hi   = localMeanA(i90, hw);

        // Stima di p1 = 1/τ_eff dalla coppia di punti estremi
        // Da A = p0·exp(p1·TOT):
        //   ln(A_hi/A_lo) = p1·(TOT_hi - TOT_lo)
        //   → p1 = ln(A_hi/A_lo) / (TOT_hi - TOT_lo)
        double ratio = std::max(A_hi, 1.0) / std::max(A_lo, 1.0);
        double dTOT  = TOT_hi - TOT_lo;
        double p1_init = (dTOT > 1e-3 && ratio > 1.0) ?
                         log(ratio) / dTOT : 0.08;

        // Stima di p0 dal punto basso: p0 = A_lo · exp(-p1·TOT_lo)
        double p0_init = A_lo * exp(-p1_init * TOT_lo);
        // Sanity: p0 deve essere > 0 e ragionevole (non più di 2×q_ref)
        if (p0_init <= 0 || p0_init > 500.0)
            p0_init = q_ref;
        // Sanity: p1 deve dare τ_eff tra 2 ns e 100 ns
        if (p1_init < 0.01 || p1_init > 0.5)
            p1_init = 0.08;

        std::cout << "  [" << name << "] init: p0 = "
                  << Form("%.1f", p0_init) << " mV, p1 = "
                  << Form("%.4f", p1_init) << " ns^-1 (τ_eff = "
                  << Form("%.1f", 1.0 / p1_init) << " ns)"
                  << std::endl;

        // ---- Step 5: fit esponenziale sul campione pulito ----
        double tot_min_c, tot_max_c, A_min_c, A_max_c;
        gr_clean->ComputeRange(tot_min_c, A_min_c, tot_max_c, A_max_c);

        TF1 *f_exp = new TF1(Form("%s_exp", name),
                             "[0]*TMath::Exp([1]*x)",
                             tot_min_c, tot_max_c);
        f_exp->SetParName(0, "p0");
        f_exp->SetParName(1, "p1");
        f_exp->SetParameters(p0_init, p1_init);

        // Vincoli: p0 > 0, p1 > 0 (monotonia crescente), con range
        // fisicamente ragionevoli per evitare che MINUIT vada in regioni
        // esotiche.
        f_exp->SetParLimits(0, 1.0,   500.0);     // p0: 1–500 mV
        f_exp->SetParLimits(1, 0.005,   0.5);      // p1: τ_eff tra 2 e 200 ns

        // Fit con χ² standard (errori non specificati → unitari, accettabile
        // perché vogliamo la curva media, non gli errori sui parametri).
        // "M" = IMPROVE: chiede a MINUIT di cercare più a fondo il minimo
        // globale (utile se ci sono minimi locali).
        gr_clean->Fit(f_exp, "Q R N M");

        p0_out      = f_exp->GetParameter(0);
        p0_err_out  = f_exp->GetParError(0);
        p1_out      = f_exp->GetParameter(1);
        p1_err_out  = f_exp->GetParError(1);
        chi2ndf_out = (f_exp->GetNDF() > 0) ?
                      f_exp->GetChisquare() / (double)f_exp->GetNDF() : -1.0;

        // ---- Step 6: aggancio la curva al TGraph originale ----
        // Così nel canvas si vedono TUTTI i punti (anche gli outlier
        // rifiutati) con la curva sovrapposta: controllo visivo immediato.
        TF1 *f_show = new TF1(Form("%s_show", name),
                              "[0]*TMath::Exp([1]*x)",
                              tot_min_c, tot_max_c);
        f_show->SetParameters(p0_out, p1_out);
        f_show->SetLineColor(kGreen + 2);
        f_show->SetLineWidth(3);
        gr->GetListOfFunctions()->Add(f_show);

        delete gr_clean;
        delete f_exp;
    };

    double p0_1, p0_1_err, p1_1, p1_1_err, chi2_1;
    double p0_2, p0_2_err, p1_2, p1_2_err, chi2_2;
    FitOneTOT(g_PMT1, "fit_TOT_PMT1", p0_1, p0_1_err, p1_1, p1_1_err, chi2_1);
    FitOneTOT(g_PMT2, "fit_TOT_PMT2", p0_2, p0_2_err, p1_2, p1_2_err, chi2_2);

    // ---- Salva nelle globali ----
    gTOT_p0_PMT[0]   = p0_1;
    gTOT_p1_PMT[0]   = p1_1;
    gTOT_Amax_PMT[0] = A_max_PMT1;
    gTOT_p0_PMT[1]   = p0_2;
    gTOT_p1_PMT[1]   = p1_2;
    gTOT_Amax_PMT[1] = A_max_PMT2;
    gTOT_calibrated  = true;

    // ---- Stampa diagnostica ----
    std::cout << "\n  RISULTATI CALIBRAZIONE TOT:" << std::endl;
    std::cout << "  ----------------------------------------" << std::endl;
    std::cout << "  PMT1: A = ("
              << Form("%.2f ± %.2f", p0_1, p0_1_err) << ") · exp(("
              << Form("%.4f ± %.4f", p1_1, p1_1_err) << ") · TOT)"
              << "   χ²/ndf=" << Form("%.2f", chi2_1)
              << "   τ_eff = " << Form("%.2f ns", 1.0 / p1_1) << std::endl;
    std::cout << "  PMT2: A = ("
              << Form("%.2f ± %.2f", p0_2, p0_2_err) << ") · exp(("
              << Form("%.4f ± %.4f", p1_2, p1_2_err) << ") · TOT)"
              << "   χ²/ndf=" << Form("%.2f", chi2_2)
              << "   τ_eff = " << Form("%.2f ns", 1.0 / p1_2) << std::endl;

    // ---- Canvas di output ----
    TCanvas *c_calib_tot = new TCanvas("c_calib_tot",
                                        "Calibrazione TOT → A", 1200, 500);
    c_calib_tot->Divide(2, 1);

    auto DrawOne = [&](TVirtualPad* pad, TGraph* gr, int color,
                       const char* tag, double p0, double p1,
                       double chi2, int nptscur) {
        pad->cd();
        gPad->SetGrid(1, 1);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.4);
        gr->SetMarkerColor(color);
        gr->Draw("AP");

        TPaveText *pt = new TPaveText(0.15, 0.62, 0.55, 0.88, "NDC");
        pt->SetFillColorAlpha(kWhite, 0.85);
        pt->SetTextFont(42);
        pt->SetTextSize(0.04);
        pt->SetTextAlign(12);
        pt->AddText(Form("%s: A = p_{0} #upoint exp(p_{1} #upoint TOT)", tag));
        pt->AddText(Form("p_{0} = %.2f mV", p0));
        pt->AddText(Form("p_{1} = %.4f ns^{-1}  (#tau_{eff} = %.2f ns)",
                         p1, 1.0/p1));
        pt->AddText(Form("#chi^{2}/ndf = %.2f", chi2));
        pt->AddText(Form("N punti = %d", nptscur));
        pt->Draw();
    };

    DrawOne(c_calib_tot->cd(1), g_PMT1, kBlue + 1, "PMT1",
            p0_1, p1_1, chi2_1, n_pts_PMT1);
    DrawOne(c_calib_tot->cd(2), g_PMT2, kRed + 1, "PMT2",
            p0_2, p1_2, chi2_2, n_pts_PMT2);

    c_calib_tot->Update();
    if (fout) {
        fout->cd();
        c_calib_tot->Write("Calibration_TOT_to_amplitude");
    }

    std::cout << "  [DONE] Calibrazione TOT completata." << std::endl;
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
    FitResult dt12;         // Risultato fit Δt₁₂ = t₁ − t₂  (metodo slew)
    FitResult dt13;         // Risultato fit Δt₁₃ = t₁ − t₃  (metodo slew)
    FitResult dt23;         // Risultato fit Δt₂₃ = t₂ − t₃  (metodo slew)
    FitResult C_fit;        // Risultato fit C = t₃ − (t₁+t₂)/2 (metodo slew)
    int       n_good;       // Numero di eventi buoni (metodo slew)

    // Risultati paralleli per il metodo TOT
    FitResult dt12_tot;     // Risultato fit Δt₁₂ (metodo TOT)
    FitResult dt13_tot;     // Risultato fit Δt₁₃ (metodo TOT)
    FitResult dt23_tot;     // Risultato fit Δt₂₃ (metodo TOT)
    FitResult C_fit_tot;    // Risultato fit C    (metodo TOT)
    int       n_good_tot;   // Numero di eventi buoni (metodo TOT)
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
// ---- Istogrammi PARALLELI per il metodo TOT ----
    // Stessi binning e range degli istogrammi slew, con suffisso "_tot"
    // per evitare conflitti di nomi nel file ROOT.
    TH1D *h_dt12_tot = new TH1D(Form("h_dt12_tot_%s", tree_name.c_str()),
                                 Form("#Delta t_{12} (TOT) @ %s;"
                                      "#Delta t_{12} = t_{1} - t_{2} [ns];Conteggi",
                                      tree_name.c_str()),
                                 NBINS_DT, DT_LO, DT_HI);

    TH1D *h_dt13_tot = new TH1D(Form("h_dt13_tot_%s", tree_name.c_str()),
                                 Form("#Delta t_{13} (TOT) @ %s;"
                                      "#Delta t_{13} = t_{1} - t_{3} [ns];Conteggi",
                                      tree_name.c_str()),
                                 NBINS_DT, DT_LO, DT_HI);

    TH1D *h_dt23_tot = new TH1D(Form("h_dt23_tot_%s", tree_name.c_str()),
                                 Form("#Delta t_{23} (TOT) @ %s;"
                                      "#Delta t_{23} = t_{2} - t_{3} [ns];Conteggi",
                                      tree_name.c_str()),
                                 NBINS_DT, DT_LO, DT_HI);

    TH1D *h_C_tot = new TH1D(Form("h_C_tot_%s", tree_name.c_str()),
                              Form("C (TOT) = t_{3} - (t_{1}+t_{2})/2 @ %s;"
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

    // Branch numero di crossing positivi per ogni canale (diagnostica
    // del nuovo algoritmo Schmitt trigger). Permette analisi del tipo:
    //   tree->Draw("n_cross[1]", "amp[1] > 440")  → distribuzione dei
    //   crossing su PMT2 per i clippati. Aiuta a validare la soglia
    //   OSC_NCROSS_MAX e a capire la natura degli scarti.
    Int_t n_cross[3];
    tree->Branch("n_cross", n_cross, "n_cross[3]/I");

    // ---- Branch del METODO TOT (parallelo allo slew) ----
    // amp_tot[k]      = ampiezza ricostruita via TOT [mV] (0 se non recuperato)
    // t_cfd_tot[k]    = tempo CFD ricalcolato sulla soglia f·A_tot [ns]
    // cfd_tot_ok[k]   = 1 se il CFD via TOT è valido per il canale k
    // tot_ref[k]      = TOT [ns] alla soglia di riferimento (diagnostica)
    // dt12_tot, dt13_tot, dt23_tot, C_tot = differenze temporali calcolate
    //                  con il "metodo TOT" — equivalenti a dt12/dt13/dt23/C
    //                  ma con i tempi dei clippati corretti via TOT invece
    //                  che via slew rate.
    // good_tot        = 1 se l'evento è utilizzabile col metodo TOT
    Float_t amp_tot[3], t_cfd_tot[3], tot_ref[3];
    Int_t   cfd_tot_ok[3];
    Float_t dt12_tot, dt13_tot, dt23_tot, C_tot;
    Int_t   good_tot;
    tree->Branch("amp_tot",    amp_tot,    "amp_tot[3]/F");
    tree->Branch("t_cfd_tot",  t_cfd_tot,  "t_cfd_tot[3]/F");
    tree->Branch("cfd_tot_ok", cfd_tot_ok, "cfd_tot_ok[3]/I");
    tree->Branch("tot_ref",    tot_ref,    "tot_ref[3]/F");
    tree->Branch("dt12_tot",   &dt12_tot,  "dt12_tot/F");
    tree->Branch("dt13_tot",   &dt13_tot,  "dt13_tot/F");
    tree->Branch("dt23_tot",   &dt23_tot,  "dt23_tot/F");
    tree->Branch("C_tot",      &C_tot,     "C_tot/F");
    tree->Branch("good_tot",   &good_tot,  "good_tot/I");
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
        // Eseguiamo IN PARALLELO entrambi i metodi di recupero (slew rate
        // e TOT). I risultati vanno in campi separati della struct
        // (amplitude_rec/t_cfd_recovered per slew; amplitude_tot/
        // t_cfd_rec_tot per TOT) e sono confrontabili evento per evento
        // a posteriori sul TTree.
        // Solo PMT1 e PMT2 (indici 0,1) vengono recuperati. PMT3 non clippa.
        EventData &e_mod = events[ev];
        for (int k = 0; k < 2; k++) {
            if (e_mod.ch[k].is_clipped) {
                RecoverClippedCFD    (e_mod.ch[k], k);   // metodo slew rate
                RecoverClippedCFD_TOT(e_mod.ch[k], k);   // metodo TOT
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

        // Buffer per il "tempo migliore TOT": parallelo a t_used[]
        Float_t t_used_tot[3];
        bool    used_ok_tot[3];

        const int kref = TOT_CALIB_THR_IDX;   // indice soglia TOT di riferimento

        for (int k = 0; k < 3; k++) {
            const ChannelData &cd = e_mod.ch[k];

            // ---- Branch del metodo SLEW (invariati) ----
            t_cfd[k]   = (Float_t)cd.t_cfd;
            amp[k]     = (Float_t)cd.amplitude;
            bl[k]      = (Float_t)cd.baseline;
            amp_rec[k] = (Float_t)cd.amplitude_rec;
            n_cross[k] = cd.n_pos_crossings;

            // ---- Branch del metodo TOT (nuovi) ----
            amp_tot[k]   = (Float_t)cd.amplitude_tot;          // 0 se non rec.
            t_cfd_tot[k] = (Float_t)(cd.cfd_recovered_tot ?
                                     cd.t_cfd_rec_tot :
                                     cd.t_cfd);                // se non clip. uso CFD
            tot_ref[k]   = (cd.tot_ok[kref]) ?
                           (Float_t)cd.tot_q[kref] : 0.0f;     // diagnostica

            if (cd.is_clipped)     any_clipped     = true;
            if (cd.is_oscillating) any_oscillating = true;

            // ---- Selezione tempo "migliore" SLEW (logica originale) ----
            if (!cd.is_clipped) {
                t_used[k]  = (Float_t)cd.t_cfd;
                used_ok[k] = cd.cfd_ok;
            } else if (cd.cfd_recovered) {
                t_used[k]       = (Float_t)cd.t_cfd_recovered;
                used_ok[k]      = true;
                ch_recovered[k] = true;
                any_recovered   = true;
            } else {
                t_used[k]  = -999.0f;
                used_ok[k] = false;
            }

            // ---- Selezione tempo "migliore" TOT (parallela) ----
            // Stessa logica ma usa t_cfd_rec_tot al posto di t_cfd_recovered
            // sui clippati. Per i non clippati, usiamo lo stesso CFD standard.
            if (!cd.is_clipped) {
                t_used_tot[k]  = (Float_t)cd.t_cfd;
                used_ok_tot[k] = cd.cfd_ok;
                cfd_tot_ok[k]  = cd.cfd_ok ? 1 : 0;
            } else if (cd.cfd_recovered_tot) {
                t_used_tot[k]  = (Float_t)cd.t_cfd_rec_tot;
                used_ok_tot[k] = true;
                cfd_tot_ok[k]  = 1;
            } else {
                t_used_tot[k]  = -999.0f;
                used_ok_tot[k] = false;
                cfd_tot_ok[k]  = 0;
            }
        }

        bool all_t_ok_tot = used_ok_tot[0] && used_ok_tot[1] && used_ok_tot[2];
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
            // (metodo slew rate: t_cfd standard o t_cfd_recovered).
            dt12  = t_used[0] - t_used[1];
            dt13  = t_used[0] - t_used[2];
            dt23  = t_used[1] - t_used[2];
            C_val = t_used[2] - (t_used[0] + t_used[1]) / 2.0;

            good_event = 1;
            result.n_good++;

            // Riempi gli istogrammi (basati sul metodo slew, come prima)
            h_dt12->Fill(dt12);
            h_dt13->Fill(dt13);
            h_dt23->Fill(dt23);
            h_C->Fill(C_val);

        } else {
            dt12 = -999; dt13 = -999; dt23 = -999; C_val = -999;
            good_event = 0;
        }

        // ---- Calcolo PARALLELO con il metodo TOT ----
        // Nota: il flag good_tot è indipendente da good_event. Un evento
        // potrebbe essere "buono" col metodo slew ma non col TOT (o viceversa)
        // se il recupero fallisce in modi diversi. Il filtro di oscillazione
        // si applica a entrambi.
        bool good_quality_tot = true;
        if (ENABLE_CFD_CUT && !all_t_ok_tot)    good_quality_tot = false;
        if (ENABLE_OSC_CUT && any_oscillating)  good_quality_tot = false;

        if (good_quality_tot) {
            dt12_tot = t_used_tot[0] - t_used_tot[1];
            dt13_tot = t_used_tot[0] - t_used_tot[2];
            dt23_tot = t_used_tot[1] - t_used_tot[2];
            C_tot    = t_used_tot[2] - (t_used_tot[0] + t_used_tot[1]) / 2.0;
            good_tot = 1;

            // Riempi gli istogrammi TOT (paralleli a quelli slew)
            h_dt12_tot->Fill(dt12_tot);
            h_dt13_tot->Fill(dt13_tot);
            h_dt23_tot->Fill(dt23_tot);
            h_C_tot->Fill(C_tot);

        } else {
            dt12_tot = -999; dt13_tot = -999; dt23_tot = -999; C_tot = -999;
            good_tot = 0;
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
              << " eventi buoni (slew)"
              << " (di cui " << n_recovered << " recuperati via slew rate)"
              << std::endl;
    std::cout << "       TOT: "
              << result.n_good_tot << "/" << events.size()
              << " eventi buoni (metodo TOT)" << std::endl;
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

    // ---- Fit degli istogrammi TOT (paralleli) ----
    result.dt12_tot  = FitGaussianCore(h_dt12_tot, 0.90,
                                       Form("fit_dt12_tot_%s", tree_name.c_str()));
    result.dt13_tot  = FitGaussianCore(h_dt13_tot, 0.90,
                                       Form("fit_dt13_tot_%s", tree_name.c_str()));
    result.dt23_tot  = FitGaussianCore(h_dt23_tot, 0.90,
                                       Form("fit_dt23_tot_%s", tree_name.c_str()));
    result.C_fit_tot = FitGaussianCore(h_C_tot, 0.90,
                                       Form("fit_C_tot_%s", tree_name.c_str()));

    // Contatore degli eventi buoni per il metodo TOT
    // (lo calcoliamo qui perché il contatore nel loop eventi non era
    // stato accumulato in result — usiamo gli entries dell'istogramma)
    result.n_good_tot = (int)h_dt12_tot->GetEntries();

    // ---- Salvataggio ----
    outfile->cd();
    tree->Write();
    h_dt12->Write();
    h_dt13->Write();
    h_dt23->Write();
    h_C->Write();
    h_dt12_tot->Write();
    h_dt13_tot->Write();
    h_dt23_tot->Write();
    h_C_tot->Write();

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
    // Metodo slew (originale)
    std::vector<double> x_pos;
    std::vector<double> dx_pos;
    std::vector<double> dt12_mean;
    std::vector<double> dt12_err;
    std::vector<double> dt12_width;
    std::vector<double> dt12_werr;
    std::vector<double> C_mean;
    std::vector<double> C_err;
    std::vector<std::string> labels;

    // Metodo TOT (parallelo)
    std::vector<double> dt12_mean_tot;
    std::vector<double> dt12_err_tot;
    std::vector<double> dt12_width_tot;
    std::vector<double> dt12_werr_tot;
    std::vector<double> C_mean_tot;
    std::vector<double> C_err_tot;

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

    // Calibrazione parallela: relazione A vs TOT alla soglia di
    // riferimento. Stessa logica di selezione (|x| ≤ CALIB_K_X_ABS_MAX,
    // solo non clippati) ma usa il TOT come variabile predittiva.
    // Permette ricostruzione dei clippati indipendente dallo slew rate
    // → cross-check fra i due metodi a livello di TTree.
    CalibrateTOT(folder, DEFAULT_POSITIONS, fout);

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
        double sx_k = ComputeSigmaX(pos.x_cm, effective_parallax_file);
        dx_pos.push_back(sx_k);
        dt12_mean.push_back(dres.dt12.center);
        dt12_err.push_back(dres.dt12.center_err);
        dt12_width.push_back(dres.dt12.width);
        dt12_werr.push_back(dres.dt12.width_err);
        C_mean.push_back(dres.C_fit.center);
        C_err.push_back(dres.C_fit.center_err);
        labels.push_back(pos.name);

        // Risultati paralleli TOT
        dt12_mean_tot.push_back(dres.dt12_tot.center);
        dt12_err_tot.push_back(dres.dt12_tot.center_err);
        dt12_width_tot.push_back(dres.dt12_tot.width);
        dt12_werr_tot.push_back(dres.dt12_tot.width_err);
        C_mean_tot.push_back(dres.C_fit_tot.center);
        C_err_tot.push_back(dres.C_fit_tot.center_err);

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
                  << ": x = " << pos.x_cm << " cm" << std::endl;
        std::cout << "   SLEW:  Δt₁₂ = " << Form("%.3f ± %.3f", dres.dt12.center, dres.dt12.center_err)
                  << " ns, σ = " << Form("%.3f", dres.dt12.width) << " ns"
                  << ", C = " << Form("%.3f ± %.3f", dres.C_fit.center, dres.C_fit.center_err)
                  << " ns  (N=" << dres.n_good << ")" << std::endl;
        std::cout << "   TOT:   Δt₁₂ = " << Form("%.3f ± %.3f", dres.dt12_tot.center, dres.dt12_tot.center_err)
                  << " ns, σ = " << Form("%.3f", dres.dt12_tot.width) << " ns"
                  << ", C = " << Form("%.3f ± %.3f", dres.C_fit_tot.center, dres.C_fit_tot.center_err)
                  << " ns  (N=" << dres.n_good_tot << ")" << std::endl;
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
// Fit lineare C(x) = C_0 + C_1·x (sui dati TOT per consistenza).
    // Dichiarati qui a livello di funzione perché servono sia nel canvas
    // di confronto Slew/TOT sia nel TTree fit_params.
    double C_lin_p0 = 0.0, C_lin_p0_err = 0.0;
    double C_lin_p1 = 0.0, C_lin_p1_err = 0.0;
    double C_lin_chi2ndf = -1.0;

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
    //  CANVAS 6: CONFRONTO Δt₁₂ vs x — SLEW vs TOT
    // ==================================================================
    // Due rette di calibrazione sovrapposte: una per il metodo slew
    // (punti blu, fit blu) e una per il metodo TOT (punti rossi, fit
    // rosso). Idealmente dovrebbero coincidere; differenze sistematiche
    // segnalano bias nel recupero dei clippati.
    {
        TCanvas *c_cmp = new TCanvas("c_dt12_compare",
                                      "#Delta t_{12} vs x: confronto Slew vs TOT",
                                      1000, 800);
        c_cmp->Divide(1, 2);

        // --- Pannello superiore: dati e fit ---
        c_cmp->cd(1);
        gPad->SetPad(0.0, 0.35, 1.0, 1.0);
        gPad->SetBottomMargin(0.02);
        gPad->SetTopMargin(0.08);
        gPad->SetGrid(1, 1);

        // Grafico SLEW (invariato, lo ridisegniamo)
        TGraphErrors *g_slew = new TGraphErrors(npts,
                                                x_pos.data(), dt12_mean.data(),
                                                dx_pos.data(), dt12_err.data());
        g_slew->SetMarkerStyle(20);
        g_slew->SetMarkerSize(1.0);
        g_slew->SetMarkerColor(kBlue + 1);
        g_slew->SetLineColor(kBlue + 1);
        g_slew->SetTitle(";Posizione x [cm];#Delta t_{12} [ns]");
        g_slew->GetXaxis()->SetLabelSize(0);

        TF1 *f_slew = new TF1("f_slew_cmp", "[0]+[1]*x", -150, 150);
        f_slew->SetLineColor(kBlue + 1);
        f_slew->SetLineWidth(2);
        g_slew->Fit(f_slew, "S R Q N");

        g_slew->Draw("AP");
        f_slew->Draw("same");

        // Grafico TOT (sovrapposto)
        TGraphErrors *g_tot_cmp = new TGraphErrors(npts,
                                                   x_pos.data(), dt12_mean_tot.data(),
                                                   dx_pos.data(), dt12_err_tot.data());
        g_tot_cmp->SetMarkerStyle(21);
        g_tot_cmp->SetMarkerSize(1.0);
        g_tot_cmp->SetMarkerColor(kRed + 1);
        g_tot_cmp->SetLineColor(kRed + 1);

        TF1 *f_tot_cmp = new TF1("f_tot_cmp", "[0]+[1]*x", -150, 150);
        f_tot_cmp->SetLineColor(kRed + 1);
        f_tot_cmp->SetLineWidth(2);
        f_tot_cmp->SetLineStyle(2);   // tratteggiato per distinguerlo
        g_tot_cmp->Fit(f_tot_cmp, "S R Q N");

        g_tot_cmp->Draw("P same");
        f_tot_cmp->Draw("same");

        // Velocità derivate
        double v_slew     = 2.0 / fabs(f_slew->GetParameter(1));
        double v_slew_err = 2.0 * f_slew->GetParError(1) /
                            (f_slew->GetParameter(1) * f_slew->GetParameter(1));
        double v_tot_cal  = 2.0 / fabs(f_tot_cmp->GetParameter(1));
        double v_tot_err  = 2.0 * f_tot_cmp->GetParError(1) /
                            (f_tot_cmp->GetParameter(1) * f_tot_cmp->GetParameter(1));

        TPaveText *pt_cmp = new TPaveText(0.15, 0.55, 0.58, 0.92, "NDC");
        pt_cmp->SetFillColorAlpha(kWhite, 0.85);
        pt_cmp->SetLineColor(kBlack);
        pt_cmp->SetTextFont(42);
        pt_cmp->SetTextSize(0.04);
        pt_cmp->SetTextAlign(12);
        pt_cmp->AddText("#color[4]{Slew rate (cerchi pieni)}");
        pt_cmp->AddText(Form("#color[4]{v_{eff} = %.2f #pm %.2f cm/ns}",
                             v_slew, v_slew_err));
        pt_cmp->AddText(Form("#color[4]{#chi^{2}/ndf = %.2f}",
                             f_slew->GetChisquare() / std::max(f_slew->GetNDF(), 1)));
        pt_cmp->AddText("#color[2]{TOT (quadrati pieni)}");
        pt_cmp->AddText(Form("#color[2]{v_{eff} = %.2f #pm %.2f cm/ns}",
                             v_tot_cal, v_tot_err));
        pt_cmp->AddText(Form("#color[2]{#chi^{2}/ndf = %.2f}",
                             f_tot_cmp->GetChisquare() / std::max(f_tot_cmp->GetNDF(), 1)));
        pt_cmp->Draw();

        // Legenda
        TLegend *leg = new TLegend(0.60, 0.75, 0.90, 0.92);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);
        leg->AddEntry(g_slew,    "Slew rate", "p");
        leg->AddEntry(g_tot_cmp, "TOT",       "p");
        leg->Draw();

        // --- Pannello inferiore: residui TOT − Slew (punto per punto) ---
        c_cmp->cd(2);
        gPad->SetPad(0.0, 0.0, 1.0, 0.35);
        gPad->SetTopMargin(0.02);
        gPad->SetBottomMargin(0.25);
        gPad->SetGrid(1, 1);

        std::vector<double> diff(npts), diff_err(npts);
        for (int i = 0; i < npts; i++) {
            diff[i]     = dt12_mean_tot[i] - dt12_mean[i];
            diff_err[i] = sqrt(dt12_err_tot[i] * dt12_err_tot[i] +
                               dt12_err[i] * dt12_err[i]);
        }
        TGraphErrors *g_diff = new TGraphErrors(npts,
                                                x_pos.data(), diff.data(),
                                                dx_pos.data(), diff_err.data());
        g_diff->SetMarkerStyle(20);
        g_diff->SetMarkerSize(1.0);
        g_diff->SetMarkerColor(kMagenta + 1);
        g_diff->SetLineColor(kMagenta + 1);
        g_diff->SetTitle(";Posizione x [cm];"
                         "#Delta t_{12}^{TOT} - #Delta t_{12}^{slew} [ns]");
        g_diff->GetXaxis()->SetTitleSize(0.08);
        g_diff->GetXaxis()->SetLabelSize(0.07);
        g_diff->GetYaxis()->SetTitleSize(0.08);
        g_diff->GetYaxis()->SetLabelSize(0.07);
        g_diff->GetYaxis()->SetTitleOffset(0.5);
        g_diff->Draw("AP");

        double xmin_d = *std::min_element(x_pos.begin(), x_pos.end()) - 10;
        double xmax_d = *std::max_element(x_pos.begin(), x_pos.end()) + 10;
        TLine *l0 = new TLine(xmin_d, 0, xmax_d, 0);
        l0->SetLineColor(kBlack);
        l0->SetLineStyle(2);
        l0->Draw("same");

        c_cmp->Update();
        fout->cd();
        c_cmp->Write("Compare_Dt12_Slew_vs_TOT");
    }

    // ==================================================================
    //  CANVAS 7: CONFRONTO RISOLUZIONE σ(Δt₁₂) vs x — SLEW vs TOT
    // ==================================================================
    {
        TCanvas *c_res_cmp = new TCanvas("c_resolution_compare",
                                          "#sigma(#Delta t_{12}) vs x: Slew vs TOT",
                                          800, 500);
        c_res_cmp->SetGrid(1, 1);

        TGraphErrors *g_w_slew = new TGraphErrors(npts,
                                                  x_pos.data(), dt12_width.data(),
                                                  dx_pos.data(), dt12_werr.data());
        g_w_slew->SetMarkerStyle(20);
        g_w_slew->SetMarkerSize(1.0);
        g_w_slew->SetMarkerColor(kBlue + 1);
        g_w_slew->SetLineColor(kBlue + 1);
        g_w_slew->SetTitle("Risoluzione temporale: confronto;"
                           "Posizione x [cm];"
                           "#sigma(#Delta t_{12}) [ns]");
        g_w_slew->Draw("AP");

        TGraphErrors *g_w_tot = new TGraphErrors(npts,
                                                 x_pos.data(), dt12_width_tot.data(),
                                                 dx_pos.data(), dt12_werr_tot.data());
        g_w_tot->SetMarkerStyle(21);
        g_w_tot->SetMarkerSize(1.0);
        g_w_tot->SetMarkerColor(kRed + 1);
        g_w_tot->SetLineColor(kRed + 1);
        g_w_tot->Draw("P same");

        TLegend *leg2 = new TLegend(0.60, 0.75, 0.90, 0.92);
        leg2->SetTextFont(42);
        leg2->SetTextSize(0.04);
        leg2->AddEntry(g_w_slew, "Slew rate", "p");
        leg2->AddEntry(g_w_tot,  "TOT",       "p");
        leg2->Draw();

        c_res_cmp->Update();
        fout->cd();
        c_res_cmp->Write("Compare_Resolution_Slew_vs_TOT");
    }

    // ==================================================================
    //  CANVAS 8: CONFRONTO C vs x — SLEW vs TOT
    // ==================================================================
    {
        TCanvas *c_C_cmp = new TCanvas("c_C_compare",
                                        "C vs x: confronto Slew vs TOT",
                                        800, 600);
        c_C_cmp->SetGrid(1, 1);

        TGraphErrors *g_C_slew = new TGraphErrors(npts,
                                                  x_pos.data(), C_mean.data(),
                                                  dx_pos.data(), C_err.data());
        g_C_slew->SetMarkerStyle(22);
        g_C_slew->SetMarkerSize(1.2);
        g_C_slew->SetMarkerColor(kBlue + 1);
        g_C_slew->SetLineColor(kBlue + 1);
        g_C_slew->SetTitle("Costante C: confronto Slew vs TOT;"
                           "Posizione x [cm];"
                           "C = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]");

        TF1 *f_C_slew = new TF1("f_C_slew", "[0]", -150, 150);
        g_C_slew->Fit(f_C_slew, "S R Q N");
        g_C_slew->Draw("AP");
        f_C_slew->SetLineColor(kBlue + 1);
        f_C_slew->SetLineWidth(2);
        f_C_slew->Draw("same");

        TGraphErrors *g_C_tot_cmp = new TGraphErrors(npts,
                                                     x_pos.data(), C_mean_tot.data(),
                                                     dx_pos.data(), C_err_tot.data());
        g_C_tot_cmp->SetMarkerStyle(23);
        g_C_tot_cmp->SetMarkerSize(1.2);
        g_C_tot_cmp->SetMarkerColor(kRed + 1);
        g_C_tot_cmp->SetLineColor(kRed + 1);

// Fit lineare C(x) = p0 + p1·x sui dati TOT
        TF1 *f_C_tot = new TF1("f_C_tot", "[0]+[1]*x", -150, 150);
        f_C_tot->SetParNames("C_{0}", "C_{1}");
        f_C_tot->SetParameters(-12.3, 0.001);
        g_C_tot_cmp->Fit(f_C_tot, "S E R Q N");
        g_C_tot_cmp->Draw("P same");
        f_C_tot->SetLineColor(kRed + 1);
        f_C_tot->SetLineWidth(2);
        f_C_tot->SetLineStyle(2);
        f_C_tot->Draw("same");

        // Salva i risultati del fit lineare nelle variabili a livello di funzione
        C_lin_p0     = f_C_tot->GetParameter(0);
        C_lin_p0_err = f_C_tot->GetParError(0);
        C_lin_p1     = f_C_tot->GetParameter(1);
        C_lin_p1_err = f_C_tot->GetParError(1);
        C_lin_chi2ndf = (f_C_tot->GetNDF() > 0) ?
                         f_C_tot->GetChisquare() / f_C_tot->GetNDF() : -1.0;

        TPaveText *pt_CC = new TPaveText(0.15, 0.68, 0.62, 0.92, "NDC");
        pt_CC->SetFillColorAlpha(kWhite, 0.85);
        pt_CC->SetTextFont(42);
        pt_CC->SetTextSize(0.035);
        pt_CC->SetTextAlign(12);
        pt_CC->AddText(Form("#color[4]{C_{slew} = %.3f #pm %.3f ns}",
                            f_C_slew->GetParameter(0), f_C_slew->GetParError(0)));
        pt_CC->AddText(Form("#color[2]{C(x) = (%.4f #pm %.4f) + (%.6f #pm %.6f) #upoint x}",
                            C_lin_p0, C_lin_p0_err, C_lin_p1, C_lin_p1_err));
        pt_CC->AddText(Form("#color[2]{#chi^{2}/ndf = %.2f}", C_lin_chi2ndf));
        pt_CC->Draw();

        TLegend *leg3 = new TLegend(0.60, 0.75, 0.90, 0.92);
        leg3->SetTextFont(42);
        leg3->SetTextSize(0.04);
        leg3->AddEntry(g_C_slew,    "Slew rate", "p");
        leg3->AddEntry(g_C_tot_cmp, "TOT",       "p");
        leg3->Draw();

        c_C_cmp->Update();
        fout->cd();
c_C_cmp->Write("Compare_C_Slew_vs_TOT");
    }

    // ==================================================================
    //  CANVAS 8b: FIT LINEARE C(x) TOT CON RESIDUI
    // ==================================================================
    // Pad superiore: dati TOT + retta C(x) = C_0 + C_1·x
    // Pad inferiore: residui normalizzati (pull)
    {
        TCanvas *c_C_lin = new TCanvas("c_C_lin",
                                        "Fit lineare C(x) TOT",
                                        900, 700);
        c_C_lin->Divide(1, 2);

        // --- Pad superiore: dati + fit ---
        c_C_lin->cd(1);
        gPad->SetPad(0.0, 0.35, 1.0, 1.0);
        gPad->SetBottomMargin(0.02);
        gPad->SetGrid(1, 1);

        TGraphErrors *g_C_lin = new TGraphErrors(npts,
                                                  x_pos.data(), C_mean_tot.data(),
                                                  dx_pos.data(), C_err_tot.data());
        g_C_lin->SetMarkerStyle(23);
        g_C_lin->SetMarkerSize(1.3);
        g_C_lin->SetMarkerColor(kRed + 1);
        g_C_lin->SetLineColor(kRed + 1);
        g_C_lin->SetTitle("");
        g_C_lin->GetYaxis()->SetTitle("C = t_{3} - #frac{t_{1}+t_{2}}{2} [ns]");
        g_C_lin->GetYaxis()->SetTitleSize(0.06);
        g_C_lin->GetYaxis()->SetLabelSize(0.05);
        g_C_lin->GetXaxis()->SetLabelSize(0);   // nascosto, lo mostra il pad sotto
        g_C_lin->Draw("AP");

        TF1 *f_C_lin_draw = new TF1("f_C_lin_draw", "[0]+[1]*x", -150, 150);
        f_C_lin_draw->SetParameters(C_lin_p0, C_lin_p1);
        f_C_lin_draw->SetLineColor(kRed + 1);
        f_C_lin_draw->SetLineWidth(2);
        f_C_lin_draw->Draw("same");

        // Banda ±1σ del fit (propagazione errori: σ_C² = σ_p0² + x²·σ_p1²)
        const int nband = 200;
        double xb[nband], yb_up[nband], yb_dn[nband];
        for (int i = 0; i < nband; i++) {
            xb[i] = -150.0 + 300.0 * i / (nband - 1);
            double c_val = C_lin_p0 + C_lin_p1 * xb[i];
            double c_err = sqrt(C_lin_p0_err * C_lin_p0_err
                              + xb[i] * xb[i] * C_lin_p1_err * C_lin_p1_err);
            yb_up[i] = c_val + c_err;
            yb_dn[i] = c_val - c_err;
        }
        TGraph *g_band_up = new TGraph(nband, xb, yb_up);
        TGraph *g_band_dn = new TGraph(nband, xb, yb_dn);
        g_band_up->SetLineColor(kRed - 7);
        g_band_up->SetLineStyle(7);
        g_band_dn->SetLineColor(kRed - 7);
        g_band_dn->SetLineStyle(7);
        g_band_up->Draw("L same");
        g_band_dn->Draw("L same");
        g_C_lin->Draw("P same");   // ridisegna punti sopra la banda

        TPaveText *pt_Clin = new TPaveText(0.15, 0.15, 0.65, 0.45, "NDC");
        pt_Clin->SetFillColorAlpha(kWhite, 0.85);
        pt_Clin->SetTextFont(42);
        pt_Clin->SetTextSize(0.045);
        pt_Clin->SetTextAlign(12);
        pt_Clin->AddText("C(x) = C_{0} + C_{1} #upoint x  (dati TOT)");
        pt_Clin->AddText(Form("C_{0} = %.4f #pm %.4f ns",
                               C_lin_p0, C_lin_p0_err));
        pt_Clin->AddText(Form("C_{1} = %.6f #pm %.6f ns/cm",
                               C_lin_p1, C_lin_p1_err));
        pt_Clin->AddText(Form("#chi^{2}/ndf = %.2f", C_lin_chi2ndf));
        pt_Clin->Draw();

        // --- Pad inferiore: residui (pull) ---
        c_C_lin->cd(2);
        gPad->SetPad(0.0, 0.0, 1.0, 0.35);
        gPad->SetTopMargin(0.02);
        gPad->SetBottomMargin(0.25);
        gPad->SetGrid(1, 1);

        std::vector<double> pull_x(npts), pull_y(npts), pull_ex(npts), pull_ey(npts);
        for (int i = 0; i < npts; i++) {
            double C_fit_i = C_lin_p0 + C_lin_p1 * x_pos[i];
            double res_i   = C_mean_tot[i] - C_fit_i;
            pull_x[i]  = x_pos[i];
            pull_ex[i] = dx_pos[i];
            pull_y[i]  = (C_err_tot[i] > 0) ? res_i / C_err_tot[i] : 0.0;
            pull_ey[i] = (C_err_tot[i] > 0) ? 1.0 : 0.0;
        }
        TGraphErrors *g_pull = new TGraphErrors(npts,
                                                 pull_x.data(), pull_y.data(),
                                                 pull_ex.data(), pull_ey.data());
        g_pull->SetMarkerStyle(20);
        g_pull->SetMarkerSize(1.0);
        g_pull->SetMarkerColor(kRed + 1);
        g_pull->SetLineColor(kRed + 1);
        g_pull->SetTitle(";Posizione x [cm];Pull (C_{mis} - C_{fit}) / #sigma");
        g_pull->GetXaxis()->SetTitleSize(0.10);
        g_pull->GetXaxis()->SetLabelSize(0.08);
        g_pull->GetYaxis()->SetTitleSize(0.09);
        g_pull->GetYaxis()->SetTitleOffset(0.5);
        g_pull->GetYaxis()->SetLabelSize(0.08);
        g_pull->GetYaxis()->SetRangeUser(-3.5, 3.5);
        g_pull->GetXaxis()->SetLimits(-150, 150);
        g_pull->Draw("AP");

        TLine *l_zero = new TLine(-150, 0, 150, 0);
        l_zero->SetLineColor(kBlack);
        l_zero->SetLineWidth(2);
        l_zero->SetLineStyle(2);
        l_zero->Draw("same");
        TLine *l_p1 = new TLine(-150, 1, 150, 1);
        l_p1->SetLineColor(kGray + 1); l_p1->SetLineStyle(3); l_p1->Draw("same");
        TLine *l_m1 = new TLine(-150, -1, 150, -1);
        l_m1->SetLineColor(kGray + 1); l_m1->SetLineStyle(3); l_m1->Draw("same");
        TLine *l_p2 = new TLine(-150, 2, 150, 2);
        l_p2->SetLineColor(kGray + 2); l_p2->SetLineStyle(3); l_p2->Draw("same");
        TLine *l_m2 = new TLine(-150, -2, 150, -2);
        l_m2->SetLineColor(kGray + 2); l_m2->SetLineStyle(3); l_m2->Draw("same");

        c_C_lin->Update();
        fout->cd();
        c_C_lin->Write("C_linear_fit_TOT");
    }

    // ==================================================================
    //  CANVAS 9 e 10: GALLERIE ISTOGRAMMI TOT
    // ==================================================================
    {
        int ncols = (npts <= 4) ? npts : (npts <= 9 ? 3 : 4);
        int nrows = (npts + ncols - 1) / ncols;

        TCanvas *c_gal_dt12_tot = new TCanvas("c_dt12_tot_gallery",
                                               "Distribuzioni #Delta t_{12} (TOT)",
                                               350 * ncols, 300 * nrows);
        c_gal_dt12_tot->Divide(ncols, nrows);
        for (int i = 0; i < npts; i++) {
            c_gal_dt12_tot->cd(i + 1);
            gPad->SetGrid(1, 1);
            TH1D *h = (TH1D*)fout->Get(Form("h_dt12_tot_%s", labels[i].c_str()));
            if (h) {
                h->SetLineColor(kRed + 1);
                h->SetFillColorAlpha(kRed - 9, 0.3);
                h->Draw();
            }
        }
        c_gal_dt12_tot->Update();
        fout->cd();
        c_gal_dt12_tot->Write("Gallery_Dt12_TOT");

        TCanvas *c_gal_C_tot = new TCanvas("c_C_tot_gallery",
                                            "Distribuzioni C (TOT)",
                                            350 * ncols, 300 * nrows);
        c_gal_C_tot->Divide(ncols, nrows);
        for (int i = 0; i < npts; i++) {
            c_gal_C_tot->cd(i + 1);
            gPad->SetGrid(1, 1);
            TH1D *h = (TH1D*)fout->Get(Form("h_C_tot_%s", labels[i].c_str()));
            if (h) {
                h->SetLineColor(kOrange + 1);
                h->SetFillColorAlpha(kOrange - 9, 0.3);
                h->Draw();
            }
        }
        c_gal_C_tot->Update();
        fout->cd();
        c_gal_C_tot->Write("Gallery_C_TOT");
    }

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
    std::cout << "  --- Metodo SLEW RATE ---" << std::endl;
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
    std::cout << "  --- Fit lineare C(x) TOT ---" << std::endl;
    std::cout << "  C(x) = " << Form("%.4f ± %.4f", C_lin_p0, C_lin_p0_err)
              << " + (" << Form("%.6f ± %.6f", C_lin_p1, C_lin_p1_err)
              << ") · x  ns" << std::endl;
    std::cout << "  χ²/ndf (C lineare):       = "
              << Form("%.2f", C_lin_chi2ndf) << std::endl;          
std::cout << "\n  Output salvato in: " << outname << std::endl;
    std::cout << "  TTree 'summary' con tutti i parametri dei fit." << std::endl;
    std::cout << "=============================================" << std::endl;

    // ==================================================================
    //  TTree "fit_params" — parametri caricabili da TOF_Analysis
    // ==================================================================
    // Salviamo in un piccolo TTree con un singolo entry tutti i parametri
    // numerici che TOF_Analysis dovrebbe usare automaticamente:
    //   - m, q della retta Δt₁₂ vs x (per ricostruire la posizione)
    //   - k_PMT1, k_PMT2 della relazione A = k·s (recupero slew rate)
    //   - p0_PMT1/2, p1_PMT1/2 della relazione A = p0·exp(p1·TOT)
    //     (recupero TOT)
    //   - A_max_PMT1/2 (limite sup di estrapolazione TOT)
    //
    // In questo modo TOF_Analysis::LoadCalibrationFromFile() può leggere
    // tutto il necessario senza che l'utente metta i numeri a mano.
    {
        fout->cd();
        TTree *fit_params = new TTree("fit_params",
                                       "Parametri retta + recupero clippati (slew + TOT)");

        // Parametri retta Δt₁₂ vs x
        Float_t fp_m,        fp_m_err;
        Float_t fp_q,        fp_q_err;
        Float_t fp_v_eff,    fp_v_eff_err;
        Float_t fp_chi2ndf;
        // Parametri recupero slew rate
        Float_t fp_k_PMT1,   fp_k_PMT1_err;
        Float_t fp_k_PMT2,   fp_k_PMT2_err;
        // Parametri recupero TOT (modello esponenziale)
        Float_t fp_p0_TOT_PMT1, fp_p0_TOT_PMT1_err;
        Float_t fp_p1_TOT_PMT1, fp_p1_TOT_PMT1_err;
        Float_t fp_p0_TOT_PMT2, fp_p0_TOT_PMT2_err;
        Float_t fp_p1_TOT_PMT2, fp_p1_TOT_PMT2_err;
        Float_t fp_Amax_TOT_PMT1, fp_Amax_TOT_PMT2;
        // Soglia di riferimento TOT [mV] e indice
        Float_t fp_TOT_thr_ref;
        Int_t   fp_TOT_calib_idx;

        // ---- Branch: parametri retta ----
        fit_params->Branch("m",          &fp_m,         "m/F");
        fit_params->Branch("m_err",      &fp_m_err,     "m_err/F");
        fit_params->Branch("q",          &fp_q,         "q/F");
        fit_params->Branch("q_err",      &fp_q_err,     "q_err/F");
        fit_params->Branch("v_eff",      &fp_v_eff,     "v_eff/F");
        fit_params->Branch("v_eff_err",  &fp_v_eff_err, "v_eff_err/F");
        fit_params->Branch("chi2ndf",    &fp_chi2ndf,   "chi2ndf/F");

// ---- Branch: parametri slew rate (modello A = k·s + q) ----
        fit_params->Branch("k_PMT1",     &fp_k_PMT1,     "k_PMT1/F");
        fit_params->Branch("k_PMT1_err", &fp_k_PMT1_err, "k_PMT1_err/F");
        fit_params->Branch("k_PMT2",     &fp_k_PMT2,     "k_PMT2/F");
        fit_params->Branch("k_PMT2_err", &fp_k_PMT2_err, "k_PMT2_err/F");
        Float_t fp_q_PMT1, fp_q_PMT1_err, fp_q_PMT2, fp_q_PMT2_err;
        fit_params->Branch("q_PMT1",     &fp_q_PMT1,     "q_PMT1/F");
        fit_params->Branch("q_PMT1_err", &fp_q_PMT1_err, "q_PMT1_err/F");
        fit_params->Branch("q_PMT2",     &fp_q_PMT2,     "q_PMT2/F");
        fit_params->Branch("q_PMT2_err", &fp_q_PMT2_err, "q_PMT2_err/F");

        // ---- Branch: parametri TOT (modello A = p0·exp(p1·TOT)) ----
        fit_params->Branch("p0_TOT_PMT1",     &fp_p0_TOT_PMT1,     "p0_TOT_PMT1/F");
        fit_params->Branch("p0_TOT_PMT1_err", &fp_p0_TOT_PMT1_err, "p0_TOT_PMT1_err/F");
        fit_params->Branch("p1_TOT_PMT1",     &fp_p1_TOT_PMT1,     "p1_TOT_PMT1/F");
        fit_params->Branch("p1_TOT_PMT1_err", &fp_p1_TOT_PMT1_err, "p1_TOT_PMT1_err/F");
        fit_params->Branch("p0_TOT_PMT2",     &fp_p0_TOT_PMT2,     "p0_TOT_PMT2/F");
        fit_params->Branch("p0_TOT_PMT2_err", &fp_p0_TOT_PMT2_err, "p0_TOT_PMT2_err/F");
        fit_params->Branch("p1_TOT_PMT2",     &fp_p1_TOT_PMT2,     "p1_TOT_PMT2/F");
        fit_params->Branch("p1_TOT_PMT2_err", &fp_p1_TOT_PMT2_err, "p1_TOT_PMT2_err/F");
        fit_params->Branch("Amax_TOT_PMT1",   &fp_Amax_TOT_PMT1,   "Amax_TOT_PMT1/F");
        fit_params->Branch("Amax_TOT_PMT2",   &fp_Amax_TOT_PMT2,   "Amax_TOT_PMT2/F");

        // ---- Branch: soglia TOT di riferimento ----
        fit_params->Branch("TOT_thr_ref",     &fp_TOT_thr_ref,     "TOT_thr_ref/F");
        fit_params->Branch("TOT_calib_idx",   &fp_TOT_calib_idx,   "TOT_calib_idx/I");

        // ---- Branch: fit lineare C(x) = C_0 + C_1·x (dati TOT) ----
        Float_t fp_C_lin_p0, fp_C_lin_p0_err;
        Float_t fp_C_lin_p1, fp_C_lin_p1_err;
        Float_t fp_C_lin_chi2ndf;
        fit_params->Branch("C_lin_p0",        &fp_C_lin_p0,        "C_lin_p0/F");
        fit_params->Branch("C_lin_p0_err",    &fp_C_lin_p0_err,    "C_lin_p0_err/F");
        fit_params->Branch("C_lin_p1",        &fp_C_lin_p1,        "C_lin_p1/F");
        fit_params->Branch("C_lin_p1_err",    &fp_C_lin_p1_err,    "C_lin_p1_err/F");
        fit_params->Branch("C_lin_chi2ndf",   &fp_C_lin_chi2ndf,   "C_lin_chi2ndf/F");

        // ---- Popolamento dei valori ----
        fp_m         = (Float_t)m_fit;
        fp_m_err     = (Float_t)m_err;
        fp_q         = (Float_t)q_fit;
        fp_q_err     = (Float_t)q_err;
        fp_v_eff     = (Float_t)v_eff;
        fp_v_eff_err = (Float_t)v_eff_err;
        fp_chi2ndf   = (Float_t)chi2ndf;

        fp_k_PMT1    = (Float_t)gK_PMT[0];
        fp_k_PMT1_err= (Float_t)gK_PMT_err[0];
        fp_k_PMT2    = (Float_t)gK_PMT[1];
        fp_k_PMT2_err= (Float_t)gK_PMT_err[1];
        fp_q_PMT1    = (Float_t)gQ_PMT[0];
        fp_q_PMT1_err= (Float_t)gQ_PMT_err[0];
        fp_q_PMT2    = (Float_t)gQ_PMT[1];
        fp_q_PMT2_err= (Float_t)gQ_PMT_err[1];

        // I parametri TOT sono nelle globali popolate da CalibrateTOT.
        // Gli errori non sono stati salvati come globali, quindi metto 0
        // (non sono critici: il recupero usa solo i valori centrali).
        fp_p0_TOT_PMT1     = (Float_t)gTOT_p0_PMT[0];
        fp_p0_TOT_PMT1_err = 0.0f;
        fp_p1_TOT_PMT1     = (Float_t)gTOT_p1_PMT[0];
        fp_p1_TOT_PMT1_err = 0.0f;
        fp_p0_TOT_PMT2     = (Float_t)gTOT_p0_PMT[1];
        fp_p0_TOT_PMT2_err = 0.0f;
        fp_p1_TOT_PMT2     = (Float_t)gTOT_p1_PMT[1];
        fp_p1_TOT_PMT2_err = 0.0f;
        fp_Amax_TOT_PMT1   = (Float_t)gTOT_Amax_PMT[0];
        fp_Amax_TOT_PMT2   = (Float_t)gTOT_Amax_PMT[1];

        fp_TOT_thr_ref     = (Float_t)TOT_THR_MV[TOT_CALIB_THR_IDX];
        fp_TOT_calib_idx   = TOT_CALIB_THR_IDX;
        fp_C_lin_p0        = (Float_t)C_lin_p0;
        fp_C_lin_p0_err    = (Float_t)C_lin_p0_err;
        fp_C_lin_p1        = (Float_t)C_lin_p1;
        fp_C_lin_p1_err    = (Float_t)C_lin_p1_err;
        fp_C_lin_chi2ndf   = (Float_t)C_lin_chi2ndf;
        fit_params->Fill();
        fit_params->Write();

        std::cout << "\n[INFO] TTree 'fit_params' salvato:" << std::endl;
        std::cout << "       m, q, v_eff (retta calibrazione)" << std::endl;
        std::cout << "       k_PMT1 = " << gK_PMT[0]
                  << ", k_PMT2 = " << gK_PMT[1] << " (recupero slew)" << std::endl;
        if (gTOT_calibrated) {
            std::cout << "       p0_TOT/p1_TOT PMT1 = " << gTOT_p0_PMT[0]
                      << " / " << gTOT_p1_PMT[0]
                      << " (τ_eff = " << 1.0/gTOT_p1_PMT[0] << " ns)" << std::endl;
            std::cout << "       p0_TOT/p1_TOT PMT2 = " << gTOT_p0_PMT[1]
                      << " / " << gTOT_p1_PMT[1]
                      << " (τ_eff = " << 1.0/gTOT_p1_PMT[1] << " ns)" << std::endl;
            std::cout << "       Soglia TOT riferimento: "
                      << TOT_THR_MV[TOT_CALIB_THR_IDX] << " mV" << std::endl;
        } else {
            std::cout << "       [WARNING] Calibrazione TOT non disponibile!"
                      << std::endl;
        }
    }

    // Chiudi il file di output
    fout->Close();

    std::cout << "\n[DONE] Calibrazione completata." << std::endl;
}
