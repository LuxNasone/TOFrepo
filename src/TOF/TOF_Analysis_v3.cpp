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
const double CFD_FRACTION = 0.35;
const double NOISE_THRESH = 5.0;

// --- Soglie di qualità per scartare eventi non fisici ---
// CLIP_V_LO: se V_min scende sotto questa soglia [mV], il segnale è saturato.
//   La DRS4 satura a circa ±500 mV. Con margine: −450 mV.
const double CLIP_V_LO = -450.0;

// CLIP_V_HI: se V_max supera questa soglia [mV], il segnale ha una componente
//   positiva anomala (oscillazione, ringing). Un impulso fisico negativo non
//   dovrebbe mai superare significativamente la baseline (~0 mV).
const double CLIP_V_HI = 100.0;

// OSC_FRACTION: se la frazione di campioni con |V − baseline| > OSC_AMP_THRESH
//   × ampiezza supera OSC_FRACTION, il segnale è oscillante (non un impulso fisico).
//   Un impulso di scintillazione è breve: tipicamente <10% dei campioni sono
//   "grandi". Un'oscillazione come l'evento 13/4341 ne ha >30%.
const double OSC_FRACTION   = 0.15;  // Soglia sulla frazione (15%)
const double OSC_AMP_THRESH = 0.20;  // Soglia sull'ampiezza (20% del picco)

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
const double X_CUT_LO = -150.0;     // [cm]  bordo sinistro fisico
const double X_CUT_HI =  150.0;     // [cm]  bordo destro fisico

// --- Geometria del setup TOF (Guida B) ---
// h: distanza VERTICALE tra il centro della barra e il centro di PMT3 [cm]
//    Manca & Selicato usano 156 cm; il vostro setup usa ~101 cm.
//    MODIFICARE in base alla misura effettiva.
double PAR_H = 101.0;

// x_PMT3: posizione orizzontale del centro di PMT3 rispetto al centro
//         della barra [cm]. 0 = allineato al centro.
//          ~139.9 cm dalla media degli x_imp, ovvero
//         circa al centro della barra nel loro riferimento (0-280).
//         Nel nostro riferimento (centro = 0): x_PMT3 ≈ 0.
//         MODIFICARE se PMT3 è disallineato.
double PAR_X_PMT3 = 0.0;

// --- Parametri di calibrazione ---
// Retta: Δt₁₂ = m · x + q   →   x = (Δt₁₂ − q) / m
// QUESTI VALORI DEVONO ESSERE AGGIORNATI DAI RISULTATI DI TOF_Calibration.cpp
double CAL_M = 0.132;     // Pendenza retta calibrazione [ns/cm]
double CAL_Q = -1.805;      // Intercetta retta calibrazione [ns]

// --- C(x_k): costante di offset per il TOF ---
// Coppie {posizione_x [cm], valore_C [ns]} dalle calibrazioni.
// La funzione GetC(x) restituisce il valore C del punto di calibrazione
// più vicino a x (funzione a scalini / piecewise constant).
// AGGIORNARE con i risultati del fit di C per ogni posizione.
struct CPoint {
    double x;    // Posizione di calibrazione [cm]
    double C;    // Valore di C = t₃ − (t₁+t₂)/2 misurato [ns]
};

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

    // ---- CONTROLLO OSCILLAZIONE ----
    // Un impulso di scintillazione è un singolo picco negativo breve (~10 ns).
    // La maggior parte della forma d'onda (200 ns) è vicina alla baseline.
    // Un'oscillazione patologica ha molti campioni con |V − bl| grande.
    // Contiamo la frazione di campioni sopra il 20% dell'ampiezza.
    double osc_thresh = OSC_AMP_THRESH * amp;
    int n_large = 0;
    for (int i = 0; i < ns; i++) {
        if (fabs(v[i] - bl) > osc_thresh) n_large++;
    }
    double frac_large = (double)n_large / ns;
    if (frac_large > OSC_FRACTION) {
        cd.is_oscillating = true;
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

double GetC(double x) {

    int n = (int)CAL_C_POINTS.size();
    if (n == 0) {
        std::cerr << "[ERRORE] Nessun punto di calibrazione C definito!" << std::endl;
        return 0.0;
    }
    if (n == 1) return CAL_C_POINTS[0].C;

    // Trova il punto di calibrazione più vicino a x
    int i_nearest = 0;
    double min_dist = fabs(x - CAL_C_POINTS[0].x);
    for (int i = 1; i < n; i++) {
        double dist = fabs(x - CAL_C_POINTS[i].x);
        if (dist < min_dist) {
            min_dist = dist;
            i_nearest = i;
        }
    }
    return CAL_C_POINTS[i_nearest].C;
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

    // --- Leggi m e q dal fit della retta (salvato nel canvas) ---
    // I parametri della retta sono nel TF1 "f_lin" associato al TGraphErrors.
    // In alternativa, l'utente li imposta manualmente.
    // Per sicurezza, non sovrascriviamo automaticamente m e q,
    // ma stampiamo un reminder.
    std::cout << "\n[REMINDER] Aggiorna manualmente CAL_M e CAL_Q in cima al file\n"
              << "           oppure usa SetCalibration(m, q) dopo il caricamento."
              << std::endl;

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

    // --- Istogrammi ---
    TH1D *h_tof    = new TH1D("h_tof",    "Tempo di volo;TOF [ns];Conteggi",
                               NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta   = new TH1D("h_beta",   "Distribuzione #beta = v/c;#beta;Conteggi",
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
    int total_good      = 0;
    int rej_no_cfd      = 0;   // Contatore: scartati per CFD mancante
    int rej_clipped     = 0;   // Contatore: scartati per segnale clippato/saturato
    int rej_oscillating = 0;   // Contatore: scartati per oscillazione patologica
    int rej_x_out       = 0;   // Contatore: scartati per x fuori dalla barra
    int rej_tof_unphys  = 0;   // Contatore: scartati per TOF < h/c (non fisico)

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
            //  TAGLIO 1: segnale non clippato su nessun canale
            // ============================================================
            // Un segnale che satura la DRS4 (come l'evento 13/4341 con
            // V oscillante da −500 a +500 mV) non ha un leading edge
            // definito → il tempo CFD è privo di senso fisico.
            {
                bool any_clipped = false;
                for (int k = 0; k < 3; k++) {
                    if (e.ch[k].is_clipped) { any_clipped = true; break; }
                }
                if (any_clipped) {
                    good = 0; dt12 = -999; x_imp = -999;
                    T_meas = -999; C_used = -999; tof = -999;
                    path_len = -999; theta = -999; vel = -999;
                    beta = -999; inv_vel = -999;
                    tree->Fill();
                    rej_clipped++;
                    continue;
                }
            }

            // ============================================================
            //  TAGLIO 2: segnale non oscillante su nessun canale
            // ============================================================
            // Un impulso fisico di scintillazione è un singolo picco
            // negativo breve. Se >15% dei campioni hanno |V−bl| > 20%
            // dell'ampiezza, il segnale è un'oscillazione patologica.
            {
                bool any_osc = false;
                for (int k = 0; k < 3; k++) {
                    if (e.ch[k].is_oscillating) { any_osc = true; break; }
                }
                if (any_osc) {
                    good = 0; dt12 = -999; x_imp = -999;
                    T_meas = -999; C_used = -999; tof = -999;
                    path_len = -999; theta = -999; vel = -999;
                    beta = -999; inv_vel = -999;
                    tree->Fill();
                    rej_oscillating++;
                    continue;
                }
            }

            // ---- PASSO 1: Δt₁₂ e posizione di impatto ----
            dt12 = t_cfd[0] - t_cfd[1];
            x_imp = (dt12 - CAL_Q) / CAL_M;

// ============================================================
            //  TAGLIO 3: x_imp dentro i limiti FISICI della barra
            // ============================================================
            // La barra ha estensione [X_CUT_LO, X_CUT_HI] = [-140, +140] cm.
            // Un evento ricostruito fuori da questo intervallo NON può
            // essere un muone reale che ha attraversato la barra: è
            // certamente un artefatto del CFD su rumore/ringing.
            if (x_imp < X_CUT_LO || x_imp > X_CUT_HI) {
                good = 0;
                T_meas = -999; C_used = -999; tof = -999;
                path_len = -999; theta = -999; vel = -999;
                beta = -999; inv_vel = -999;
                tree->Fill();
                rej_x_out++;
                continue;
            }

// ---- PASSO 2: tempo di volo ----
            T_meas = t_cfd[2] - (t_cfd[0] + t_cfd[1]) / 2.0;
            C_used = GetC(x_imp);
            tof = T_meas - C_used;

            // ============================================================
            //  TAGLIO 4: TOF non inferiore al minimo fisicamente possibile
            // ============================================================
            // Per qualunque traiettoria, la distanza percorsa è
            //     ℓ = √[(x_imp - x_PMT3)² + h²] ≥ h
            // e per una particella materiale (β ≤ 1):
            //     TOF = ℓ / (β·c) ≥ ℓ/c ≥ h/c
            // Quindi h/c è un LIMITE INFERIORE assoluto, valido per
            // qualunque x_imp: eventi con TOF < h/c violano la causalità
            // per qualunque traiettoria possibile e sono certamente
            // artefatti (re-trigger, doppio impulso, accidentali...).
            //
            // NB: questo è un cut LASCO. Eventi con h/c < TOF < ℓ/c (cioè
            //     1 < β < ℓ/h) vengono CONSERVATI nel TTree (good=0 ma
            //     tof,beta valorizzati): sono fisicamente non leciti
            //     ma diagnosticamente interessanti, in particolare per
            //     studiare l'eventuale "cluster superluminale" che P&R
            //     e M&S non hanno indagato a fondo.
            if (tof < tof_min_phys) {
                good = 0;
                // Conservo dt12, x_imp, T_meas, C_used, tof (sono valori
                // CALCOLATI, anche se TOF non è fisico). Azzero solo i
                // campi che non hanno senso senza un TOF lecito.
                path_len = -999; theta = -999; vel = -999;
                beta    = -999; inv_vel = -999;
                tree->Fill();
                rej_tof_unphys++;
                continue;
            }

            // ---- PASSO 3: distanza percorsa ----
            double dx = x_imp - PAR_X_PMT3;
            path_len = sqrt(dx * dx + PAR_H * PAR_H);

            // ---- PASSO 4: angolo ----
            theta = atan2(dx, PAR_H);

            // ---- PASSO 5: velocità e β ----
            if (tof > 0.1) {
                vel     = path_len / tof;
                beta    = vel / C_LIGHT;
                inv_vel = tof / path_len;
            } else {
                vel     = -999;
                beta    = -999;
                inv_vel = -999;
            }

            good = 1;
            total_good++;

            // Riempi istogrammi (solo eventi che passano TUTTI i tagli)
            h_T_meas->Fill(T_meas);
            h_x_imp->Fill(x_imp);

            if (tof > 0.1) {
                h_tof->Fill(tof);
                h_path->Fill(path_len);
                h_theta->Fill(theta);
                h2_tof_x->Fill(x_imp, tof);

                if (beta > 0 && beta < 3.0) {
                    h_beta->Fill(beta);
                    h2_beta_x->Fill(x_imp, beta);
                }
                if (inv_vel > 0 && inv_vel < 0.2) {
                    h_invv->Fill(inv_vel);
                }
            }

            tree->Fill();
        }
    }

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  RIEPILOGO TAGLI DI QUALITA'                " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Eventi totali:        " << total_events << std::endl;
    std::cout << "  Scartati (no CFD):    " << rej_no_cfd
              << " (" << (total_events>0 ? 100.0*rej_no_cfd/total_events : 0)
              << "%)" << std::endl;
    std::cout << "  Scartati (clipped):   " << rej_clipped
              << " (" << (total_events>0 ? 100.0*rej_clipped/total_events : 0)
              << "%)" << std::endl;
    std::cout << "  Scartati (oscill.):   " << rej_oscillating
              << " (" << (total_events>0 ? 100.0*rej_oscillating/total_events : 0)
              << "%)" << std::endl;
    std::cout << "  Scartati (x fuori):   " << rej_x_out
              << " (" << (total_events>0 ? 100.0*rej_x_out/total_events : 0)
              << "%)" << std::endl;
    std::cout << "  Scartati (TOF<h/c):   " << rej_tof_unphys
              << " (" << (total_events>0 ? 100.0*rej_tof_unphys/total_events : 0)
              << "%)" << std::endl;
    std::cout << "  ─────────────────────────────────────" << std::endl;
    std::cout << "  Eventi buoni:         " << total_good
              << " (" << (total_events>0 ? 100.0*total_good/total_events : 0)
              << "%)" << std::endl;
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

    fout->Close();

    // --- Riepilogo finale ---
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  ANALISI TOF COMPLETATA                      " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  File processati: " << file_list.size() << std::endl;
    std::cout << "  Eventi totali:   " << total_events << std::endl;
    std::cout << "  Eventi buoni:    " << total_good << std::endl;
    std::cout << "  Parametri usati:" << std::endl;
    std::cout << "    h       = " << PAR_H << " cm" << std::endl;
    std::cout << "    x_PMT3  = " << PAR_X_PMT3 << " cm" << std::endl;
    std::cout << "    m (cal) = " << CAL_M << " ns/cm" << std::endl;
    std::cout << "    q (cal) = " << CAL_Q << " ns" << std::endl;
    std::cout << "  Output: " << outname << std::endl;
    std::cout << "=============================================" << std::endl;
}
