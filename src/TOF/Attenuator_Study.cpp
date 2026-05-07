// ==========================================================================
//  Attenuator_Study.cpp — Studio dell'effetto dell'attenuatore -6 dB
//                        sui segnali PMT1 e PMT2 in ingresso alla DRS4
// ==========================================================================
//
//  SCOPO:
//    Confrontare le caratteristiche dei segnali PMT1 (CH1) e PMT2 (CH2)
//    in due dataset acquisiti nella stessa identica configurazione (PMT3 in
//    x_3 = 0), uno CON attenuatore -6 dB inserito tra il PMT e l'ingresso
//    DRS4, l'altro SENZA attenuatore. Tutto il resto (HV, soglie, posizioni,
//    impostazioni DRS, logica di trigger) è identico tra i due dataset.
//
//    L'attenuatore -6 dB riduce idealmente l'ampiezza in tensione di un
//    fattore 10^(-6/20) ≈ 0.5012 (≈ 1/2). Se l'attenuatore fosse perfetto:
//      - V_peak              → V_peak / 2
//      - SNR                 → SNR / 2 (rumore DRS invariato in V)
//      - rise time 20-80%    → INVARIATO (forma del fronte)
//      - slew rate 20-80%    → ridotto di 1/2 (proporzionale a V/t)
//      - t_CFD               → INVARIATO (la frazione 20% è scale-invariant)
//      - Δt_12, t_3-(t_1+t_2)/2 → INVARIATI (sono differenze di tempi)
//
//    Le deviazioni dall'idealità ci dicono se l'attenuatore introduce
//    bias temporali, distorsioni di forma, o cambia la risoluzione
//    temporale dei singoli canali.
//
//  GRANDEZZE CONFRONTATE (per ognuna: distribuzioni ATT vs NOATT sovrapposte)
//    1. Ampiezza CH1, CH2  (rispetto alla baseline, NON picco-picco)
//    2. SNR        CH1, CH2  (= ampiezza / RMS_baseline)
//    3. Rise time 20-80%  CH1, CH2
//    4. Slew rate 20-80%  CH1, CH2 (sul fronte di discesa, segnali negativi)
//    5. Tempo CFD 20%      CH1, CH2 (assoluti, scala DRS)
//    6. Δt_12 = t_1 - t_2
//    7. t_3 - (t_1+t_2)/2
//
//  UTILIZZO:
//    root -l 'Attenuator_Study.cpp("/path/with_att.xml", "/path/without_att.xml")'
//
//    Oppure dalla console ROOT:
//      .L Attenuator_Study.cpp
//      Attenuator_Study("file_with_att.xml", "file_without_att.xml");
//
//  OUTPUT:
//    File ROOT (default: Attenuator_Study_output.root) contenente:
//      - 2 TTree event-by-event (tree_with, tree_without)
//      - 24 istogrammi (12 grandezze × 2 dataset)
//      - 12 canvas di confronto (overlay ATT vs NOATT)
//      - 1 TTree summary con le statistiche di confronto
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
#include <iomanip>
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
const int MAX_CHANNELS = 4;      // Canali della DRS4 Eval Board

// --- Parametri di analisi della forma d'onda ---
const int    NBL_SAMPLES    = 50;     // Campioni per il calcolo della baseline
const double CFD_FRACTION   = 0.20;   // Frazione CFD primaria (richiesta: 20%)
const double RT_LOW_FRAC    = 0.20;   // Soglia bassa per rise time (20%)
const double RT_HIGH_FRAC   = 0.80;   // Soglia alta per rise time  (80%)
const double NOISE_THRESH   = 5.0;    // Soglia minima ampiezza [mV] per pulse

// --- Soglie di qualità per scartare eventi non fisici ---
// Stesse soglie del codice di calibrazione TOF_Calibration_v5.cpp:
//   - is_clipped     ← V_min < CLIP_V_LO  (saturazione DRS4)
//   - is_oscillating ← V_max > OSC_VMAX_THRESH (escursione positiva anomala)
//
// ATTENZIONE FISICA: nel dataset SENZA attenuatore l'ampiezza assoluta è
// circa il doppio rispetto a quello CON attenuatore. Quindi è normale che
// la frazione di eventi clippati sia significativamente diversa tra i due
// dataset. Il programma stampa il breakdown in modo che l'effetto sia ben
// quantificato. Se il dataset NOATT ha molti clipping, la statistica
// "buona" diminuisce ma il confronto resta valido (su statistiche grandi).

const double CLIP_V_LO       = -450.0;   // [mV]
const double OSC_VMAX_THRESH =  100.0;   // [mV]

const bool ENABLE_CFD_CUT  = true;
const bool ENABLE_CLIP_CUT = true;
const bool ENABLE_OSC_CUT  = true;

// --- Parametri di binning istogrammi ---
// Range scelti per essere adatti tipicamente a entrambi i dataset.
// Se i tuoi dati hanno range diversi, modifica qui prima di rilanciare.
//
// Le grandezze "ampiezza-dipendenti" (amp, SNR, slew_rate) hanno range
// generosi per accomodare sia il caso ATT (~ metà del NOATT) sia NOATT.

// Ampiezza [mV]
const int    NB_AMP   = 100;
const double AMP_LO   =   0.0;
const double AMP_HI   = 500.0;

// SNR (= amp / RMS_baseline)
const int    NB_SNR   = 100;
const double SNR_LO   =   0.0;
const double SNR_HI   = 250.0;

// Rise time 20-80 [ns]
const int    NB_RT    = 100;
const double RT_LO    =   0.0;
const double RT_HI    =  10.0;

// Slew rate 20-80 [mV/ns]  (in modulo)
const int    NB_SR    = 100;
const double SR_LO    =   0.0;
const double SR_HI    = 250.0;

// Tempo CFD assoluto (scala DRS) [ns]
// Range largo per accomodare la finestra DRS (1024 cell × 0.2 ns @ 5 GS/s = 200 ns)
// e tenere conto di possibili shift introdotti dall'attenuatore.
const int    NB_TCFD  = 250;
const double TCFD_LO  =   0.0;
const double TCFD_HI = 250.0;

// Δt_12 = t_1 - t_2 [ns]
const int    NB_DT12  = 200;
const double DT12_LO  = -25.0;
const double DT12_HI  =  25.0;

// t_3 - (t_1+t_2)/2 [ns]
const int    NB_C     = 200;
const double C_LO     = -50.0;
const double C_HI     =  50.0;


// ==========================================================================
//  SEZIONE 2: STRUTTURE DATI
// ==========================================================================

/// Dati di un singolo canale. Estende la versione di TOF_Calibration_v5.cpp
/// aggiungendo le grandezze necessarie a questo studio (rise time, slew
/// rate, SNR, e i due tempi di crossing al 20% e 80%).
struct ChannelData {
    int    nsamples;                  // Numero di campioni (tipicamente 1024)
    float  time[MAX_SAMPLES];         // Tempi dei campioni [ns]
    float  voltage[MAX_SAMPLES];      // Tensione dei campioni [mV]

    // Grandezze base estratte da AnalyzeChannel():
    double baseline;       // Media tensione nei primi NBL_SAMPLES [mV]
    double baseline_rms;   // Deviazione standard della baseline [mV]
    double amplitude;      // Ampiezza: baseline - V_min [mV] (positiva)
    double snr;            // amp / baseline_rms (adimensionale)
    double v_min;          // Tensione minima (picco negativo) [mV]
    double v_max;          // Tensione massima [mV]
    double t_min;          // Tempo del campione con V minima [ns]

    // Grandezze di timing:
    double t_cfd;          // Tempo CFD a CFD_FRACTION (20%) [ns]
    double t_cfd_low;      // Tempo crossing al 20% (per rise time) [ns]
    double t_cfd_high;     // Tempo crossing al 80% (per rise time) [ns]
    double rise_time;      // Rise time 20-80% [ns]
    double slew_rate;      // |dV/dt| medio sul fronte 20-80% [mV/ns]

    // Flag di qualità:
    bool   has_pulse;
    bool   cfd_ok;          // CFD primario al 20% trovato
    bool   rise_ok;         // entrambi i crossing 20% e 80% trovati
    bool   is_clipped;
    bool   is_oscillating;
};

/// Dati completi di un evento DRS4 (header + canali). Identica al codice
/// di calibrazione.
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

/// Statistiche descrittive di una distribuzione (per ogni grandezza, per
/// ognuno dei due dataset). Usato per costruire la tabella di confronto.
struct DistStats {
    double mean;        // media aritmetica
    double mean_err;    // errore standard sulla media: RMS / sqrt(N)
    double rms;         // deviazione standard (campionaria, N-1)
    double rms_err;     // errore sulla RMS: RMS / sqrt(2(N-1))
    double median;      // mediana (quantile 50%)
    double q16;         // quantile 16% (-1 sigma per gaussiana)
    double q84;         // quantile 84% (+1 sigma per gaussiana)
    int    n;           // numero di entries

    // Per le grandezze temporali, anche risultato del fit gaussiano
    double fit_mu;      // centro del fit gaussiano sul core 90%
    double fit_mu_err;
    double fit_sigma;   // larghezza del fit gaussiano
    double fit_sigma_err;
    double fit_chi2_ndf;
    bool   fit_done;
};


// ==========================================================================
//  SEZIONE 3: ANALISI DELLA FORMA D'ONDA — VERSIONE ESTESA
// ==========================================================================

/// Funzione helper: trova l'istante (interpolato linearmente) in cui la
/// tensione attraversa una soglia v_thr scendendo (perché segnali negativi),
/// cercando in avanti dall'inizio del leading edge fino al campione del
/// minimo.
///
/// PARAMETRI:
///   t          — array dei tempi
///   v          — array delle tensioni
///   i_start    — indice di partenza della ricerca
///   i_end      — indice di fine ricerca (incluso, esempio: indice del minimo)
///   v_thr      — soglia di crossing [mV]
///   t_cross_out (output) — tempo di crossing interpolato
///
/// RETURN:
///   true se ha trovato un crossing valido, false altrimenti.
///
/// Cerca il primo intervallo [i, i+1] in cui v[i] > v_thr e v[i+1] <= v_thr,
/// poi interpola linearmente:
///   t_cross = t[i] + (v_thr - v[i]) / (v[i+1] - v[i]) * (t[i+1] - t[i])
static bool FindFallingCrossing(const float* t, const float* v,
                                int i_start, int i_end,
                                double v_thr, double &t_cross_out)
{
    for (int i = i_start; i < i_end; ++i) {
        if (v[i] > v_thr && v[i+1] <= v_thr) {
            double dv = (double)v[i+1] - v[i];
            if (fabs(dv) > 1e-9) {
                t_cross_out = t[i] + (v_thr - v[i]) / dv * (t[i+1] - t[i]);
                return true;
            }
        }
    }
    return false;
}

/// AnalyzeChannel(): versione estesa rispetto a TOF_Calibration_v5.cpp.
///
/// Calcola tutte le grandezze base (baseline, ampiezza, V_min, V_max, t_min,
/// flag di clipping/oscillazione) e in più:
///   - t_cfd_low  : istante crossing al 20% (RT_LOW_FRAC)
///   - t_cfd_high : istante crossing al 80% (RT_HIGH_FRAC)
///   - rise_time  : t_cfd_high - t_cfd_low [ns]
///   - slew_rate  : (RT_HIGH_FRAC - RT_LOW_FRAC) * amp / rise_time [mV/ns]
///                  (in modulo, sempre positivo)
///   - t_cfd      : istante crossing al CFD_FRACTION (20%) — coincide con
///                  t_cfd_low qui, ma è separato per chiarezza concettuale
///                  (in caso si volesse cambiare CFD_FRACTION rispetto a
///                  RT_LOW_FRAC)
///   - snr        : amp / baseline_rms

static void AnalyzeChannel(ChannelData &cd) {

    // --- Inizializzazione a valori "non calcolato" ---
    cd.has_pulse      = false;
    cd.cfd_ok         = false;
    cd.rise_ok        = false;
    cd.is_clipped     = false;
    cd.is_oscillating = false;
    cd.baseline       = 0.0;
    cd.baseline_rms   = 0.0;
    cd.amplitude      = 0.0;
    cd.snr            = 0.0;
    cd.v_min          = 0.0;
    cd.v_max          = 0.0;
    cd.t_min          = 0.0;
    cd.t_cfd          = -999.0;
    cd.t_cfd_low      = -999.0;
    cd.t_cfd_high     = -999.0;
    cd.rise_time      = -999.0;
    cd.slew_rate      = -999.0;

    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return;

    float *t = cd.time;
    float *v = cd.voltage;

    // ---- PASSO 1: BASELINE ----
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < NBL_SAMPLES; ++i) {
        sum  += v[i];
        sum2 += (double)v[i] * v[i];
    }
    double bl     = sum / NBL_SAMPLES;
    double bl_rms = sqrt(fabs(sum2 / NBL_SAMPLES - bl * bl));
    cd.baseline     = bl;
    cd.baseline_rms = bl_rms;

    // ---- PASSO 2: MIN / MAX GLOBALI ----
    double vmin = v[0], vmax = v[0];
    int    imin = 0;
    for (int i = 1; i < ns; ++i) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
        if (v[i] > vmax)   vmax = v[i];
    }
    cd.v_min = vmin;
    cd.v_max = vmax;
    cd.t_min = t[imin];

    // ---- PASSO 3: CLIPPING ----
    if (vmin < CLIP_V_LO) cd.is_clipped = true;

    // ---- PASSO 4: AMPIEZZA E SNR ----
    double amp = bl - vmin;            // segnale negativo → amp positiva
    cd.amplitude = amp;
    if (bl_rms > 1e-9) cd.snr = amp / bl_rms;

    if (amp < NOISE_THRESH) return;
    cd.has_pulse = true;

    // ---- PASSO 5: OSCILLAZIONE ----
    if (vmax > OSC_VMAX_THRESH) cd.is_oscillating = true;

    // ---- PASSO 6: CROSSING AL 20% (CFD primario) E AL 80% (per rise time)
    //
    // Soglie:
    //   v_thr_low  = bl - 0.20 * amp   (vicino alla baseline)
    //   v_thr_high = bl - 0.80 * amp   (vicino al picco)
    //
    // Sul leading edge (V che scende), v_thr_low viene ATTRAVERSATO PRIMA di
    // v_thr_high. Quindi t_low < t_high e rise_time = t_high - t_low > 0.
    //
    // Range di ricerca: da NBL_SAMPLES (fine baseline) a imin (campione del
    // minimo). Questo evita falsi crossing su strutture parassite del tail.

    double v_thr_low  = bl - RT_LOW_FRAC  * amp;   // 20%
    double v_thr_high = bl - RT_HIGH_FRAC * amp;   // 80%

    double t_low = -999.0, t_high = -999.0;
    bool ok_low  = FindFallingCrossing(t, v, NBL_SAMPLES, imin, v_thr_low,  t_low);
    bool ok_high = FindFallingCrossing(t, v, NBL_SAMPLES, imin, v_thr_high, t_high);

    cd.t_cfd_low  = ok_low  ? t_low  : -999.0;
    cd.t_cfd_high = ok_high ? t_high : -999.0;

    // ---- PASSO 7: TEMPO CFD PRIMARIO (al CFD_FRACTION) ----
    //
    // Se CFD_FRACTION coincide con RT_LOW_FRAC (caso tipico, entrambi 20%),
    // riutilizziamo t_low. Altrimenti calcoliamo un terzo crossing.
    if (fabs(CFD_FRACTION - RT_LOW_FRAC) < 1e-6) {
        cd.t_cfd  = t_low;
        cd.cfd_ok = ok_low;
    } else {
        double v_thr_cfd = bl - CFD_FRACTION * amp;
        double t_cfd_x = -999.0;
        bool ok_cfd = FindFallingCrossing(t, v, NBL_SAMPLES, imin, v_thr_cfd, t_cfd_x);
        cd.t_cfd  = ok_cfd ? t_cfd_x : -999.0;
        cd.cfd_ok = ok_cfd;
    }

    // ---- PASSO 8: RISE TIME E SLEW RATE ----
    //
    // Definizione operativa:
    //   rise_time = t(80%) - t(20%)
    //   slew_rate = |dV/dt| medio = (V_thr_low - V_thr_high) / (t_high - t_low)
    //             = (RT_HIGH_FRAC - RT_LOW_FRAC) * amp / rise_time
    //   (V_thr_low - V_thr_high = (bl - 0.20·amp) - (bl - 0.80·amp) = 0.60·amp)
    //
    // Tutto in valore assoluto.

    if (ok_low && ok_high && (t_high > t_low)) {
        cd.rise_time = t_high - t_low;
        cd.slew_rate = (RT_HIGH_FRAC - RT_LOW_FRAC) * amp / cd.rise_time;
        cd.rise_ok   = true;
    }
}


// ==========================================================================
//  SEZIONE 4: PARSING XML DRS4
// ==========================================================================
//
//  Parser identico a quello di TOF_Calibration_v5.cpp (macchina a stati
//  finiti IDLE → IN_EVENT → IN_BOARD → IN_CHANNEL). Riportato qui per
//  rendere il file self-contained.

static int ParseXML(const char* filename, std::vector<EventData> &events) {

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "[ERRORE] Impossibile aprire: " << filename << std::endl;
        return -1;
    }
    std::cout << "[INFO] Parsing " << filename << " ..." << std::flush;

    enum State { IDLE, IN_EVENT, IN_BOARD, IN_CHANNEL };
    State state = IDLE;

    EventData current_evt;
    int ch_idx       = -1;
    int sample_count = 0;
    int scaler_idx   = 0;
    std::string line;
    char buf[512];

    while (std::getline(infile, line)) {

        if (!line.empty() && line.back() == '\r') line.pop_back();

        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';

        // --- IDLE: attesa <Event> ---
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

        // --- IN_EVENT: header (Serial, Time, Board) ---
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
                for (int i = 0; i < current_evt.nchannels; ++i) {
                    AnalyzeChannel(current_evt.ch[i]);
                }
                events.push_back(current_evt);
                state = IDLE;
            }
            continue;
        }

        // --- IN_BOARD: trigger cell, scaler, apertura canali ---
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

        // --- IN_CHANNEL: campioni <Data>t,v</Data> ---
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
//  SEZIONE 5: FIT GAUSSIANO SUL CORE 90% (per le grandezze temporali)
// ==========================================================================
//
//  Riusato dal codice TOF_Calibration_v5.cpp. Restituisce una FitResult
//  semplificata (mu, sigma, errori, chi2/ndf, fit_ok).

struct FitOutcome {
    double mu;
    double mu_err;
    double sigma;
    double sigma_err;
    double chi2_ndf;
    bool   fit_ok;
};

/// EstimateHistFWHM(): stima robusta della FWHM direttamente dall'istogramma
/// (insensibile alle code), interpolando linearmente tra bin adiacenti.
static double EstimateHistFWHM(TH1D* h) {

    int max_bin = h->GetMaximumBin();
    double half_max = h->GetBinContent(max_bin) / 2.0;
    int nbins = h->GetNbinsX();

    double x_left = h->GetBinCenter(1);
    for (int i = max_bin; i > 1; --i) {
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

    double x_right = h->GetBinCenter(nbins);
    for (int i = max_bin; i < nbins; ++i) {
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

/// FitGaussianCore(): fit gaussiano ITERATIVO sul core della distribuzione.
///
/// MOTIVAZIONE:
///   La versione precedente faceva un singolo fit a range fisso usando i
///   quantili 5%-95% dell'istogramma. Questo approccio fallisce quando la
///   distribuzione è molto piccata (σ << range Q5-Q95) e ha code non
///   gaussiane (time walk residuo, re-trigger, eventi su CFD ambiguo): il
///   fit insegue le spalle e σ esce sovrastimata di un fattore 2-4.
///
/// METODO ITERATIVO (standard in HEP per la risoluzione temporale):
///   1) Stima iniziale: μ₀ = bin del massimo, σ₀ = FWHM_istogramma / 2.355
///   2) Per ogni iterazione:
///         range = [μᵢ − k·σᵢ , μᵢ + k·σᵢ]
///         fit gaussiano in quel range  →  μᵢ₊₁, σᵢ₊₁
///   3) Convergenza: |Δμ| < 0.01·σ  E  |Δσ| < 0.01·σ
///
///   Il primo passo usa k = 3 (ampio, per "agganciare" il picco), poi k = 2
///   (95% degli eventi gaussiani — concentrato sul core, esclude le code).
///
/// FORMA FUNZIONALE: stessa di prima ([0]·exp(-((x-[1])/[2])²/2)).
///
/// PARAMETRI:
///   h              — istogramma da fittare
///   core_fraction  — usato per calcolare il k della prima iterazione:
///                    k = √2·erfinv(core_fraction). Per 0.90 → k ≈ 1.645.
///                    Mantenuto per compatibilità con la firma vecchia.
///                    Le iterazioni successive forzano k = 2.
///   fit_name       — nome della TF1 (per evitare conflitti)

static FitOutcome FitGaussianCore(TH1D* h,
                                  double core_fraction = 0.90,
                                  const char* fit_name = "gaus_fit") {

    // ---- Parametri dell'algoritmo iterativo ----
    const int    MAX_ITER     = 6;     // iterazioni massime
    const double TOL_REL      = 0.01;  // tolleranza convergenza (in unità di σ)
    const double K_FIRST_ITER = 3.0;   // k della prima iterazione (ampio)
    const double K_REFINE     = 2.0;   // k delle iterazioni successive (stretto)

    // (core_fraction non è più usato direttamente per definire il range, ma
    //  viene mantenuto come argomento per non rompere le chiamate esistenti.
    //  Se in futuro vorrai farne uso, puoi mappare: k = √2·erfinv(core_frac).)
    (void)core_fraction;

    FitOutcome out;
    out.fit_ok = false;
    out.mu = 0; out.mu_err = 0;
    out.sigma = 0; out.sigma_err = 0;
    out.chi2_ndf = -1;

    int nentries = (int)h->GetEntries();
    if (nentries < 30) {
        std::cerr << "[WARNING] Troppi pochi eventi (" << nentries
                  << ") per il fit di " << h->GetName() << std::endl;
        out.mu        = h->GetMean();
        out.mu_err    = h->GetMeanError();
        out.sigma     = h->GetStdDev();
        out.sigma_err = h->GetStdDevError();
        return out;
    }

    // ---- Stima iniziale ROBUSTA ----
    // μ₀ = posizione del massimo, A₀ = altezza del massimo,
    // σ₀ = FWHM stimata robustamente / 2.355.
    int    imax     = h->GetMaximumBin();
    double A0       = h->GetBinContent(imax);
    double mu_i     = h->GetBinCenter(imax);
    double fwhm0    = EstimateHistFWHM(h);
    double sig_i    = (fwhm0 > 0) ? fwhm0 / 2.355 : 0.5 * h->GetStdDev();
    double bin_w    = h->GetBinWidth(1);

    // Sanity floor: σ_init non può essere più piccola di mezza bin width
    if (sig_i < 0.5 * bin_w) sig_i = 0.5 * bin_w;

    // Costruzione UNA SOLA VOLTA della TF1 (la riusiamo cambiando range).
    // Range placeholder: viene aggiornato dentro il loop con SetRange().
    TF1 *fG = new TF1(fit_name,
        "[0]*TMath::Exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",
        h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    fG->SetParName(0, "A_peak");
    fG->SetParName(1, "mu");
    fG->SetParName(2, "sigma");

    // Stato finale (riempito dentro il loop)
    double mu_final = mu_i, mu_err_final = 0;
    double sig_final = sig_i, sig_err_final = 0;
    double chi2_final = -1;
    bool   converged = false;
    int    last_status = -1;

    // ---- LOOP ITERATIVO ----
    for (int it = 0; it < MAX_ITER; ++it) {

        // Scelta del fattore k: ampio alla prima, stretto dopo
        double k = (it == 0) ? K_FIRST_ITER : K_REFINE;
        double x_lo = mu_i - k * sig_i;
        double x_hi = mu_i + k * sig_i;

        // Aggiorna il range della TF1 (NON ricreiamo la TF1 ogni volta)
        fG->SetRange(x_lo, x_hi);

        // Inizializzazione parametri (uso A₀ del MASSIMO dell'istogramma,
        // che è una stima esatta dell'altezza del picco se la gaussiana
        // è ben centrata sul massimo).
        fG->SetParameters(A0, mu_i, sig_i);

        // Limiti morbidi:
        //   A   in [0, 100·A₀]
        //   μ   nel range [x_lo, x_hi]
        //   σ   in [bin_w/2 , (x_hi - x_lo)]   ← evita σ patologiche
        fG->SetParLimits(0, 0.0, 100.0 * A0);
        fG->SetParLimits(1, x_lo, x_hi);
        fG->SetParLimits(2, 0.5 * bin_w, x_hi - x_lo);

        // ---- Esecuzione del fit ----
        // Opzioni: "Q R 0"
        //   Q : silenzia Minuit
        //   R : usa il range della TF1 (cioè il [x_lo, x_hi] che abbiamo appena impostato)
        //   0 : non disegnare automaticamente sul canvas corrente (cosmetico)
        // NB: NON usiamo "L". Il chi-square è più stabile quando i bin
        //     centrali del fit sono ben popolati (caso tipico per le
        //     distribuzioni temporali con O(100) eventi/bin nel core).
        last_status = h->Fit(fG, "Q R 0");

        // Estrai i nuovi parametri
        double mu_new  = fG->GetParameter(1);
        double sig_new = fabs(fG->GetParameter(2));   // sicurezza segno

        // Salva sempre l'ultimo valore (anche se non converge)
        mu_final      = mu_new;
        mu_err_final  = fG->GetParError(1);
        sig_final     = sig_new;
        sig_err_final = fG->GetParError(2);
        chi2_final    = (fG->GetNDF() > 0) ?
                        fG->GetChisquare() / fG->GetNDF() : -1.0;

        // Test di convergenza: ENTRAMBI Δμ e Δσ piccoli rispetto a σ
        double d_mu  = fabs(mu_new  - mu_i);
        double d_sig = fabs(sig_new - sig_i);
        bool   ok    = (d_mu  < TOL_REL * sig_new) &&
                       (d_sig < TOL_REL * sig_new);

        // Aggiorna le stime per la prossima iterazione
        mu_i  = mu_new;
        sig_i = sig_new;

        if (ok && it >= 1) {            // almeno 2 iterazioni prima di dichiarare convergenza
            converged = true;
            break;
        }

        // Sanity: se σ è scappata fuori da limiti ragionevoli, abortisci
        // (può succedere se l'istogramma è patologico o la stima iniziale
        // era totalmente fuori dal picco).
        if (sig_i <= 0 || sig_i > 0.5 * (h->GetXaxis()->GetXmax() -
                                          h->GetXaxis()->GetXmin())) {
            std::cerr << "[WARNING] Fit iterativo: σ patologica per "
                      << h->GetName() << " all'iterazione " << it
                      << " (σ = " << sig_i << "). Esco." << std::endl;
            break;
        }
    }

    // ---- Compilazione del risultato ----
    out.mu        = mu_final;
    out.mu_err    = mu_err_final;
    out.sigma     = sig_final;
    out.sigma_err = sig_err_final;
    out.chi2_ndf  = chi2_final;
    out.fit_ok    = converged && (last_status == 0);

    if (!converged) {
        std::cerr << "[WARNING] Fit iterativo NON convergente per "
                  << h->GetName() << " dopo " << MAX_ITER
                  << " iterazioni. Uso ultima stima: μ=" << mu_final
                  << " σ=" << sig_final << "." << std::endl;
    }

    return out;
}
// ==========================================================================
//  SEZIONE 6: STATISTICHE DESCRITTIVE DA UN ISTOGRAMMA
// ==========================================================================
//
//  Riempie una DistStats con: mean, RMS, mediana, percentili 16/84,
//  + (opzionale) il fit gaussiano sul core 90% per le grandezze temporali.

static DistStats ExtractStats(TH1D* h, bool do_gaussian_fit = false,
                              const char* fit_name = "fit") {
    DistStats s;
    s.n         = (int)h->GetEntries();
    s.mean      = h->GetMean();
    s.mean_err  = h->GetMeanError();
    s.rms       = h->GetStdDev();
    s.rms_err   = h->GetStdDevError();

    // Quantili 16, 50, 84
    double xq[3] = { 0.16, 0.50, 0.84 };
    double yq[3] = { 0.0, 0.0, 0.0 };
    h->GetQuantiles(3, yq, xq);
    s.q16    = yq[0];
    s.median = yq[1];
    s.q84    = yq[2];

    // Fit gaussiano (solo per le grandezze temporali)
    s.fit_done       = false;
    s.fit_mu         = 0;
    s.fit_mu_err     = 0;
    s.fit_sigma      = 0;
    s.fit_sigma_err  = 0;
    s.fit_chi2_ndf   = -1;

    if (do_gaussian_fit && s.n >= 30) {
        FitOutcome out = FitGaussianCore(h, 0.90, fit_name);
        s.fit_mu        = out.mu;
        s.fit_mu_err    = out.mu_err;
        s.fit_sigma     = out.sigma;
        s.fit_sigma_err = out.sigma_err;
        s.fit_chi2_ndf  = out.chi2_ndf;
        s.fit_done      = out.fit_ok;
    }
    return s;
}


// ==========================================================================
//  SEZIONE 7: ELABORAZIONE DI UN SINGOLO DATASET (RUN)
// ==========================================================================
//
//  ProcessRun(): per un dato file XML, parsa, applica i tagli di qualità,
//  riempie un TTree event-by-event e i 12 istogrammi richiesti, e calcola
//  le statistiche descrittive di ogni grandezza.
//
//  Mappatura canali (come nel codice di calibrazione):
//    CH1 (idx 0) → PMT1
//    CH2 (idx 1) → PMT2
//    CH3 (idx 2) → PMT3 (riferimento)
//
//  Logica di scarto:
//    - "good_12":   evento OK su CH1 e CH2 (per le grandezze di canale e Δt₁₂)
//    - "good_123":  evento OK su CH1, CH2 e CH3 (per t₃-(t₁+t₂)/2)
//
//  Le grandezze di singolo canale (amp, SNR, rise time, slew, t_cfd) sono
//  riempite ognuna con la quality del PROPRIO canale (un evento può essere
//  buono per CH1 e cattivo per CH2, per esempio).

/// Contenitore dei risultati della run completa.
struct RunResults {
    // 12 istogrammi di output
    TH1D* h_amp[2];        // [0] = CH1, [1] = CH2
    TH1D* h_snr[2];
    TH1D* h_rt[2];         // rise time 20-80
    TH1D* h_sr[2];         // slew rate 20-80
    TH1D* h_tcfd[2];       // tempo CFD assoluto al 20%
    TH1D* h_dt12;          // t_1 - t_2
    TH1D* h_C;             // t_3 - (t_1+t_2)/2

    // Statistiche descrittive
    DistStats st_amp[2];
    DistStats st_snr[2];
    DistStats st_rt[2];
    DistStats st_sr[2];
    DistStats st_tcfd[2];
    DistStats st_dt12;
    DistStats st_C;

    // Contatori di qualità (utili per il riepilogo finale)
    int n_total;          // eventi parsati
    int n_no_cfd_ch1;     // CFD non trovato su CH1
    int n_no_cfd_ch2;
    int n_no_cfd_ch3;
    int n_clip_ch1;       // saturazione su CH1
    int n_clip_ch2;
    int n_clip_ch3;
    int n_osc_ch1;        // oscillazione su CH1
    int n_osc_ch2;
    int n_osc_ch3;
    int n_good_ch1;       // CH1 buono per le sue grandezze di canale
    int n_good_ch2;       // CH2 buono per le sue grandezze
    int n_good_12;        // CH1 e CH2 buoni (per Δt₁₂)
    int n_good_123;       // CH1, CH2, CH3 buoni (per t_3-(t_1+t_2)/2)
};

/// ProcessRun(): legge file XML, riempie TTree e istogrammi, fa i fit.
///
/// PARAMETRI:
///   xml_path  — path al file XML della DRS4
///   tag       — tag breve che identifica il dataset ("with" o "without")
///   outfile   — file ROOT di output dove salvare TTree e istogrammi
///
/// RETURN: RunResults popolata.

static RunResults ProcessRun(const char* xml_path,
                             const std::string &tag,
                             TFile* outfile) {

    RunResults res;
    res.n_total = 0;
    res.n_no_cfd_ch1 = res.n_no_cfd_ch2 = res.n_no_cfd_ch3 = 0;
    res.n_clip_ch1   = res.n_clip_ch2   = res.n_clip_ch3   = 0;
    res.n_osc_ch1    = res.n_osc_ch2    = res.n_osc_ch3    = 0;
    res.n_good_ch1   = res.n_good_ch2   = 0;
    res.n_good_12    = res.n_good_123   = 0;

    // ---- Parsing del file XML ----
    std::vector<EventData> events;
    int n_parsed = ParseXML(xml_path, events);
    if (n_parsed <= 0) {
        std::cerr << "[ERRORE] Nessun evento letto da " << xml_path << std::endl;
        // Inizializza i puntatori a nullptr per evitare segfault
        for (int k = 0; k < 2; ++k) {
            res.h_amp[k] = res.h_snr[k] = res.h_rt[k] =
            res.h_sr[k]  = res.h_tcfd[k] = nullptr;
        }
        res.h_dt12 = res.h_C = nullptr;
        return res;
    }
    res.n_total = (int)events.size();

    // ---- Creazione degli istogrammi ----
    // Stesso binning per ATT e NOATT (cruciale per il confronto).
    // I nomi includono il tag per evitare conflitti nel file ROOT.
    outfile->cd();
    const char* chn_label[2] = { "ch1", "ch2" };

    for (int k = 0; k < 2; ++k) {
        res.h_amp[k]  = new TH1D(Form("h_amp_%s_%s",  chn_label[k], tag.c_str()),
                                 Form("Ampiezza %s [%s];Ampiezza [mV];Conteggi",
                                      chn_label[k], tag.c_str()),
                                 NB_AMP, AMP_LO, AMP_HI);
        res.h_snr[k]  = new TH1D(Form("h_snr_%s_%s",  chn_label[k], tag.c_str()),
                                 Form("SNR %s [%s];SNR;Conteggi",
                                      chn_label[k], tag.c_str()),
                                 NB_SNR, SNR_LO, SNR_HI);
        res.h_rt[k]   = new TH1D(Form("h_rt_%s_%s",   chn_label[k], tag.c_str()),
                                 Form("Rise time 20-80%% %s [%s];Rise time [ns];Conteggi",
                                      chn_label[k], tag.c_str()),
                                 NB_RT, RT_LO, RT_HI);
        res.h_sr[k]   = new TH1D(Form("h_sr_%s_%s",   chn_label[k], tag.c_str()),
                                 Form("Slew rate 20-80%% %s [%s];|Slew rate| [mV/ns];Conteggi",
                                      chn_label[k], tag.c_str()),
                                 NB_SR, SR_LO, SR_HI);
        res.h_tcfd[k] = new TH1D(Form("h_tcfd_%s_%s", chn_label[k], tag.c_str()),
                                 Form("t_{CFD,20%%} %s [%s];t_{CFD} [ns];Conteggi",
                                      chn_label[k], tag.c_str()),
                                 NB_TCFD, TCFD_LO, TCFD_HI);
    }

    res.h_dt12 = new TH1D(Form("h_dt12_%s", tag.c_str()),
                          Form("#Delta t_{12} = t_{1} - t_{2} [%s];#Delta t_{12} [ns];Conteggi",
                               tag.c_str()),
                          NB_DT12, DT12_LO, DT12_HI);
    res.h_C    = new TH1D(Form("h_C_%s", tag.c_str()),
                          Form("t_{3} - (t_{1}+t_{2})/2 [%s];"
                               "t_{3} - (t_{1}+t_{2})/2 [ns];Conteggi",
                               tag.c_str()),
                          NB_C, C_LO, C_HI);

    // ---- Creazione del TTree event-by-event ----
    // Salviamo TUTTE le grandezze utili: questo permette di rifare l'analisi
    // (cambiare CFD, range, tagli) senza dover ri-parsare l'XML.
    TTree* tree = new TTree(Form("tree_%s", tag.c_str()),
                            Form("Dataset %s (event-by-event)", tag.c_str()));

    // Branch per le grandezze per canale
    Float_t amp_v[3], snr_v[3], rt_v[3], sr_v[3], tcfd_v[3];
    Float_t bl_v[3], bl_rms_v[3], vmin_v[3];
    Int_t   clip_v[3], osc_v[3], cfd_ok_v[3], rise_ok_v[3];

    tree->Branch("amp",      amp_v,     "amp[3]/F");      // [mV]
    tree->Branch("snr",      snr_v,     "snr[3]/F");      // adimensionale
    tree->Branch("rise_time", rt_v,     "rise_time[3]/F"); // [ns]
    tree->Branch("slew_rate", sr_v,     "slew_rate[3]/F"); // [mV/ns]
    tree->Branch("t_cfd",    tcfd_v,    "t_cfd[3]/F");    // [ns]
    tree->Branch("baseline", bl_v,      "baseline[3]/F"); // [mV]
    tree->Branch("baseline_rms", bl_rms_v, "baseline_rms[3]/F"); // [mV]
    tree->Branch("vmin",     vmin_v,    "vmin[3]/F");     // [mV]
    tree->Branch("clip",     clip_v,    "clip[3]/I");
    tree->Branch("osc",      osc_v,     "osc[3]/I");
    tree->Branch("cfd_ok",   cfd_ok_v,  "cfd_ok[3]/I");
    tree->Branch("rise_ok",  rise_ok_v, "rise_ok[3]/I");

    // Differenze temporali derivate
    Float_t dt12, t3_avg12;
    tree->Branch("dt12",     &dt12,     "dt12/F");        // t_1 - t_2 [ns]
    tree->Branch("t3_avg12", &t3_avg12, "t3_avg12/F");    // t_3 - (t_1+t_2)/2 [ns]

    // Flag aggregati
    Int_t good_12, good_123;
    tree->Branch("good_12",  &good_12,  "good_12/I");
    tree->Branch("good_123", &good_123, "good_123/I");

    // ---- Loop sugli eventi ----
    for (size_t ev = 0; ev < events.size(); ++ev) {

        EventData &e = events[ev];

        // Inizializza tutto a "non valido"
        for (int k = 0; k < 3; ++k) {
            amp_v[k] = snr_v[k] = rt_v[k] = sr_v[k] = -999;
            tcfd_v[k] = bl_v[k] = bl_rms_v[k] = vmin_v[k] = -999;
            clip_v[k] = osc_v[k] = cfd_ok_v[k] = rise_ok_v[k] = 0;
        }
        dt12 = t3_avg12 = -999;
        good_12 = good_123 = 0;

        if (e.nchannels < 2) {
            tree->Fill();
            continue;
        }

        // Quanti canali abbiamo? (almeno 2, idealmente 3)
        int nch = std::min(3, e.nchannels);

        // Estrai grandezze e flag per ogni canale
        for (int k = 0; k < nch; ++k) {
            amp_v[k]      = (Float_t)e.ch[k].amplitude;
            snr_v[k]      = (Float_t)e.ch[k].snr;
            rt_v[k]       = (Float_t)e.ch[k].rise_time;
            sr_v[k]       = (Float_t)e.ch[k].slew_rate;
            tcfd_v[k]     = (Float_t)e.ch[k].t_cfd;
            bl_v[k]       = (Float_t)e.ch[k].baseline;
            bl_rms_v[k]   = (Float_t)e.ch[k].baseline_rms;
            vmin_v[k]     = (Float_t)e.ch[k].v_min;
            clip_v[k]     = e.ch[k].is_clipped     ? 1 : 0;
            osc_v[k]      = e.ch[k].is_oscillating ? 1 : 0;
            cfd_ok_v[k]   = e.ch[k].cfd_ok         ? 1 : 0;
            rise_ok_v[k]  = e.ch[k].rise_ok        ? 1 : 0;
        }

        // Aggiorna i contatori (per ogni canale, conta le ragioni di scarto)
        for (int k = 0; k < nch; ++k) {
            if (!e.ch[k].cfd_ok) {
                if (k == 0) res.n_no_cfd_ch1++;
                if (k == 1) res.n_no_cfd_ch2++;
                if (k == 2) res.n_no_cfd_ch3++;
            }
            if (e.ch[k].is_clipped) {
                if (k == 0) res.n_clip_ch1++;
                if (k == 1) res.n_clip_ch2++;
                if (k == 2) res.n_clip_ch3++;
            }
            if (e.ch[k].is_oscillating) {
                if (k == 0) res.n_osc_ch1++;
                if (k == 1) res.n_osc_ch2++;
                if (k == 2) res.n_osc_ch3++;
            }
        }

        // Definizione di canale "buono":
        //  - cfd_ok (se ENABLE_CFD_CUT)
        //  - !is_clipped (se ENABLE_CLIP_CUT)
        //  - !is_oscillating (se ENABLE_OSC_CUT)
        // L'evento può essere "buono" per un canale e non per un altro.
        bool ch_good[3] = { false, false, false };
        for (int k = 0; k < nch; ++k) {
            ch_good[k] = true;
            if (ENABLE_CFD_CUT  && !e.ch[k].cfd_ok)        ch_good[k] = false;
            if (ENABLE_CLIP_CUT &&  e.ch[k].is_clipped)    ch_good[k] = false;
            if (ENABLE_OSC_CUT  &&  e.ch[k].is_oscillating) ch_good[k] = false;
        }

        // ---- Riempimento istogrammi per le grandezze di canale ----
        // Per ognuna delle grandezze di canale, riempiamo SOLO se quel
        // canale è buono. Inoltre per rise time e slew rate richiediamo
        // ANCHE rise_ok (entrambi i crossing 20% e 80% trovati).
        for (int k = 0; k < 2; ++k) {     // solo CH1 e CH2 nei plot
            if (!ch_good[k]) continue;

            res.h_amp[k]->Fill(e.ch[k].amplitude);
            res.h_snr[k]->Fill(e.ch[k].snr);
            res.h_tcfd[k]->Fill(e.ch[k].t_cfd);

            if (e.ch[k].rise_ok) {
                res.h_rt[k]->Fill(e.ch[k].rise_time);
                res.h_sr[k]->Fill(e.ch[k].slew_rate);
            }
        }

        // Conteggio per-canale dei "buoni"
        if (ch_good[0]) res.n_good_ch1++;
        if (ch_good[1]) res.n_good_ch2++;

        // ---- Differenze temporali ----
        // Δt_12 richiede CH1 e CH2 buoni
        if (ch_good[0] && ch_good[1]) {
            dt12 = (Float_t)(e.ch[0].t_cfd - e.ch[1].t_cfd);
            res.h_dt12->Fill(dt12);
            good_12 = 1;
            res.n_good_12++;
        }

        // t_3 - (t_1+t_2)/2 richiede tutti e 3 i canali buoni
        if (nch >= 3 && ch_good[0] && ch_good[1] && ch_good[2]) {
            t3_avg12 = (Float_t)(e.ch[2].t_cfd - 0.5*(e.ch[0].t_cfd + e.ch[1].t_cfd));
            res.h_C->Fill(t3_avg12);
            good_123 = 1;
            res.n_good_123++;
        }

        tree->Fill();
    }

    // ---- Stampa diagnostica ----
    std::cout << "[INFO] Run '" << tag << "' completata: "
              << res.n_total << " eventi totali" << std::endl;
    std::cout << "       Buoni: CH1=" << res.n_good_ch1
              << "  CH2=" << res.n_good_ch2
              << "  CH1&CH2=" << res.n_good_12
              << "  CH1&CH2&CH3=" << res.n_good_123 << std::endl;
    std::cout << "       Scartati per CFD non trovato: "
              << res.n_no_cfd_ch1 << "/" << res.n_no_cfd_ch2 << "/"
              << res.n_no_cfd_ch3 << " (CH1/CH2/CH3)" << std::endl;
    std::cout << "       Scartati per clipping:        "
              << res.n_clip_ch1 << "/" << res.n_clip_ch2 << "/"
              << res.n_clip_ch3 << " (CH1/CH2/CH3)" << std::endl;
    std::cout << "       Scartati per oscillazione:    "
              << res.n_osc_ch1 << "/" << res.n_osc_ch2 << "/"
              << res.n_osc_ch3 << " (CH1/CH2/CH3)" << std::endl;

    // ---- Estrazione delle statistiche descrittive ----
    // Per le grandezze "ampiezza-like" (amp, SNR, rise time, slew rate):
    // statistiche grezze (mean, RMS, mediana, percentili). NON faccio fit
    // perché queste distribuzioni non sono gaussiane (Landau, Cauchy, etc).
    //
    // Per le grandezze temporali (t_cfd, Δt_12, t_3-(t_1+t_2)/2): faccio
    // anche il fit gaussiano sul core 90% per estrarre μ ± σ.
    for (int k = 0; k < 2; ++k) {
        res.st_amp[k]  = ExtractStats(res.h_amp[k],  false);
        res.st_snr[k]  = ExtractStats(res.h_snr[k],  false);
        res.st_rt[k]   = ExtractStats(res.h_rt[k],   false);
        res.st_sr[k]   = ExtractStats(res.h_sr[k],   false);
        res.st_tcfd[k] = ExtractStats(res.h_tcfd[k], true,
                                       Form("fit_tcfd_%s_%s",
                                            chn_label[k], tag.c_str()));
    }
    res.st_dt12 = ExtractStats(res.h_dt12, true,
                               Form("fit_dt12_%s", tag.c_str()));
    res.st_C    = ExtractStats(res.h_C,    true,
                               Form("fit_C_%s",    tag.c_str()));

    // ---- Salvataggio ----
    outfile->cd();
    tree->Write();
    for (int k = 0; k < 2; ++k) {
        res.h_amp[k]->Write();
        res.h_snr[k]->Write();
        res.h_rt[k]->Write();
        res.h_sr[k]->Write();
        res.h_tcfd[k]->Write();
    }
    res.h_dt12->Write();
    res.h_C->Write();

    return res;
}


// ==========================================================================
//  SEZIONE 8: CONFRONTO E PLOT DI OVERLAY
// ==========================================================================
//
//  Per ognuna delle 12 distribuzioni si genera un canvas con due pad:
//    - sinistra: istogrammi sovrapposti in conteggi (statistica grezza)
//    - destra:   istogrammi sovrapposti normalizzati (forma)
//  Più una legenda con le statistiche e una box con KS-test e
//  rapporto delle medie.

/// MakeOverlayCanvas(): crea un canvas TCanvas con due pad confrontando
/// l'istogramma "with attenuator" (rosso) con "without" (blu).
/// Restituisce il puntatore al canvas (già scritto su file).
///
/// PARAMETRI:
///   h_with        — istogramma con attenuatore
///   h_without     — istogramma senza attenuatore
///   st_with       — statistiche dataset ATT
///   st_without    — statistiche dataset NOATT
///   canvas_name   — nome del canvas (es. "c_amp_ch1")
///   canvas_title  — titolo del canvas
///   x_label       — etichetta dell'asse X (con unità)
///   show_fit      — se true mostra il risultato del fit gaussiano
///   outfile       — file di output dove salvare il canvas
static TCanvas* MakeOverlayCanvas(TH1D* h_with, TH1D* h_without,
                                  const DistStats &st_with,
                                  const DistStats &st_without,
                                  const char* canvas_name,
                                  const char* canvas_title,
                                  const char* x_label,
                                  bool show_fit,
                                  TFile* outfile)
{
    if (!h_with || !h_without) return nullptr;

    TCanvas* c = new TCanvas(canvas_name, canvas_title, 1200, 500);
    c->Divide(2, 1);

    // ------------------------------------------------------------------
    // PAD 1 — SOVRAPPOSIZIONE IN CONTEGGI (statistica grezza)
    // ------------------------------------------------------------------
    c->cd(1);
    gPad->SetGrid(1, 1);

    // Stile: linee colorate, fill semi-trasparente
    h_with   ->SetLineColor(kRed + 1);
    h_with   ->SetLineWidth(2);
    h_with   ->SetFillColorAlpha(kRed - 7, 0.3);

    h_without->SetLineColor(kBlue + 1);
    h_without->SetLineWidth(2);
    h_without->SetFillColorAlpha(kBlue - 7, 0.3);

    // Determina il massimo per scalare l'asse Y correttamente
    double ymax = std::max(h_with->GetMaximum(), h_without->GetMaximum());
    h_with->GetYaxis()->SetRangeUser(0, 1.2 * ymax);
    h_with->GetXaxis()->SetTitle(x_label);
    h_with->SetTitle(Form("%s — Conteggi", canvas_title));

    h_with   ->Draw("HIST");
    h_without->Draw("HIST SAME");

    // Legenda con statistiche
    TLegend* leg1 = new TLegend(0.55, 0.60, 0.93, 0.90);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.035);
    leg1->SetFillColorAlpha(kWhite, 0.85);
    leg1->AddEntry(h_with,
        Form("ATT  N=%d  #mu=%.3g  RMS=%.3g",
             st_with.n, st_with.mean, st_with.rms), "lf");
    leg1->AddEntry(h_without,
        Form("NOATT  N=%d  #mu=%.3g  RMS=%.3g",
             st_without.n, st_without.mean, st_without.rms), "lf");
    leg1->Draw();

    // ------------------------------------------------------------------
    // PAD 2 — SOVRAPPOSIZIONE NORMALIZZATA (forma della distribuzione)
    // ------------------------------------------------------------------
    c->cd(2);
    gPad->SetGrid(1, 1);

    // Cloni normalizzati (così non alteriamo gli originali nel file)
    TH1D* h_with_n    = (TH1D*)h_with   ->Clone(Form("%s_norm", h_with   ->GetName()));
    TH1D* h_without_n = (TH1D*)h_without->Clone(Form("%s_norm", h_without->GetName()));

    if (h_with_n   ->Integral() > 0) h_with_n   ->Scale(1.0 / h_with_n   ->Integral("width"));
    if (h_without_n->Integral() > 0) h_without_n->Scale(1.0 / h_without_n->Integral("width"));

    // Re-imposta lo stile (lo Scale resetta alcune cose)
    h_with_n   ->SetLineColor(kRed + 1);
    h_with_n   ->SetLineWidth(2);
    h_with_n   ->SetFillColorAlpha(kRed - 7, 0.3);
    h_without_n->SetLineColor(kBlue + 1);
    h_without_n->SetLineWidth(2);
    h_without_n->SetFillColorAlpha(kBlue - 7, 0.3);

    double ymax_n = std::max(h_with_n->GetMaximum(), h_without_n->GetMaximum());
    h_with_n->GetYaxis()->SetRangeUser(0, 1.2 * ymax_n);
    h_with_n->GetYaxis()->SetTitle("PDF (normalizzata)");
    h_with_n->GetXaxis()->SetTitle(x_label);
    h_with_n->SetTitle(Form("%s — Normalizzata", canvas_title));

    h_with_n   ->Draw("HIST");
    h_without_n->Draw("HIST SAME");

    // ------------------------------------------------------------------
    // BOX CON CONFRONTO QUANTITATIVO
    // ------------------------------------------------------------------
    // KS-test a due campioni: TH1::KolmogorovTest restituisce la
    // probabilità (p-value) che le due distribuzioni siano consistenti.
    // p < 0.01 → distribuzioni significativamente diverse.
    double ks_p = h_with->KolmogorovTest(h_without);

    // Rapporto delle medie con propagazione degli errori:
    //   R = mean_A / mean_B
    //   σ_R/R = sqrt((σ_A/mean_A)² + (σ_B/mean_B)²)
    double R = 0, R_err = 0;
    if (fabs(st_without.mean) > 1e-12) {
        R = st_with.mean / st_without.mean;
        double rel_err_with    = (fabs(st_with.mean) > 1e-12) ?
                                 st_with.mean_err / fabs(st_with.mean) : 0;
        double rel_err_without = st_without.mean_err / fabs(st_without.mean);
        R_err = fabs(R) * sqrt(rel_err_with*rel_err_with +
                               rel_err_without*rel_err_without);
    }

    // Differenza delle medie (per le grandezze temporali, atteso ~0):
    double dMean = st_with.mean - st_without.mean;
    double dMean_err = sqrt(st_with.mean_err*st_with.mean_err +
                            st_without.mean_err*st_without.mean_err);

    TPaveText* pt = new TPaveText(0.55, 0.55, 0.93, 0.92, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.85);
    pt->SetLineColor(kBlack);
    pt->SetTextFont(42);
    pt->SetTextSize(0.033);
    pt->SetTextAlign(12);
    // pt->AddText(Form("KS p-value = %.3g", ks_p));
    pt->AddText(Form("R = #mu_{ATT}/#mu_{NOATT} = %.4f #pm %.4f", R, R_err));
    // pt->AddText(Form("#Delta = #mu_{ATT} - #mu_{NOATT} = %.4g #pm %.4g",
    //                  dMean, dMean_err));
    if (show_fit && st_with.fit_done && st_without.fit_done) {
        pt->AddText("");
        pt->AddText(Form("Fit ATT:    #mu=%.4g  #sigma=%.4g",
                         st_with.fit_mu, st_with.fit_sigma));
        pt->AddText(Form("Fit NOATT: #mu=%.4g  #sigma=%.4g",
                         st_without.fit_mu, st_without.fit_sigma));
        pt->AddText(Form("Rapporto #sigma: %.3f",
                         (st_without.fit_sigma > 1e-12) ?
                         st_with.fit_sigma / st_without.fit_sigma : 0.0));
    }
    pt->Draw();

    // Legenda compatta sul pad normalizzato
    // TLegend* leg2 = new TLegend(0.13, 0.78, 0.50, 0.92);
    // leg2->SetTextFont(42);
    // leg2->SetTextSize(0.035);
    // leg2->SetFillColorAlpha(kWhite, 0.85);
    // leg2->AddEntry(h_with_n,    "ATT (rosso)",    "lf");
    // leg2->AddEntry(h_without_n, "NOATT (blu)", "lf");
    // leg2->Draw();

    c->Update();
    outfile->cd();
    c->Write();

    return c;
}


// ==========================================================================
//  SEZIONE 9: TABELLA RIASSUNTIVA
// ==========================================================================
//
//  Stampa a video e salva nel file ROOT una tabella di confronto numerica
//  con tutte le grandezze, i loro valori medi e i rapporti ATT/NOATT.

/// Stampa una riga della tabella con formato pulito.
static void PrintTableRow(const std::string &label,
                          const std::string &unit,
                          const DistStats &a,        // ATT
                          const DistStats &b,        // NOATT
                          double expected_ratio = 0.0,  // 0 = non stampare
                          bool   show_fit = false)
{
    // Rapporto e differenza con propagazione errori
    double R = 0, R_err = 0;
    if (fabs(b.mean) > 1e-12) {
        R = a.mean / b.mean;
        double r1 = (fabs(a.mean) > 1e-12) ? a.mean_err / fabs(a.mean) : 0;
        double r2 = b.mean_err / fabs(b.mean);
        R_err = fabs(R) * sqrt(r1*r1 + r2*r2);
    }
    double dMean = a.mean - b.mean;
    double dMean_err = sqrt(a.mean_err*a.mean_err + b.mean_err*b.mean_err);

    std::cout << "  " << std::left << std::setw(28) << (label + " [" + unit + "]")
              << std::right
              << "  ATT: " << std::setw(10) << std::fixed << std::setprecision(4) << a.mean
              << " #" << std::setw(2) << a.n
              << "  NOATT: " << std::setw(10) << std::fixed << std::setprecision(4) << b.mean
              << " #" << std::setw(2) << b.n
              << "  R="    << std::setw(8) << std::fixed << std::setprecision(4) << R
              << " ±" << std::setw(7) << std::fixed << std::setprecision(4) << R_err;
    if (expected_ratio > 0) {
        std::cout << " (atteso " << std::setprecision(2) << expected_ratio << ")";
    }
    std::cout << std::endl;

    if (show_fit && a.fit_done && b.fit_done) {
        std::cout << "    fit core 90%:  ATT μ=" << std::fixed << std::setprecision(4)
                  << a.fit_mu << "±" << a.fit_mu_err
                  << " σ=" << a.fit_sigma << "±" << a.fit_sigma_err
                  << "  |  NOATT μ=" << b.fit_mu << "±" << b.fit_mu_err
                  << " σ=" << b.fit_sigma << "±" << b.fit_sigma_err << std::endl;
        std::cout << "    Δμ_fit = " << (a.fit_mu - b.fit_mu)
                  << " ± " << sqrt(a.fit_mu_err*a.fit_mu_err + b.fit_mu_err*b.fit_mu_err)
                  << "   σ_ratio = "
                  << ((b.fit_sigma > 1e-12) ? a.fit_sigma/b.fit_sigma : 0.0) << std::endl;
    }
}

/// Stampa l'intera tabella di confronto e salva un TTree summary nel file.
static void PrintAndSaveSummary(const RunResults &A,    // ATT
                                const RunResults &B,    // NOATT
                                TFile* outfile)
{
    std::cout << "\n=================================================================" << std::endl;
    std::cout << "  TABELLA RIASSUNTIVA: confronto ATT vs NOATT" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "  Eventi totali:  ATT = " << A.n_total
              << "   NOATT = " << B.n_total << std::endl;
    std::cout << "  Buoni (CH1&CH2): ATT = " << A.n_good_12
              << "   NOATT = " << B.n_good_12 << std::endl;
    std::cout << "  Buoni (tutti):   ATT = " << A.n_good_123
              << "   NOATT = " << B.n_good_123 << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;

    // Grandezze di canale
    std::cout << "\n[Grandezze di canale: rapporto atteso ATT/NOATT ≈ 0.5012 per amp/SNR/slew, ≈ 1.0 per rise time]\n";
    PrintTableRow("Ampiezza CH1",    "mV",     A.st_amp[0],  B.st_amp[0],  0.5012);
    PrintTableRow("Ampiezza CH2",    "mV",     A.st_amp[1],  B.st_amp[1],  0.5012);
    PrintTableRow("SNR CH1",         "-",      A.st_snr[0],  B.st_snr[0],  0.5012);
    PrintTableRow("SNR CH2",         "-",      A.st_snr[1],  B.st_snr[1],  0.5012);
    PrintTableRow("Rise time CH1",   "ns",     A.st_rt[0],   B.st_rt[0],   1.0);
    PrintTableRow("Rise time CH2",   "ns",     A.st_rt[1],   B.st_rt[1],   1.0);
    PrintTableRow("Slew rate CH1",   "mV/ns",  A.st_sr[0],   B.st_sr[0],   0.5012);
    PrintTableRow("Slew rate CH2",   "mV/ns",  A.st_sr[1],   B.st_sr[1],   0.5012);

    // Grandezze temporali (con fit)
    std::cout << "\n[Grandezze temporali: Δμ atteso ≈ 0; il rapporto σ informa su come l'attenuatore cambia il jitter]\n";
    PrintTableRow("t_CFD CH1",       "ns",     A.st_tcfd[0], B.st_tcfd[0], 0.0, true);
    PrintTableRow("t_CFD CH2",       "ns",     A.st_tcfd[1], B.st_tcfd[1], 0.0, true);
    PrintTableRow("Δt_12",           "ns",     A.st_dt12,    B.st_dt12,    0.0, true);
    PrintTableRow("t_3-(t_1+t_2)/2", "ns",     A.st_C,       B.st_C,       0.0, true);
    std::cout << "=================================================================\n" << std::endl;

    // ---- Salvataggio del summary in un TTree con una sola entry ----
    outfile->cd();
    TTree* summary = new TTree("summary",
                                "Riepilogo confronto ATT vs NOATT");

    // Per ogni grandezza: mean ATT, mean NOATT, rms ATT, rms NOATT,
    // ratio, ratio_err, dMean, dMean_err.
    // Per i tempi: anche fit_mu e fit_sigma.

    // Variabili scalari (creo un blocco unico)
    Float_t v[200];
    Char_t  name[64];
    int idx = 0;

    auto add_quantity = [&](const char* base,
                             const DistStats &a, const DistStats &b,
                             bool with_fit) {
        // Indici di partenza per questa quantità
        int i0 = idx;
        v[idx++] = (Float_t)a.mean;       v[idx++] = (Float_t)a.mean_err;
        v[idx++] = (Float_t)a.rms;        v[idx++] = (Float_t)a.rms_err;
        v[idx++] = (Float_t)a.median;     v[idx++] = (Float_t)a.q16;       v[idx++] = (Float_t)a.q84;
        v[idx++] = (Float_t)a.n;
        v[idx++] = (Float_t)b.mean;       v[idx++] = (Float_t)b.mean_err;
        v[idx++] = (Float_t)b.rms;        v[idx++] = (Float_t)b.rms_err;
        v[idx++] = (Float_t)b.median;     v[idx++] = (Float_t)b.q16;       v[idx++] = (Float_t)b.q84;
        v[idx++] = (Float_t)b.n;
        if (with_fit) {
            v[idx++] = (Float_t)a.fit_mu;     v[idx++] = (Float_t)a.fit_mu_err;
            v[idx++] = (Float_t)a.fit_sigma;  v[idx++] = (Float_t)a.fit_sigma_err;
            v[idx++] = (Float_t)b.fit_mu;     v[idx++] = (Float_t)b.fit_mu_err;
            v[idx++] = (Float_t)b.fit_sigma;  v[idx++] = (Float_t)b.fit_sigma_err;
        }
        // Crea i branch corrispondenti
        const char* fields_basic[] = {
            "att_mean","att_mean_err","att_rms","att_rms_err",
            "att_median","att_q16","att_q84","att_n",
            "noatt_mean","noatt_mean_err","noatt_rms","noatt_rms_err",
            "noatt_median","noatt_q16","noatt_q84","noatt_n"
        };
        const char* fields_fit[] = {
            "att_fit_mu","att_fit_mu_err","att_fit_sigma","att_fit_sigma_err",
            "noatt_fit_mu","noatt_fit_mu_err","noatt_fit_sigma","noatt_fit_sigma_err"
        };
        int nb = 16;
        int nf = with_fit ? 8 : 0;
        for (int k = 0; k < nb; ++k) {
            std::string bname = std::string(base) + "_" + fields_basic[k];
            summary->Branch(bname.c_str(), &v[i0 + k], (bname + "/F").c_str());
        }
        for (int k = 0; k < nf; ++k) {
            std::string bname = std::string(base) + "_" + fields_fit[k];
            summary->Branch(bname.c_str(), &v[i0 + nb + k], (bname + "/F").c_str());
        }
    };

    add_quantity("amp_ch1",  A.st_amp[0],  B.st_amp[0],  false);
    add_quantity("amp_ch2",  A.st_amp[1],  B.st_amp[1],  false);
    add_quantity("snr_ch1",  A.st_snr[0],  B.st_snr[0],  false);
    add_quantity("snr_ch2",  A.st_snr[1],  B.st_snr[1],  false);
    add_quantity("rt_ch1",   A.st_rt[0],   B.st_rt[0],   false);
    add_quantity("rt_ch2",   A.st_rt[1],   B.st_rt[1],   false);
    add_quantity("sr_ch1",   A.st_sr[0],   B.st_sr[0],   false);
    add_quantity("sr_ch2",   A.st_sr[1],   B.st_sr[1],   false);
    add_quantity("tcfd_ch1", A.st_tcfd[0], B.st_tcfd[0], true);
    add_quantity("tcfd_ch2", A.st_tcfd[1], B.st_tcfd[1], true);
    add_quantity("dt12",     A.st_dt12,    B.st_dt12,    true);
    add_quantity("C",        A.st_C,       B.st_C,       true);

    summary->Fill();
    summary->Write();
}


// ==========================================================================
//  SEZIONE 10: FUNZIONE PRINCIPALE
// ==========================================================================

/// Attenuator_Study(): funzione principale.
///
/// PARAMETRI:
///   xml_with    — path al file XML acquisito CON attenuatore -6 dB
///   xml_without — path al file XML acquisito SENZA attenuatore
///   outname     — nome del file ROOT di output
///                 (default: "Attenuator_Study_output.root")
///
/// UTILIZZO:
///   root -l 'Attenuator_Study.cpp("/path/with_att.xml", "/path/without_att.xml")'
///
///   Oppure dalla console ROOT:
///     .L Attenuator_Study.cpp
///     Attenuator_Study("with.xml", "without.xml");

void Attenuator_Study(const char* xml_with,
                      const char* xml_without,
                      const char* outname = "Attenuator_Study_output.root")
{
    std::cout << "================================================================" << std::endl;
    std::cout << "  Attenuator_Study — confronto segnali con/senza attenuatore" << std::endl;
    std::cout << "  CON attenuatore (-6 dB): " << xml_with << std::endl;
    std::cout << "  SENZA attenuatore:        " << xml_without << std::endl;
    std::cout << "  Output:                   " << outname << std::endl;
    std::cout << "  Frazione CFD:             " << CFD_FRACTION
              << "  (rise time: " << RT_LOW_FRAC << "-" << RT_HIGH_FRAC << ")" << std::endl;
    std::cout << "================================================================" << std::endl;

    // Setup stile grafico
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetTitleSize(0.045, "xy");
    gStyle->SetTitleOffset(1.1, "y");

    // Apri il file di output ROOT
    TFile* fout = new TFile(outname, "RECREATE");
    if (!fout->IsOpen()) {
        std::cerr << "[ERRORE] Impossibile creare " << outname << std::endl;
        return;
    }

    // ==================================================================
    //  PROCESSING DEI DUE DATASET
    // ==================================================================
    std::cout << "\n--- Processing dataset CON attenuatore ---" << std::endl;
    RunResults A = ProcessRun(xml_with,    "with",    fout);

    std::cout << "\n--- Processing dataset SENZA attenuatore ---" << std::endl;
    RunResults B = ProcessRun(xml_without, "without", fout);

    if (!A.h_amp[0] || !B.h_amp[0]) {
        std::cerr << "[ERRORE] Uno dei dataset non ha prodotto istogrammi validi."
                  << std::endl;
        fout->Close();
        return;
    }

    // ==================================================================
    //  CANVAS DI CONFRONTO (overlay ATT vs NOATT)
    // ==================================================================
    std::cout << "\n--- Generazione canvas di confronto ---" << std::endl;

    // Per ogni grandezza, creiamo un canvas con due pad (counts + normalizzato)
    // e una box con KS-test, rapporto medie, e (per le grandezze temporali)
    // i parametri del fit gaussiano.
    const char* chn_label[2] = { "CH1", "CH2" };

    for (int k = 0; k < 2; ++k) {
        MakeOverlayCanvas(A.h_amp[k], B.h_amp[k], A.st_amp[k], B.st_amp[k],
                          Form("c_amp_ch%d", k+1),
                          Form("Ampiezza %s (rispetto a baseline)", chn_label[k]),
                          "Ampiezza [mV]",
                          false, fout);

        MakeOverlayCanvas(A.h_snr[k], B.h_snr[k], A.st_snr[k], B.st_snr[k],
                          Form("c_snr_ch%d", k+1),
                          Form("SNR %s (= ampiezza / RMS_{baseline})", chn_label[k]),
                          "SNR",
                          false, fout);

        MakeOverlayCanvas(A.h_rt[k], B.h_rt[k], A.st_rt[k], B.st_rt[k],
                          Form("c_rt_ch%d", k+1),
                          Form("Rise time 20-80%% %s", chn_label[k]),
                          "Rise time 20-80% [ns]",
                          false, fout);

        MakeOverlayCanvas(A.h_sr[k], B.h_sr[k], A.st_sr[k], B.st_sr[k],
                          Form("c_sr_ch%d", k+1),
                          Form("|Slew rate| 20-80%% %s (fronte di discesa)",
                               chn_label[k]),
                          "|Slew rate| [mV/ns]",
                          false, fout);

        MakeOverlayCanvas(A.h_tcfd[k], B.h_tcfd[k], A.st_tcfd[k], B.st_tcfd[k],
                          Form("c_tcfd_ch%d", k+1),
                          Form("Tempo CFD 20%% %s (scala DRS)", chn_label[k]),
                          "t_{CFD} [ns]",
                          true, fout);
    }

    MakeOverlayCanvas(A.h_dt12, B.h_dt12, A.st_dt12, B.st_dt12,
                      "c_dt12",
                      "#Delta t_{12} = t_{1} - t_{2}",
                      "#Delta t_{12} [ns]",
                      true, fout);

    MakeOverlayCanvas(A.h_C, B.h_C, A.st_C, B.st_C,
                      "c_t3_avg12",
                      "t_{3} - (t_{1}+t_{2})/2",
                      "t_{3} - (t_{1}+t_{2})/2 [ns]",
                      true, fout);

    // ==================================================================
    //  TABELLA RIASSUNTIVA + TTREE summary
    // ==================================================================
    PrintAndSaveSummary(A, B, fout);

    // ==================================================================
    //  CHIUSURA
    // ==================================================================
    fout->Close();
    std::cout << "[DONE] Studio dell'attenuatore completato." << std::endl;
    std::cout << "       Output: " << outname << std::endl;
    std::cout << "       Per visualizzare: root -l " << outname
              << " -> TBrowser b" << std::endl;
}
