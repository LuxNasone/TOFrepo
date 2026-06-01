// ==========================================================================
//  TOF_Analysis_v9.cpp — Analisi del Tempo di Volo e ricostruzione beta, theta
// ==========================================================================
//
//  SCOPO:
//    A partire dai parametri di calibrazione (retta Dt12 vs x e funzione C(x))
//    e dai dati acquisiti con PMT3 in posizione Guida B (sotto la barra),
//    ricostruire evento per evento:
//      - la posizione di impatto x sulla barra
//      - il tempo di volo (TOF)
//      - la distanza percorsa l
//      - l'angolo di incidenza theta
//      - la velocita' v e beta = v/c
//
//  FISICA:
//    Per ogni evento con tripla coincidenza:
//
//    1) POSIZIONE DI IMPATTO sulla barra:
//       Dalla retta di calibrazione  Dt12 = m*x + q   si inverte:
//         x_imp = (Dt12 - q) / m
//
//    2) TEMPO DI VOLO:
//         T_meas = t3 - (t1 + t2)/2
//         TOF = T_meas - C(x_imp)
//       dove C(x_imp) e' calcolata con il modello quadratico
//       C(x) = C_poly_p0 + C_poly_p1*x + C_poly_p2*x^2 (PRINCIPALE),
//       oppure con il modello lineare o interpolazione (fallback).
//       Una pipeline parallela "Cconst" usa C = <C> costante per stimare
//       la sistematica introdotta dalla scelta del modello C(x).
//
//    3) DISTANZA PERCORSA (Pitagora):
//         l = sqrt[ (x_imp - x_PMT3)^2 + h^2 ]
//
//    4) ANGOLO DI INCIDENZA rispetto alla verticale:
//         theta = arctan[ (x_imp - x_PMT3) / h ]
//
//    5) VELOCITA' e beta:
//         v = l / TOF       [cm/ns]
//         beta = v / c      (c = 29.9792 cm/ns)
//         1/v = TOF / l     [ns/cm]
//
// ==========================================================================
//  VERSIONE 9 — ALLINEAMENTO AL PARADIGMA TOT-POLINOMIALE DELLA v13
// ==========================================================================
//
//  Questa versione e' il contraltare di TOF_Calibration_v13.cpp. Cambiamenti
//  strutturali rispetto alla v8:
//
//    - RIMOSSO il recupero clippati via "slew rate" (RecoverClippedCFD, tutte
//      le costanti SLEW_*, le globali gK_PMT/gQ_PMT, i campi slew_* di
//      ChannelData, il fit lineare del fronte in AnalyzeChannel).
//    - RIMOSSO il modello esponenziale A = q + p0*(exp(p1*TOT)-1): il recupero
//      dei clippati usa ora ESCLUSIVAMENTE il polinomio cubico A(TOT) prodotto
//      dalla calibrazione v13, caricato da fit_params.
//    - CFD_FRACTION portata a 0.15 (era 0.20) per coerenza ESATTA con la
//      calibrazione: la retta Dt12(x) e la funzione C(x) dipendono dalla
//      frazione CFD, quindi calibrazione e analisi devono usare lo stesso valore.
//    - AnalyzeChannel allineato a quello di v13 (clipping = solo V_min sotto
//      CLIP_V_LO; oscillazione = solo Schmitt trigger; TOT multi-soglia).
//
//  PIPELINE DI ANALISI: due pipeline parallele e indipendenti.
//    - hybrid_tot : METODO FISICO ufficiale. CFD classico sui canali non
//                   clippati + ricostruzione TOT polinomiale sui soli clippati.
//                   Popola i branch "principali" del TTree (senza suffisso).
//    - noclip     : DIAGNOSTICO. Usa solo i canali non clippati (nessun
//                   recupero); un evento con un canale clippato viene scartato.
//                   E' il riferimento "pulito" privo di ricostruzione.
//                   Popola i branch con suffisso "_noclip".
//
//  CONTRATTO CON LA CALIBRAZIONE (file ROOT prodotto da TOF_Calibration_v13):
//    - TTree "summary"    -> branch 'x' e 'C_hybrid_mu' per i punti C(x_k)
//    - TTree "fit_params" -> m, q, v_eff, C_poly_p0/p1/p2 (quadratico, principale),
//      C_const_val (costante, confronto), C_lin_p0/p1 (lineare, legacy), e i
//      coefficienti del polinomio cubico A(TOT) per PMT1/2/3.
//    Se i branch del polinomio mancano, LoadCalibrationFromFile() emette un
//    ERRORE: in v13 il polinomio e' l'unico modello, non c'e' fallback.
//
//  UTILIZZO:
//    root -l 'TOF_Analysis_v9.cpp("file1.xml,file2.xml", "out.root", "cal_v13.root")'
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
#include <TH2D.h>
#include <TF1.h>
#include <TFitResult.h>  // necessario per fr->IsValid() con ACLiC
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TAxis.h>
#include <TSystem.h>     // gSystem->mkdir() per creare la cartella di cache
#include <TParameter.h>  // marcatore di versione del formato della cache
#include <TDirectory.h>  // salvataggio/ripristino della directory ROOT corrente
#include <TMultiGraph.h> // grafici comparativi con piu' serie

// --- RooFit: fit non parametrico RooKeysPdf (KDE) per le distribuzioni ---
//  In modalita' interpretata (.L) la libreria va caricata esplicitamente:
//  R__LOAD_LIBRARY viene onorata da Cling al caricamento della macro. Se nel
//  proprio setup non bastasse, eseguire una volta gSystem->Load("libRooFit")
//  prima di .L. RooFit serve solo al metodo RooKeysPdf (vedi FitRooKeys()).
#if defined(__CLING__)
R__LOAD_LIBRARY(libRooFit)
#endif
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "RooCmdArg.h"
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
#include <functional>    // std::function per i selettori dei canvas comparativi
#include <sys/stat.h>    // stat() per leggere il tempo di modifica dei file


// ==========================================================================
//  SEZIONE 1: PARAMETRI CONFIGURABILI
// ==========================================================================

// --- Costanti fisiche ---
const double C_LIGHT = 29.9792458;    // Velocita' della luce [cm/ns]

// --- Parametri hardware DRS4 (identici alla calibrazione v13) ---
const int    MAX_SAMPLES  = 1024;
const int    MAX_CHANNELS = 4;
const int    NBL_SAMPLES  = 50;
// FRAZIONE CFD: DEVE essere identica a quella di TOF_Calibration_v13.cpp
// (0.15). La retta Dt12(x) e la funzione C(x) sono state ricavate con questa
// frazione: usarne una diversa in analisi introdurrebbe un bias sistematico.
const double CFD_FRACTION = 0.15;
const double NOISE_THRESH = 5.0;

// CLIP_V_LO: se V_min scende sotto questa soglia [mV] il segnale e' saturato.
const double CLIP_V_LO = -499.0;

// ============================================================
//  RILEVAMENTO OSCILLAZIONE: discriminatore a isteresi (Schmitt trigger)
// ============================================================
//  Unico criterio di oscillazione (allineato a v13): si contano i crossing
//  in salita di OSC_V_HIGH dopo il picco negativo, con isteresi su OSC_V_LOW.
//  >= OSC_NCROSS_MAX crossing -> oscillazione patologica -> evento scartato.
const double OSC_V_HIGH     = 50.0;   // [mV] soglia alta crossing
const double OSC_V_LOW      = 25.0;   // [mV] soglia bassa isteresi
const int    OSC_NCROSS_MAX = 3;      // crossing massimi ammessi

// X_CUT_LO, X_CUT_HI: limiti su x_imp ricostruito. La barra BC408 e' lunga
// 280 cm (da -140 a +140 nel sistema centrato); un x_imp fuori da questo
// intervallo (con un piccolo margine per la risoluzione finita) indica che
// il modello Dt12 = m*x + q non e' applicabile a quell'evento.
const double X_CUT_LO = -145.0;     // [cm]
const double X_CUT_HI =  145.0;     // [cm]

// --- Geometria del setup TOF (Guida B) ---
// h: distanza VERTICALE tra il centro della barra e il centro di PMT3 [cm].
double PAR_H = 101.0;
// x_PMT3: posizione orizzontale del centro di PMT3 rispetto al centro barra [cm].
double PAR_X_PMT3 = 0.0;

// --- Parametri della retta di calibrazione ---
// Dt12 = m*x + q  ->  x = (Dt12 - q) / m. Sovrascritti da LoadCalibrationFromFile().
double CAL_M = 0.132;     // pendenza [ns/cm]
double CAL_Q = -1.805;    // intercetta [ns]

// ============================================================
//  PARAMETRI TIME-OVER-THRESHOLD (TOT) — recupero clippati
// ============================================================
//  Identici alla calibrazione v13. Il segnale viene caratterizzato dalla
//  durata sopra soglia (TOT), da cui si stima l'ampiezza con il polinomio
//  cubico A(TOT) caricato da fit_params.
const int    NTOT_THR        = 5;
const double TOT_THR_MV[NTOT_THR] = { 50.0, 100.0, 150.0, 200.0, 300.0 };

// Indice (0-based) della soglia di riferimento per il calcolo del TOT usato
// nel recupero. Default 0 (= 50 mV); puo' essere AGGIORNATO dal valore letto
// nel file di calibrazione (branch TOT_calib_idx di fit_params) in modo da
// restare automaticamente coerente con la calibrazione.
int          TOT_CALIB_THR_IDX = 0;

// Margine sopra |CLIP_V_LO| per non calcolare il TOT a soglie troppo vicine
// al plateau saturato (il fall crossing sarebbe distorto).
const double TOT_CLIP_MARGIN  = 5.0;   // [mV]

// Numero minimo di campioni tra rise e fall per accettare il TOT.
const int    TOT_NMIN_SAMPLES = 3;

// Cap di estrapolazione: A_rec accettata solo se <= TOT_AREC_MAX_FACTOR * Amax.
const double TOT_AREC_MAX_FACTOR = 3.0;

// Grado del polinomio A(TOT) e dimensione degli array statici.
// TOT_POLY_DEGREE puo' essere aggiornato dal valore letto in fit_params.
int          TOT_POLY_DEGREE  = 3;
const int    TOT_POLY_MAX_DEG = 5;

// ============================================================
//  GLOBALI DEL MODELLO A(TOT) — caricate da LoadCalibrationFromFile()
// ============================================================
//  A(TOT) = sum_k gTOT_poly_PMT[ch][k] * TOT^k  (k = 0..TOT_POLY_DEGREE).
//  gTOT_poly_loaded[ch] = true se la calibrazione e' disponibile per quel
//  canale. Indici 0,1,2 = PMT1,2,3 (PMT3 incluso: entra in C = t3-(t1+t2)/2).
static bool   gTOT_poly_loaded    [MAX_CHANNELS]                       = {false};
static double gTOT_poly_PMT       [MAX_CHANNELS][TOT_POLY_MAX_DEG + 1]  = {{0.0}};
static double gTOT_poly_err_PMT   [MAX_CHANNELS][TOT_POLY_MAX_DEG + 1]  = {{0.0}};
static double gTOT_poly_chi2ndf_PMT[MAX_CHANNELS]                       = {0.0};
static double gTOT_Amax_PMT       [MAX_CHANNELS]                        = {0.0};
static int    gTOT_poly_degree    = 3;   // grado effettivo letto dal file

// ============================================================
//  VETO OFFLINE SUL PMT4 (CH4 della DRS)
// ============================================================
//  Il PMT4 e' il PMT del veto, sotto le lastre di piombo. Logica hardware:
//    START = (1 & 2 & 3) & NOT(4)
//  Se il PMT4 NON ha visto il muone, il muone si presume fermato nel piombo
//  (candidato al decadimento). Il discriminatore NIM ha pero' una soglia fissa
//  e una larghezza finita: alcuni eventi possono PASSARE il veto hardware pur
//  avendo il muone transitato. Avendo la waveform completa del PMT4 sul CH4
//  della DRS4, applichiamo un VETO OFFLINE con soglia piu' stringente e
//  finestra temporale stretta in coincidenza col muone.
const int    VETO_CH_INDEX   = 3;     // indice del canale PMT4
double VETO_N_SIGMA    = 5.0;         // soglia in sigma del rumore di baseline
double VETO_AMP_FLOOR  = 5.0;         // [mV] soglia assoluta minima
double VETO_T_WIN_LO   = -20.0;       // [ns] finestra prima di t_mu
double VETO_T_WIN_HI   =  20.0;       // [ns] finestra dopo  t_mu
int    VETO_NSAMP_MIN  = 3;           // campioni consecutivi sotto soglia minimi
const double VETO_SIGMA_FLOOR = 0.3;  // [mV] floor su sigma_baseline

// ============================================================
//  FUNZIONE C(x): modello lineare + interpolazione di backup
// ============================================================
//  Punto di calibrazione {x_k, C_k}: usato come backup se il modello lineare
//  non e' disponibile. La via principale e' il modello lineare C(x)=C0+C1*x.
struct CPoint {
    double x;    // posizione di calibrazione [cm]
    double C;    // valore di C = t3 - (t1+t2)/2 misurato [ns]
};

// Modello lineare C(x) = gC_lin_p0 + gC_lin_p1*x, caricato da fit_params.
// Se gC_lin_loaded e' true, GetC() usa questo modello.
static double gC_lin_p0     = 0.0;
static double gC_lin_p0_err = 0.0;
static double gC_lin_p1     = 0.0;
static double gC_lin_p1_err = 0.0;
static bool   gC_lin_loaded = false;

// Modello quadratico C(x) = gC_poly_p0 + gC_poly_p1*x + gC_poly_p2*x^2.
// MODELLO PRINCIPALE per C(x): cattura la curvatura parabolica dovuta al
// time-walk di PMT3. Priorita': se caricato, GetC() usa questo modello;
// altrimenti cade sul lineare; altrimenti sull'interpolazione di backup.
static double gC_poly_p0     = 0.0;
static double gC_poly_p0_err = 0.0;
static double gC_poly_p1     = 0.0;
static double gC_poly_p1_err = 0.0;
static double gC_poly_p2     = 0.0;
static double gC_poly_p2_err = 0.0;
static bool   gC_poly_loaded = false;

// Modello costante C(x) = gC_const: media pesata dei C(x_k).
// Usato dalla pipeline di confronto per stimare la sistematica del modello C(x).
static double gC_const_val     = 0.0;
static double gC_const_err     = 0.0;
static bool   gC_const_loaded  = false;

// Parametri della calibrazione corretta per rise time (Pietro&Rick), se
// presenti nel file prodotto da TOF_Calibration_v14.cpp. La pipeline rise_corr
// li usa per ricostruire x_imp e C(x) con la stessa convenzione della calibrazione.
static double CAL_M_CORR = 0.0;
static double CAL_Q_CORR = 0.0;
static bool   gRiseCorrCal_loaded = false;

static double gC_corr_poly_p0     = 0.0;
static double gC_corr_poly_p1     = 0.0;
static double gC_corr_poly_p2     = 0.0;
static bool   gC_corr_poly_loaded = false;

// Parametri della calibrazione RISTRETTA al range lineare (corrected), usati
// dalla pipeline Lin_Range. Caricati da LoadCalibrationFromFile() dai branch
// m_corr_lin/q_corr_lin e C_const_corr_lin (modello 2b: C costante ristretto).
// Sono distinti dai corrected globali: la retta e' fittata solo su x in
// [-84,+70] cm e C e' una costante stimata nello stesso range.
static double CAL_M_LIN = 0.0;          // pendenza retta ristretta [ns/cm]
static double CAL_Q_LIN = 0.0;          // intercetta retta ristretta [ns]
static bool   gLinRangeCal_loaded = false;

// C costante ristretto (pipeline 2b, quella scelta): C(x) = gC_const_lin_val
// per ogni x nel range lineare.
static double gC_const_lin_val     = 0.0;
static double gC_const_lin_err     = 0.0;
static bool   gC_const_lin_loaded  = false;

// C lineare ristretto (pipeline 2a, gia' caricato ma non usato finche' non si
// attiva la 2a): C(x) = gC_lin_corr_p0 + gC_lin_corr_p1 * x nel range lineare.
static double gC_lin_corr_p0       = 0.0;
static double gC_lin_corr_p1       = 0.0;
static bool   gC_lin_corr_loaded   = false;

// Punti C(x_k) di backup (interpolazione lineare). Riempiti da
// LoadCalibrationFromFile() leggendo 'x' e 'C_hybrid_mu' dal TTree summary.
std::vector<CPoint> CAL_C_POINTS = {
    {-130.0, -13.10}, {-112.0, -13.10}, { -84.0, -13.10},
    { -56.0, -13.10}, { -28.0, -13.10}, {   0.0, -13.10},
    {  28.0, -13.10}, {  56.0, -13.10}, {  84.0, -13.10},
    { 112.0, -13.10}, { 130.0, -13.10}
};

// --- Parametri di binning istogrammi ---
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

// Finestra di posizione in cui il rise time dei due PMT della barra e' ben
// approssimato da una relazione lineare. Usata dalla terza pipeline richiesta:
// tempi corretti per rise time + restrizione al range lineare.
const double RISE_LINEAR_X_LO = -84.0;  // [cm]
const double RISE_LINEAR_X_HI =  70.0;  // [cm]

// ============================================================
//  TAGLIO OUTLIER SUL RISE TIME
// ============================================================
//  Lo scatterplot rise time T90-T10 vs x mostra, su PMT1/2/3, eventi sparsi a
//  rise time anomalo (decine di ns, fino a ~140 ns) ben staccati dalla banda
//  fisica (~pochi ns). Sono forme d'onda rumorose/oscillanti che sfuggono allo
//  Schmitt trigger: il loro fronte e' mal definito, quindi anche il CFD lo e'.
//  Introdurrebbero code spurie in TOF/beta. Li rimuoviamo con una soglia rigida
//  applicata a TUTTE le pipeline a monte della ricostruzione.
//
//  CONVENZIONE: si scarta l'evento se ALMENO UN canale ha un rise time VALIDO
//  (rise_time_ok == true) e maggiore di RISE_TIME_MAX. Un canale con rise time
//  NON valido (T90 non trovato, p.es. su clippati) non provoca lo scarto: il
//  rise time semplicemente non e' misurabile, non e' "anomalo".
const double RISE_TIME_MAX = 10.0;   // [ns]

// ============================================================
//  PROFILO RISE TIME MEDIO vs POSIZIONE
// ============================================================
//  Per studiare quantitativamente la dipendenza del rise time dalla posizione
//  d'impatto, si divide la barra in bin equispaziati e si calcola il rise time
//  MEDIO per bin (con errore standard della media). Range e passo come da
//  specifica: [-130,+130] cm a passi di 10 cm -> 26 bin.
const double RISEPROF_X_LO   = -130.0;  // [cm] estremo sinistro del primo bin
const double RISEPROF_X_HI   =  130.0;  // [cm] estremo destro dell'ultimo bin
const double RISEPROF_BIN_W  =   10.0;  // [cm] larghezza di un bin
// Numero minimo di eventi in un bin per disegnarne media ed errore: sotto
// questa soglia la media e l'errore standard non sono statisticamente
// significativi e il punto verrebbe dominato dal rumore.
const int    RISEPROF_N_MIN  =    5;

// ============================================================
//  PARAMETRI DELLA CACHE DI PARSING (XML -> ROOT)
// ============================================================
//  Il parsing testuale degli XML DRS4 e' la fase piu' lenta della macro.
//  Per ogni file XML viene creato un ROOT di cache nella sottocartella
//  PARSE_CACHE_SUBDIR, accanto al dataset originale. La cache viene riusata
//  solo se esiste ed e' almeno recente quanto l'XML sorgente.
//
//  Scelta intenzionale: nella cache salvo SOLO le waveform grezze e gli header.
//  Le grandezze derivate da AnalyzeChannel() vengono ricalcolate a ogni lettura:
//  cosi' cambi futuri a CFD_FRACTION, TOT o rise-time non richiedono di
//  cancellare manualmente tutte le cache gia' generate.
const bool  USE_PARSE_CACHE    = true;            // false = disattiva la cache
const bool  FORCE_REPARSE      = false;           // true = ignora la cache
const char* PARSE_CACHE_SUBDIR = "parsed_cache";  // sottocartella cache
const int   PARSE_CACHE_FORMAT = 2;               // aumenta per invalidare
//  v2: aggiunto il campo 'timestamp' (<Time> XML) alla cache. Necessario per
//  Lead_analysis(), che usa il timestamp DRS per il matching col FIFO. Le cache
//  v1 (prive di timestamp) vengono riconosciute incompatibili e rigenerate.


// ==========================================================================
//  SEZIONE 1b: PARAMETRI E STRUTTURE DEL FILTRO FIFO (DE10-Nano)
// ==========================================================================
//
//  CONTESTO FISICO
//  ---------------
//  Nel run con i PIOMBI la DRS4 viene triggerata dal segnale di START
//      START = 1 & 2 & 3 & NOT(4)
//  cioe' un muone che attraversa la barra (PMT1,2,3) ma NON il PMT4 sotto il
//  piombo: candidato a essersi fermato nel piombo. La sola tripla coincidenza
//  pero' non garantisce che il muone si sia FERMATO e DECADUTO; per quello
//  serve vedere l'elettrone di decadimento (mu -> e + 2nu) che riattraversa il
//  PMT4 entro la finestra caratteristica del decadimento (tau_mu ~ 2.2 us).
//
//  In parallelo alla DRS, la DE10-Nano (FPGA) registra una FIFO di eventi su
//  due canali logici:
//      CH0 (codifica col1 = 2^0 = 1) = START  : la STESSA coincidenza
//                                               1 & 2 & 3 & NOT(4) del trigger DRS
//      CH1 (codifica col1 = 2^1 = 2) = STOP   : 4 & GATE, dove GATE e' una
//                                               finestra (Dual-timer) aperta sul
//                                               fronte di discesa dello START e
//                                               lunga ~15 us. Lo STOP e' quindi
//                                               l'elettrone di decadimento visto
//                                               dal PMT4 entro la finestra.
//
//  Una COPPIA DI DECADIMENTO VALIDA e' uno START seguito da uno STOP entro la
//  GATE. Questa selezione e' INTERNA al FIFO (non coinvolge la DRS). Solo dopo
//  averla costruita confrontiamo il timestamp dello START di ogni coppia valida
//  con i timestamp di acquisizione della DRS: se coincidono entro una finestra
//  (dominata dalla risoluzione al ms del timestamp DRS) l'evento DRS e' un vero
//  candidato a decadimento nel piombo e va analizzato.
//
//  NOTA SUL MODELLO DI MATCHING (correzione rispetto al fork rudimentale):
//  lo START del FIFO e il trigger della DRS sono LO STESSO impulso fisico
//  arrivato a due schede diverse. L'ancora del match e' percio' lo START della
//  coppia valida, NON un generico "canale 1 vicino seguito da canale 2".
//
//  FORMATO DEL FILE FIFO (verificato sui dati reali FIFOread_*.txt)
//  ----------------------------------------------------------------
//  Ogni riga ha due colonne intere separate da TAB/spazi:
//      col1 = tipo evento : 1 = START (CH0), 2 = STOP (CH1)
//      col2 = timestamp   : 0 .. 2^30-1  (cicli di clock da 5 ns nel buffer)
//  Le righe di RESET hanno:
//      col1 = 2^31              (sentinella fissa)
//      col2 = 2^31 + k          (k = 0 al 1o reset, 1 al 2o, ...) contatore reset
//  Il contatore di clock si azzera ad ogni reset, cioe' ogni 2^30 cicli =
//  2^30 * 5 ns = 5.36870912 s (durata di un "buffer").
//
//  Il file inizia tipicamente a META' del primo buffer (buffer 0): le righe
//  prima del PRIMO reset appartengono a un buffer parziale. Per la ricerca
//  delle coppie START-STOP scartiamo questo buffer parziale (le coppie a cavallo
//  dell'inizio file sarebbero ambigue), MA NON perdiamo il conteggio del tempo:
//  nell'interpretazione A l'orario del nome file e' l'istante in cui il contatore
//  FPGA valeva 0 (inizio buffer 0), quindi il buffer 0 dura comunque 5.369 s e il
//  primo reset cade a orario_nome_file + 5.369 s. Cosi' i timestamp assoluti
//  restano ancorati all'orologio reale.
//
//  NOME FILE FIFO
//  --------------
//  "FIFOread_YYYYMMDD-HHMMSS.txt": YYYYMMDD = data, HHMMSS = ora di avvio
//  acquisizione FPGA (istante in cui il clock vale 0). L'orologio DRS e'
//  avanti di FIFO_TZ_OFFSET_S secondi (default 2 h) rispetto a questo.

// --- Periodo del clock FPGA e durata del buffer ---
// Un ciclo di clock dura 5 ns; un buffer sono 2^30 cicli.
static const double   FIFO_CLK_PERIOD_S   = 5.0e-9;          // [s] durata di 1 ciclo
static const long long FIFO_BUF_SIZE_CLK  = (1LL << 30);     // 2^30 cicli per buffer
static const double   FIFO_BUF_DURATION_S = 5.36870912;      // [s] = 2^30 * 5 ns

// --- Sentinella di reset: col1 di una riga di reset ---
static const long long FIFO_RESET_SENTINEL = (1LL << 31);    // 2^31 = 2147483648

// --- Codifica dei canali logici nella col1 ---
static const int FIFO_CH_START = 1;   // CH0 = 2^0 -> START (1 & 2 & 3 & !4)
static const int FIFO_CH_STOP  = 2;   // CH1 = 2^1 -> STOP  (4 & GATE)

// --- Finestra della coppia di decadimento START -> STOP [cicli di clock] ---
// Deve riprodurre la GATE hardware (~15 us) con un piccolo dead-time minimo per
// escludere il rimbalzo immediato dello START stesso. 15 us / 5 ns = 3000 cicli.
// I valori sono modificabili per studiare la sensibilita' della selezione.
//   min: 20 cicli  = 100 ns  (oltre il dead-time / ringing immediato)
//   max: 3000 cicli = 15 us  (lunghezza della GATE)
static long long FIFO_DECAY_DT_MIN = 20LL;     // [clock] 100 ns
static long long FIFO_DECAY_DT_MAX = 3000LL;   // [clock] 15 us

// --- Finestra di matching DRS <-> coppia FIFO [secondi] ---
// Dominata dalla risoluzione del timestamp DRS, che e' al MILLISECONDO
// (campo <Time> = "YYYY/MM/DD HH:MM:SS.sss"). Una finestra di +-100 ms e' un
// punto di partenza prudente; va STRETTA progressivamente (50 ms, 20 ms, ...)
// per trovare il valore minimo che non perde match veri (vedi studio di tuning).
// E' un parametro a runtime: lo si passa a Lead_analysis() o lo si imposta con
// SetFifoMatchWindow().
static double FIFO_MATCH_WINDOW_S = 10;     // [s] +-100 ms (default)

// --- Offset di fuso orario DRS - FPGA [secondi] ---
// T_fpga = T_drs - FIFO_TZ_OFFSET_S. Default 2 h (CEST). Parametro a runtime.
static int FIFO_TZ_OFFSET_S = 2 * 3600;        // [s] 7200

// --- Offset fine di origine clock FPGA [secondi] ---
// Correzione residua se l'orario del nome file non coincide ESATTAMENTE con
// ts=0 del buffer 0 (p.es. latenza di apertura file). T_fpga_riga =
// orario_nome_file + FIFO_T0_OFFSET_S + ts_abs_clock * 5 ns. Si determina
// empiricamente dal primo match riuscito; default 0 (interpretazione A pura).
static double FIFO_T0_OFFSET_S = 0.0;          // [s]


// ==========================================================================
//  STRUTTURE DATI DEL FILTRO FIFO
// ==========================================================================

/// Una coppia di decadimento valida ricostruita dal FIFO: START seguito da uno
/// STOP entro la finestra di GATE. E' l'unita' su cui si fa il matching col DRS.
struct FifoDecayPair {
    long long ts_start_clk;   // timestamp assoluto dello START [cicli di clock]
    long long ts_stop_clk;    // timestamp assoluto dello STOP  [cicli di clock]
    long long dt_clk;         // ts_stop - ts_start [cicli] (= tempo di decadimento)
    double    t_start_abs_s;  // istante assoluto dello START [s dall'inizio giorno, ref. FPGA]
    bool      used;           // true se gia' associata a un evento DRS (consumo 1-a-1)
};

/// Contenitore del contenuto di UN file FIFO, dopo il caricamento e la
/// costruzione delle coppie di decadimento. Una di queste strutture per ogni
/// .txt elaborato (accoppiato posizionalmente al rispettivo .xml).
struct FifoData {
    // Orario di avvio acquisizione FPGA, ricavato dal nome file:
    //   day_index   = giorni dall'epoch (per confronti multi-giorno robusti)
    //   start_s_day = secondi dall'inizio del giorno di avvio (HH*3600+MM*60+SS)
    double    start_s_day;        // [s] HH:MM:SS del nome file
    long long start_day_index;    // giorno assoluto (vedi DaysFromCivil) dell'avvio

    int       n_resets;           // numero di reset (= buffer completati) trovati
    long long n_lines;            // righe totali lette dal file

    // Coppie di decadimento valide, ordinate per ts_start_clk crescente.
    std::vector<FifoDecayPair> pairs;
};

// ==========================================================================
//  SEZIONE 2: STRUTTURE DATI (identiche a TOF_Calibration_v13.cpp)
// ==========================================================================

/// Dati di un singolo canale di un singolo evento.
struct ChannelData {
    int    nsamples;
    float  time[MAX_SAMPLES];
    float  voltage[MAX_SAMPLES];

    // ---- Grandezze estratte da AnalyzeChannel() ----
    double baseline, baseline_rms, amplitude, v_min, v_max, t_min, t_cfd;
    bool   has_pulse, cfd_ok;
    bool   is_clipped;      // true se V_min < CLIP_V_LO (saturazione DRS4)
    bool   is_oscillating;  // true se n_pos_crossings >= OSC_NCROSS_MAX
    int    n_pos_crossings; // crossing Schmitt in salita dopo il picco

    // ---- Time-Over-Threshold multi-soglia ----
    double tot_q   [NTOT_THR];   // durata sopra TOT_THR_MV[k] [ns]
    double tot_rise[NTOT_THR];   // istante del crossing di salita [ns]
    double tot_fall[NTOT_THR];   // istante del crossing di discesa [ns]
    bool   tot_ok  [NTOT_THR];   // true se il TOT alla soglia k e' valido

    // ---- Recupero clippati via TOT polinomiale (metodo hybrid_tot) ----
    // Popolati SOLO sui canali clippati da RecoverClippedCFD_TOT().
    double amplitude_tot;     // ampiezza ricostruita A(TOT) [mV]
    double t_cfd_rec_tot;     // tempo CFD ricalcolato sulla soglia f*A(TOT) [ns]
    bool   cfd_recovered_tot; // true se il recupero del clippato e' riuscito

    // ---- Rise time e crossing 10%/30% per correzione Pietro&Rick ----
    // rise_time = T90 - T10 sul fronte di salita del segnale negativo.
    // t_10 e t_30 servono per stimare dT = t_30 - t_10, cioe' il rise time
    // effettivo da sottrarre al CFD nella pipeline rise_corr.
    double rise_time;     // T90 - T10 [ns] (-999 se non valido)
    bool   rise_time_ok;  // true se T10 e T90 sono stati trovati
    double t_10;          // crossing al 10% dell'ampiezza [ns]
    double t_30;          // crossing al 30% dell'ampiezza [ns]
    bool   t10_ok;        // true se t_10 e' valido
    bool   t30_ok;        // true se t_30 e' valido
};

/// Dati completi di un evento DRS4.
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

/// VetoResult: risultati dell'analisi di veto sul CH4 (PMT4).
struct VetoResult {
    bool   has_signal;       // true se rilevato segnale -> evento vetato
    double amplitude;        // (baseline - vmin_in_window) [mV], >= 0
    double t_peak;           // tempo del minimo nella finestra [ns]
    double baseline;         // baseline media stimata su NBL_SAMPLES [mV]
    double baseline_rms;     // rms baseline (dopo floor) [mV]
    double significance;     // amplitude / baseline_rms  (~ "N sigma")
    double v_thr_used;       // soglia in tensione effettivamente applicata [mV]
    int    n_samples_below;  // campioni consecutivi sotto soglia attorno al minimo
    bool   window_ok;        // true se la finestra temporale era esplorabile
};


// ==========================================================================
//  SEZIONE 3: ANALISI DELLA FORMA D'ONDA (identica a TOF_Calibration_v13)
// ==========================================================================

/// AnalyzeChannel(): estrae da un canale baseline, ampiezza, tempo CFD, i flag
/// di qualita' (is_clipped, is_oscillating) e il Time-Over-Threshold multi-soglia.
/// Implementazione IDENTICA a quella della calibrazione v13: questo garantisce
/// che il significato operativo di t_cfd, is_clipped, tot_q sia lo stesso in
/// calibrazione e in analisi.
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

    for (int kth = 0; kth < NTOT_THR; kth++) {
        cd.tot_q   [kth] = 0.0;
        cd.tot_rise[kth] = 0.0;
        cd.tot_fall[kth] = 0.0;
        cd.tot_ok  [kth] = false;
    }
    cd.amplitude_tot     = 0.0;
    cd.t_cfd_rec_tot     = -999.0;
    cd.cfd_recovered_tot = false;
    cd.rise_time         = -999.0;
    cd.rise_time_ok      = false;
    cd.t_10              = -999.0;
    cd.t_30              = -999.0;
    cd.t10_ok            = false;
    cd.t30_ok            = false;

    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return;

    float *t = cd.time;
    float *v = cd.voltage;

    // ---- PASSO 1: BASELINE (media + RMS dei primi NBL_SAMPLES campioni) ----
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < NBL_SAMPLES; i++) {
        sum  += v[i];
        sum2 += (double)v[i] * v[i];
    }
    double bl     = sum / NBL_SAMPLES;
    double bl_rms = sqrt(fabs(sum2 / NBL_SAMPLES - bl * bl));
    cd.baseline     = bl;
    cd.baseline_rms = bl_rms;

    // ---- PASSO 2: MINIMO E MASSIMO GLOBALI ----
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
    if (vmin < CLIP_V_LO) cd.is_clipped = true;

    // ---- PASSO 4: AMPIEZZA E SOGLIA DI RUMORE ----
    double amp = bl - vmin;
    cd.amplitude = amp;
    if (amp < NOISE_THRESH) return;
    cd.has_pulse = true;

    // ---- PASSO 5: RILEVAMENTO OSCILLAZIONE (Schmitt trigger con isteresi) ----
    {
        int  n_cross = 0;
        bool above   = false;
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

    // ---- PASSO 6b: RISE TIME T90-T10 + CROSSING T10/T30 ----
    // I tre crossing vengono calcolati sul fronte di salita del segnale
    // negativo. T90-T10 e' la diagnostica del rise time; t_10 e t_30 sono
    // salvati separatamente per la correzione Pietro&Rick:
    //   dT_n = t_30,n - t_10,n
    //   t_corr,n = t_cfd,n - dT_n
    // La correzione viene poi applicata solo a eventi non clippati, perche'
    // sui clippati l'ampiezza misurata e' sottostimata e le frazioni 10/30/90%
    // non rappresentano piu' le frazioni della vera ampiezza.
    {
        double v_10 = bl - 0.10 * amp;
        double v_30 = bl - 0.30 * amp;
        double v_90 = bl - 0.90 * amp;

        auto findLeadingCrossing = [&](double v_level, double &t_cross) -> bool {
            t_cross = -999.0;
            for (int i = NBL_SAMPLES; i < imin; i++) {
                if (v[i] > v_level && v[i + 1] <= v_level) {
                    double dv = (double)v[i + 1] - v[i];
                    if (fabs(dv) > 1e-6) {
                        t_cross = t[i] + (v_level - v[i]) / dv * (t[i + 1] - t[i]);
                        return true;
                    }
                    break;
                }
            }
            return false;
        };

        double t10 = -999.0, t30 = -999.0, t90 = -999.0;
        bool ok10 = findLeadingCrossing(v_10, t10);
        bool ok30 = findLeadingCrossing(v_30, t30);
        bool ok90 = findLeadingCrossing(v_90, t90);

        if (ok10) {
            cd.t_10   = t10;
            cd.t10_ok = true;
        }
        if (ok10 && ok30 && t30 > t10) {
            cd.t_30   = t30;
            cd.t30_ok = true;
        }
        if (ok10 && ok90 && t90 > t10) {
            cd.rise_time    = t90 - t10;
            cd.rise_time_ok = true;
        }
    }

    // ---- PASSO 7: TIME-OVER-THRESHOLD MULTI-SOGLIA ----
    {
        const double q_max_safe = fabs(CLIP_V_LO) - TOT_CLIP_MARGIN;
        double sample_period = (ns > 1)
                             ? ((double)t[ns - 1] - (double)t[0]) / (ns - 1)
                             : 0.2;   // fallback 0.2 ns = 5 GS/s
        double dt_min_tot = TOT_NMIN_SAMPLES * sample_period;

        for (int kth = 0; kth < NTOT_THR; kth++) {
            double q = TOT_THR_MV[kth];
            if (q > q_max_safe) continue;
            if (amp < q)        continue;

            double v_thr_tot = bl - q;

            // crossing sul fronte di salita
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

            // crossing sul fronte di discesa
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
            if (tot_val < dt_min_tot) continue;

            cd.tot_rise[kth] = t_rise;
            cd.tot_fall[kth] = t_fall;
            cd.tot_q   [kth] = tot_val;
            cd.tot_ok  [kth] = true;
        }
    }
}


// ==========================================================================
//  SEZIONE 3b: ANALISI DEL CANALE DI VETO (CH4 = PMT4)
// ==========================================================================
//
//  AnalyzeVeto(): esamina la waveform del PMT4 e decide se e' presente un
//  segnale vero (non rumore) IN COINCIDENZA con il passaggio del muone.
//  Funzione invariata rispetto alla v8.
//
//  ALGORITMO (4 passi):
//   1) BASELINE E RUMORE sui primi NBL_SAMPLES campioni (regione pre-trigger).
//   2) SOGLIA: A_thr = max(VETO_N_SIGMA * bl_rms_eff, VETO_AMP_FLOOR);
//      v_thr = bl - A_thr (negativa, per impulsi negativi).
//   3) RICERCA DEL MINIMO nella finestra [t_mu + WIN_LO, t_mu + WIN_HI].
//   4) VALIDAZIONE: il candidato e' "segnale" se vmin_in_win < v_thr e se
//      almeno VETO_NSAMP_MIN campioni consecutivi attorno al minimo sono
//      sotto v_thr (durata fisica del picco).
//
//  PARAMETRI: cd = ChannelData del PMT4; t_mu = tempo stimato di passaggio
//  del muone [ns] (tipicamente (t_cfd_PMT1 + t_cfd_PMT2)/2).
VetoResult AnalyzeVeto(const ChannelData &cd, double t_mu) {

    VetoResult r;
    r.has_signal      = false;
    r.amplitude       = 0.0;
    r.t_peak          = -999.0;
    r.baseline        = 0.0;
    r.baseline_rms    = 0.0;
    r.significance    = 0.0;
    r.v_thr_used      = 0.0;
    r.n_samples_below = 0;
    r.window_ok       = false;

    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return r;   // troppo corto: non analizzabile

    const float *t = cd.time;
    const float *v = cd.voltage;

    // --- (1) baseline e rms sui primi NBL_SAMPLES ---
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < NBL_SAMPLES; i++) {
        sum  += (double)v[i];
        sum2 += (double)v[i] * (double)v[i];
    }
    double bl     = sum / NBL_SAMPLES;
    double bl_rms = sqrt(fabs(sum2 / NBL_SAMPLES - bl * bl));
    r.baseline     = bl;
    r.baseline_rms = bl_rms;

    double bl_rms_eff = std::max(bl_rms, VETO_SIGMA_FLOOR);

    // --- (2) soglia di ampiezza e di tensione ---
    double A_thr = std::max(VETO_N_SIGMA * bl_rms_eff, VETO_AMP_FLOOR);
    double v_thr = bl - A_thr;
    r.v_thr_used = v_thr;

    // --- (3) ricerca del minimo nella finestra temporale ---
    double t_lo = t_mu + VETO_T_WIN_LO;
    double t_hi = t_mu + VETO_T_WIN_HI;

    double vmin_in_win = bl;
    double tmin_in_win = t_mu;
    int    imin_in_win = -1;
    bool   any_in_win  = false;

    for (int i = 0; i < ns; i++) {
        double ti = (double)t[i];
        if (ti < t_lo || ti > t_hi) continue;
        any_in_win = true;
        if ((double)v[i] < vmin_in_win) {
            vmin_in_win = (double)v[i];
            tmin_in_win = ti;
            imin_in_win = i;
        }
    }
    r.window_ok = any_in_win;
    if (!any_in_win || imin_in_win < 0) return r;

    r.amplitude    = bl - vmin_in_win;
    r.t_peak       = tmin_in_win;
    r.significance = (bl_rms_eff > 1e-9) ? r.amplitude / bl_rms_eff : 0.0;

    // --- (4) validazione del candidato ---
    if (vmin_in_win >= v_thr) return r;   // il minimo e' sopra soglia -> rumore

    // durata: campioni CONSECUTIVI sotto soglia attorno al minimo
    int n_consec = 1;
    for (int i = imin_in_win - 1; i >= 0; i--) {
        if ((double)v[i] < v_thr) n_consec++;
        else break;
    }
    for (int i = imin_in_win + 1; i < ns; i++) {
        if ((double)v[i] < v_thr) n_consec++;
        else break;
    }
    r.n_samples_below = n_consec;

    if (n_consec >= VETO_NSAMP_MIN) r.has_signal = true;
    return r;
}


// ==========================================================================
//  SEZIONE 3c: STRUTTURA PipelineHistos E FUNZIONI DI SERVIZIO
// ==========================================================================
//
//  MOTIVAZIONE DEL REFACTOR
//  ------------------------
//  L'analisi confronta 4 strategie di ricostruzione che condividono ESATTAMENTE
//  lo stesso insieme di osservabili fisiche. Invece di dichiarare, riempire e
//  disegnare a mano un set di istogrammi per ciascuna pipeline (codice
//  duplicato 4 volte, fragile e difficile da mantenere), raccogliamo tutti gli
//  istogrammi di UNA pipeline in una struttura PipelineHistos e operiamo su di
//  essa con funzioni generiche. Aggiungere/rimuovere una pipeline diventa una
//  riga, e i 4 canvas avranno per costruzione gli STESSI grafici.
//
//  Le 4 pipeline comparabili (stesso set di grafici) sono:
//    "hybrid"   -> Hybrid_tot : segnali validi + ricostruzione TOT dei clippati
//    "noclip"   -> No_clip    : solo segnali non clippati (nessun recupero)
//    "risecorr" -> Rise_correct : validi + TOT + correzione rise time
//    "linrange" -> Lin_Range  : come risecorr, ristretta alla zona lineare di x
//
//  La pipeline "Cconst" (sistematica del modello C(x)) NON e' una delle 4
//  comparabili: resta gestita a parte, fuori da questa struttura e dai canvas
//  di confronto.

/// PipelineHistos: contenitore di tutti gli istogrammi di UNA pipeline.
///
///  Le 9 osservabili richieste per ogni pipeline:
///    beta        : distribuzione di beta = v/c                 (h_beta)
///    beta_norm   : stessa di beta, ma normalizzata ad area (PDF) (h_beta_norm)
///    tof         : tempo di volo                                (h_tof)
///    invv        : 1/v                                          (h_invv)
///    theta       : angolo di incidenza                          (h_theta)
///    x_imp       : posizione di impatto                         (h_ximp)
///    path        : distanza percorsa l                          (h_path)
///    tof_vs_x    : TOF vs x_imp (2D)                            (h2_tof_x)
///    beta_vs_x   : beta vs x_imp (2D)                           (h2_beta_x)
///
///  Tutti i puntatori sono inizializzati a nullptr; vengono allocati da
///  BookPipeline(). La struttura NON possiede gli istogrammi nel senso di
///  doverli deallocare: in ROOT, creando un TH1 con un TFile aperto in
///  RECREATE, l'ownership e' del file, che li distrugge alla chiusura. Per
///  questo NON definiamo un distruttore che faccia delete (causerebbe un
///  double-free alla chiusura del TFile).
struct PipelineHistos {
    std::string tag;     // suffisso dei nomi ROOT (es. "hybrid")
    std::string label;   // etichetta leggibile nei titoli (es. "Hybrid_tot")

    TH1D* h_beta      = nullptr;
    TH1D* h_beta_norm = nullptr;
    TH1D* h_tof       = nullptr;
    TH1D* h_invv      = nullptr;
    TH1D* h_theta     = nullptr;
    TH1D* h_ximp      = nullptr;
    TH1D* h_path      = nullptr;
    TH1D* h_Tmeas     = nullptr;   // T_meas: utile come diagnostica, non comparato
    TH2D* h2_tof_x    = nullptr;
    TH2D* h2_beta_x   = nullptr;
};

/// BookPipeline(): alloca tutti gli istogrammi di una pipeline con nomi e
/// titoli derivati dal tag/label. I parametri di binning sono le costanti
/// globali (NBINS_*, *_LO, *_HI), cosi' i 4 set sono identici per costruzione
/// e i confronti sovrapposti sono leciti bin a bin.
///
///   ph    : struttura da popolare (passata per riferimento)
///   tag   : suffisso dei nomi ROOT (univoco per pipeline)
///   label : etichetta leggibile, compare tra parentesi nei titoli
void BookPipeline(PipelineHistos& ph,
                  const std::string& tag,
                  const std::string& label) {
    ph.tag   = tag;
    ph.label = label;
    const char* t = tag.c_str();
    const char* L = label.c_str();

    ph.h_beta = new TH1D(Form("h_beta_%s", t),
        Form("Distribuzione #beta = v/c (%s);#beta;Conteggi", L),
        NBINS_BETA, BETA_LO, BETA_HI);

    ph.h_beta_norm = new TH1D(Form("h_beta_norm_%s", t),
        Form("Distribuzione #beta normalizzata (%s);#beta;PDF", L),
        NBINS_BETA, BETA_LO, BETA_HI);

    ph.h_tof = new TH1D(Form("h_tof_%s", t),
        Form("Tempo di volo (%s);TOF [ns];Conteggi", L),
        NBINS_TOF, TOF_LO, TOF_HI);

    ph.h_invv = new TH1D(Form("h_invv_%s", t),
        Form("Distribuzione 1/v (%s);1/v [ns/cm];Conteggi", L),
        NBINS_INVV, INVV_LO, INVV_HI);

    ph.h_theta = new TH1D(Form("h_theta_%s", t),
        Form("Distribuzione angolare (%s);#theta [rad];Conteggi", L),
        NBINS_THETA, THETA_LO, THETA_HI);

    ph.h_ximp = new TH1D(Form("h_ximp_%s", t),
        Form("Posizione di impatto (%s);x_{imp} [cm];Conteggi", L),
        NBINS_X, X_LO, X_HI);

    ph.h_path = new TH1D(Form("h_path_%s", t),
        Form("Distanza percorsa (%s);l [cm];Conteggi", L),
        100, 90, 250);

    ph.h_Tmeas = new TH1D(Form("h_Tmeas_%s", t),
        Form("T_{meas} = t_{3}-(t_{1}+t_{2})/2 (%s);T_{meas} [ns];Conteggi", L),
        150, -5, 25);

    ph.h2_tof_x = new TH2D(Form("h2_tof_x_%s", t),
        Form("TOF vs posizione (%s);x_{imp} [cm];TOF [ns]", L),
        70, X_LO, X_HI, 75, TOF_LO, TOF_HI);

    ph.h2_beta_x = new TH2D(Form("h2_beta_x_%s", t),
        Form("#beta vs posizione (%s);x_{imp} [cm];#beta", L),
        70, X_LO, X_HI, 100, BETA_LO, BETA_HI);
}

/// FillPipeline(): riempie gli istogrammi di una pipeline da una SINGOLA
/// osservazione gia' ricostruita e validata a monte (good == 1, tof > 0.1).
/// Replica esattamente la logica di riempimento usata finora per hybrid_tot,
/// inclusi i tagli di sicurezza sui range di beta e 1/v.
///
///  La separazione tra "evento buono" (deciso dal chiamante) e "riempimento"
///  (questa funzione) mantiene la logica fisica nel loop principale e qui solo
///  la meccanica di Fill: la funzione NON applica tagli fisici, solo i clamp
///  di range istogramma identici a quelli preesistenti.
void FillPipeline(PipelineHistos& ph,
                  double tof, double beta, double invv,
                  double theta, double x_imp, double path, double T_meas) {
    ph.h_Tmeas->Fill(T_meas);
    ph.h_ximp->Fill(x_imp);
    if (tof > 0.1) {
        ph.h_tof->Fill(tof);
        ph.h_path->Fill(path);
        ph.h_theta->Fill(theta);
        ph.h2_tof_x->Fill(x_imp, tof);
        if (beta > 0.0 && beta < 3.0) {
            ph.h_beta->Fill(beta);
            ph.h_beta_norm->Fill(beta);
            ph.h2_beta_x->Fill(x_imp, beta);
        }
        if (invv > 0.0 && invv < 0.2) ph.h_invv->Fill(invv);
    }
}

/// NormalizeBetaPDF(): trasforma h_beta_norm in una densita' di probabilita'
/// dividendo per l'integrale "con larghezza di bin". Dopo questa operazione
/// l'asse Y e' in unita' di [1/beta] e pipeline con statistica diversa sono
/// direttamente confrontabili. h_beta (per le statistiche) NON viene toccato.
void NormalizeBetaPDF(PipelineHistos& ph) {
    if (ph.h_beta_norm && ph.h_beta_norm->Integral() > 0.0) {
        double nf = ph.h_beta_norm->Integral("width");
        ph.h_beta_norm->Scale(1.0 / nf);
        ph.h_beta_norm->GetYaxis()->SetTitle("PDF [1/unit. #beta]");
        std::cout << "[INFO] h_beta_norm_" << ph.tag
                  << " normalizzata (fattore = " << nf << ")." << std::endl;
    } else {
        std::cout << "[WARNING] h_beta_norm_" << ph.tag
                  << " vuota: normalizzazione saltata." << std::endl;
    }
}

/// WritePipeline(): scrive nel file ROOT corrente tutti gli istogrammi della
/// pipeline. Da chiamare con il TFile di output gia' come directory corrente.
void WritePipeline(PipelineHistos& ph) {
    if (ph.h_beta)      ph.h_beta->Write();
    if (ph.h_beta_norm) ph.h_beta_norm->Write();
    if (ph.h_tof)       ph.h_tof->Write();
    if (ph.h_invv)      ph.h_invv->Write();
    if (ph.h_theta)     ph.h_theta->Write();
    if (ph.h_ximp)      ph.h_ximp->Write();
    if (ph.h_path)      ph.h_path->Write();
    if (ph.h_Tmeas)     ph.h_Tmeas->Write();
    if (ph.h2_tof_x)    ph.h2_tof_x->Write();
    if (ph.h2_beta_x)   ph.h2_beta_x->Write();
}

/// DrawPipelineSummary(): canvas riassuntivo 3x3 di una singola pipeline, con
/// le 9 osservabili richieste. Viene scritto nel file ROOT con un nome che
/// include il tag, cosi' le 4 pipeline producono 4 canvas omologhi.
void DrawPipelineSummary(PipelineHistos& ph) {
    TCanvas* c = new TCanvas(Form("c_summary_%s", ph.tag.c_str()),
        Form("Riassunto %s", ph.label.c_str()), 1500, 1300);
    c->Divide(3, 3);

    int col = kBlue + 1;   // colore neutro per i 1D del riassunto

    auto draw1D = [&](int pad, TH1D* h) {
        c->cd(pad);
        gPad->SetGrid();
        if (!h) return;
        h->SetLineColor(col);
        h->SetLineWidth(2);
        h->Draw("HIST");
    };

    draw1D(1, ph.h_beta);
    draw1D(2, ph.h_beta_norm);
    draw1D(3, ph.h_tof);
    draw1D(4, ph.h_invv);
    draw1D(5, ph.h_theta);
    draw1D(6, ph.h_ximp);
    draw1D(7, ph.h_path);

    c->cd(8); gPad->SetGrid();
    if (ph.h2_tof_x) ph.h2_tof_x->Draw("COLZ");

    c->cd(9); gPad->SetGrid();
    if (ph.h2_beta_x) {
        ph.h2_beta_x->Draw("COLZ");
        TLine* l = new TLine(X_LO, 1.0, X_HI, 1.0);
        l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw("same");
    }

    c->Write(Form("Canvas_Summary_%s", ph.tag.c_str()));
}

/// DrawPipelineComparison(): sovrappone la STESSA osservabile delle 4 pipeline
/// in un unico pad, con legenda. E' il mattone dei 6 canvas comparativi.
///
///   phs[]    : array delle 4 pipeline
///   nphs     : numero di pipeline (4)
///   selector : funzione che, data una PipelineHistos, restituisce il TH1D
///              dell'osservabile da confrontare (es. [](PipelineHistos&p){return p.h_beta;})
///   cname    : nome ROOT del canvas
///   ctitle   : titolo del canvas (anche titolo del primo istogramma disegnato)
///   draw_beta1_line : se true disegna la linea verticale a beta = 1
void DrawPipelineComparison(PipelineHistos phs[], int nphs,
                            std::function<TH1D*(PipelineHistos&)> selector,
                            const char* cname, const char* ctitle,
                            bool draw_beta1_line = false) {
    // Palette e stili distinti per le 4 pipeline.
    const int cols  [4] = { kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1 };
    const int styles[4] = { 1, 2, 1, 2 };

    TCanvas* c = new TCanvas(cname, ctitle, 900, 650);
    c->SetGrid();

    // Massimo globale per non tagliare nessuna curva.
    double ymax = 0.0;
    for (int i = 0; i < nphs; i++) {
        TH1D* h = selector(phs[i]);
        if (h) ymax = std::max(ymax, h->GetMaximum());
    }

    TLegend* leg = new TLegend(0.60, 0.66, 0.89, 0.89);
    leg->SetTextFont(42);
    leg->SetTextSize(0.032);

    bool first = true;
    for (int i = 0; i < nphs; i++) {
        TH1D* h = selector(phs[i]);
        if (!h) continue;
        h->SetLineColor(cols[i % 4]);
        h->SetLineWidth(2);
        h->SetLineStyle(styles[i % 4]);
        if (first) {
            h->SetTitle(ctitle);
            h->SetMaximum((ymax > 0.0) ? 1.18 * ymax : 1.0);
            h->Draw("HIST");
            first = false;
        } else {
            h->Draw("HIST same");
        }
        leg->AddEntry(h, phs[i].label.c_str(), "l");
    }

if (draw_beta1_line && ymax > 0.0) {
        TLine* l = new TLine(1.0, 0.0, 1.0, 1.18 * ymax);
        l->SetLineColor(kBlack); l->SetLineStyle(2); l->SetLineWidth(2);
        l->Draw("same");
    }
    leg->Draw();
    c->Write(cname);
}

/// BuildRiseProfile(): costruisce il profilo del rise time medio vs posizione
/// per un singolo PMT, a partire dalle coppie (x_imp, rise_time) raccolte nel
/// loop eventi (campione Hybrid_tot).
///
///  ALGORITMO
///   1. Si definiscono nbin = (X_HI - X_LO)/BIN_W bin equispaziati in x.
///   2. Per ogni evento, si individua il bin di appartenenza dalla x_imp e si
///      accumulano i momenti del rise time: N, somma (S1) e somma dei quadrati
///      (S2). Questi tre accumulatori bastano per media e varianza (formule a
///      un passaggio, niente secondo loop sui dati).
///   3. Per ogni bin con N >= N_MIN:
///         media   <tau>   = S1 / N
///         var. camp. s^2  = (S2 - N*<tau>^2) / (N - 1)   [stimatore non distorto]
///         err. media      = s / sqrt(N)                  [errore standard]
///      Il punto e' posto al CENTRO del bin; la barra d'errore orizzontale e'
///      la semi-larghezza del bin (BIN_W/2), quella verticale e' l'errore della
///      media. I bin con N < N_MIN vengono saltati (non disegnati).
///
///  PARAMETRI
///   vx, vy   : vettori paralleli (x_imp [cm], rise_time [ns]) del PMT.
///   x_lo,x_hi,bin_w,n_min : binning e taglio di statistica minima.
///   name,title,color : identita' grafica del TGraphErrors restituito.
///
///  RITORNA un TGraphErrors* (vuoto se nessun bin supera N_MIN). L'ownership
///  e' del chiamante (va disegnato e/o scritto su file).
TGraphErrors* BuildRiseProfile(const std::vector<double>& vx,
                               const std::vector<double>& vy,
                               double x_lo, double x_hi,
                               double bin_w, int n_min,
                               const char* name, const char* title,
                               int color) {
    int nbin = (int)std::lround((x_hi - x_lo) / bin_w);
    if (nbin < 1) nbin = 1;

    std::vector<long>   N (nbin, 0);
    std::vector<double> S1(nbin, 0.0);   // somma dei rise time del bin
    std::vector<double> S2(nbin, 0.0);   // somma dei quadrati dei rise time

    // --- Accumulo dei momenti per bin ---
    for (size_t i = 0; i < vx.size() && i < vy.size(); i++) {
        double x = vx[i];
        double y = vy[i];
        if (x < x_lo || x >= x_hi) continue;          // fuori range -> ignorato
        int ib = (int)((x - x_lo) / bin_w);           // indice di bin
        if (ib < 0 || ib >= nbin) continue;           // guardia
        N [ib] += 1;
        S1[ib] += y;
        S2[ib] += y * y;
    }

    // --- Calcolo media/errore e riempimento del grafico ---
    TGraphErrors* g = new TGraphErrors();
    g->SetName(name);
    g->SetTitle(title);
    int ip = 0;
    for (int ib = 0; ib < nbin; ib++) {
        if (N[ib] < n_min) continue;
        double n     = (double)N[ib];
        double mean  = S1[ib] / n;
        double var   = (n > 1.0) ? (S2[ib] - n * mean * mean) / (n - 1.0) : 0.0;
        if (var < 0.0) var = 0.0;                     // sicurezza numerica
        double sdev  = std::sqrt(var);
        double emean = (n > 0.0) ? sdev / std::sqrt(n) : 0.0;

        double x_center = x_lo + (ib + 0.5) * bin_w;
        g->SetPoint(ip, x_center, mean);
        g->SetPointError(ip, bin_w / 2.0, emean);     // ex = semi-bin, ey = err media
        ip++;
    }

g->SetMarkerStyle(20);
    g->SetMarkerSize(0.9);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
    return g;
}


// ==========================================================================
//  SEZIONE 3d: FUNZIONI DI FIT (copiate da TOF_Calibration_v14, autocontenute)
// ==========================================================================
//
//  Per fittare le distribuzioni fisiche (beta, TOF, 1/v) riutilizziamo gli
//  stessi strumenti di fit validati in calibrazione, copiati qui per mantenere
//  la macro autocontenuta:
//    - FitResult         : contenitore standard di un risultato di fit;
//    - EstimateHistFWHM  : stima della FWHM dell'istogramma (init dei fit);
//    - FitGaussianCore   : gaussiana ristretta al core (iterativa);
//    - FitLorentzian     : Lorentziana/Cauchy (diagnostica code, non usata
//                          nel flusso principale ma mantenuta per completezza);
//    - RiseTimeEMG       : funzione EMG (gaussiana + coda esponenziale destra);
//    - FitRiseTimeEMG    : fit EMG con fallback gaussiano.
//  Le firme e la logica sono IDENTICHE a quelle della calibrazione: cosi' i
//  valori sono confrontabili e il metodo e' gia' documentato/validato.

/// FitResult: contenitore di un risultato di fit a singolo picco.
struct FitResult {
    double center;       // centro/posizione caratteristica del picco
    double center_err;   // errore sul centro
    double width;        // larghezza (sigma per la gaussiana)
    double width_err;    // errore sulla larghezza
    double chi2_ndf;     // chi2/ndf (diagnostica; non e' la quantita' minimizzata
                         // nei fit di log-likelihood)
    int    nentries;     // numero di eventi nell'istogramma
    bool   fit_ok;       // true se il fit e' convergito con risultati sensati
};

/// EstimateHistFWHM(): stima della FWHM dell'istogramma dal massimo e dai due
/// attraversamenti della meta'-altezza (interpolati linearmente). Usata solo
/// per inizializzare i fit; non e' una misura finale.
double EstimateHistFWHM(TH1D* h) {
    int    max_bin  = h->GetMaximumBin();
    double max_val  = h->GetBinContent(max_bin);
    double half_max = max_val / 2.0;
    int    nbins    = h->GetNbinsX();

    // Bordo sinistro (primo bin sotto meta'-altezza scendendo dal massimo).
    double x_left = h->GetBinCenter(1);
    for (int i = max_bin; i >= 1; i--) {
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

    // Bordo destro (primo bin sotto meta'-altezza salendo dal massimo).
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

/// FitGaussianCore(): gaussiana ristretta al "core" (default 90% centrale),
/// con seconda iterazione in [mu +/- 2 sigma] per stabilizzare il centro su
/// distribuzioni con code asimmetriche. Restituisce mu come center e sigma
/// come width. Log-likelihood Poissoniana ("L").
FitResult FitGaussianCore(TH1D* h,
                          double core_fraction = 0.90,
                          const char* fit_name = "gaus_fit") {
    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h->GetEntries();

    if (res.nentries < 30) {
        res.center     = h->GetMean();
        res.center_err = h->GetMeanError();
        res.width      = h->GetStdDev();
        res.width_err  = h->GetStdDevError();
        res.chi2_ndf   = -1.0;
        return res;
    }

    double q_lo = 0.5 - core_fraction / 2.0;
    double q_hi = 0.5 + core_fraction / 2.0;
    const int nq = 2;
    double xq[nq] = { q_lo, q_hi };
    double yq[nq] = { 0.0, 0.0 };
    h->GetQuantiles(nq, yq, xq);
    double x_lo = yq[0];
    double x_hi = yq[1];

    if (x_hi - x_lo < 5.0 * h->GetBinWidth(1)) {
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

    int fit_status = h->Fit(fG, "L S R Q 0");

    if (fit_status == 0 && fG->GetParameter(2) > 0) {
        double mu_1    = fG->GetParameter(1);
        double sig_1   = fabs(fG->GetParameter(2));
        double fit2_lo = std::max(mu_1 - 2.0 * sig_1, x_lo);
        double fit2_hi = std::min(mu_1 + 2.0 * sig_1, x_hi);
        if (fit2_hi - fit2_lo > 5.0 * h->GetBinWidth(1)) {
            fG->SetRange(fit2_lo, fit2_hi);
            fG->SetParLimits(1, fit2_lo, fit2_hi);
            int    imax_loc = h->GetMaximumBin();
            double hmax_loc = h->GetBinContent(imax_loc);
            fG->SetParameter(0, hmax_loc);
            fG->SetParameter(1, h->GetBinCenter(imax_loc));
            fG->SetParameter(2, sig_1);
            fit_status = h->Fit(fG, "L S R Q 0");
        }
    }

    res.center     = fG->GetParameter(1);
    res.center_err = fG->GetParError(1);
    res.width      = fabs(fG->GetParameter(2));
    res.width_err  = fG->GetParError(2);
    res.chi2_ndf   = (fG->GetNDF() > 0) ? fG->GetChisquare() / fG->GetNDF() : -1.0;
    res.fit_ok     = (fit_status == 0);
    return res;
}

/// FitLorentzian(): Lorentziana/Cauchy con parametrizzazione ad altezza di
/// picco e log-likelihood Poissoniana. Mantenuta per completezza/diagnostica
/// (code ~1/x^2); non e' nel flusso principale dell'analisi.
FitResult FitLorentzian(TH1D* h, const char* fit_name = "bw_fit") {
    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h->GetEntries();

    if (res.nentries < 30) {
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

    int fit_status = h->Fit(fL, "L S R Q 0");

    res.center     = fL->GetParameter(1);
    res.center_err = fL->GetParError(1);
    res.width      = fL->GetParameter(2);
    res.width_err  = fL->GetParError(2);
    res.chi2_ndf   = (fL->GetNDF() > 0) ? fL->GetChisquare() / fL->GetNDF() : -1.0;
    res.fit_ok     = (fit_status == 0);
    return res;
}

/// RiseTimeEMG(): Exponentially-Modified Gaussian (gaussiana + coda esponenziale
/// DESTRA), con ramo asintotico stabile per x << mu. Parametri:
///   par[0]=A (scala), par[1]=mu, par[2]=sigma, par[3]=tau.
/// Media = mu + tau; la moda (picco) e' in (mu, mu+tau).
double RiseTimeEMG(double* xx, double* par) {
    double x   = xx[0];
    double A   = par[0];
    double mu  = par[1];
    double sig = par[2];
    double tau = par[3];
    if (sig <= 0.0 || tau <= 0.0) return 0.0;

    double r = sig / tau;
    double u = (x - mu) / sig;
    double z = (r - u) / TMath::Sqrt2();

    if (z > 6.0) {
        double g = TMath::Exp(-0.5 * u * u);
        return A * g / (z * TMath::Sqrt(TMath::Pi()));
    }
    double arg_exp = 0.5 * r * r - (x - mu) / tau;
    return A * TMath::Exp(arg_exp) * TMath::Erfc(z);
}

/// FitRiseTimeEMG(): fit EMG (coda destra). Restituisce come center la MODA
/// della curva (picco), come width il sigma gaussiano. Fallback su
/// FitGaussianCore con < 50 eventi o se l'EMG non converge.
///   fit_name   : nome univoco della TF1.
///   line_color : colore curva (>=0 per disegnarla colorata, -1 default).
///   attach     : se false usa "N" (non disegna/non aggancia all'istogramma).
FitResult FitRiseTimeEMG(TH1D* h,
                         const char* fit_name = "emg_fit",
                         int line_color = -1,
                         bool attach = true) {
    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h->GetEntries();

    if (res.nentries < 50) {
        return FitGaussianCore(h, 0.90, fit_name);
    }

    int    imax     = h->GetMaximumBin();
    double peak_hgt = h->GetBinContent(imax);
    double peak_pos = h->GetBinCenter(imax);
    double rms      = h->GetStdDev();
    double xlo      = h->GetXaxis()->GetXmin();
    double xhi      = h->GetXaxis()->GetXmax();
    double bw       = h->GetBinWidth(1);

    double sig0 = std::max(0.4 * rms, 2.0 * bw);
    double tau0 = std::max(0.5 * rms, 2.0 * bw);
    double mu0  = peak_pos - 0.5 * tau0;

    TF1* fEMG = new TF1(fit_name, RiseTimeEMG, xlo, xhi, 4);
    fEMG->SetParNames("A_scale", "mu", "sigma", "tau");
    fEMG->SetParameters(1.0, mu0, sig0, tau0);
    double shape_at_peak = fEMG->Eval(peak_pos);
    double A0 = (shape_at_peak > 1e-12) ? peak_hgt / shape_at_peak : peak_hgt;
    fEMG->SetParameter(0, A0);

    fEMG->SetParLimits(0, 1e-3 * peak_hgt, 1e4 * peak_hgt);
    fEMG->SetParLimits(1, xlo, xhi);
    fEMG->SetParLimits(2, bw, 0.7 * (xhi - xlo));
    fEMG->SetParLimits(3, 0.05 * bw, xhi - xlo);

    if (line_color >= 0) {
        fEMG->SetLineColor(line_color);
        fEMG->SetLineWidth(2);
    }

    TString opt = attach ? "L S R Q" : "L S R Q N";
    int fit_status = h->Fit(fEMG, opt.Data());

    bool sane = (fit_status == 0)
             && (fEMG->GetParameter(2) > 0.0)
             && (fEMG->GetParameter(3) > 0.0);

    if (!sane) {
        std::cerr << "[WARNING] EMG non convergente per " << h->GetName()
                  << " -> fallback FitGaussianCore." << std::endl;
        if (attach && h->GetListOfFunctions())
            h->GetListOfFunctions()->Remove(fEMG);
        delete fEMG;
        return FitGaussianCore(h, 0.90, Form("%s_gfb", fit_name));
    }

    double mode = fEMG->GetMaximumX(xlo, xhi);

res.center     = mode;
    res.center_err = fEMG->GetParError(1);
    res.width      = fabs(fEMG->GetParameter(2));
    res.width_err  = fEMG->GetParError(2);
    res.chi2_ndf   = (fEMG->GetNDF() > 0)
                   ? fEMG->GetChisquare() / fEMG->GetNDF() : -1.0;
    res.fit_ok     = true;
    return res;
}


// ==========================================================================
//  SEZIONE 3e: FITTER DEDICATI ALLE OSSERVABILI (beta, TOF, 1/v)
// ==========================================================================
//
//  Tre metodi indipendenti e confrontabili per estrarre il "valore tipico" e
//  l'incertezza di una distribuzione fisica asimmetrica:
//    [1] FitGaussianPeak : gaussiana ristretta a una finestra attorno al PICCO
//                          (ignora le code, valore robusto ma "ingenuo");
//    [2] FitObservableEMG: EMG (gaussiana + coda esponenziale), con riflessione
//                          per le distribuzioni a coda SINISTRA (beta);
//    [3] FitRooKeys      : KDE non parametrico (RooKeysPdf), moda + mediana,
//                          controllo modello-indipendente.
//
//  FISICA: per beta la coda e' a SINISTRA (TOF sovrastimato -> beta basso),
//  per TOF e 1/v la coda e' a DESTRA (stessa popolazione, variabile non
//  invertita). La moda dell'EMG e' il valore "tipico" non distorto dalla coda.

/// FitGaussianPeak(): gaussiana ristretta a [x_peak - k*sigma0, x_peak + k*sigma0],
/// dove x_peak e' il centro del bin massimo e sigma0 = FWHM/2.355 stimata
/// dall'istogramma. A differenza di FitGaussianCore (quantili 5-95%), la
/// finestra e' ANCORATA AL PICCO e simmetrica: cattura il "core" del picco
/// escludendo le code asimmetriche.
///
///   k_sigma : semiampiezza della finestra in unita' di sigma0 (default 1.5).
/// Restituisce mu (center) e sigma (width). Log-likelihood Poissoniana.
FitResult FitGaussianPeak(TH1D* h,
                          double k_sigma = 1.5,
                          const char* fit_name = "gpeak_fit") {
    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h->GetEntries();

    if (res.nentries < 30) {
        res.center     = h->GetMean();
        res.center_err = h->GetMeanError();
        res.width      = h->GetStdDev();
        res.width_err  = h->GetStdDevError();
        res.chi2_ndf   = -1.0;
        return res;
    }

    int    imax     = h->GetMaximumBin();
    double x_peak   = h->GetBinCenter(imax);
    double peak_hgt = h->GetBinContent(imax);
    double fwhm     = EstimateHistFWHM(h);
    double sigma0   = (fwhm > 0.0) ? fwhm / 2.355 : 2.0 * h->GetBinWidth(1);

    double x_lo = std::max(x_peak - k_sigma * sigma0, h->GetXaxis()->GetXmin());
    double x_hi = std::min(x_peak + k_sigma * sigma0, h->GetXaxis()->GetXmax());
    // Garanzia di un minimo di bin nel range.
    if (x_hi - x_lo < 5.0 * h->GetBinWidth(1)) {
        x_lo = std::max(x_peak - 5.0 * h->GetBinWidth(1), h->GetXaxis()->GetXmin());
        x_hi = std::min(x_peak + 5.0 * h->GetBinWidth(1), h->GetXaxis()->GetXmax());
    }

    TF1 *fG = new TF1(fit_name,
        "[0]*TMath::Exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",
        x_lo, x_hi);
    fG->SetParName(0, "A_peak");
    fG->SetParName(1, "mu");
    fG->SetParName(2, "sigma");
    fG->SetParameters(peak_hgt, x_peak, sigma0);
    fG->SetParLimits(0, 0.1, 100.0 * peak_hgt);
    fG->SetParLimits(1, x_lo, x_hi);
    fG->SetParLimits(2, h->GetBinWidth(1), x_hi - x_lo);

    int fit_status = h->Fit(fG, "L S R Q 0");

    res.center     = fG->GetParameter(1);
    res.center_err = fG->GetParError(1);
    res.width      = fabs(fG->GetParameter(2));
    res.width_err  = fG->GetParError(2);
    res.chi2_ndf   = (fG->GetNDF() > 0) ? fG->GetChisquare() / fG->GetNDF() : -1.0;
    res.fit_ok     = (fit_status == 0);
    return res;
}

/// FitObservableEMG(): fit EMG di un'osservabile, con gestione della coda.
///
///  left_tail = false (TOF, 1/v): coda a destra. Si chiama direttamente l'EMG
///    dritta (FitRiseTimeEMG). Per il TOF, il troncamento a sinistra (taglio
///    causale tof > h/c) cade sul fianco sinistro ripido: l'EMG modella picco
///    + coda destra senza problemi, e la MODA restituita non e' influenzata dal
///    troncamento (e' un massimo locale interno).
///
///  left_tail = true (beta): coda a SINISTRA. Strategia "rifletti-fitta-rifletti":
///    1. si sceglie un pivot P (estremo destro del range dell'istogramma);
///    2. si costruisce un istogramma ausiliario in x' = P - x (la coda sinistra
///       diventa destra: l'EMG dritta diventa applicabile);
///    3. si fitta l'EMG sull'istogramma riflesso;
///    4. si riflette indietro: moda_x = P - moda_x'. La sigma e' invariante per
///       riflessione; l'errore sul centro idem.
///
///  L'istogramma riflesso ausiliario viene creato, fittato e poi cancellato:
///  non viene scritto su file (e' solo uno strumento di calcolo).
///
///   h_obs     : istogramma dell'osservabile.
///   left_tail : true se la coda fisica e' a sinistra (beta).
///   fit_name  : nome univoco della TF1 (deve essere distinto per ogni chiamata).
FitResult FitObservableEMG(TH1D* h_obs,
                           bool left_tail,
                           const char* fit_name,
                           double& pivot_used) {
    pivot_used = 0.0;   // significativo solo per left_tail = true
    if (!left_tail) {
        // Coda a destra: EMG diretta. attach=false: non disegniamo qui (il
        // disegno e' gestito nel canvas dedicato di FitObservable).
        return FitRiseTimeEMG(h_obs, fit_name, -1, /*attach=*/false);
    }

// ---- Coda a sinistra: riflessione attorno al pivot P ----
    //  SCELTA DEL PIVOT: NON l'estremo destro dell'asse (sprecherebbe il range
    //  e lascerebbe il picco riflesso lontano dall'origine, con un init pessimo
    //  per l'EMG), ma il bordo destro della REGIONE POPOLATA: ultimo bin con
    //  conteggi + 2 larghezze di bin di margine. Cosi' dopo la riflessione il
    //  picco resta vicino all'origine e la coda (ora a destra) e' ben campionata.
    double xlo  = h_obs->GetXaxis()->GetXmin();
    int    nb   = h_obs->GetNbinsX();
    double bw   = h_obs->GetBinWidth(1);

    int last_nonempty = 1;
    for (int i = nb; i >= 1; i--) {
        if (h_obs->GetBinContent(i) > 0.0) { last_nonempty = i; break; }
    }
    double P = h_obs->GetBinCenter(last_nonempty) + 2.0 * bw;   // pivot adattivo
    pivot_used = P;     // comunicato al chiamante per il ridisegno della curva
    // Nuovo asse riflesso x' = P - x, crescente da 0 (= P - P) fino a (P - xlo).
    // Numero di bin proporzionale alla regione effettivamente coperta, con la
    // stessa larghezza di bin dell'originale (preserva la risoluzione).
    double new_hi = P - xlo;
    if (new_hi <= 0.0) new_hi = (h_obs->GetXaxis()->GetXmax() - xlo);  // guardia
    int nb_ref = std::max(10, (int)std::lround(new_hi / bw));
    TH1D* h_ref = new TH1D(Form("%s_reflected", h_obs->GetName()),
                           "reflected;x';conteggi", nb_ref, 0.0, new_hi);
    for (int i = 1; i <= nb; i++) {
        double xc = h_obs->GetBinCenter(i);
        double c  = h_obs->GetBinContent(i);
        if (c <= 0.0) continue;
        double xp = P - xc;
        if (xp < 0.0 || xp > new_hi) continue;   // fuori dal nuovo asse
        h_ref->Fill(xp, c);          // specchia il contenuto del bin
    }

    // Fit EMG dritta sull'istogramma riflesso.
    FitResult r_ref = FitRiseTimeEMG(h_ref, fit_name, -1, /*attach=*/false);

    // Riflessione inversa dei risultati: il centro torna in x = P - center'.
    FitResult res = r_ref;
    res.center = P - r_ref.center;   // moda nel sistema originale
    // width (sigma), center_err, width_err, chi2_ndf, nentries, fit_ok: invariati.

    delete h_ref;
    return res;
}

/// FitRooKeys(): stima non parametrica della densita' con RooKeysPdf (KDE), da
/// cui estrae MODA (argmax su griglia fine) e MEDIANA (dalla CDF cumulata).
///
///  RooKeysPdf costruisce una densita' liscia come somma di kernel gaussiani
///  centrati sui dati, con larghezza adattiva (opzione MirrorBoth per ridurre
///  il bias ai bordi). Non assume alcuna forma funzionale: e' il controllo
///  modello-indipendente richiesto (stesso approccio di Pietro&Rick).
///
///  STRATEGIA DI CALCOLO (senza dipendere da metodi interni di RooFit):
///    - si valuta la PDF su una griglia fine di NGRID punti nel range;
///    - MODA = ascissa del massimo della griglia;
///    - CDF discreta = somma cumulata (regola dei trapezi); MEDIANA = ascissa a
///      cui la CDF normalizzata raggiunge 0.5 (interpolata linearmente);
///    - errore sulla moda ~ FWHM_KDE / (2.355 * sqrt(N)) (proxy gaussiano
///      dell'errore standard del centro: la FWHM della KDE diviso 2.355 stima
///      una sigma efficace, /sqrt(N) la converte in errore della media).
///
///  RITORNA un FitResult con center = moda, width = sigma efficace (FWHM/2.355),
///  e in piu' la mediana viene passata via parametro di output 'median_out'.
///  Se RooFit fallisce o ci sono pochi eventi, fit_ok = false (fallback gestito
///  dal chiamante).
///
///   h_obs      : istogramma dell'osservabile (usato come RooDataHist).
///   var_lo,var_hi : range della variabile.
///   var_label  : nome della variabile (per RooRealVar).
///   median_out : (output) mediana stimata dalla KDE.
FitResult FitRooKeys(TH1D* h_obs,
                     double var_lo, double var_hi,
                     const char* var_label,
                     double& median_out,
                     double& median_err_out,
                     TGraph** kde_curve_out = nullptr) {
    FitResult res;
    res.fit_ok   = false;
    res.nentries = (int)h_obs->GetEntries();
    median_out     = -999.0;
    median_err_out = 0.0;

    if (res.nentries < 50) {
        res.center     = h_obs->GetMean();
        res.center_err = h_obs->GetMeanError();
        res.width      = h_obs->GetStdDev();
        res.width_err  = h_obs->GetStdDevError();
        res.chi2_ndf   = -1.0;
        return res;
    }

// --- Costruzione della KDE da RooDataSet ---
    //  RooKeysPdf richiede dati NON binnati (RooDataSet), non un RooDataHist.
    //  Ricostruiamo quindi un RooDataSet dai bin dell'istogramma: per ogni bin
    //  con contenuto c (arrotondato all'intero piu' vicino) inseriamo c entry
    //  poste al CENTRO del bin. Con i bin fini usati qui la discretizzazione al
    //  centro del bin e' un'approssimazione trascurabile della forma vera.
    //  Silenziamo i messaggi di RooFit per non inondare lo stdout (molti fit).
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    RooRealVar x(var_label, var_label, var_lo, var_hi);
    RooDataSet ds(Form("ds_%s", h_obs->GetName()),
                  "RooDataSet da istogramma", RooArgSet(x));
    for (int i = 1; i <= h_obs->GetNbinsX(); i++) {
        double xc = h_obs->GetBinCenter(i);
        if (xc < var_lo || xc > var_hi) continue;        // fuori dominio di x
        long   c  = (long)std::lround(h_obs->GetBinContent(i));
        if (c <= 0) continue;
        x.setVal(xc);
        for (long j = 0; j < c; j++) ds.add(RooArgSet(x));
    }
    if (ds.numEntries() < 50) {
        // Troppo pochi dati ricostruiti: fallback statistico semplice.
        res.center     = h_obs->GetMean();
        res.center_err = h_obs->GetMeanError();
        res.width      = h_obs->GetStdDev();
        res.width_err  = h_obs->GetStdDevError();
        res.chi2_ndf   = -1.0;
        return res;
    }
    // RooKeysPdf: MirrorBoth riduce il bias ai bordi; rho=2 regola la larghezza
    // dei kernel (gia' adattiva sui dati).
    RooKeysPdf kde(Form("kde_%s", h_obs->GetName()), "KDE",
                   x, ds, RooKeysPdf::MirrorBoth, 2.0);

    // --- Campionamento della PDF su griglia fine ---
    //  Per evitare del tutto il vincolo di RooFit recente su getVal() (che
    //  vieta RooArgSet temporanei nel set di normalizzazione), NON valutiamo la
    //  pdf punto-punto: chiediamo a RooFit un istogramma gia' campionato della
    //  densita' tramite createHistogram(). Internamente RooFit gestisce la
    //  normalizzazione in modo corretto. Poi lavoriamo sui bin di questo TH1.
    const int NGRID = 2000;
    x.setBins(NGRID);
    TH1* h_kde = kde.createHistogram(Form("hkde_%s", h_obs->GetName()),
                                     x, RooFit::Binning(NGRID, var_lo, var_hi));
    if (!h_kde || h_kde->GetEntries() == 0 || h_kde->Integral() <= 0.0) {
        if (h_kde) delete h_kde;
        res.center     = h_obs->GetMean();
        res.center_err = h_obs->GetMeanError();
        res.width      = h_obs->GetStdDev();
        res.width_err  = h_obs->GetStdDevError();
        res.chi2_ndf   = -1.0;
        return res;
    }

    std::vector<double> gx(NGRID), gp(NGRID);
    double pmax = -1.0; int imax = 0;
    for (int i = 0; i < NGRID; i++) {
        // createHistogram restituisce bin 1..NGRID; centro e contenuto del bin.
        double xi = h_kde->GetBinCenter(i + 1);
        double pi = h_kde->GetBinContent(i + 1);
        gx[i] = xi; gp[i] = pi;
        if (pi > pmax) { pmax = pi; imax = i; }
    }
    double dxg = (var_hi - var_lo) / (NGRID - 1);
    if (pmax <= 0.0) { delete h_kde; return res; }   // KDE degenere

    // --- MODA = ascissa del massimo ---
    double mode = gx[imax];

    // --- CDF cumulata (trapezi) e MEDIANA ---
    std::vector<double> cdf(NGRID, 0.0);
    for (int i = 1; i < NGRID; i++)
        cdf[i] = cdf[i - 1] + 0.5 * (gp[i] + gp[i - 1]) * dxg;
    double norm = cdf[NGRID - 1];
    double median = mode;
    if (norm > 0.0) {
        for (int i = 1; i < NGRID; i++) {
            double c0 = cdf[i - 1] / norm;
            double c1 = cdf[i]     / norm;
            if (c0 <= 0.5 && c1 >= 0.5) {
                double frac = (fabs(c1 - c0) > 1e-12) ? (0.5 - c0) / (c1 - c0) : 0.0;
                median = gx[i - 1] + frac * (gx[i] - gx[i - 1]);
                break;
            }
        }
    }
    median_out = median;

    // --- FWHM della KDE -> sigma efficace -> errore sulla moda ---
    double half = 0.5 * pmax;
    // bordo sinistro
    double xL = gx[0];
    for (int i = imax; i >= 1; i--) {
        if (gp[i] < half) {
            double frac = (fabs(gp[i + 1] - gp[i]) > 1e-12)
                        ? (half - gp[i]) / (gp[i + 1] - gp[i]) : 0.0;
            xL = gx[i] + frac * (gx[i + 1] - gx[i]);
            break;
        }
    }
    // bordo destro
    double xR = gx[NGRID - 1];
    for (int i = imax; i < NGRID - 1; i++) {
        if (gp[i] < half) {
            double frac = (fabs(gp[i - 1] - gp[i]) > 1e-12)
                        ? (half - gp[i]) / (gp[i - 1] - gp[i]) : 0.0;
            xR = gx[i] + frac * (gx[i - 1] - gx[i]);
            break;
        }
    }
    double fwhm_kde   = (xR > xL) ? (xR - xL) : h_obs->GetStdDev() * 2.355;
    double sigma_eff  = fwhm_kde / 2.355;
    double n          = (double)res.nentries;
    double mode_err   = (n > 0.0) ? sigma_eff / std::sqrt(n) : 0.0;
    // Errore standard della MEDIANA per un campione di n eventi:
    //   sigma_mediana ~ sqrt(pi/2) * sigma / sqrt(n) ~ 1.2533 * sigma/sqrt(n).
    //  (vale asintoticamente per distribuzione ~gaussiana; qui usiamo sigma_eff
    //  della KDE come stima di sigma). E' ~25% maggiore dell'errore della media.
    median_err_out = (n > 0.0) ? 1.2533 * sigma_eff / std::sqrt(n) : 0.0;

res.center     = mode;
    res.center_err = mode_err;
    res.width      = sigma_eff;
    res.width_err  = 0.0;            // non stimato per la KDE
    res.chi2_ndf   = -1.0;          // non applicabile alla KDE
    res.fit_ok     = true;

    // --- Curva KDE sovrapponibile all'istogramma fisico (se richiesta) ---
    //  La KDE e' una densita' (area = 1). Per sovrapporla all'istogramma in
    //  conteggi va riscalata per (N_entries * bin_width): cosi' l'area sotto la
    //  curva eguaglia l'area dell'istogramma e il confronto e' visivamente
    //  corretto. La normalizzazione gp e' gia' tale che sum(gp)*dxg ~ 1, quindi
    //  moltiplichiamo ogni punto per N*binwidth.
    if (kde_curve_out) {
        double Nentries = (double)h_obs->GetEntries();
        double binw     = h_obs->GetBinWidth(1);
        // Normalizzazione esatta: ricalcolo l'integrale della griglia (trapezi)
        // per non dipendere da assunzioni, e scalo a N*binw.
        double area_grid = 0.0;
        for (int i = 1; i < NGRID; i++)
            area_grid += 0.5 * (gp[i] + gp[i - 1]) * dxg;
        double scale = (area_grid > 0.0) ? (Nentries * binw / area_grid) : 0.0;

        TGraph* gk = new TGraph(NGRID);
        for (int i = 0; i < NGRID; i++) gk->SetPoint(i, gx[i], gp[i] * scale);
        gk->SetName(Form("gkde_%s", h_obs->GetName()));
        *kde_curve_out = gk;
    }

    delete h_kde;                    // l'istogramma di servizio non serve piu'
    return res;
}


// ==========================================================================
//  SEZIONE 3f: ORCHESTRATORE DEI FIT SULLE OSSERVABILI
// ==========================================================================
//
//  FitObservable() applica i 3 metodi (gaussiano-picco, EMG, RooKeys) a un
//  singolo istogramma, disegna un canvas con istogramma + curve + pannello dei
//  risultati, e raccoglie tutto in un ObsFitResult per la tabella riepilogativa.

/// ObsFitResult: risultati dei 3 metodi di fit su una (pipeline, osservabile).
struct ObsFitResult {
    std::string pipeline;    // etichetta pipeline (es. "Hybrid_tot")
    std::string observable;  // etichetta osservabile (es. "beta")
    std::string unit;        // unita' di misura (es. "", "ns", "ns/cm")

    FitResult gpeak;         // gaussiana ristretta al picco
    FitResult emg;           // EMG (moda come center)
    FitResult rkeys;         // RooKeysPdf (moda come center)
    double    rkeys_median;      // mediana dalla KDE
    double    rkeys_median_err;  // errore standard della mediana

    bool valid = false;      // true se almeno un metodo e' andato a buon fine
};
/// FmtVE(): formatta "valore +/- errore" scegliendo i decimali in modo che
/// l'errore mostri ~2 cifre significative. Risolve il caso di osservabili
/// piccole (es. 1/v ~ 0.035, errore ~ 3e-5) che con "%.4f" darebbero "0.0000".
///   ndig_min : decimali minimi (default 4), per uniformita' con le altre colonne.
/// Restituisce una std::string pronta da stampare/inserire in un TLatex.
std::string FmtVE(double v, double e, int ndig_min = 4) {
    // Numero di decimali: tanti quanti servono perche' l'errore (se > 0) abbia
    // almeno 2 cifre significative; mai meno di ndig_min.
    int nd = ndig_min;
    if (e > 0.0) {
        // posizione della prima cifra significativa dell'errore:
        // floor(log10(e)) = -k  =>  servono k+1 decimali per 2 cifre signif.
        int k = (int)std::floor(std::log10(e));
        int need = -k + 1;            // 2 cifre significative
        if (need > nd) nd = need;
        if (nd > 8) nd = 8;           // tetto di sicurezza
    }
    char buf[96];
    snprintf(buf, sizeof(buf), "%.*f #pm %.*f", nd, v, nd, e);
    return std::string(buf);
}

/// FmtV(): come FmtVE ma per un solo valore (es. mediana senza errore esplicito).
std::string FmtV(double v, int ndig = 4) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.*f", ndig, v);
    return std::string(buf);
}

/// FitObservable(): esegue i 3 fit su h, disegna il canvas e ritorna i risultati.
///
///   h          : istogramma dell'osservabile (NON viene modificato nei contenuti).
///   var_lo,var_hi : range della variabile (per RooKeys e per le curve).
///   var_label  : nome variabile per RooFit (univoco: include pipeline+osservabile).
///   pipe_label : etichetta pipeline (per titolo/tabella).
///   obs_label  : etichetta osservabile (per titolo/tabella).
///   unit       : unita' di misura (per la tabella).
///   left_tail  : true se la coda fisica e' a sinistra (beta) -> EMG riflessa.
///   ref_value  : valore di riferimento fisico da marcare con linea (es. 1.0
///                per beta, 1/c per 1/v); passare un NaN per non disegnarla.
///   fout       : file ROOT corrente (il canvas viene scritto qui).
ObsFitResult FitObservable(TH1D* h,
                           double var_lo, double var_hi,
                           const char* var_label,
                           const char* pipe_label,
                           const char* obs_label,
                           const char* unit,
                           bool left_tail,
                           double ref_value,
                           TFile* fout) {
    ObsFitResult out;
    out.pipeline   = pipe_label;
    out.observable = obs_label;
    out.unit       = unit;

// Tag univoco per i nomi delle TF1/oggetti ROOT di questa combinazione.
    std::string tag = std::string(obs_label) + "_" + pipe_label;
    // Sanitizza i caratteri non ammessi nei nomi ROOT: lo SLASH '/' e' un
    // separatore di directory (romperebbe il salvataggio del canvas e il
    // ritrovamento delle TF1 via FindObject), lo spazio idem. Li mappiamo:
    //   '/' -> 'p'  (es. "1/v" -> "1pv"),   ' ' -> '_'.
    for (auto& ch : tag) {
        if (ch == ' ') ch = '_';
        else if (ch == '/') ch = 'p';
    }

    if (h->GetEntries() < 30) {
        std::cout << "[WARNING] FitObservable: pochi eventi ("
                  << h->GetEntries() << ") per " << tag
                  << ", fit saltati." << std::endl;
        return out;   // valid resta false
    }

    // ---- [1] Gaussiana ristretta al picco ----
    std::string n_gp = "gpeak_" + tag;
    out.gpeak = FitGaussianPeak(h, 1.5, n_gp.c_str());

// ---- [2] EMG (riflessa per coda sinistra) ----
    std::string n_emg = "emg_" + tag;
    double emg_pivot = 0.0;
    out.emg = FitObservableEMG(h, left_tail, n_emg.c_str(), emg_pivot);

// ---- [3] RooKeysPdf (moda + mediana + curva KDE per il disegno) ----
    TGraph* kde_curve = nullptr;
    out.rkeys = FitRooKeys(h, var_lo, var_hi, var_label, out.rkeys_median,
                           out.rkeys_median_err, &kde_curve);

    out.valid = out.gpeak.fit_ok || out.emg.fit_ok || out.rkeys.fit_ok;

    // ======================================================================
    //  CANVAS: istogramma + curve dei 3 metodi + pannello risultati
    // ======================================================================
    if (fout) fout->cd();
    TCanvas* c = new TCanvas(Form("c_fit_%s", tag.c_str()),
        Form("Fit %s (%s)", obs_label, pipe_label), 900, 650);
    c->SetGrid();

    h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    h->SetTitle(Form("Fit %s (%s);%s%s;Conteggi",
                     obs_label, pipe_label, obs_label,
                     (strlen(unit) > 0) ? Form(" [%s]", unit) : ""));
    h->Draw("HIST");
    double hmax = h->GetMaximum();

    // ---- Curva gaussiana-picco (verde), ridisegnata dal suo range ----
    if (out.gpeak.fit_ok) {
        TF1* fgp = (TF1*)gROOT->FindObject(n_gp.c_str());
        if (fgp) {
            fgp->SetLineColor(kGreen + 2);
            fgp->SetLineWidth(2);
            fgp->SetLineStyle(1);
            fgp->Draw("same");
        }
    }

    // ---- Curva EMG (rosso) ----
    //  Per coda destra (TOF/1v) la TF1 vive nel sistema originale: la disegno
    //  direttamente. Per coda sinistra (beta) la TF1 vive nel sistema riflesso
    //  x' = P - x: la ricostruisco come TGraph valutandola e riflettendo le
    //  ascisse, cosi' appare correttamente nel sistema fisico.
    if (out.emg.fit_ok) {
        if (!left_tail) {
            TF1* femg = (TF1*)gROOT->FindObject(n_emg.c_str());
            if (femg) {
                femg->SetLineColor(kRed + 1);
                femg->SetLineWidth(2);
                femg->Draw("same");
            }
        } else {
            // EMG riflessa: la TF1 (nome n_emg) e' definita in x' su [0, P-xlo].
            TF1* femg = (TF1*)gROOT->FindObject(n_emg.c_str());
            if (femg) {
                double P = emg_pivot;   // pivot ADATTIVO usato in FitObservableEMG
                const int NG = 500;
                std::vector<double> gx(NG), gy(NG);
                double xpmin = femg->GetXmin();
                double xpmax = femg->GetXmax();
                for (int i = 0; i < NG; i++) {
                    double xp = xpmin + (xpmax - xpmin) * i / (NG - 1);
                    double yv = femg->Eval(xp);
                    gx[i] = P - xp;       // riflessione inversa dell'ascissa
                    gy[i] = yv;
                }
                TGraph* gemg = new TGraph(NG, gx.data(), gy.data());
                gemg->SetName(Form("gemg_%s", tag.c_str()));
                gemg->SetLineColor(kRed + 1);
                gemg->SetLineWidth(2);
                gemg->Draw("L same");
            }
        }
    }

// ---- Curva KDE RooKeysPdf (blu) + moda/mediana ----
    if (out.rkeys.fit_ok && kde_curve) {
        kde_curve->SetLineColor(kBlue + 1);
        kde_curve->SetLineWidth(2);
        kde_curve->SetLineStyle(1);
        kde_curve->Draw("L same");   // curva KDE sovrapposta all'istogramma
    }
    // ---- Moda e mediana RooKeys (blu pieno / blu tratteggiato) ----
    if (out.rkeys.fit_ok) {
        TLine* l_mode = new TLine(out.rkeys.center, 0.0, out.rkeys.center, hmax * 0.95);
        l_mode->SetLineColor(kBlue + 1);
        l_mode->SetLineWidth(2);
        l_mode->Draw("same");
        if (out.rkeys_median > -900.0) {
            TLine* l_med = new TLine(out.rkeys_median, 0.0,
                                     out.rkeys_median, hmax * 0.95);
            l_med->SetLineColor(kBlue + 1);
            l_med->SetLineStyle(2);
            l_med->SetLineWidth(2);
            l_med->Draw("same");
        }
    }

    // ---- Linea di riferimento fisico (es. beta = 1, 1/v = 1/c) ----
    if (ref_value == ref_value /* non-NaN */ &&
        ref_value >= var_lo && ref_value <= var_hi) {
        TLine* l_ref = new TLine(ref_value, 0.0, ref_value, hmax * 0.95);
        l_ref->SetLineColor(kBlack);
        l_ref->SetLineStyle(3);
        l_ref->SetLineWidth(2);
        l_ref->Draw("same");
    }

    // ---- Pannello dei risultati ----
    TPaveText* pv = new TPaveText(0.58, 0.50, 0.89, 0.89, "NDC");
    pv->SetFillColor(0);
    pv->SetBorderSize(1);
    pv->SetTextAlign(12);
    pv->SetTextFont(42);
    pv->SetTextSize(0.028);
    pv->AddText(Form("Entries = %lld", (Long64_t)h->GetEntries()));
    pv->AddText("");
    pv->AddText("#color[418]{#bf{Gauss-picco}}");
    if (out.gpeak.fit_ok)
        pv->AddText(Form("#mu = %s", FmtVE(out.gpeak.center, out.gpeak.center_err).c_str()));
    else
        pv->AddText("non convergente");
    pv->AddText("#color[633]{#bf{EMG (moda)}}");
    if (out.emg.fit_ok)
        pv->AddText(Form("moda = %s", FmtVE(out.emg.center, out.emg.center_err).c_str()));
    else
        pv->AddText("non convergente");
    pv->AddText("#color[600]{#bf{RooKeysPdf}}");
    if (out.rkeys.fit_ok) {
        pv->AddText(Form("moda = %s", FmtVE(out.rkeys.center, out.rkeys.center_err).c_str()));
        pv->AddText(Form("mediana = %s",
                         FmtVE(out.rkeys_median, out.rkeys_median_err).c_str()));
    } else {
        pv->AddText("non disponibile");
    }
    pv->Draw();

    c->Write(Form("Canvas_fit_%s", tag.c_str()));
    return out;
}


// ==========================================================================
//  SEZIONE 4: MODELLO A(TOT) E RECUPERO DEI SEGNALI CLIPPATI
// ==========================================================================
//
//  In v9 il recupero dei clippati usa ESCLUSIVAMENTE il polinomio cubico
//  A(TOT) calibrato dalla v13 e caricato da fit_params: nessun modello
//  slew-rate, nessun modello esponenziale.

/// EvalAmplitudeFromTOT(): valuta A(TOT) in mV usando il polinomio cubico
/// caricato per il canale ch_index, con lo schema di Horner.
///
/// RITORNA: ampiezza ricostruita [mV], oppure -1.0 se la calibrazione non e'
///   disponibile per quel canale o se il TOT non e' fisico (tot <= 0).
double EvalAmplitudeFromTOT(double tot, int ch_index) {
    if (ch_index < 0 || ch_index >= MAX_CHANNELS) return -1.0;
    if (!gTOT_poly_loaded[ch_index])              return -1.0;
    if (tot <= 0.0)                               return -1.0;

    const double *coeff = gTOT_poly_PMT[ch_index];
    double A = coeff[gTOT_poly_degree];
    for (int k = gTOT_poly_degree - 1; k >= 0; --k) {
        A = A * tot + coeff[k];
    }
    return A;
}

/// TOTCalibrationAvailable(): true se la calibrazione polinomiale A(TOT) e'
/// stata caricata e utilizzabile per il canale richiesto.
bool TOTCalibrationAvailable(int ch_index) {
    if (ch_index < 0 || ch_index >= MAX_CHANNELS) return false;
    if (!gTOT_poly_loaded[ch_index])              return false;
    if (gTOT_poly_degree < 1)                     return false;
    if (gTOT_Amax_PMT[ch_index] <= 0.0)           return false;
    return true;
}

/// RecoverClippedCFD_TOT(): unica funzione di recupero dei segnali clippati.
///   Per un canale CLIPPATO stima l'ampiezza vera dal TOT alla soglia di
///   riferimento (polinomio cubico) e ricalcola il tempo CFD su f*A_rec.
///
/// ALGORITMO (identico alla v13):
///   1. se il canale non e' clippato -> non fa nulla (si usa il CFD standard);
///   2. se il canale e' oscillante -> fallisce;
///   3. se il TOT alla soglia di riferimento non e' valido -> fallisce;
///   4. A_rec = EvalAmplitudeFromTOT(tot_ref, ch_index); deve essere > 0;
///   5. cap di estrapolazione: A_rec <= TOT_AREC_MAX_FACTOR * Amax;
///   6. soglia CFD: V_cfd = baseline - CFD_FRACTION * A_rec;
///   7. V_cfd deve cadere SOPRA il plateau saturato (V_cfd > CLIP_V_LO +
///      TOT_CLIP_MARGIN), altrimenti il crossing non esiste -> fallisce;
///   8. ricerca del crossing sul fronte di salita con interpolazione lineare.
///
/// EFFETTO: popola cd.amplitude_tot, cd.t_cfd_rec_tot, cd.cfd_recovered_tot.
void RecoverClippedCFD_TOT(ChannelData &cd, int ch_index) {

    if (!cd.is_clipped)                           return;
    if (cd.is_oscillating)                        return;
    if (ch_index < 0 || ch_index >= MAX_CHANNELS) return;
    if (!TOTCalibrationAvailable(ch_index))       return;

    const int kref = TOT_CALIB_THR_IDX;
    if (kref < 0 || kref >= NTOT_THR) return;
    if (!cd.tot_ok[kref])             return;

    double tot_ref = cd.tot_q[kref];

    // Stima dell'ampiezza ricostruita (polinomio cubico).
    double A_rec = EvalAmplitudeFromTOT(tot_ref, ch_index);
    if (A_rec <= 0.0) return;
    cd.amplitude_tot = A_rec;

    // Cap di estrapolazione (anti-outlier).
    double A_cap = TOT_AREC_MAX_FACTOR * gTOT_Amax_PMT[ch_index];
    if (A_cap > 0.0 && A_rec > A_cap) return;

    // Soglia CFD ricostruita.
    double v_thr = cd.baseline - CFD_FRACTION * A_rec;

    // VERIFICA CRITICA: la soglia CFD deve cadere SOPRA il plateau di clipping.
    if (v_thr <= CLIP_V_LO + TOT_CLIP_MARGIN) return;

    // Ricerca del crossing sul fronte di salita.
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


// ==========================================================================
//  SEZIONE 5: PARSING XML DRS4 (identico alla calibrazione)
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
        strncpy(buf, line.c_str(), sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';

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
//  SEZIONE 5b: CACHE DI PARSING (XML -> ROOT)
// ==========================================================================

/// GetFileMTime(): ritorna il tempo di ultima modifica del file [secondi epoch],
/// oppure -1 se il file non esiste o non e' accessibile.
static long GetFileMTime(const char* path) {
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    return (long)st.st_mtime;
}

/// DatasetNameFromPath(): ricava un nome stabile dal path XML.
/// Esempio: "/dati/run_01.xml" -> "run_01". Il nome viene usato solo per
/// nominare il ROOT di cache nella sottocartella parsed_cache.
static std::string DatasetNameFromPath(const std::string& xml_path) {
    size_t slash = xml_path.find_last_of("/\\");
    std::string base = (slash == std::string::npos)
                     ? xml_path
                     : xml_path.substr(slash + 1);
    size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    return base.empty() ? std::string("dataset") : base;
}

/// WriteDatasetCache(): salva in ROOT le waveform grezze di un dataset.
///
/// La cache NON contiene t_cfd, TOT, flag di qualita' o grandezze fisiche:
/// contiene solo cio' che viene letto dall'XML. Questo mantiene la cache
/// compatibile con modifiche future agli algoritmi di ricostruzione.
bool WriteDatasetCache(const char* cache_path,
                       const std::vector<EventData>& events) {

    TDirectory* save = gDirectory;   // non perdere il file ROOT di output aperto

    TFile* fc = new TFile(cache_path, "RECREATE");
    if (!fc || !fc->IsOpen()) {
        if (save) save->cd();
        delete fc;
        return false;
    }
    fc->cd();

    TParameter<int> fmt("PARSE_CACHE_FORMAT", PARSE_CACHE_FORMAT);
    fmt.Write();

    TTree* t = new TTree("events_cache", "Cache forme d'onda DRS4 parsate");

Int_t serial, trigger_cell, board_serial, nchannels;
    Int_t channel_ids[MAX_CHANNELS], scaler[MAX_CHANNELS], nsamp[MAX_CHANNELS];
    static Float_t wf_t[MAX_CHANNELS][MAX_SAMPLES];
    static Float_t wf_v[MAX_CHANNELS][MAX_SAMPLES];
    // Timestamp <Time> dell'XML salvato come array di char a lunghezza fissa.
    // Il formato "YYYY/MM/DD HH:MM:SS.sss" e' lungo 23 caratteri; 64 e' ampio.
    char ts_str[64];

    t->Branch("serial",       &serial,       "serial/I");
    t->Branch("trigger_cell", &trigger_cell, "trigger_cell/I");
    t->Branch("board_serial", &board_serial, "board_serial/I");
    t->Branch("nchannels",    &nchannels,    "nchannels/I");
    t->Branch("channel_ids",  channel_ids,   Form("channel_ids[%d]/I", MAX_CHANNELS));
    t->Branch("scaler",       scaler,        Form("scaler[%d]/I",      MAX_CHANNELS));
    t->Branch("nsamp",        nsamp,         Form("nsamp[%d]/I",       MAX_CHANNELS));
    t->Branch("wf_t",         wf_t,          Form("wf_t[%d][%d]/F", MAX_CHANNELS, MAX_SAMPLES));
    t->Branch("wf_v",         wf_v,          Form("wf_v[%d][%d]/F", MAX_CHANNELS, MAX_SAMPLES));
    t->Branch("ts_str",       ts_str,        "ts_str[64]/C");   // <Time> dell'XML

for (const auto& e : events) {
        serial       = e.serial;
        trigger_cell = e.trigger_cell;
        board_serial = e.board_serial;
        nchannels    = e.nchannels;

        // Copia sicura del timestamp nell'array a lunghezza fissa (troncamento
        // difensivo + terminatore nullo garantito).
        strncpy(ts_str, e.timestamp.c_str(), sizeof(ts_str) - 1);
        ts_str[sizeof(ts_str) - 1] = '\0';

        for (int k = 0; k < MAX_CHANNELS; k++) {
            channel_ids[k] = (k < e.nchannels) ? e.channel_ids[k] : 0;
            scaler[k]      = e.scaler[k];

            int ns = (k < e.nchannels) ? e.ch[k].nsamples : 0;
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

    if (save) save->cd();
    return true;
}

/// ReadDatasetCache(): rilegge un dataset dalla cache ROOT.
///
/// Dopo aver ricostruito EventData, riesegue AnalyzeChannel() su ogni canale:
/// la cache velocizza il parsing XML, ma non congela la logica di analisi.
/// Ritorna il numero di eventi caricati, oppure -1 se la cache non e' valida.
int ReadDatasetCache(const char* cache_path, std::vector<EventData>& events) {

    TDirectory* save = gDirectory;

    TFile* fc = new TFile(cache_path, "READ");
    if (!fc || !fc->IsOpen()) {
        if (save) save->cd();
        delete fc;
        return -1;
    }

    TParameter<int>* fmt = (TParameter<int>*)fc->Get("PARSE_CACHE_FORMAT");
    if (!fmt || fmt->GetVal() != PARSE_CACHE_FORMAT) {
        std::cerr << "[CACHE] Formato incompatibile in " << cache_path
                  << " -> riparsing." << std::endl;
        fc->Close();
        delete fc;
        if (save) save->cd();
        return -1;
    }

    TTree* t = (TTree*)fc->Get("events_cache");
    if (!t) {
        fc->Close();
        delete fc;
        if (save) save->cd();
        return -1;
    }

Int_t serial, trigger_cell, board_serial, nchannels;
    Int_t channel_ids[MAX_CHANNELS], scaler[MAX_CHANNELS], nsamp[MAX_CHANNELS];
    static Float_t wf_t[MAX_CHANNELS][MAX_SAMPLES];
    static Float_t wf_v[MAX_CHANNELS][MAX_SAMPLES];
    char ts_str[64] = {0};   // timestamp <Time> (cache v2+)

    t->SetBranchAddress("serial",       &serial);
    t->SetBranchAddress("trigger_cell", &trigger_cell);
    t->SetBranchAddress("board_serial", &board_serial);
    t->SetBranchAddress("nchannels",    &nchannels);
    t->SetBranchAddress("channel_ids",  channel_ids);
    t->SetBranchAddress("scaler",       scaler);
    t->SetBranchAddress("nsamp",        nsamp);
    t->SetBranchAddress("wf_t",         wf_t);
    t->SetBranchAddress("wf_v",         wf_v);
    // Il branch ts_str esiste solo nelle cache v2+. Il controllo di
    // PARSE_CACHE_FORMAT a monte garantisce che qui la cache sia gia' v2,
    // ma collego il branch in modo difensivo solo se presente.
    if (t->GetBranch("ts_str")) t->SetBranchAddress("ts_str", ts_str);

    Long64_t nentries = t->GetEntries();
    events.clear();
    events.reserve((size_t)nentries);

    for (Long64_t ie = 0; ie < nentries; ie++) {
        t->GetEntry(ie);

EventData e = EventData();
        e.serial       = serial;
        e.trigger_cell = trigger_cell;
        e.board_serial = board_serial;
        e.nchannels    = std::max(0, std::min((int)nchannels, MAX_CHANNELS));
        ts_str[sizeof(ts_str) - 1] = '\0';   // sicurezza: stringa terminata
        e.timestamp    = std::string(ts_str);

        for (int k = 0; k < MAX_CHANNELS; k++) {
            e.channel_ids[k] = channel_ids[k];
            e.scaler[k]      = scaler[k];

            int ns = std::max(0, std::min((int)nsamp[k], MAX_SAMPLES));
            e.ch[k].nsamples = ns;
            for (int i = 0; i < ns; i++) {
                e.ch[k].time[i]    = wf_t[k][i];
                e.ch[k].voltage[i] = wf_v[k][i];
            }
        }

        for (int k = 0; k < e.nchannels; k++) AnalyzeChannel(e.ch[k]);
        events.push_back(e);
    }

    fc->Close();
    delete fc;
    if (save) save->cd();
    return (int)events.size();
}

/// LoadOrParseDataset(): punto unico per ottenere gli eventi di un dataset.
///
/// Se la cache e' attiva e aggiornata, carica il ROOT gia' parsato. Altrimenti
/// legge l'XML e poi scrive/aggiorna la cache per le esecuzioni successive.
int LoadOrParseDataset(const char* xml_path, std::vector<EventData>& events) {

    std::string xmls(xml_path);
    size_t slash = xmls.find_last_of("/\\");
    std::string dir = (slash == std::string::npos)
                    ? std::string(".")
                    : xmls.substr(0, slash);
    std::string name       = DatasetNameFromPath(xmls);
    std::string cache_dir  = dir + "/" + PARSE_CACHE_SUBDIR;
    std::string cache_path = cache_dir + "/" + name + ".root";

    if (USE_PARSE_CACHE && !FORCE_REPARSE) {
        long t_xml   = GetFileMTime(xml_path);
        long t_cache = GetFileMTime(cache_path.c_str());

        if (t_xml >= 0 && t_cache >= 0 && t_cache >= t_xml) {
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

    int n = ParseXML(xml_path, events);
    if (n <= 0) return n;

    if (USE_PARSE_CACHE) {
        gSystem->mkdir(cache_dir.c_str(), kTRUE);
        if (WriteDatasetCache(cache_path.c_str(), events)) {
            std::cout << "[CACHE] Salvato " << cache_path << std::endl;
        } else {
            std::cerr << "[CACHE] Impossibile scrivere " << cache_path
                      << " (parsing comunque OK)." << std::endl;
        }
    }

return n;
}


// ==========================================================================
//  SEZIONE 5c: FUNZIONI AUSILIARIE DEL FILTRO FIFO (DE10-Nano)
// ==========================================================================

/// DaysFromCivil(): converte una data civile (anno, mese, giorno) nel numero di
/// giorni trascorsi dall'epoch 1970-01-01 (algoritmo di Howard Hinnant, esatto
/// per qualunque data del calendario gregoriano proltettico). Serve a confrontare
/// in modo robusto timestamp che cadono in giorni diversi, senza dipendere dal
/// fuso orario del sistema o da <ctime>.
///
/// Riferimento: http://howardhinnant.github.io/date_algorithms.html
static long long DaysFromCivil(int y, unsigned m, unsigned d) {
    // Sposta l'inizio dell'anno a marzo: febbraio (e il suo 29) finisce in coda.
    y -= (m <= 2);
    const int      era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = (unsigned)(y - era * 400);              // anno nell'era [0,399]
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1; // giorno nell'anno [0,365]
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;  // giorno nell'era [0,146096]
    return (long long)era * 146097 + (long long)doe - 719468;    // -719468 = offset al 1970-01-01
}


/// ParseFifoStartTime(): estrae dal NOME del file FIFO la data e l'ora di avvio
/// acquisizione FPGA, e popola fd.start_s_day e fd.start_day_index.
///
/// Formato atteso: "...FIFOread_YYYYMMDD-HHMMSS.txt" (eventuali directory nel
/// path vengono ignorate: si usa solo il basename).
///
/// Nell'interpretazione A questo istante e' quello in cui il contatore di clock
/// FPGA vale 0 (inizio del buffer 0).
///
/// Ritorna true in caso di successo, false in caso di errore di formato.
static bool ParseFifoStartTime(const std::string &fifo_path, FifoData &fd) {

    // --- Estrai il basename (dopo l'ultima '/' o '\') ---
    size_t slash = fifo_path.find_last_of("/\\");
    std::string fname = (slash == std::string::npos)
                        ? fifo_path
                        : fifo_path.substr(slash + 1);

    // --- Localizza il pattern "FIFOread_" ---
    const std::string tag = "FIFOread_";
    size_t pos = fname.find(tag);
    if (pos == std::string::npos) {
        std::cerr << "[ERRORE] ParseFifoStartTime: pattern 'FIFOread_' assente in: "
                  << fname << std::endl;
        return false;
    }

    // Dopo "FIFOread_" servono 8 cifre (data) + '-' + 6 cifre (ora) = 15 caratteri.
    size_t d0 = pos + tag.size();
    if (fname.size() < d0 + 15) {
        std::cerr << "[ERRORE] ParseFifoStartTime: nome file troppo corto: "
                  << fname << std::endl;
        return false;
    }

    std::string block = fname.substr(d0, 15);   // "YYYYMMDD-HHMMSS"

    // --- Verifica del formato: 8 cifre, '-', 6 cifre ---
    bool fmt_ok = true;
    for (int i = 0; i < 15 && fmt_ok; i++) {
        if (i == 8) { if (block[i] != '-')            fmt_ok = false; }
        else        { if (!isdigit((unsigned char)block[i])) fmt_ok = false; }
    }
    if (!fmt_ok) {
        std::cerr << "[ERRORE] ParseFifoStartTime: formato data/ora non valido: '"
                  << block << "'" << std::endl;
        return false;
    }

    // --- Estrazione dei campi ---
    int YYYY = std::stoi(block.substr(0, 4));
    int MM   = std::stoi(block.substr(4, 2));
    int DD   = std::stoi(block.substr(6, 2));
    int hh   = std::stoi(block.substr(9, 2));
    int mm   = std::stoi(block.substr(11, 2));
    int ss   = std::stoi(block.substr(13, 2));

    fd.start_day_index = DaysFromCivil(YYYY, (unsigned)MM, (unsigned)DD);
    fd.start_s_day     = (double)(hh * 3600 + mm * 60 + ss);

    std::cout << "[FIFO] Avvio acquisizione FPGA: "
              << YYYY << "-" << MM << "-" << DD << " "
              << hh << ":" << mm << ":" << ss
              << "  (day_index = " << fd.start_day_index
              << ", " << fd.start_s_day << " s dall'inizio del giorno)"
              << std::endl;
    return true;
}


/// ParseDRSTimestampAbs(): converte il campo <Time> dell'XML DRS
/// ("YYYY/MM/DD HH:MM:SS.sss") in un istante assoluto espresso come
/// (numero di secondi reali) misurato sullo STESSO riferimento del FPGA.
///
/// Per gestire correttamente acquisizioni che attraversano la mezzanotte e/o
/// dataset multipli su piu' giorni, costruiamo un tempo assoluto continuo:
///     t_abs_s = (day_index - start_day_index) * 86400 + (s nel giorno) - tz_offset
/// dove start_day_index e' il giorno di avvio del FIFO accoppiato. In questo modo
/// il risultato e' direttamente confrontabile con t_start_abs_s delle coppie FIFO
/// (che usano lo stesso start_day_index come origine dei giorni).
///
/// drs_abs_s in uscita: secondi dall'inizio del giorno di avvio FPGA, gia'
/// corretti per il fuso orario (riferimento FPGA). Ritorna false se il parsing
/// fallisce.
static bool ParseDRSTimestampAbs(const std::string &ts,
                                 long long start_day_index,
                                 double &drs_abs_s) {
    // Formato: "2026/05/21 18:21:16.512"
    //           0123456789...
    int YYYY, MM, DD, hh, mm;
    double sec = 0.0;

    // Leggiamo data e ora; i secondi (con eventuale parte decimale) come double.
    int matched = sscanf(ts.c_str(), "%d/%d/%d %d:%d:%lf",
                         &YYYY, &MM, &DD, &hh, &mm, &sec);
    if (matched < 6) {
        std::cerr << "[ERRORE] ParseDRSTimestampAbs: parsing fallito per: '"
                  << ts << "'" << std::endl;
        return false;
    }

    long long day_index = DaysFromCivil(YYYY, (unsigned)MM, (unsigned)DD);
    double s_in_day = (double)(hh * 3600 + mm * 60) + sec;

    // Tempo assoluto sul riferimento FPGA, origine = inizio giorno di avvio FIFO.
    drs_abs_s = (double)(day_index - start_day_index) * 86400.0
              + s_in_day
              - (double)FIFO_TZ_OFFSET_S;

    return true;
}


/// LoadFifoAndBuildPairs(): legge un file FIFO, ricostruisce i timestamp
/// assoluti in cicli di clock (interpretazione A: nessuna perdita di tempo) e
/// costruisce la lista delle COPPIE DI DECADIMENTO valide (START seguito da STOP
/// entro la finestra di GATE). Popola fd.pairs (ordinate per ts_start crescente).
///
/// ALGORITMO
///   1. Apre il file. Le righe prima del PRIMO reset appartengono al buffer 0
///      parziale: NON costruiamo coppie con esse (sarebbero ambigue), ma il
///      conteggio dei buffer parte comunque da 0 cosi' i timestamp restano
///      ancorati all'orario del nome file.
///   2. Ad ogni reset (col1 == 2^31) incrementa il contatore di buffer.
///   3. Per ogni riga dati valida (col1 = 1 o 2) calcola
///         ts_abs = n_buffer_corrente * 2^30 + col2
///      e accumula START e STOP in sequenza temporale.
///   4. Scorre la sequenza: per ogni START cerca il PRIMO STOP successivo con
///         dt = ts_stop - ts_start in [FIFO_DECAY_DT_MIN, FIFO_DECAY_DT_MAX].
///      Se trovato, registra la coppia; l'eventuale STOP "consumato" non viene
///      riusato per START successivi.
///
/// PARAMETRO discard_buffer0: se true (default) si scarta il buffer 0 parziale
/// per la COSTRUZIONE delle coppie. Restituisce il numero di coppie trovate,
/// oppure -1 in caso di errore di apertura/contenuto.
static int LoadFifoAndBuildPairs(const std::string &fifo_path, FifoData &fd,
                                 bool discard_buffer0 = true) {

    std::ifstream fin(fifo_path);
    if (!fin.is_open()) {
        std::cerr << "[ERRORE] LoadFifoAndBuildPairs: impossibile aprire '"
                  << fifo_path << "'" << std::endl;
        return -1;
    }
    std::cout << "[FIFO] Lettura ed elaborazione: " << fifo_path << " ..."
              << std::flush;

    fd.n_resets = 0;
    fd.n_lines  = 0;
    fd.pairs.clear();

    // Sequenza temporale di tutti gli eventi (START/STOP) con ts assoluto.
    // chan: FIFO_CH_START o FIFO_CH_STOP ; ts: cicli assoluti dall'avvio.
    std::vector<int>       seq_chan;
    std::vector<long long> seq_ts;
    seq_chan.reserve(1 << 16);
    seq_ts.reserve(1 << 16);

    long long col1, col2;
    long long buf_index    = 0;       // indice del buffer corrente (0-based)
    bool      first_reset  = false;   // true dopo aver visto il primo reset

    while (fin >> col1 >> col2) {
        fd.n_lines++;

        // ---- Riga di RESET: col1 == 2^31 ; col2 == 2^31 + k ----
        if (col1 == FIFO_RESET_SENTINEL) {
            // Il reset segna la FINE del buffer corrente e l'inizio del successivo.
            buf_index++;
            fd.n_resets++;
            first_reset = true;
            continue;
        }

        // ---- Righe del buffer 0 parziale (prima del primo reset) ----
        // Nell'interpretazione A il loro tempo NON e' perso (buf_index resta 0,
        // i ts sono gia' corretti), ma per la costruzione delle coppie le
        // saltiamo se discard_buffer0 e' attivo.
        if (discard_buffer0 && !first_reset) continue;

        // ---- Righe dati normali: col1 = 1 (START) o 2 (STOP) ----
        if (col1 != FIFO_CH_START && col1 != FIFO_CH_STOP) continue; // riga anomala

        // Timestamp assoluto in cicli di clock dall'avvio acquisizione.
        long long ts_abs = buf_index * FIFO_BUF_SIZE_CLK + col2;

        seq_chan.push_back((int)col1);
        seq_ts.push_back(ts_abs);
    }
    fin.close();

    // ----------------------------------------------------------------------
    //  COSTRUZIONE DELLE COPPIE DI DECADIMENTO START -> STOP
    // ----------------------------------------------------------------------
    // La sequenza e' gia' in ordine temporale (il file e' scritto in ordine).
    // Per ogni START cerco in avanti il primo STOP con dt nella finestra GATE.
    // Uno STOP gia' usato da una coppia non viene riassegnato a START successivi.
    std::vector<char> stop_used(seq_ts.size(), 0);

    for (size_t i = 0; i < seq_ts.size(); i++) {
        if (seq_chan[i] != FIFO_CH_START) continue;

        long long ts_start = seq_ts[i];
        // Cerco il primo STOP valido dopo lo START.
        for (size_t j = i + 1; j < seq_ts.size(); j++) {
            long long dt = seq_ts[j] - ts_start;
            if (dt > FIFO_DECAY_DT_MAX) break;        // oltre la GATE: nessun STOP utile
            if (seq_chan[j] != FIFO_CH_STOP) continue;
            if (stop_used[j])                continue;
            if (dt < FIFO_DECAY_DT_MIN)      continue; // dentro il dead-time: ignora

            // Coppia valida trovata.
            FifoDecayPair p;
            p.ts_start_clk = ts_start;
            p.ts_stop_clk  = seq_ts[j];
            p.dt_clk       = dt;
            // Istante assoluto dello START [s, riferimento FPGA, origine = inizio
            // giorno di avvio]: orario del nome file + offset fine + ts*5ns.
            p.t_start_abs_s = fd.start_s_day
                            + FIFO_T0_OFFSET_S
                            + (double)ts_start * FIFO_CLK_PERIOD_S;
            p.used = false;

            fd.pairs.push_back(p);
            stop_used[j] = 1;
            break;   // un solo STOP per START
        }
    }

    // Le coppie sono gia' in ordine di ts_start (gli START si incontrano in
    // ordine crescente nella sequenza), ma ordiniamo per sicurezza: il matching
    // col DRS usa lower_bound e richiede l'ordinamento.
    std::sort(fd.pairs.begin(), fd.pairs.end(),
              [](const FifoDecayPair &a, const FifoDecayPair &b) {
                  return a.ts_start_clk < b.ts_start_clk;
              });

    std::cout << " " << fd.n_lines << " righe, " << fd.n_resets << " reset, "
              << seq_ts.size() << " eventi (post buffer0), "
              << fd.pairs.size() << " coppie di decadimento." << std::endl;

    return (int)fd.pairs.size();
}


/// MatchDRStoFifo(): cerca, fra le coppie di decadimento del FIFO, una il cui
/// START coincida col timestamp di acquisizione dell'evento DRS entro la finestra
/// FIFO_MATCH_WINDOW_S. Implementa il consumo 1-a-1: una coppia gia' usata non
/// puo' validare un secondo evento DRS.
///
/// ARGOMENTI
///   fd          : dati FIFO del file accoppiato (pairs deve essere ordinato)
///   drs_abs_s   : istante dell'evento DRS [s, riferimento FPGA, stessa origine]
///   matched_dt  : (out) |t_DRS - t_START| della coppia scelta [s] (diagnostica)
///
/// RITORNA l'indice della coppia associata in fd.pairs, oppure -1 se nessuna
/// coppia non ancora usata cade nella finestra.
///
/// Strategia: con lower_bound trovo la prima coppia con t_start >= (drs - win),
/// poi scorro finche' t_start <= (drs + win) e scelgo la PIU' VICINA non usata.
static int MatchDRStoFifo(FifoData &fd, double drs_abs_s, double &matched_dt) {

    matched_dt = -1.0;
    if (fd.pairs.empty()) return -1;

    const double win = FIFO_MATCH_WINDOW_S;
    const double lo  = drs_abs_s - win;
    const double hi  = drs_abs_s + win;

    // lower_bound sul campo t_start_abs_s (le coppie sono ordinate per ts_start,
    // che e' monotono in t_start_abs_s essendo una trasformazione affine
    // crescente). Cerco la prima coppia con t_start_abs_s >= lo.
    auto it = std::lower_bound(fd.pairs.begin(), fd.pairs.end(), lo,
                  [](const FifoDecayPair &p, double v) {
                      return p.t_start_abs_s < v;
                  });

    int    best_idx  = -1;
    double best_dist = win;   // accetto solo distanze <= win

    for (; it != fd.pairs.end() && it->t_start_abs_s <= hi; ++it) {
        if (it->used) continue;
        double dist = std::fabs(it->t_start_abs_s - drs_abs_s);
        if (dist <= best_dist) {
            best_dist = dist;
            best_idx  = (int)(it - fd.pairs.begin());
        }
    }

    if (best_idx >= 0) {
        fd.pairs[best_idx].used = true;   // consumo 1-a-1
        matched_dt = best_dist;
    }
    return best_idx;
}


/// SetFifoMatchWindow(): imposta a runtime la semi-finestra di matching
/// DRS<->FIFO [secondi], senza ricompilare. Utile per lo studio di tuning.
void SetFifoMatchWindow(double win_s) {
    FIFO_MATCH_WINDOW_S = win_s;
    std::cout << "[FIFO] Finestra di matching DRS<->FIFO impostata a +-"
              << FIFO_MATCH_WINDOW_S * 1e3 << " ms" << std::endl;
}

/// SetFifoTiming(): imposta a runtime l'offset di fuso orario [s], l'offset fine
/// di origine clock [s] e la finestra della coppia di decadimento [cicli].
/// Tutti i parametri sono opzionali: passare un valore negativo a dt_min/dt_max
/// li lascia invariati.
void SetFifoTiming(int tz_offset_s = -1,
                   double t0_offset_s = 0.0,
                   long long decay_dt_min = -1,
                   long long decay_dt_max = -1) {
    if (tz_offset_s >= 0)   FIFO_TZ_OFFSET_S = tz_offset_s;
    FIFO_T0_OFFSET_S = t0_offset_s;
    if (decay_dt_min >= 0)  FIFO_DECAY_DT_MIN = decay_dt_min;
    if (decay_dt_max >= 0)  FIFO_DECAY_DT_MAX = decay_dt_max;
    std::cout << "[FIFO] tz_offset = " << FIFO_TZ_OFFSET_S << " s, "
              << "t0_offset = " << FIFO_T0_OFFSET_S << " s, "
              << "GATE = [" << FIFO_DECAY_DT_MIN << ", " << FIFO_DECAY_DT_MAX
              << "] clock = [" << FIFO_DECAY_DT_MIN * FIFO_CLK_PERIOD_S * 1e6
              << ", " << FIFO_DECAY_DT_MAX * FIFO_CLK_PERIOD_S * 1e6
              << "] us" << std::endl;
}



/// ComputeCFDTimeAtFraction(): ricalcola il tempo CFD a una frazione arbitraria.
/// Usa baseline, ampiezza e waveform gia' disponibili in ChannelData; la logica
/// di crossing e interpolazione e' la stessa di AnalyzeChannel(), ma il valore
/// della frazione CFD e' passato come argomento.
double ComputeCFDTimeAtFraction(const ChannelData& cd, double frac, bool& ok) {
    ok = false;
    if (!cd.has_pulse) return -999.0;
    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return -999.0;

    const float* t = cd.time;
    const float* v = cd.voltage;
    double bl  = cd.baseline;
    double amp = cd.amplitude;
    if (amp < NOISE_THRESH) return -999.0;

    int imin = 0;
    double vmin = v[0];
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
    }

    double v_thr = bl - frac * amp;
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


// ==========================================================================
//  SEZIONE 6: FUNZIONE C(x)
// ==========================================================================

/// GetC(): restituisce la costante di offset C per la posizione x sulla barra.
///
/// GERARCHIA DEI MODELLI:
///   1) Quadratico (gC_poly_loaded): C(x) = p0 + p1*x + p2*x^2
///      Modello PRINCIPALE, cattura la curvatura parabolica.
///   2) Lineare (gC_lin_loaded):     C(x) = p0 + p1*x
///      Fallback se il quadratico non e' disponibile.
///   3) Interpolazione dei punti C(x_k) di CAL_C_POINTS (con clamp ai bordi).
///      Ultimo fallback se nessun modello analitico e' caricato.
double GetC(double x) {

    // Via principale: modello quadratico caricato dalla calibrazione.
    if (gC_poly_loaded) {
        return gC_poly_p0 + gC_poly_p1 * x + gC_poly_p2 * x * x;
    }

    // Via secondaria: modello lineare.
    if (gC_lin_loaded) {
        return gC_lin_p0 + gC_lin_p1 * x;
    }

    // Via di backup: interpolazione lineare tra i punti C(x_k).
    int n = (int)CAL_C_POINTS.size();
    if (n == 0) {
        std::cerr << "[ERRORE] Nessun punto di calibrazione C definito!" << std::endl;
        return 0.0;
    }
    if (n == 1) return CAL_C_POINTS[0].C;

    // Clamp ai bordi: estrapolazione costante.
    if (x <= CAL_C_POINTS[0].x)     return CAL_C_POINTS[0].C;
    if (x >= CAL_C_POINTS[n - 1].x) return CAL_C_POINTS[n - 1].C;

    // Trova l'intervallo [k, k+1] che contiene x (CAL_C_POINTS ordinato per x).
    int k = 0;
    for (int i = 0; i < n - 1; i++) {
        if (x >= CAL_C_POINTS[i].x && x <= CAL_C_POINTS[i + 1].x) { k = i; break; }
    }
    double x_lo = CAL_C_POINTS[k].x,   x_hi = CAL_C_POINTS[k + 1].x;
    double C_lo = CAL_C_POINTS[k].C,   C_hi = CAL_C_POINTS[k + 1].C;
    double dx = x_hi - x_lo;
    if (fabs(dx) < 1e-9) return 0.5 * (C_lo + C_hi);
    double frac = (x - x_lo) / dx;
    return C_lo + frac * (C_hi - C_lo);
}

/// GetC_const(): restituisce il valore COSTANTE di C (media pesata).
/// Usato dalla pipeline di confronto "Cconst" per quantificare l'impatto
/// della modellizzazione di C(x) sulla distribuzione di beta.
double GetC_const(double /*x*/) {
    return gC_const_val;
}

/// GetC_const_lin(): C COSTANTE nel range lineare ristretto (pipeline 2b).
/// Restituisce il valore costante stimato dalla calibrazione sui soli punti
/// x in [-84,+70] cm con il metodo corrected. Indipendente da x per costruzione.
/// Se non caricato (branch assente), ripiega sul C costante globale GetC_const().
double GetC_const_lin(double /*x*/) {
    if (gC_const_lin_loaded) return gC_const_lin_val;
    return gC_const_val;   // fallback: C costante globale
}

/// GetC_lin_corr(): C LINEARE nel range ristretto (pipeline 2a, non ancora
/// attiva). Pronta per un futuro uso: C(x) = p0 + p1*x con i parametri ristretti.
double GetC_lin_corr(double x) {
    if (gC_lin_corr_loaded) return gC_lin_corr_p0 + gC_lin_corr_p1 * x;
    return GetC_const_lin(x);   // fallback ragionevole nel range
}

/// GetC_corr(): modello C(x) coerente con la calibrazione rise-time corrected.
/// Se il file di calibrazione non contiene C_poly_corr_*, usa GetC() come
/// fallback esplicito per mantenere la macro eseguibile con file legacy.
double GetC_corr(double x) {
    if (gC_corr_poly_loaded) {
        return gC_corr_poly_p0 + gC_corr_poly_p1 * x + gC_corr_poly_p2 * x * x;
    }
    return GetC(x);
}


// ==========================================================================
//  SEZIONE 7: CARICAMENTO DELLA CALIBRAZIONE DA FILE ROOT
// ==========================================================================
//
//  LoadCalibrationFromFile(): legge il file ROOT prodotto da
//  TOF_Calibration_v13.cpp e popola TUTTE le globali di calibrazione.
//
//  COSA LEGGE (contratto con la v13):
//    Dal TTree "summary":
//      - 'x'             -> posizione di calibrazione x_k
//      - 'C_hybrid_mu'   -> valore di C misurato col metodo hybrid_tot
//      Le coppie (x_k, C_k) popolano CAL_C_POINTS, usato come BACKUP da GetC().
//    Dal TTree "fit_params" (un solo entry):
//      - m, m_err, q, q_err, v_eff, v_eff_err, chi2ndf -> retta Dt12(x)
//      - C_lin_p0/p1 (+ err, chi2ndf)                  -> modello lineare C(x)
//      - TOT_thr_ref, TOT_calib_idx, poly_degree       -> metadati A(TOT)
//      - Amax_TOT_PMT1/2/3                             -> cap di estrapolazione
//      - poly_TOT_PMT1/2/3[6] (+ err, chi2ndf)         -> polinomio cubico A(TOT)
//
//  COMPORTAMENTO IN CASO DI BRANCH MANCANTI:
//    In v13 il polinomio cubico A(TOT) e' l'UNICO modello di ricostruzione dei
//    clippati: non esiste piu' un modello esponenziale o slew-rate di fallback.
//    Pertanto, se i branch poly_TOT_* NON sono presenti nel file, la funzione
//    emette un ERRORE esplicito e lascia gTOT_poly_loaded[] a false: il
//    recupero dei clippati sara' disabilitato e gli eventi con un canale
//    clippato verranno scartati dalla pipeline hybrid_tot. NON c'e' un
//    fallback silenzioso a un altro modello.
//
//  RITORNA: true se il file e' stato aperto e il TTree fit_params letto;
//           false se il file non e' apribile o fit_params manca.

bool LoadCalibrationFromFile(const char* cal_filename) {

    TFile *fcal = TFile::Open(cal_filename, "READ");
    if (!fcal || fcal->IsZombie()) {
        std::cerr << "[ERRORE] Non riesco ad aprire il file di calibrazione: "
                  << cal_filename << std::endl;
        std::cerr << "         Verranno usati i parametri hardcoded "
                  << "(CAL_M, CAL_Q, CAL_C_POINTS) e il recupero clippati "
                  << "sara' DISABILITATO." << std::endl;
        return false;
    }

    // ----------------------------------------------------------------------
    //  (1) Punti C(x_k) dal TTree "summary"  (branch espliciti hybrid_tot)
    // ----------------------------------------------------------------------
    TTree *summary = (TTree*)fcal->Get("summary");
    if (summary) {
        // I branch del summary della v13 hanno nomi ESPLICITI per metodo:
        // la posizione e' 'x', il C del metodo fisico e' 'C_hybrid_mu'.
        if (summary->GetBranch("x") && summary->GetBranch("C_hybrid_mu")) {
            Float_t s_x = 0.0f, s_C_hy_mu = 0.0f;
            summary->SetBranchAddress("x",           &s_x);
            summary->SetBranchAddress("C_hybrid_mu", &s_C_hy_mu);

            CAL_C_POINTS.clear();
            int nentries = (int)summary->GetEntries();
            for (int i = 0; i < nentries; i++) {
                summary->GetEntry(i);
                CAL_C_POINTS.push_back({ (double)s_x, (double)s_C_hy_mu });
            }
            // Ordina per x crescente: precondizione di GetC() (via di backup).
            std::sort(CAL_C_POINTS.begin(), CAL_C_POINTS.end(),
                      [](const CPoint &a, const CPoint &b) { return a.x < b.x; });

            std::cout << "[INFO] Caricati " << CAL_C_POINTS.size()
                      << " punti C(x) (hybrid_tot) dal TTree summary." << std::endl;
            for (auto &p : CAL_C_POINTS)
                std::cout << "       x = " << p.x << " cm, C = " << p.C << " ns"
                          << std::endl;
        } else {
            std::cerr << "[WARNING] Il TTree 'summary' non ha i branch attesi "
                      << "'x' / 'C_hybrid_mu'. Punti C(x) di backup invariati."
                      << std::endl;
        }
    } else {
        std::cerr << "[WARNING] TTree 'summary' non trovato in " << cal_filename
                  << ". Punti C(x) di backup invariati." << std::endl;
    }

    // ----------------------------------------------------------------------
    //  (2) Parametri scalari e polinomio A(TOT) dal TTree "fit_params"
    // ----------------------------------------------------------------------
    TTree *fp = (TTree*)fcal->Get("fit_params");
    if (!fp || fp->GetEntries() <= 0) {
        std::cerr << "[ERRORE] TTree 'fit_params' assente o vuoto in "
                  << cal_filename << "." << std::endl;
        std::cerr << "         Impossibile caricare retta, C(x) e polinomio "
                  << "A(TOT): si usano i parametri hardcoded e il recupero "
                  << "clippati resta DISABILITATO." << std::endl;
        fcal->Close();
        return false;
    }

    // ROOT mantiene gli indirizzi passati a SetBranchAddress() finche' non
    // vengono resettati. Qui sotto leggiamo fit_params in piu' blocchi con
    // variabili locali: senza ResetBranchAddresses() un GetEntry successivo
    // potrebbe scrivere in variabili gia' uscite dallo scope, corrompendo i
    // parametri di calibrazione e quindi il TOF ricostruito.
    auto ResetFitParamBranches = [&]() {
        fp->ResetBranchAddresses();
    };

    // --- (2a) Retta di calibrazione Dt12 = m*x + q ---
    {
        Float_t fp_m = 0, fp_m_err = 0, fp_q = 0, fp_q_err = 0;
        Float_t fp_v_eff = 0, fp_v_eff_err = 0, fp_chi2ndf = 0;
        fp->SetBranchAddress("m",         &fp_m);
        fp->SetBranchAddress("m_err",     &fp_m_err);
        fp->SetBranchAddress("q",         &fp_q);
        fp->SetBranchAddress("q_err",     &fp_q_err);
        fp->SetBranchAddress("v_eff",     &fp_v_eff);
        fp->SetBranchAddress("v_eff_err", &fp_v_eff_err);
        fp->SetBranchAddress("chi2ndf",   &fp_chi2ndf);
        fp->GetEntry(0);

        CAL_M = (double)fp_m;
        CAL_Q = (double)fp_q;

        std::cout << "[INFO] Retta di calibrazione caricata da fit_params:"
                  << std::endl;
        std::cout << "       m     = " << CAL_M << " +/- " << fp_m_err
                  << " ns/cm   (chi2/ndf = " << fp_chi2ndf << ")" << std::endl;
        std::cout << "       q     = " << CAL_Q << " +/- " << fp_q_err
                  << " ns" << std::endl;
        std::cout << "       v_eff = " << fp_v_eff << " +/- " << fp_v_eff_err
                  << " cm/ns" << std::endl;
        ResetFitParamBranches();
    }

    // --- (2a') Retta corretta per rise time: Dt12_corr = m_corr*x + q_corr ---
    if (fp->GetBranch("m_corr") && fp->GetBranch("q_corr")) {
        Float_t fp_m = 0, fp_m_err = 0, fp_q = 0, fp_q_err = 0;
        Float_t fp_v_eff = 0, fp_v_eff_err = 0, fp_chi2ndf = 0;
        fp->SetBranchAddress("m_corr",         &fp_m);
        fp->SetBranchAddress("m_corr_err",     &fp_m_err);
        fp->SetBranchAddress("q_corr",         &fp_q);
        fp->SetBranchAddress("q_corr_err",     &fp_q_err);
        if (fp->GetBranch("v_eff_corr"))
            fp->SetBranchAddress("v_eff_corr", &fp_v_eff);
        if (fp->GetBranch("v_eff_corr_err"))
            fp->SetBranchAddress("v_eff_corr_err", &fp_v_eff_err);
        if (fp->GetBranch("chi2ndf_corr"))
            fp->SetBranchAddress("chi2ndf_corr", &fp_chi2ndf);
        fp->GetEntry(0);

        CAL_M_CORR = (double)fp_m;
        CAL_Q_CORR = (double)fp_q;
        gRiseCorrCal_loaded = true;

        std::cout << "[INFO] Retta rise-time corrected caricata:" << std::endl;
        std::cout << "       m_corr = " << CAL_M_CORR << " +/- " << fp_m_err
                  << " ns/cm   (chi2/ndf = " << fp_chi2ndf << ")" << std::endl;
        std::cout << "       q_corr = " << CAL_Q_CORR << " +/- " << fp_q_err
                  << " ns" << std::endl;
        if (fp_v_eff != 0.0f)
            std::cout << "       v_eff_corr = " << fp_v_eff << " +/- "
                      << fp_v_eff_err << " cm/ns" << std::endl;
        ResetFitParamBranches();
} else {
        gRiseCorrCal_loaded = false;
        std::cout << "[INFO] Branch m_corr/q_corr assenti: la pipeline "
                  << "rise_corr usera' m/q classici come fallback." << std::endl;
    }

    // --- (2a'') Retta RISTRETTA al range lineare: Dt12_corr_lin = m_corr_lin*x + q_corr_lin ---
    //  Usata dalla pipeline Lin_Range (ricostruzione indipendente). Caricata
    //  solo se la calibrazione ha salvato i branch ristretti (linrange_valid).
    if (fp->GetBranch("m_corr_lin") && fp->GetBranch("q_corr_lin")) {
        Float_t fp_ml = 0, fp_ml_err = 0, fp_ql = 0, fp_ql_err = 0;
        Float_t fp_vl = 0, fp_vl_err = 0, fp_chil = 0;
        Int_t   fp_lrvalid = 0;
        fp->SetBranchAddress("m_corr_lin",         &fp_ml);
        fp->SetBranchAddress("q_corr_lin",         &fp_ql);
        if (fp->GetBranch("m_corr_lin_err"))
            fp->SetBranchAddress("m_corr_lin_err", &fp_ml_err);
        if (fp->GetBranch("q_corr_lin_err"))
            fp->SetBranchAddress("q_corr_lin_err", &fp_ql_err);
        if (fp->GetBranch("v_eff_corr_lin"))
            fp->SetBranchAddress("v_eff_corr_lin", &fp_vl);
        if (fp->GetBranch("v_eff_corr_lin_err"))
            fp->SetBranchAddress("v_eff_corr_lin_err", &fp_vl_err);
        if (fp->GetBranch("chi2ndf_corr_lin"))
            fp->SetBranchAddress("chi2ndf_corr_lin", &fp_chil);
        if (fp->GetBranch("linrange_valid"))
            fp->SetBranchAddress("linrange_valid", &fp_lrvalid);
        fp->GetEntry(0);

        CAL_M_LIN = (double)fp_ml;
        CAL_Q_LIN = (double)fp_ql;
        // Considera la calibrazione ristretta valida solo se il flag e' 1
        // (oppure assente, nel qual caso ci fidiamo della presenza dei branch).
        gLinRangeCal_loaded = (fp->GetBranch("linrange_valid"))
                            ? (fp_lrvalid == 1) : true;

        if (gLinRangeCal_loaded) {
            std::cout << "[INFO] Retta RISTRETTA (range lineare) caricata:" << std::endl;
            std::cout << "       m_corr_lin = " << CAL_M_LIN << " +/- " << fp_ml_err
                      << " ns/cm   (chi2/ndf = " << fp_chil << ")" << std::endl;
            std::cout << "       q_corr_lin = " << CAL_Q_LIN << " +/- " << fp_ql_err
                      << " ns" << std::endl;
            if (fp_vl != 0.0f)
                std::cout << "       v_eff_corr_lin = " << fp_vl << " +/- "
                          << fp_vl_err << " cm/ns" << std::endl;
        } else {
            std::cout << "[INFO] linrange_valid = 0: calibrazione ristretta "
                << "non valida, Lin_Range usa il fallback." << std::endl;        }
        ResetFitParamBranches();
    } else {
        gLinRangeCal_loaded = false;
        std::cout << "[INFO] Branch m_corr_lin/q_corr_lin assenti: la pipeline "
                  << "Lin_Range usera' la retta corrected globale come fallback."
                  << std::endl;
    }

    // --- (2a''') C costante ristretto (pipeline 2b) e C lineare ristretto (2a) ---
    if (fp->GetBranch("C_const_corr_lin")) {
        Float_t fp_Ccl = 0, fp_Ccl_err = 0;
        fp->SetBranchAddress("C_const_corr_lin", &fp_Ccl);
        if (fp->GetBranch("C_const_corr_lin_err"))
            fp->SetBranchAddress("C_const_corr_lin_err", &fp_Ccl_err);
        fp->GetEntry(0);
        gC_const_lin_val    = (double)fp_Ccl;
        gC_const_lin_err    = (double)fp_Ccl_err;
        gC_const_lin_loaded = true;
        std::cout << "[INFO] C costante ristretto (Lin_Range, 2b): <C>_lin = "
                  << gC_const_lin_val << " +/- " << gC_const_lin_err << " ns"
                  << std::endl;
        ResetFitParamBranches();
    } else {
        gC_const_lin_loaded = false;
        std::cout << "[INFO] Branch C_const_corr_lin assente: Lin_Range usera' "
                  << "il C costante globale come fallback." << std::endl;
    }

    if (fp->GetBranch("C_corr_lin_p0") && fp->GetBranch("C_corr_lin_p1")) {
        Float_t fp_p0 = 0, fp_p1 = 0;
        fp->SetBranchAddress("C_corr_lin_p0", &fp_p0);
        fp->SetBranchAddress("C_corr_lin_p1", &fp_p1);
        fp->GetEntry(0);
        gC_lin_corr_p0     = (double)fp_p0;
        gC_lin_corr_p1     = (double)fp_p1;
        gC_lin_corr_loaded = true;
        std::cout << "[INFO] C lineare ristretto (Lin_Range, 2a pronto): C(x) = "
                  << gC_lin_corr_p0 << " + (" << gC_lin_corr_p1 << ")*x  ns"
                  << std::endl;
        ResetFitParamBranches();
    } else {
        gC_lin_corr_loaded = false;
    }

    // --- (2b) Modello lineare C(x) = C_lin_p0 + C_lin_p1*x ---
    if (fp->GetBranch("C_lin_p0") && fp->GetBranch("C_lin_p1")) {
        Float_t fp_Cp0 = 0, fp_Cp0_err = 0, fp_Cp1 = 0, fp_Cp1_err = 0;
        Float_t fp_Cchi2 = 0;
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
        std::cout << "       C(x) = " << gC_lin_p0 << " + (" << gC_lin_p1
                  << ") * x  ns   (chi2/ndf = " << fp_Cchi2 << ")" << std::endl;
        std::cout << "       -> C usata da GetC(): MODELLO LINEARE" << std::endl;
        ResetFitParamBranches();
    } else {
        gC_lin_loaded = false;
        std::cerr << "[WARNING] Branch C_lin_* assenti in fit_params. "
                  << "GetC() user&agrave; l'INTERPOLAZIONE dei punti C(x_k)."
                  << std::endl;
    }

    // --- (2b') Modello quadratico C(x) = C_poly_p0 + C_poly_p1*x + C_poly_p2*x^2 ---
    // MODELLO PRINCIPALE per C(x): priorita' > lineare nella gerarchia di GetC().
    if (fp->GetBranch("C_poly_p0") && fp->GetBranch("C_poly_p1")
                                    && fp->GetBranch("C_poly_p2")) {
        Float_t fp_Cp0 = 0, fp_Cp0_err = 0, fp_Cp1 = 0, fp_Cp1_err = 0;
        Float_t fp_Cp2 = 0, fp_Cp2_err = 0, fp_Cchi2 = 0;
        fp->SetBranchAddress("C_poly_p0",      &fp_Cp0);
        fp->SetBranchAddress("C_poly_p0_err",  &fp_Cp0_err);
        fp->SetBranchAddress("C_poly_p1",      &fp_Cp1);
        fp->SetBranchAddress("C_poly_p1_err",  &fp_Cp1_err);
        fp->SetBranchAddress("C_poly_p2",      &fp_Cp2);
        fp->SetBranchAddress("C_poly_p2_err",  &fp_Cp2_err);
        fp->SetBranchAddress("C_poly_chi2ndf", &fp_Cchi2);
        fp->GetEntry(0);

        gC_poly_p0     = (double)fp_Cp0;
        gC_poly_p0_err = (double)fp_Cp0_err;
        gC_poly_p1     = (double)fp_Cp1;
        gC_poly_p1_err = (double)fp_Cp1_err;
        gC_poly_p2     = (double)fp_Cp2;
        gC_poly_p2_err = (double)fp_Cp2_err;
        gC_poly_loaded = true;

        std::cout << "[INFO] Modello quadratico C(x) caricato (PRINCIPALE):"
                  << std::endl;
        std::cout << "       C(x) = " << gC_poly_p0 << " + (" << gC_poly_p1
                  << ")*x + (" << gC_poly_p2 << ")*x^2  ns"
                  << "   (chi2/ndf = " << fp_Cchi2 << ")" << std::endl;
        std::cout << "       -> C usata da GetC(): MODELLO QUADRATICO" << std::endl;
        ResetFitParamBranches();
    } else {
        gC_poly_loaded = false;
        std::cout << "[INFO] Branch C_poly_* assenti in fit_params. "
                  << "GetC() usera' il modello lineare o l'interpolazione."
                  << std::endl;
    }

    // --- (2b'') Modello quadratico C_corr(x) per tempi rise-time corrected ---
    if (fp->GetBranch("C_poly_corr_p0") && fp->GetBranch("C_poly_corr_p1")
                                        && fp->GetBranch("C_poly_corr_p2")) {
        Float_t fp_Cp0 = 0, fp_Cp1 = 0, fp_Cp2 = 0, fp_Cchi2 = 0;
        fp->SetBranchAddress("C_poly_corr_p0", &fp_Cp0);
        fp->SetBranchAddress("C_poly_corr_p1", &fp_Cp1);
        fp->SetBranchAddress("C_poly_corr_p2", &fp_Cp2);
        if (fp->GetBranch("C_poly_corr_chi2ndf"))
            fp->SetBranchAddress("C_poly_corr_chi2ndf", &fp_Cchi2);
        fp->GetEntry(0);

        gC_corr_poly_p0     = (double)fp_Cp0;
        gC_corr_poly_p1     = (double)fp_Cp1;
        gC_corr_poly_p2     = (double)fp_Cp2;
        gC_corr_poly_loaded = true;

        std::cout << "[INFO] Modello C_corr(x) caricato (rise-time corrected):"
                  << std::endl;
        std::cout << "       C_corr(x) = " << gC_corr_poly_p0 << " + ("
                  << gC_corr_poly_p1 << ")*x + (" << gC_corr_poly_p2
                  << ")*x^2  ns   (chi2/ndf = " << fp_Cchi2 << ")"
                  << std::endl;
        ResetFitParamBranches();
    } else {
        gC_corr_poly_loaded = false;
        std::cout << "[INFO] Branch C_poly_corr_* assenti: la pipeline "
                  << "rise_corr usera' C(x) classico come fallback." << std::endl;
    }

    // --- (2b'') Modello costante C = <C> (media pesata) ---
    if (fp->GetBranch("C_const_val")) {
        Float_t fp_Cval = 0, fp_Cerr = 0, fp_Cchi2 = 0;
        fp->SetBranchAddress("C_const_val",     &fp_Cval);
        fp->SetBranchAddress("C_const_err",     &fp_Cerr);
        fp->SetBranchAddress("C_const_chi2ndf", &fp_Cchi2);
        fp->GetEntry(0);

        gC_const_val    = (double)fp_Cval;
        gC_const_err    = (double)fp_Cerr;
        gC_const_loaded = true;

        std::cout << "[INFO] Modello costante C caricato (confronto):"
                  << std::endl;
        std::cout << "       <C> = " << gC_const_val << " +/- " << gC_const_err
                  << " ns   (chi2/ndf = " << fp_Cchi2 << ")" << std::endl;
        ResetFitParamBranches();
    } else {
        gC_const_loaded = false;
        std::cout << "[INFO] Branch C_const_val assente: pipeline Cconst disabilitata."
                  << std::endl;
    }

    // --- (2c) Metadati del modello A(TOT) ---
    {
        Float_t fp_TOT_thr_ref = 0;
        Int_t   fp_TOT_calib_idx = 0, fp_poly_degree = 0;
        if (fp->GetBranch("TOT_thr_ref")) {
            fp->SetBranchAddress("TOT_thr_ref", &fp_TOT_thr_ref);
        }
        if (fp->GetBranch("TOT_calib_idx")) {
            fp->SetBranchAddress("TOT_calib_idx", &fp_TOT_calib_idx);
        }
        if (fp->GetBranch("poly_degree")) {
            fp->SetBranchAddress("poly_degree", &fp_poly_degree);
        }
        fp->GetEntry(0);

        // Aggiorna l'indice della soglia di riferimento per il TOT, se valido.
        if (fp->GetBranch("TOT_calib_idx") &&
            fp_TOT_calib_idx >= 0 && fp_TOT_calib_idx < NTOT_THR) {
            TOT_CALIB_THR_IDX = (int)fp_TOT_calib_idx;
        }
        // Aggiorna il grado del polinomio, se valido.
        if (fp->GetBranch("poly_degree") &&
            fp_poly_degree >= 1 && fp_poly_degree <= TOT_POLY_MAX_DEG) {
            gTOT_poly_degree = (int)fp_poly_degree;
            TOT_POLY_DEGREE  = (int)fp_poly_degree;
        }

        std::cout << "[INFO] Metadati A(TOT): soglia di riferimento indice "
                  << TOT_CALIB_THR_IDX << " (= " << TOT_THR_MV[TOT_CALIB_THR_IDX]
                  << " mV";
        if (fp->GetBranch("TOT_thr_ref"))
            std::cout << ", file = " << fp_TOT_thr_ref << " mV";
        std::cout << "), grado polinomio = " << gTOT_poly_degree << std::endl;
        ResetFitParamBranches();
    }

    // --- (2d) Coefficienti del polinomio A(TOT) per PMT1, PMT2, PMT3 ---
    // I branch sono array di dimensione fissa TOT_POLY_MAX_DEG+1 = 6.
    // Se mancano, il recupero clippati resta DISABILITATO (nessun fallback).
    {
        const char* poly_names[3] = { "poly_TOT_PMT1",
                                      "poly_TOT_PMT2",
                                      "poly_TOT_PMT3" };
        const char* poly_err_names[3] = { "poly_TOT_err_PMT1",
                                          "poly_TOT_err_PMT2",
                                          "poly_TOT_err_PMT3" };
        const char* poly_chi2_names[3] = { "poly_TOT_chi2ndf_PMT1",
                                           "poly_TOT_chi2ndf_PMT2",
                                           "poly_TOT_chi2ndf_PMT3" };
        const char* amax_names[3] = { "Amax_TOT_PMT1",
                                      "Amax_TOT_PMT2",
                                      "Amax_TOT_PMT3" };

        bool all_poly_present = true;
        for (int k = 0; k < 3; k++) {
            if (!fp->GetBranch(poly_names[k])) all_poly_present = false;
        }

        if (!all_poly_present) {
            std::cerr << "[ERRORE] Branch del polinomio A(TOT) "
                      << "(poly_TOT_PMT1/2/3) ASSENTI in fit_params." << std::endl;
            std::cerr << "         Il file di calibrazione non e' compatibile "
                      << "con TOF_Analysis_v9 (atteso output di "
                      << "TOF_Calibration_v13)." << std::endl;
            std::cerr << "         Recupero dei clippati DISABILITATO: gli "
                      << "eventi con un canale clippato saranno scartati "
                      << "dalla pipeline hybrid_tot." << std::endl;
            for (int k = 0; k < MAX_CHANNELS; k++) gTOT_poly_loaded[k] = false;
        } else {
            // Buffer temporanei per la lettura degli array.
            Float_t buf_poly [TOT_POLY_MAX_DEG + 1];
            Float_t buf_err  [TOT_POLY_MAX_DEG + 1];
            Float_t buf_chi2 = 0;
            Float_t buf_amax = 0;

            for (int k = 0; k < 3; k++) {
                ResetFitParamBranches();
                for (int j = 0; j <= TOT_POLY_MAX_DEG; j++) {
                    buf_poly[j] = 0.0f;
                    buf_err[j]  = 0.0f;
                }
                buf_chi2 = 0.0f;
                buf_amax = 0.0f;

                fp->SetBranchAddress(poly_names[k], buf_poly);
                if (fp->GetBranch(poly_err_names[k]))
                    fp->SetBranchAddress(poly_err_names[k], buf_err);
                if (fp->GetBranch(poly_chi2_names[k]))
                    fp->SetBranchAddress(poly_chi2_names[k], &buf_chi2);
                if (fp->GetBranch(amax_names[k]))
                    fp->SetBranchAddress(amax_names[k], &buf_amax);

                fp->GetEntry(0);

                // Copia nelle globali (double).
                for (int j = 0; j <= TOT_POLY_MAX_DEG; j++) {
                    gTOT_poly_PMT[k][j]     = (double)buf_poly[j];
                    gTOT_poly_err_PMT[k][j] = (double)buf_err[j];
                }
                gTOT_poly_chi2ndf_PMT[k] = (double)buf_chi2;
                gTOT_Amax_PMT[k]         = (double)buf_amax;

                // La calibrazione del canale e' "caricata" solo se ha un Amax
                // sensato (> 0): senza Amax non si puo' applicare il cap di
                // estrapolazione e il polinomio non e' affidabile.
                gTOT_poly_loaded[k] = (gTOT_Amax_PMT[k] > 0.0);
            }
            ResetFitParamBranches();
            // Il quarto canale (PMT4 = veto) non ha calibrazione A(TOT).
            if (MAX_CHANNELS > 3) gTOT_poly_loaded[3] = false;

            std::cout << "[INFO] Polinomio A(TOT) caricato (grado "
                      << gTOT_poly_degree << "):" << std::endl;
            const char* pmt_lbl[3] = { "PMT1", "PMT2", "PMT3" };
            for (int k = 0; k < 3; k++) {
                std::cout << "       " << pmt_lbl[k] << ": ";
                if (gTOT_poly_loaded[k]) {
                    std::cout << "A(TOT) = ";
                    for (int j = 0; j <= gTOT_poly_degree; j++) {
                        std::cout << gTOT_poly_PMT[k][j];
                        if (j > 0) std::cout << "*TOT^" << j;
                        if (j < gTOT_poly_degree) std::cout << " + ";
                    }
                    std::cout << "  (chi2/ndf = " << gTOT_poly_chi2ndf_PMT[k]
                              << ", A_max = " << gTOT_Amax_PMT[k] << " mV)"
                              << std::endl;
                } else {
                    std::cout << "calibrazione NON disponibile "
                              << "(A_max non valido)" << std::endl;
                }
            }
            bool any_loaded = gTOT_poly_loaded[0] || gTOT_poly_loaded[1] ||
                              gTOT_poly_loaded[2];
            std::cout << "       Recupero clippati (TOT polinomiale): "
                      << (any_loaded ? "ABILITATO" : "DISABILITATO")
                      << std::endl;
        }
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
    std::cout << "       v_eff = " << 2.0 / fabs(m) << " cm/ns" << std::endl;
}

/// SetGeometry(): imposta h e x_PMT3 da terminale ROOT senza ricompilare.
void SetGeometry(double h, double x_pmt3 = 0.0) {
    PAR_H = h;
    PAR_X_PMT3 = x_pmt3;
    std::cout << "[INFO] Geometria aggiornata: h = " << h
              << " cm, x_PMT3 = " << x_pmt3 << " cm" << std::endl;
}

/// SetVetoParams(): imposta i parametri del veto offline PMT4 da terminale
/// ROOT senza ricompilare (utile per studi di sensitivita' della
/// distribuzione di beta alla soglia di veto).
void SetVetoParams(double nsigma     = 5.0,
                   double amp_floor  = 5.0,
                   double t_win_lo   = -20.0,
                   double t_win_hi   =  20.0,
                   int    nsamp_min  = 3) {
    VETO_N_SIGMA   = nsigma;
    VETO_AMP_FLOOR = amp_floor;
    VETO_T_WIN_LO  = t_win_lo;
    VETO_T_WIN_HI  = t_win_hi;
    VETO_NSAMP_MIN = nsamp_min;
    std::cout << "[INFO] Parametri veto PMT4 aggiornati:" << std::endl;
    std::cout << "       N_sigma   = " << VETO_N_SIGMA   << std::endl;
    std::cout << "       AMP_FLOOR = " << VETO_AMP_FLOOR << " mV" << std::endl;
    std::cout << "       T_WIN     = [" << VETO_T_WIN_LO << ", "
              << VETO_T_WIN_HI << "] ns" << std::endl;
    std::cout << "       NSAMP_MIN = " << VETO_NSAMP_MIN << std::endl;
}


// ==========================================================================
//  SEZIONE 8: FUNZIONE PRINCIPALE — ANALISI TOF
// ==========================================================================
//
//  TOF_Analysis(): processa uno o piu' file XML di dati TOF (PMT3 in Guida B)
//  e ricostruisce, evento per evento e con DUE pipeline parallele:
//
//    PIPELINE PRINCIPALE  -> hybrid_tot
//      CFD classico sui canali non clippati; ricostruzione TOT polinomiale
//      sui soli canali clippati. Popola i branch SENZA suffisso del TTree
//      tof_data e gli istogrammi principali. E' il risultato fisico ufficiale.
//
//    PIPELINE DIAGNOSTICA -> noclip
//      Usa esclusivamente i canali non clippati: se anche un solo canale e'
//      clippato l'evento viene scartato (nessun recupero). Popola i branch
//      con suffisso "_noclip" e gli istogrammi "_noclip". Serve da riferimento
//      "pulito" per quantificare quanto la ricostruzione TOT sposta i risultati.
//
//  Per gli eventi SENZA canali clippati le due pipeline producono numeri
//  IDENTICI (entrambe usano t_cfd). Le differenze emergono solo dove la
//  pipeline hybrid_tot recupera un canale che la noclip scarterebbe.
//
//  ARGOMENTI:
//    xml_files — percorsi dei file XML separati da virgola ("a.xml,b.xml")
//    outname   — nome del file ROOT di output
//    cal_file  — file ROOT di calibrazione prodotto da TOF_Calibration_v13
//                (se vuoto, si usano i parametri hardcoded e il recupero
//                 clippati resta disabilitato)
//
//  UTILIZZO:
//    root -l 'TOF_Analysis_v9.cpp("run1.xml,run2.xml","out.root","cal_v13.root")'

void TOF_Analysis(const char* xml_files,
                  const char* outname  = "TOF_Analysis_v9_output.root",
                  const char* cal_file = "") {

    gROOT->SetBatch(kTRUE);   // disattiva le finestre grafiche interattive

    std::cout << "=============================================" << std::endl;
    std::cout << "  TOF Analysis v9 (recupero clippati: TOT cubico)" << std::endl;
    std::cout << "=============================================" << std::endl;

    // --- Caricamento della calibrazione (retta, C(x), polinomio A(TOT)) ---
    if (strlen(cal_file) > 0) {
        std::cout << "[INFO] Caricamento calibrazione da: " << cal_file
                  << std::endl;
        LoadCalibrationFromFile(cal_file);
        } else if (!gC_poly_loaded && !gC_lin_loaded) {
        // Il warning viene emesso SOLO se le globali di calibrazione non sono
        // state caricate neppure da una chiamata precedente a
        // LoadCalibrationFromFile(). Se l'utente ha gia' chiamato la funzione
        // manualmente prima di TOF_Analysis(), le globali sono gia' popolate
        // e il warning sarebbe fuorviante.
        std::cerr << "[WARNING] Nessun file di calibrazione fornito e nessuna "
                  << "calibrazione precedentemente caricata: si usano "
                  << "i parametri hardcoded e il recupero clippati e' "
                  << "DISABILITATO." << std::endl;
    }

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "  Geometria:    h = " << PAR_H << " cm, x_PMT3 = "
              << PAR_X_PMT3 << " cm" << std::endl;
    std::cout << "  Retta cal.:   m = " << CAL_M << " ns/cm, q = "
              << CAL_Q << " ns  (v_eff = " << 2.0 / fabs(CAL_M)
              << " cm/ns)" << std::endl;
    std::cout << "  C(x):         "
              << (gC_poly_loaded  ? "modello QUADRATICO caricato (PRINCIPALE)"
               :  gC_lin_loaded   ? "modello lineare caricato (fallback)"
               :                    "interpolazione punti di backup")
              << std::endl;
    if (gC_const_loaded)
        std::cout << "  C costante:   <C> = " << gC_const_val
                  << " ns (pipeline di confronto)" << std::endl;
    std::cout << "=============================================" << std::endl;

    // --- Parsing della lista di file XML (separati da virgola) ---
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
    std::cout << "[INFO] " << file_list.size() << " file da processare."
              << std::endl;

    // --- Setup dello stile grafico ---
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetTitleSize(0.045, "xy");

    // --- Apertura del file di output ROOT ---
    TFile *fout = new TFile(outname, "RECREATE");
    if (!fout->IsOpen()) {
        std::cerr << "[ERRORE] Impossibile creare " << outname << std::endl;
        return;
    }

    // ======================================================================
    //  TTree DI OUTPUT "tof_data" — un entry per evento
    // ======================================================================
    TTree *tree = new TTree("tof_data", "Dati TOF evento per evento");

    // ---- BRANCH PRINCIPALI: pipeline hybrid_tot (senza suffisso) ----
    Float_t t_cfd[3], amp[3];
    tree->Branch("t_cfd", t_cfd, "t_cfd[3]/F");   // tempi usati (CFD o ricostr.) [ns]
    tree->Branch("amp",   amp,   "amp[3]/F");     // ampiezze (vere o ricostr.)  [mV]

    Float_t dt12, x_imp, T_meas, C_used, tof;
    Float_t path_len, theta, vel, beta, inv_vel;
    Int_t   good;
    tree->Branch("dt12",     &dt12,     "dt12/F");      // Dt12 = t1 - t2 [ns]
    tree->Branch("x_imp",    &x_imp,    "x_imp/F");     // posizione di impatto [cm]
    tree->Branch("T_meas",   &T_meas,   "T_meas/F");    // t3 - (t1+t2)/2 [ns]
    tree->Branch("C_used",   &C_used,   "C_used/F");    // C(x_imp) applicato [ns]
    tree->Branch("tof",      &tof,      "tof/F");       // TOF = T_meas - C [ns]
    tree->Branch("path_len", &path_len, "path_len/F");  // distanza percorsa [cm]
    tree->Branch("theta",    &theta,    "theta/F");     // angolo di incidenza [rad]
    tree->Branch("vel",      &vel,      "vel/F");       // velocita' [cm/ns]
    tree->Branch("beta",     &beta,     "beta/F");      // beta = v/c
    tree->Branch("inv_vel",  &inv_vel,  "inv_vel/F");   // 1/v [ns/cm]
    tree->Branch("good",     &good,     "good/I");      // 1 = evento buono (hybrid)

    // ---- BRANCH PARALLELI: pipeline noclip (suffisso "_noclip") ----
    Float_t dt12_noclip, x_imp_noclip, T_meas_noclip, C_used_noclip, tof_noclip;
    Float_t path_len_noclip, theta_noclip, vel_noclip, beta_noclip, inv_vel_noclip;
    Int_t   good_noclip;
    tree->Branch("dt12_noclip",     &dt12_noclip,     "dt12_noclip/F");
    tree->Branch("x_imp_noclip",    &x_imp_noclip,    "x_imp_noclip/F");
    tree->Branch("T_meas_noclip",   &T_meas_noclip,   "T_meas_noclip/F");
    tree->Branch("C_used_noclip",   &C_used_noclip,   "C_used_noclip/F");
    tree->Branch("tof_noclip",      &tof_noclip,      "tof_noclip/F");
    tree->Branch("path_len_noclip", &path_len_noclip, "path_len_noclip/F");
    tree->Branch("theta_noclip",    &theta_noclip,    "theta_noclip/F");
    tree->Branch("vel_noclip",      &vel_noclip,      "vel_noclip/F");
    tree->Branch("beta_noclip",     &beta_noclip,     "beta_noclip/F");
    tree->Branch("inv_vel_noclip",  &inv_vel_noclip,  "inv_vel_noclip/F");
    tree->Branch("good_noclip",     &good_noclip,     "good_noclip/I");

    // ---- BRANCH PARALLELI: pipeline Cconst (C costante — confronto) ----
    // Identica alla pipeline hybrid_tot ma usa GetC_const() al posto di GetC().
    // La posizione x_imp e' la stessa (dipende solo dalla retta Dt12(x)), ma
    // il TOF differisce per delta_C = C_poly(x) - <C>. Confrontare beta_Cconst
    // con beta (hybrid_tot) quantifica la sistematica del modello C(x).
    Float_t tof_Cconst, C_used_Cconst, beta_Cconst, inv_vel_Cconst, vel_Cconst;
    Int_t   good_Cconst;
    tree->Branch("tof_Cconst",     &tof_Cconst,     "tof_Cconst/F");
    tree->Branch("C_used_Cconst",  &C_used_Cconst,  "C_used_Cconst/F");
    tree->Branch("beta_Cconst",    &beta_Cconst,    "beta_Cconst/F");
    tree->Branch("inv_vel_Cconst", &inv_vel_Cconst, "inv_vel_Cconst/F");
    tree->Branch("vel_Cconst",     &vel_Cconst,     "vel_Cconst/F");
    tree->Branch("good_Cconst",    &good_Cconst,    "good_Cconst/I");

    // ---- BRANCH PARALLELI: pipeline rise_corr (correzione Pietro&Rick) ----
    Float_t dt12_risecorr, x_imp_risecorr, T_meas_risecorr, C_used_risecorr;
    Float_t tof_risecorr, path_len_risecorr, theta_risecorr;
    Float_t vel_risecorr, beta_risecorr, inv_vel_risecorr;
    Int_t   good_risecorr, good_risecorr_linrange;
    tree->Branch("dt12_risecorr",     &dt12_risecorr,     "dt12_risecorr/F");
    tree->Branch("x_imp_risecorr",    &x_imp_risecorr,    "x_imp_risecorr/F");
    tree->Branch("T_meas_risecorr",   &T_meas_risecorr,   "T_meas_risecorr/F");
    tree->Branch("C_used_risecorr",   &C_used_risecorr,   "C_used_risecorr/F");
    tree->Branch("tof_risecorr",      &tof_risecorr,      "tof_risecorr/F");
    tree->Branch("path_len_risecorr", &path_len_risecorr, "path_len_risecorr/F");
    tree->Branch("theta_risecorr",    &theta_risecorr,    "theta_risecorr/F");
    tree->Branch("vel_risecorr",      &vel_risecorr,      "vel_risecorr/F");
    tree->Branch("beta_risecorr",     &beta_risecorr,     "beta_risecorr/F");
    tree->Branch("inv_vel_risecorr",  &inv_vel_risecorr,  "inv_vel_risecorr/F");
    tree->Branch("good_risecorr",     &good_risecorr,     "good_risecorr/I");
    tree->Branch("good_risecorr_linrange",
                 &good_risecorr_linrange, "good_risecorr_linrange/I");

    // ---- BRANCH DIAGNOSTICI ----
    //   clipped_chan[k] = 1 se il canale k era clippato
    //   recov_tot[k]    = 1 se il canale k e' stato recuperato via TOT
    //   file_idx        = indice del file XML sorgente
    Int_t clipped_chan[3], recov_tot[3], file_idx;
    Float_t rise_time[3], t10_frac[3], t30_frac[3];
    Int_t   rise_time_ok[3], t1030_ok[3];
    tree->Branch("clipped_chan", clipped_chan, "clipped_chan[3]/I");
    tree->Branch("recov_tot",    recov_tot,    "recov_tot[3]/I");
    tree->Branch("file_idx",     &file_idx,    "file_idx/I");
    tree->Branch("rise_time",    rise_time,    "rise_time[3]/F");
    tree->Branch("t10_frac",     t10_frac,     "t10_frac[3]/F");
    tree->Branch("t30_frac",     t30_frac,     "t30_frac[3]/F");
    tree->Branch("rise_time_ok", rise_time_ok, "rise_time_ok[3]/I");
    tree->Branch("t1030_ok",     t1030_ok,     "t1030_ok[3]/I");

    // ---- BRANCH DIAGNOSTICI DEL VETO OFFLINE PMT4 ----
    // Salvati per OGNI evento (anche scartati): permettono il tuning a
    // posteriori della soglia di veto con tree->Draw().
    Float_t veto_amp, veto_t_peak, veto_sig, veto_thr;
    Float_t veto_baseline, veto_baseline_rms;
    Int_t   veto_nsamp, veto_fired, veto_avail;
    tree->Branch("veto_amp",          &veto_amp,          "veto_amp/F");
    tree->Branch("veto_t_peak",       &veto_t_peak,       "veto_t_peak/F");
    tree->Branch("veto_sig",          &veto_sig,          "veto_sig/F");
    tree->Branch("veto_thr",          &veto_thr,          "veto_thr/F");
    tree->Branch("veto_baseline",     &veto_baseline,     "veto_baseline/F");
    tree->Branch("veto_baseline_rms", &veto_baseline_rms, "veto_baseline_rms/F");
    tree->Branch("veto_nsamp",        &veto_nsamp,        "veto_nsamp/I");
    tree->Branch("veto_fired",        &veto_fired,        "veto_fired/I");
    tree->Branch("veto_avail",        &veto_avail,        "veto_avail/I");

  // ======================================================================
    //  ISTOGRAMMI — 4 PIPELINE COMPARABILI via PipelineHistos
    // ======================================================================
    //  Le 4 pipeline richieste, ciascuna con lo STESSO set di osservabili
    //  (vedi struct PipelineHistos, Sezione 3c). L'ordine dell'array fissa
    //  l'ordine di colori/legenda nei canvas comparativi:
    //    [0] Hybrid_tot   : validi + ricostruzione TOT dei clippati
    //    [1] No_clip      : solo non clippati (nessun recupero)
    //    [2] Rise_correct : validi + TOT + correzione rise time
    //    [3] Lin_Range    : come Rise_correct, ristretta a x in [-84,+70] cm
    enum PipeIdx { PIPE_HYBRID = 0, PIPE_NOCLIP, PIPE_RISECORR, PIPE_LINRANGE,
                   N_PIPES };
    PipelineHistos pipes[N_PIPES];
    BookPipeline(pipes[PIPE_HYBRID],   "hybrid",   "Hybrid_tot");
    BookPipeline(pipes[PIPE_NOCLIP],   "noclip",   "No_clip");
    BookPipeline(pipes[PIPE_RISECORR], "risecorr", "Rise_correct");
    BookPipeline(pipes[PIPE_LINRANGE], "linrange", "Lin_Range");

    // ---- Pipeline accessoria Cconst (sistematica del modello C(x)) ----
    //  NON e' una delle 4 comparabili: la teniamo a parte per la stima della
    //  sistematica di C(x), come da impostazione concordata. Riusa comunque la
    //  stessa struttura per uniformita' di booking/fill/write.
    PipelineHistos pipe_cconst;
    BookPipeline(pipe_cconst, "Cconst", "C costante");

    // ---- Istogrammi diagnostici del veto PMT4 ----
    TH1D *h_veto_amp_all = new TH1D("h_veto_amp_all",
        "Ampiezza nel CH4 (finestra di veto) - TUTTI gli eventi;A_{4} [mV];Conteggi",
        300, 0, 60);
    TH1D *h_veto_sig_all = new TH1D("h_veto_sig_all",
        "Significativit#grave{a} CH4 (finestra di veto);A_{4} / #sigma_{bl};Conteggi",
        300, 0, 60);
    TH1D *h_veto_baseline = new TH1D("h_veto_baseline",
        "Baseline media del CH4;baseline [mV];Conteggi", 200, -10, 10);
    TH1D *h_veto_baseline_rms = new TH1D("h_veto_baseline_rms",
        "RMS baseline del CH4;#sigma_{bl} [mV];Conteggi", 200, 0, 5);
    TH1D *h_veto_t_peak = new TH1D("h_veto_t_peak",
        "Tempo del minimo CH4 nella finestra;t_{peak} - t_{#mu} [ns];Conteggi",
        100, VETO_T_WIN_LO - 2, VETO_T_WIN_HI + 2);

    // ----------------------------------------------------------------------
    //  LIMITE FISICO SUL TOF
    // ----------------------------------------------------------------------
    // tof_min_phys = h/c e' il TOF minimo compatibile con la causalita'
    // (beta <= 1) per qualunque traiettoria. Definito qui (non come const
    // globale) perche' PAR_H puo' essere cambiato a runtime via SetGeometry().
    const double tof_min_phys = PAR_H / C_LIGHT;
    std::cout << "  TOF minimo fisico (h/c) = " << tof_min_phys
              << " ns  (h = " << PAR_H << " cm)" << std::endl;
    std::cout << "=============================================" << std::endl;

    // ======================================================================
    //  LOOP SUI FILE XML E SUGLI EVENTI
    // ======================================================================
int total_events = 0;
    // Contatori della pipeline hybrid_tot (principale)
    int total_good      = 0;
    int rej_no_cfd      = 0;   // < 3 canali o CFD non valido
    int rej_oscillating = 0;   // canale oscillante (vale per entrambe le pipeline)
    int rej_rise_outlier = 0;  // rise time valido > RISE_TIME_MAX su un canale
    int rej_clipped     = 0;   // canale clippato non recuperato
    int rej_x_out       = 0;   // x_imp fuori dai limiti della barra
    int rej_tof_unphys  = 0;   // TOF < h/c (non causale)
    // Contatori della pipeline noclip (diagnostica)
    int total_good_noclip     = 0;
    int rej_clipped_noclip    = 0;   // canale clippato (noclip non recupera mai)
    int rej_x_out_noclip      = 0;
    int rej_tof_unphys_noclip = 0;
    // Contatori della pipeline Cconst (confronto con C costante)
    int total_good_Cconst     = 0;
    int rej_tof_unphys_Cconst = 0;
    // Contatori della pipeline rise_corr (Pietro&Rick)
    int total_good_risecorr       = 0;
    int total_good_risecorr_lin   = 0;
    int rej_risecorr_inputs       = 0;  // evento non idoneo alla correzione
    int rej_x_out_risecorr        = 0;
    int rej_tof_unphys_risecorr   = 0;
    int rej_linrange_risecorr     = 0;  // buono rise_corr ma fuori [-84,+70]
    // Contatori del veto offline PMT4
    int rej_veto         = 0;   // eventi scartati perche' vetati dal CH4
    int n_veto_available = 0;   // eventi in cui il CH4 era cablato

    // Scatterplot rise time T90-T10 vs posizione ricostruita, uno per PMT.
    std::vector<double> rt_x[3], rt_y[3];

    for (size_t fi = 0; fi < file_list.size(); fi++) {

        std::vector<EventData> events;
        int n_parsed = LoadOrParseDataset(file_list[fi].c_str(), events);
        if (n_parsed <= 0) {
            std::cerr << "[WARNING] Nessun evento in " << file_list[fi]
                      << std::endl;
            continue;
        }

        // ------------------------------------------------------------------
        //  LOOP SUGLI EVENTI DEL FILE
        // ------------------------------------------------------------------
        for (size_t ev = 0; ev < events.size(); ev++) {

            EventData &e = events[ev];
            total_events++;
            file_idx = (Int_t)fi;

            // ---- Default "non valido" di TUTTI i branch dell'evento ----
            for (int k = 0; k < 3; k++) {
                t_cfd[k]        = -999.0f;
                amp[k]          = 0.0f;
                clipped_chan[k] = 0;
                recov_tot[k]    = 0;
                rise_time[k]    = -999.0f;
                t10_frac[k]     = -999.0f;
                t30_frac[k]     = -999.0f;
                rise_time_ok[k] = 0;
                t1030_ok[k]     = 0;
            }
            dt12 = x_imp = T_meas = C_used = tof = -999.0f;
            path_len = theta = vel = beta = inv_vel = -999.0f;
            good = 0;
            dt12_noclip = x_imp_noclip = T_meas_noclip = -999.0f;
            C_used_noclip = tof_noclip = -999.0f;
            path_len_noclip = theta_noclip = -999.0f;
            vel_noclip = beta_noclip = inv_vel_noclip = -999.0f;
            good_noclip = 0;
            tof_Cconst = C_used_Cconst = -999.0f;
            beta_Cconst = inv_vel_Cconst = vel_Cconst = -999.0f;
            good_Cconst = 0;
            dt12_risecorr = x_imp_risecorr = T_meas_risecorr = -999.0f;
            C_used_risecorr = tof_risecorr = -999.0f;
            path_len_risecorr = theta_risecorr = -999.0f;
            vel_risecorr = beta_risecorr = inv_vel_risecorr = -999.0f;
            good_risecorr = 0;
            good_risecorr_linrange = 0;
            veto_amp = veto_t_peak = veto_sig = veto_thr = -999.0f;
            veto_baseline = veto_baseline_rms = -999.0f;
            veto_nsamp = 0; veto_fired = 0; veto_avail = 0;

            // ============================================================
            //  TAGLIO 0: almeno 3 canali con impulso e CFD valido
            // ============================================================
            if (e.nchannels < 3) {
                tree->Fill();
                rej_no_cfd++;
                continue;
            }
            bool all_cfd_ok = true;
            for (int k = 0; k < 3; k++) {
                t_cfd[k] = (Float_t)e.ch[k].t_cfd;
                amp[k]   = (Float_t)e.ch[k].amplitude;
                rise_time[k]    = e.ch[k].rise_time_ok ? (Float_t)e.ch[k].rise_time : -999.0f;
                t10_frac[k]     = e.ch[k].t10_ok       ? (Float_t)e.ch[k].t_10      : -999.0f;
                t30_frac[k]     = e.ch[k].t30_ok       ? (Float_t)e.ch[k].t_30      : -999.0f;
                rise_time_ok[k] = e.ch[k].rise_time_ok ? 1 : 0;
                t1030_ok[k]     = (e.ch[k].t10_ok && e.ch[k].t30_ok) ? 1 : 0;
                if (!e.ch[k].has_pulse || !e.ch[k].cfd_ok) all_cfd_ok = false;
            }
            if (!all_cfd_ok) {
                tree->Fill();
                rej_no_cfd++;
                continue;
            }

            // ============================================================
            //  VETO OFFLINE PMT4 — calcolato SEMPRE, taglio SE disponibile
            // ============================================================
            // t_mu = media dei CFD di PMT1 e PMT2: i due PMT della barra
            // hanno timing affidabile e la media e' ~indipendente dalla
            // posizione di impatto (la luce si propaga simmetricamente verso
            // i due capi). A questo punto cfd_ok e' vero su PMT1 e PMT2.
            if (e.nchannels >= 4) {
                veto_avail = 1;
                n_veto_available++;

                double t_mu = 0.5 * (e.ch[0].t_cfd + e.ch[1].t_cfd);
                VetoResult vr = AnalyzeVeto(e.ch[VETO_CH_INDEX], t_mu);

                veto_amp          = (Float_t)vr.amplitude;
                veto_t_peak       = (Float_t)vr.t_peak;
                veto_sig          = (Float_t)vr.significance;
                veto_thr          = (Float_t)vr.v_thr_used;
                veto_baseline     = (Float_t)vr.baseline;
                veto_baseline_rms = (Float_t)vr.baseline_rms;
                veto_nsamp        = (Int_t)vr.n_samples_below;
                veto_fired        = vr.has_signal ? 1 : 0;

                // Istogrammi diagnostici del veto: riempiti per TUTTI gli
                // eventi con PMT4 disponibile (vetati o no), per vedere la
                // separazione rumore/segnale.
                if (vr.window_ok) {
                    h_veto_amp_all->Fill(vr.amplitude);
                    h_veto_sig_all->Fill(vr.significance);
                    h_veto_baseline->Fill(vr.baseline);
                    h_veto_baseline_rms->Fill(vr.baseline_rms);
                    h_veto_t_peak->Fill(vr.t_peak - t_mu);
                }

                // TAGLIO VETO: se il CH4 ha visto un segnale, il muone NON si
                // e' fermato nel piombo -> evento scartato da ENTRAMBE le pipeline.
                if (vr.has_signal) {
                    tree->Fill();
                    rej_veto++;
                    continue;
                }
            }
            // (se e.nchannels < 4: veto_avail = 0 e si prosegue senza veto)

            // ============================================================
            //  TAGLIO OSCILLAZIONE (vale per ENTRAMBE le pipeline)
            // ============================================================
            // Un'oscillazione e' una patologia della forma d'onda a monte
            // della ricostruzione: l'evento e' inutilizzabile da qualunque
            // metodo, indipendentemente dal recupero TOT.
bool any_osc = false;
            for (int k = 0; k < 3; k++) {
                if (e.ch[k].is_oscillating) { any_osc = true; break; }
            }
            if (any_osc) {
                tree->Fill();
                rej_oscillating++;
                continue;
            }

            // ============================================================
            //  TAGLIO OUTLIER RISE TIME (vale per TUTTE le pipeline)
            // ============================================================
            // Scarta l'evento se un qualsiasi canale ha rise time VALIDO e
            // > RISE_TIME_MAX (= 10 ns). Vedi commento alla costante: rimuove
            // forme d'onda rumorose con fronte mal definito che inquinerebbero
            // le distribuzioni di TOF/beta. Posto a monte del recupero TOT, agisce
            // identico su hybrid_tot, noclip, rise_corr, lin_range e Cconst.
            bool rise_outlier = false;
            for (int k = 0; k < 3; k++) {
                if (e.ch[k].rise_time_ok && e.ch[k].rise_time > RISE_TIME_MAX) {
                    rise_outlier = true;
                    break;
                }
            }
            if (rise_outlier) {
                tree->Fill();
                rej_rise_outlier++;
                continue;
            }

            // ============================================================
            //  RECUPERO DEI CLIPPATI VIA TOT POLINOMIALE
            // ============================================================
            // Eseguito sui canali clippati di PMT1, PMT2 e PMT3. PMT3 entra
            // in C = t3 - (t1+t2)/2: se clippato e recuperabile, va recuperato
            // anch'esso (a differenza della v8, dove PMT3 non era mai recuperato).
            for (int k = 0; k < 3; k++) {
                if (e.ch[k].is_clipped) RecoverClippedCFD_TOT(e.ch[k], k);
            }
            // Popolamento dei branch diagnostici (stato finale del recupero).
            for (int k = 0; k < 3; k++) {
                clipped_chan[k] = e.ch[k].is_clipped        ? 1 : 0;
                recov_tot[k]    = e.ch[k].cfd_recovered_tot ? 1 : 0;
            }

            // ============================================================
            //  LAMBDA RunPipeline: applica i tagli e calcola la cinematica
            // ============================================================
            // Riceve i 3 tempi gia' pronti (t_use) e i 3 flag "canale
            // utilizzabile" (ok_use). Restituisce uno status:
            //    0 = evento buono
            //   -1 = un canale non utilizzabile (clippato non recuperato,
            //        oppure tempo mancante)
            //   -2 = x_imp fuori dai limiti della barra
            //   -3 = TOF non causale (< h/c)
            // Gli output cinematici sono passati per riferimento.
            // c_mode: scelta del modello di C(x) da sottrarre.
            //   0 = GetC(x)          (C quadratico/lineare standard, hybrid/noclip)
            //   1 = GetC_corr(x)     (C quadratico corrected, Rise_correct)
            //   2 = GetC_const_lin(x)(C COSTANTE ristretto, Lin_Range — pipeline 2b)
            //   3 = GetC_lin_corr(x) (C lineare ristretto, Lin_Range 2a — pronto)
            auto RunPipeline = [&](const Float_t t_use[3],
                                   const bool    ok_use[3],
                                   double cal_m, double cal_q,
                                   int c_mode,
                                   Float_t &dt12_o,  Float_t &x_imp_o,
                                   Float_t &T_meas_o, Float_t &C_used_o,
                                   Float_t &tof_o,   Float_t &path_o,
                                   Float_t &theta_o, Float_t &vel_o,
                                   Float_t &beta_o,  Float_t &invv_o) -> int {

                // Inizializzazione output a sentinella.
                dt12_o = x_imp_o = T_meas_o = C_used_o = tof_o = -999.0f;
                path_o = theta_o = vel_o = beta_o = invv_o = -999.0f;

                // ---- TAGLIO 1: tutti e 3 i canali devono essere utilizzabili ----
                if (!ok_use[0] || !ok_use[1] || !ok_use[2]) return -1;

                // ---- PASSO 1: Dt12 e posizione di impatto ----
                dt12_o  = t_use[0] - t_use[1];
                x_imp_o = (Float_t)(((double)dt12_o - cal_q) / cal_m);

                // ---- TAGLIO 2: x_imp dentro la barra ----
                if (x_imp_o < X_CUT_LO || x_imp_o > X_CUT_HI) return -2;

                // ---- PASSO 2: tempo di volo ----
                T_meas_o = t_use[2] - (t_use[0] + t_use[1]) / 2.0f;
                double Cval;
                switch (c_mode) {
                    case 1:  Cval = GetC_corr((double)x_imp_o);      break;
                    case 2:  Cval = GetC_const_lin((double)x_imp_o); break;
                    case 3:  Cval = GetC_lin_corr((double)x_imp_o);  break;
                    default: Cval = GetC((double)x_imp_o);           break;
                }
                C_used_o = (Float_t)Cval;
                tof_o    = T_meas_o - C_used_o;

                // ---- TAGLIO 3: TOF fisicamente lecito (causalita') ----
                if (tof_o < tof_min_phys) return -3;

                // ---- PASSO 3: distanza percorsa (Pitagora) ----
                double dxp = (double)x_imp_o - PAR_X_PMT3;
                path_o = (Float_t)sqrt(dxp * dxp + PAR_H * PAR_H);

                // ---- PASSO 4: angolo di incidenza ----
                theta_o = (Float_t)atan2(dxp, PAR_H);

                // ---- PASSO 5: velocita', beta, 1/v ----
                if (tof_o > 0.1f) {
                    vel_o  = (Float_t)((double)path_o / (double)tof_o);
                    beta_o = (Float_t)((double)vel_o / C_LIGHT);
                    invv_o = (Float_t)((double)tof_o / (double)path_o);
                }
                return 0;
            };

            // ============================================================
            //  PIPELINE PRINCIPALE — hybrid_tot
            // ============================================================
            // Per ciascun canale: se NON clippato si usa t_cfd (gia' valido
            // dopo il TAGLIO 0); se clippato si usa t_cfd_rec_tot e il canale
            // e' utilizzabile solo se il recupero TOT e' riuscito.
            Float_t t_hy[3];
            bool    ok_hy[3];
            for (int k = 0; k < 3; k++) {
                if (!e.ch[k].is_clipped) {
                    t_hy[k]  = (Float_t)e.ch[k].t_cfd;
                    ok_hy[k] = true;
                    amp[k]   = (Float_t)e.ch[k].amplitude;
                } else if (e.ch[k].cfd_recovered_tot) {
                    t_hy[k]  = (Float_t)e.ch[k].t_cfd_rec_tot;
                    ok_hy[k] = true;
                    amp[k]   = (Float_t)e.ch[k].amplitude_tot;
                } else {
                    t_hy[k]  = -999.0f;
                    ok_hy[k] = false;
                    amp[k]   = (Float_t)e.ch[k].amplitude;
                }
                t_cfd[k] = t_hy[k];
            }
            int status_hy = RunPipeline(t_hy, ok_hy,
                                        CAL_M, CAL_Q, /*c_mode=*/0,
                                        dt12, x_imp, T_meas, C_used, tof,
                                        path_len, theta, vel, beta, inv_vel);
            switch (status_hy) {
                case  0: good = 1; total_good++;      break;
                case -1: good = 0; rej_clipped++;     break;
                case -2: good = 0; rej_x_out++;       break;
                case -3: good = 0; rej_tof_unphys++;  break;
                default: good = 0;                    break;
            }

            // ============================================================
            //  PIPELINE DIAGNOSTICA — noclip
            // ============================================================
            // Usa solo i canali NON clippati: un canale clippato rende
            // l'evento non utilizzabile (nessun recupero applicato).
            Float_t t_nc[3];
            bool    ok_nc[3];
            for (int k = 0; k < 3; k++) {
                if (!e.ch[k].is_clipped) {
                    t_nc[k]  = (Float_t)e.ch[k].t_cfd;
                    ok_nc[k] = true;
                } else {
                    t_nc[k]  = -999.0f;
                    ok_nc[k] = false;
                }
            }
            int status_nc = RunPipeline(t_nc, ok_nc,
                                        CAL_M, CAL_Q, /*c_mode=*/0,
                                        dt12_noclip, x_imp_noclip, T_meas_noclip,
                                        C_used_noclip, tof_noclip, path_len_noclip,
                                        theta_noclip, vel_noclip, beta_noclip,
                                        inv_vel_noclip);
            switch (status_nc) {
                case  0: good_noclip = 1; total_good_noclip++;     break;
                case -1: good_noclip = 0; rej_clipped_noclip++;    break;
                case -2: good_noclip = 0; rej_x_out_noclip++;      break;
                case -3: good_noclip = 0; rej_tof_unphys_noclip++; break;
                default: good_noclip = 0;                          break;
            }

            // ============================================================
            //  PIPELINE DI CONFRONTO — Cconst (C costante)
            // ============================================================
            // Usa gli stessi tempi e la stessa x_imp della pipeline hybrid_tot.
            // L'unica differenza e' la sottrazione di C: GetC_const() al posto
            // di GetC(). Questo isola l'effetto della modellizzazione C(x) sui
            // risultati di beta, permettendo una stima diretta della sistematica.
            if (good == 1 && gC_const_loaded) {
                C_used_Cconst = (Float_t)GetC_const((double)x_imp);
                tof_Cconst    = T_meas - C_used_Cconst;

                if (tof_Cconst >= tof_min_phys) {
                    double dxp = (double)x_imp - PAR_X_PMT3;
                    double l   = sqrt(dxp * dxp + PAR_H * PAR_H);
                    if (tof_Cconst > 0.1f) {
                        vel_Cconst     = (Float_t)(l / (double)tof_Cconst);
                        beta_Cconst    = (Float_t)((double)vel_Cconst / C_LIGHT);
                        inv_vel_Cconst = (Float_t)((double)tof_Cconst / l);
                    }
                    good_Cconst = 1;
                    total_good_Cconst++;
                } else {
                    good_Cconst = 0;
                    rej_tof_unphys_Cconst++;
                }
            } else {
                good_Cconst = 0;
            }

            // ============================================================
            //  PIPELINE RISE_CORR — correzione Pietro&Rick
            // ============================================================
            // Si parte dal campione hybrid_tot buono, ma la correzione viene
            // applicata solo se i tre segnali non sono clippati e hanno t_10/t_30.
            // Per ciascun canale:
            //   dT_n = t_30,n - t_10,n
            //   t_corr,n = t_cfd,n - dT_n
            // Poi x_imp, C(x), TOF, beta e 1/v vengono ricalcolati con i
            // parametri corrected letti dalla calibrazione, se disponibili.
// Variabili locali per la ricostruzione INDIPENDENTE di Lin_Range
            // (non sono branch del TTree: servono solo a riempire la pipeline
            // PIPE_LINRANGE). Inizializzate a sentinella a ogni evento.
            Float_t dt12_lin = -999.0f, x_imp_lin = -999.0f, T_meas_lin = -999.0f;
            Float_t C_used_lin = -999.0f, tof_lin = -999.0f, path_len_lin = -999.0f;
            Float_t theta_lin = -999.0f, vel_lin = -999.0f, beta_lin = -999.0f;
            Float_t inv_vel_lin = -999.0f;

            {
                bool no_clip_any = !e.ch[0].is_clipped
                                && !e.ch[1].is_clipped
                                && !e.ch[2].is_clipped;
                bool t1030_all_ok = e.ch[0].t10_ok && e.ch[0].t30_ok
                                  && e.ch[1].t10_ok && e.ch[1].t30_ok
                                  && e.ch[2].t10_ok && e.ch[2].t30_ok;

                if (good == 1 && no_clip_any && t1030_all_ok) {
                    // Tempi rise-corrected, COMUNI a Rise_correct e Lin_Range:
                    //   t_corr,n = t_cfd,n - (t_30,n - t_10,n)
                    Float_t t_rc[3];
                    bool    ok_rc[3] = { true, true, true };
                    for (int k = 0; k < 3; k++) {
                        double dT = e.ch[k].t_30 - e.ch[k].t_10;
                        t_rc[k] = (Float_t)(e.ch[k].t_cfd - dT);
                    }

                    // ---- Rise_correct: retta corrected GLOBALE + C quadratico corrected ----
                    double m_use = gRiseCorrCal_loaded ? CAL_M_CORR : CAL_M;
                    double q_use = gRiseCorrCal_loaded ? CAL_Q_CORR : CAL_Q;
                    int status_rc = RunPipeline(t_rc, ok_rc,
                                                m_use, q_use, /*c_mode=*/1,
                                                dt12_risecorr, x_imp_risecorr,
                                                T_meas_risecorr, C_used_risecorr,
                                                tof_risecorr, path_len_risecorr,
                                                theta_risecorr, vel_risecorr,
                                                beta_risecorr, inv_vel_risecorr);

                    switch (status_rc) {
                        case  0: good_risecorr = 1; total_good_risecorr++; break;
                        case -2: rej_x_out_risecorr++;                     break;
                        case -3: rej_tof_unphys_risecorr++;                break;
                        default: rej_risecorr_inputs++;                    break;
                    }

                    // ---- Lin_Range: ricostruzione INDIPENDENTE ----
                    //  Stessi tempi rise-corrected, ma con:
                    //    - retta RISTRETTA al range lineare (m_corr_lin, q_corr_lin);
                    //    - C COSTANTE ristretto (c_mode = 2 -> GetC_const_lin);
                    //    - taglio x_imp nel range lineare, applicato al NUOVO
                    //      x_imp_lin ricostruito con la retta ristretta.
                    //  Fallback: se la calibrazione ristretta non e' disponibile,
                    //  usa la retta corrected globale (gLinRangeCal_loaded == false).
                    double m_lin_use = gLinRangeCal_loaded ? CAL_M_LIN : m_use;
                    double q_lin_use = gLinRangeCal_loaded ? CAL_Q_LIN : q_use;
                    int status_lin = RunPipeline(t_rc, ok_rc,
                                                 m_lin_use, q_lin_use, /*c_mode=*/2,
                                                 dt12_lin, x_imp_lin, T_meas_lin,
                                                 C_used_lin, tof_lin, path_len_lin,
                                                 theta_lin, vel_lin, beta_lin,
                                                 inv_vel_lin);

                    if (status_lin == 0) {
                        // Taglio range lineare sul NUOVO x_imp_lin.
                        if (x_imp_lin >= RISE_LINEAR_X_LO &&
                            x_imp_lin <= RISE_LINEAR_X_HI) {
                            good_risecorr_linrange = 1;
                            total_good_risecorr_lin++;
                        } else {
                            rej_linrange_risecorr++;
                        }
                    } else if (status_lin == -2) {
                        rej_x_out_risecorr++;
                    } else if (status_lin == -3) {
                        rej_tof_unphys_risecorr++;
                    }
                } else {
                    rej_risecorr_inputs++;
                }
            }

// ============================================================
            //  RIEMPIMENTO ISTOGRAMMI — via FillPipeline (4 pipeline + Cconst)
            // ============================================================
            //  La logica fisica (good*, valori ricostruiti) e' invariata: qui
            //  smistiamo solo i valori nelle rispettive PipelineHistos. Ogni
            //  pipeline riceve adesso il SET COMPLETO di osservabili (beta,
            //  beta_norm, tof, 1/v, theta, x_imp, path, 2D), cosi' i 4 canvas
            //  riassuntivi sono per costruzione omologhi.

            // ---- Hybrid_tot ----
            //  Lo scatterplot rise time vs x viene popolato qui (solo hybrid):
            //  e' la base del grafico diagnostico del rise time vs posizione.
            if (good == 1) {
                for (int k = 0; k < 3; k++) {
                    if (e.ch[k].rise_time_ok) {
                        rt_x[k].push_back((double)x_imp);
                        rt_y[k].push_back(e.ch[k].rise_time);
                    }
                }
                FillPipeline(pipes[PIPE_HYBRID],
                             tof, beta, inv_vel, theta, x_imp, path_len, T_meas);
            }

            // ---- No_clip ----
            if (good_noclip == 1) {
                FillPipeline(pipes[PIPE_NOCLIP],
                             tof_noclip, beta_noclip, inv_vel_noclip,
                             theta_noclip, x_imp_noclip, path_len_noclip,
                             T_meas_noclip);
            }

            // ---- Rise_correct ----
            if (good_risecorr == 1) {
                FillPipeline(pipes[PIPE_RISECORR],
                             tof_risecorr, beta_risecorr, inv_vel_risecorr,
                             theta_risecorr, x_imp_risecorr, path_len_risecorr,
                             T_meas_risecorr);
            }

            // ---- Lin_Range (ricostruzione INDIPENDENTE) ----
            //  Ora Lin_Range NON e' piu' un sottoinsieme di Rise_correct: usa la
            //  retta ristretta (m_corr_lin, q_corr_lin) e il C costante ristretto
            //  (GetC_const_lin). Riempiamo con le grandezze "_lin" calcolate dalla
            //  sua ricostruzione indipendente, NON con quelle "_risecorr".
            if (good_risecorr_linrange == 1) {
                FillPipeline(pipes[PIPE_LINRANGE],
                             tof_lin, beta_lin, inv_vel_lin,
                             theta_lin, x_imp_lin, path_len_lin,
                             T_meas_lin);
            }

            // ---- Cconst (accessoria, fuori dai 4 comparativi) ----
            //  x_imp e' la stessa della pipeline hybrid (la posizione dipende
            //  solo dalla retta Dt12(x)); cambia solo C -> TOF/beta/1v.
            if (good_Cconst == 1) {
                FillPipeline(pipe_cconst,
                             tof_Cconst, beta_Cconst, inv_vel_Cconst,
                             theta, x_imp, path_len, T_meas);
            }

            tree->Fill();
        }   // fine loop eventi
    }       // fine loop file

    // ======================================================================
    //  RIEPILOGO DEI TAGLI DI QUALITA'
    // ======================================================================
    auto pct = [&](int n) -> double {
        return (total_events > 0) ? 100.0 * n / total_events : 0.0;
    };

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  RIEPILOGO TAGLI DI QUALITA'                 " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Eventi totali:           " << total_events << std::endl;
    std::cout << "  Scartati (no CFD):       " << rej_no_cfd
              << " (" << pct(rej_no_cfd) << "%)" << std::endl;
std::cout << "  Scartati (oscillazione): " << rej_oscillating
              << " (" << pct(rej_oscillating) << "%)" << std::endl;
    std::cout << "  Scartati (rise>" << RISE_TIME_MAX << "ns): " << rej_rise_outlier
              << " (" << pct(rej_rise_outlier) << "%)" << std::endl;
    std::cout << "  Scartati (VETO PMT4):    " << rej_veto
              << " (" << pct(rej_veto) << "%)   [su " << n_veto_available
              << " eventi con CH4 disponibile]" << std::endl;
    std::cout << "  --------- PIPELINE hybrid_tot (principale) ---------"
              << std::endl;
    std::cout << "  Scartati (clipped non recup.): " << rej_clipped
              << " (" << pct(rej_clipped) << "%)" << std::endl;
    std::cout << "  Scartati (x fuori barra):      " << rej_x_out
              << " (" << pct(rej_x_out) << "%)" << std::endl;
    std::cout << "  Scartati (TOF < h/c):          " << rej_tof_unphys
              << " (" << pct(rej_tof_unphys) << "%)" << std::endl;
    std::cout << "  Eventi buoni (hybrid_tot):     " << total_good
              << " (" << pct(total_good) << "%)" << std::endl;
    std::cout << "  --------- PIPELINE noclip (diagnostica) ------------"
              << std::endl;
    std::cout << "  Scartati (clipped):            " << rej_clipped_noclip
              << " (" << pct(rej_clipped_noclip) << "%)" << std::endl;
    std::cout << "  Scartati (x fuori barra):      " << rej_x_out_noclip
              << " (" << pct(rej_x_out_noclip) << "%)" << std::endl;
    std::cout << "  Scartati (TOF < h/c):          " << rej_tof_unphys_noclip
              << " (" << pct(rej_tof_unphys_noclip) << "%)" << std::endl;
    std::cout << "  Eventi buoni (noclip):         " << total_good_noclip
              << " (" << pct(total_good_noclip) << "%)" << std::endl;
    std::cout << "  --------- PIPELINE Cconst (C costante — confronto) ----"
              << std::endl;
    std::cout << "  Scartati (TOF < h/c):          " << rej_tof_unphys_Cconst
              << " (" << pct(rej_tof_unphys_Cconst) << "%)" << std::endl;
    std::cout << "  Eventi buoni (Cconst):         " << total_good_Cconst
              << " (" << pct(total_good_Cconst) << "%)" << std::endl;
    std::cout << "  --------- PIPELINE rise_corr (Pietro&Rick) --------"
              << std::endl;
    std::cout << "  Scartati (input non idonei):    " << rej_risecorr_inputs
              << " (" << pct(rej_risecorr_inputs) << "%)" << std::endl;
    std::cout << "  Scartati (x fuori barra):       " << rej_x_out_risecorr
              << " (" << pct(rej_x_out_risecorr) << "%)" << std::endl;
    std::cout << "  Scartati (TOF < h/c):           " << rej_tof_unphys_risecorr
              << " (" << pct(rej_tof_unphys_risecorr) << "%)" << std::endl;
    std::cout << "  Eventi buoni (rise_corr):       " << total_good_risecorr
              << " (" << pct(total_good_risecorr) << "%)" << std::endl;
    std::cout << "  Fuori range lineare rise_corr:  " << rej_linrange_risecorr
              << " (" << pct(rej_linrange_risecorr) << "%)" << std::endl;
    std::cout << "  Eventi buoni rise_corr range:   " << total_good_risecorr_lin
              << " (" << pct(total_good_risecorr_lin) << "%)" << std::endl;
    std::cout << "=============================================" << std::endl;

    // ======================================================================
    //  STATISTICHE DELLA DISTRIBUZIONE beta (pipeline hybrid_tot)
    // ======================================================================
    // Tre statistiche complementari sull'istogramma h_beta (solo eventi
    // good == 1, con beta nel range [BETA_LO, BETA_HI]):
    //   [1] grezza   : media e RMS dei contenuti dei bin (sensibile alle code);
    //   [2] fit gauss: fit gaussiano del picco (valore "tipico" robusto);
    //   [3] robusta  : mediana e MAD (indipendenti dalla gaussianita').
// Alias locale: le statistiche beta sono calcolate sulla pipeline Hybrid_tot.
    // Agganciamo il vecchio nome h_beta al nuovo istogramma della struct, cosi'
    // il resto del blocco statistico resta invariato.
    TH1D* h_beta = pipes[PIPE_HYBRID].h_beta;

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  STATISTICHE DISTRIBUZIONE beta = v/c  (Hybrid_tot)" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Eventi nell'istogramma h_beta: "
              << (long)h_beta->GetEntries() << std::endl;

    // ---- [1] Statistica grezza dell'istogramma ----
    double beta_mean_raw = h_beta->GetMean();
    double beta_rms_raw  = h_beta->GetStdDev();
    double beta_mean_err = (h_beta->GetEntries() > 1)
                         ? beta_rms_raw / sqrt(h_beta->GetEntries()) : 0.0;
    std::cout << "\n  [1] Statistica grezza (istogramma):" << std::endl;
    std::cout << "      <beta>   = " << beta_mean_raw
              << " +/- " << beta_mean_err << std::endl;
    std::cout << "      sigma    = " << beta_rms_raw << std::endl;

    // ---- [2] Fit gaussiano del picco ----
    double beta_mean_fit = -1, beta_sigma_fit = -1;
    double beta_mean_fit_err = 0, beta_sigma_fit_err = 0;
    double beta_chi2_ndf = -1;
    if (h_beta->GetEntries() > 50) {
        int    bin_peak = h_beta->GetMaximumBin();
        double x_peak   = h_beta->GetBinCenter(bin_peak);
        double fit_lo   = std::max(x_peak - 0.3, BETA_LO + 0.01);
        double fit_hi   = std::min(x_peak + 0.3, BETA_HI - 0.01);

        TF1 *fgaus = new TF1("fgaus_beta", "gaus", fit_lo, fit_hi);
        fgaus->SetParameters(h_beta->GetMaximum(), x_peak, 0.1);
        TFitResultPtr fr = h_beta->Fit(fgaus, "QRS0");
        if (fr.Get() && fr->IsValid()) {
            beta_mean_fit      = fgaus->GetParameter(1);
            beta_sigma_fit     = fabs(fgaus->GetParameter(2));
            beta_mean_fit_err  = fgaus->GetParError(1);
            beta_sigma_fit_err = fgaus->GetParError(2);
            int ndf = fgaus->GetNDF();
            beta_chi2_ndf = (ndf > 0) ? fgaus->GetChisquare() / ndf : -1;
            std::cout << "\n  [2] Fit gaussiano del picco "
                      << "(range [" << fit_lo << ", " << fit_hi << "]):"
                      << std::endl;
            std::cout << "      mu       = " << beta_mean_fit
                      << " +/- " << beta_mean_fit_err << std::endl;
            std::cout << "      sigma    = " << beta_sigma_fit
                      << " +/- " << beta_sigma_fit_err << std::endl;
            std::cout << "      chi2/ndf = " << beta_chi2_ndf
                      << "  (ndf = " << ndf << ")" << std::endl;
        } else {
            std::cout << "\n  [2] Fit gaussiano: NON CONVERGENTE." << std::endl;
        }
    } else {
        std::cout << "\n  [2] Fit gaussiano: SALTATO (poche entries: "
                  << h_beta->GetEntries() << ")" << std::endl;
    }

    // ---- [3] Statistica robusta: mediana e MAD ----
    double beta_median = -1, beta_mad = -1;
    if (h_beta->GetEntries() > 0) {
        double qprob[1] = { 0.5 };
        double qval[1]  = { 0.0 };
        h_beta->GetQuantiles(1, qval, qprob);
        beta_median = qval[0];

        // MAD = mediana di |beta - mediana|: costruita bin per bin.
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
        double sigma_from_mad = (beta_mad > 0) ? 1.4826 * beta_mad : -1;
        std::cout << "\n  [3] Statistica robusta:" << std::endl;
        std::cout << "      mediana    = " << beta_median << std::endl;
        std::cout << "      MAD        = " << beta_mad << std::endl;
        std::cout << "      1.4826*MAD = " << sigma_from_mad
                  << "  (~ sigma se distribuzione gaussiana)" << std::endl;
    }
    std::cout << "=============================================" << std::endl;

// ======================================================================
    //  NORMALIZZAZIONE DEGLI ISTOGRAMMI beta -> PDF (area unitaria)
    // ======================================================================
    //  Trasformiamo h_beta_norm di ciascuna pipeline in una densita' di
    //  probabilita' (vedi NormalizeBetaPDF, Sezione 3c). Dopo questa operazione
    //  l'asse Y e' in [1/beta] e le distribuzioni sono confrontabili anche con
    //  statistiche diverse. h_beta (usato per le statistiche) NON viene toccato.
    for (int p = 0; p < N_PIPES; p++) NormalizeBetaPDF(pipes[p]);
    NormalizeBetaPDF(pipe_cconst);

   // ======================================================================
    //  SALVATAGGIO E CANVAS
    // ======================================================================
    fout->cd();
    tree->Write();

    // ----------------------------------------------------------------------
    //  CANVAS 2: distribuzione beta (Hybrid_tot) con fit gaussiano e statistiche
    // ----------------------------------------------------------------------
    //  Mantenuto come pannello "ricco" della pipeline principale: riporta il
    //  fit gaussiano del picco e le tre statistiche (grezza, fit, robusta).
    //  Le altre pipeline hanno il loro riassunto in DrawPipelineSummary.
    TCanvas *c2 = new TCanvas("c_beta", "Distribuzione #beta (Hybrid_tot)", 900, 650);
    c2->SetGrid();
    h_beta->SetLineColor(kRed + 1);
    h_beta->SetLineWidth(2);
    h_beta->Draw();

    TF1 *fgaus_show = (TF1*)gROOT->FindObject("fgaus_beta");
    if (fgaus_show && beta_mean_fit > 0) {
        fgaus_show->SetLineColor(kBlue + 2);
        fgaus_show->SetLineWidth(2);
        fgaus_show->Draw("same");
    }
    TLine *l_beta1 = new TLine(1.0, 0, 1.0, h_beta->GetMaximum() * 0.95);
    l_beta1->SetLineColor(kBlack);
    l_beta1->SetLineStyle(2);
    l_beta1->SetLineWidth(2);
    l_beta1->Draw("same");

    TPaveText *pt = new TPaveText(0.60, 0.55, 0.89, 0.88, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.030);
    pt->AddText(Form("Entries = %lld", (Long64_t)h_beta->GetEntries()));
    pt->AddText("");
    pt->AddText("#bf{Statistica grezza}");
    pt->AddText(Form("#LT#beta#GT = %.4f #pm %.4f", beta_mean_raw, beta_mean_err));
    pt->AddText(Form("#sigma(#beta) = %.4f", beta_rms_raw));
    if (beta_mean_fit > 0) {
        pt->AddText("");
        pt->AddText("#bf{Fit gaussiano}");
        pt->AddText(Form("#mu = %.4f #pm %.4f", beta_mean_fit, beta_mean_fit_err));
        pt->AddText(Form("#sigma = %.4f #pm %.4f", beta_sigma_fit, beta_sigma_fit_err));
        pt->AddText(Form("#chi^{2}/ndf = %.2f", beta_chi2_ndf));
    }
    if (beta_median > 0) {
        pt->AddText("");
        pt->AddText("#bf{Robusta}");
        pt->AddText(Form("mediana = %.4f", beta_median));
        pt->AddText(Form("MAD = %.4f", beta_mad));
    }
    pt->Draw();
    c2->Write("Canvas_beta_Hybrid");

    // ----------------------------------------------------------------------
    //  CANVAS RIASSUNTIVI PER PIPELINE (3x3, una per pipeline)
    // ----------------------------------------------------------------------
    //  Un pannello omologo per ciascuna delle 4 pipeline + Cconst: stessa
    //  disposizione, stesse 9 osservabili. E' il "diagnostico interno" di ogni
    //  metodo. Nomi ROOT: Canvas_Summary_<tag>.
    for (int p = 0; p < N_PIPES; p++) DrawPipelineSummary(pipes[p]);
    DrawPipelineSummary(pipe_cconst);

    // ----------------------------------------------------------------------
    //  6 CANVAS COMPARATIVI tra le 4 pipeline (con legenda)
    // ----------------------------------------------------------------------
    //  Ogni canvas sovrappone la STESSA osservabile delle 4 pipeline. I
    //  selettori estraggono dalla struct l'istogramma giusto; DrawPipelineComparison
    //  pensa a colori/stili/legenda/scala. Le 6 osservabili richieste:
    //    beta, beta_norm (PDF), TOF, 1/v, x_imp, distanza percorsa.
    DrawPipelineComparison(pipes, N_PIPES,
        [](PipelineHistos& p){ return p.h_beta; },
        "Canvas_Compare_beta",
        "#beta: confronto 4 metodi;#beta;Conteggi",
        /*draw_beta1_line=*/true);

    DrawPipelineComparison(pipes, N_PIPES,
        [](PipelineHistos& p){ return p.h_beta_norm; },
        "Canvas_Compare_beta_PDF",
        "PDF #beta: confronto 4 metodi;#beta;PDF",
        /*draw_beta1_line=*/true);

    DrawPipelineComparison(pipes, N_PIPES,
        [](PipelineHistos& p){ return p.h_tof; },
        "Canvas_Compare_TOF",
        "TOF: confronto 4 metodi;TOF [ns];Conteggi");

    DrawPipelineComparison(pipes, N_PIPES,
        [](PipelineHistos& p){ return p.h_invv; },
        "Canvas_Compare_invv",
        "1/v: confronto 4 metodi;1/v [ns/cm];Conteggi");

    DrawPipelineComparison(pipes, N_PIPES,
        [](PipelineHistos& p){ return p.h_ximp; },
        "Canvas_Compare_ximp",
        "x_{imp}: confronto 4 metodi;x_{imp} [cm];Conteggi");

    DrawPipelineComparison(pipes, N_PIPES,
        [](PipelineHistos& p){ return p.h_path; },
        "Canvas_Compare_path",
        "Distanza percorsa: confronto 4 metodi;l [cm];Conteggi");

    // ----------------------------------------------------------------------
    //  CANVAS ACCESSORIO: sistematica del modello C(x) (Hybrid vs Cconst vs noclip)
    // ----------------------------------------------------------------------
    //  Mantiene la diagnostica sulla sistematica di C(x): confronta la pipeline
    //  Hybrid_tot (C quadratico) con Cconst (C costante) e No_clip su TOF, beta,
    //  1/v. Se Hybrid e Cconst coincidono, la sistematica di C(x) e' trascurabile.
    {
        TCanvas *c_Ccmp = new TCanvas("c_compare_Cmodels",
                                      "Sistematica C(x): poly vs costante", 1200, 800);
        c_Ccmp->Divide(2, 2);

        auto drawTriple = [&](int pad, TH1D* h_poly, TH1D* h_cc, TH1D* h_nc,
                              const char* tit) {
            c_Ccmp->cd(pad);
            gPad->SetGrid();
            h_poly->SetLineColor(kRed + 1);   h_poly->SetLineWidth(2);
            h_cc->SetLineColor(kOrange + 2);  h_cc->SetLineWidth(2);
            h_cc->SetLineStyle(7);
            h_nc->SetLineColor(kBlue + 1);    h_nc->SetLineWidth(2);
            h_nc->SetLineStyle(2);
            double ymax = std::max({h_poly->GetMaximum(),
                                    h_cc->GetMaximum(),
                                    h_nc->GetMaximum()});
            h_poly->SetTitle(tit);
            h_poly->SetMaximum(1.15 * ymax);
            h_poly->Draw("HIST");
            h_cc->Draw("HIST same");
            h_nc->Draw("HIST same");
            TLegend *lg = new TLegend(0.55, 0.68, 0.89, 0.89);
            lg->SetTextFont(42); lg->SetTextSize(0.035);
            lg->AddEntry(h_poly, "C(x) quadratico (Hybrid)", "l");
            lg->AddEntry(h_cc,   "C costante",                "l");
            lg->AddEntry(h_nc,   "No_clip",                   "l");
            lg->Draw();
        };
        drawTriple(1, pipes[PIPE_HYBRID].h_tof,  pipe_cconst.h_tof,
                      pipes[PIPE_NOCLIP].h_tof,
                   "TOF: C(x) poly vs costante;TOF [ns];Conteggi");
        drawTriple(2, pipes[PIPE_HYBRID].h_beta, pipe_cconst.h_beta,
                      pipes[PIPE_NOCLIP].h_beta,
                   "#beta: C(x) poly vs costante;#beta;Conteggi");
        drawTriple(3, pipes[PIPE_HYBRID].h_invv, pipe_cconst.h_invv,
                      pipes[PIPE_NOCLIP].h_invv,
                   "1/v: C(x) poly vs costante;1/v [ns/cm];Conteggi");

        c_Ccmp->cd(4);
        gPad->SetGrid();
        pipe_cconst.h2_beta_x->Draw("COLZ");
        TLine *l_b1_cc = new TLine(X_LO, 1.0, X_HI, 1.0);
        l_b1_cc->SetLineColor(kRed); l_b1_cc->SetLineStyle(2);
        l_b1_cc->Draw("same");

        c_Ccmp->Write("Canvas_Compare_Cmodels");
    }

    // ---- Canvas 10d: scatterplot rise time T90-T10 vs posizione ricostruita ----
    {
        const char* pmt_lab[3] = { "PMT1", "PMT2", "PMT3" };
        const int   pmt_col[3] = { kRed + 1, kBlue + 1, kGreen + 2 };
        TCanvas *c_rt = new TCanvas("c_risetime_vs_x",
            "Rise time T90-T10 vs posizione ricostruita", 1200, 420);
        c_rt->Divide(3, 1);

        for (int k = 0; k < 3; k++) {
            c_rt->cd(k + 1);
            gPad->SetGrid();
            if (!rt_x[k].empty()) {
                TGraph *g_rt = new TGraph((int)rt_x[k].size(),
                                          rt_x[k].data(), rt_y[k].data());
                g_rt->SetName(Form("g_risetime_vs_x_%s", pmt_lab[k]));
                g_rt->SetTitle(Form("Rise time T_{90}-T_{10} vs x - %s;"
                                    "x_{imp} [cm];Rise time [ns]", pmt_lab[k]));
                g_rt->SetMarkerStyle(20);
                g_rt->SetMarkerSize(0.35);
                g_rt->SetMarkerColor(pmt_col[k]);
                g_rt->SetLineColor(pmt_col[k]);
                g_rt->Draw("AP");
                g_rt->Write();
} else {
                TLatex lat;
                lat.SetNDC();
                lat.SetTextFont(42);
                lat.SetTextSize(0.045);
                lat.DrawLatex(0.20, 0.50, Form("%s: nessun rise time valido", pmt_lab[k]));
            }
        }
        c_rt->Write("Canvas_RiseTime_vs_x");
    }

    // ---- Canvas 10e: PROFILO rise time medio vs posizione (bin da 10 cm) ----
    //  Per ciascun PMT, il rise time medio per bin di 10 cm in x_imp, con
    //  errore standard della media (barra verticale) e semi-larghezza del bin
    //  (barra orizzontale). Costruito dagli stessi (x_imp, rise_time) raccolti
    //  per lo scatterplot del campione Hybrid_tot (rt_x/rt_y).
    //
    //  Lettura fisica attesa: andamento monotono per PMT1 e PMT2 (di pendenza
    //  opposta, perche' i due PMT stanno ai capi opposti della barra) e profilo
    //  ~piatto per PMT3 (distanza dalla barra ~costante). La pendenza quantifica
    //  il rallentamento del fronte con il cammino ottico nella barra.
    {
        const char* pmt_lab[3] = { "PMT1", "PMT2", "PMT3" };
        const int   pmt_col[3] = { kRed + 1, kBlue + 1, kGreen + 2 };

        TCanvas *c_rtp = new TCanvas("c_risetime_profile_vs_x",
            "Profilo rise time medio vs posizione", 1200, 420);
        c_rtp->Divide(3, 1);

        for (int k = 0; k < 3; k++) {
            c_rtp->cd(k + 1);
            gPad->SetGrid();

            TGraphErrors* g_prof = BuildRiseProfile(
                rt_x[k], rt_y[k],
                RISEPROF_X_LO, RISEPROF_X_HI, RISEPROF_BIN_W, RISEPROF_N_MIN,
                Form("g_riseprofile_%s", pmt_lab[k]),
                Form("Rise time medio vs x - %s;x_{imp} [cm];#LT#tau_{rise}#GT [ns]",
                     pmt_lab[k]),
                pmt_col[k]);

            if (g_prof->GetN() > 0) {
                // "A" assi, "P" punti, "E" barre d'errore.
                g_prof->Draw("APE");
                g_prof->Write();
            } else {
                TLatex lat;
                lat.SetNDC();
                lat.SetTextFont(42);
                lat.SetTextSize(0.045);
                lat.DrawLatex(0.18, 0.50,
                    Form("%s: nessun bin con N >= %d", pmt_lab[k], RISEPROF_N_MIN));
                delete g_prof;
            }
        }
        c_rtp->Write("Canvas_RiseTime_Profile_vs_x");
    }

    // ---- Canvas 10: diagnostica del veto PMT4 (2x2) ----
    TCanvas *c_veto = new TCanvas("c_veto", "Diagnostica veto PMT4", 1200, 800);
    c_veto->Divide(2, 2);
    c_veto->cd(1); gPad->SetGrid(); gPad->SetLogy();
    h_veto_amp_all->SetLineColor(kViolet + 1);
    h_veto_amp_all->SetLineWidth(2);
    h_veto_amp_all->Draw("HIST");
    {
        TLine *l_thr = new TLine(VETO_AMP_FLOOR, 0.5,
                                 VETO_AMP_FLOOR, h_veto_amp_all->GetMaximum());
        l_thr->SetLineColor(kRed); l_thr->SetLineStyle(2); l_thr->SetLineWidth(2);
        l_thr->Draw("same");
    }
    c_veto->cd(2); gPad->SetGrid(); gPad->SetLogy();
    h_veto_sig_all->SetLineColor(kViolet + 1);
    h_veto_sig_all->SetLineWidth(2);
    h_veto_sig_all->Draw("HIST");
    {
        TLine *l_nsig = new TLine(VETO_N_SIGMA, 0.5,
                                  VETO_N_SIGMA, h_veto_sig_all->GetMaximum());
        l_nsig->SetLineColor(kRed); l_nsig->SetLineStyle(2); l_nsig->SetLineWidth(2);
        l_nsig->Draw("same");
    }
    c_veto->cd(3); gPad->SetGrid();
    h_veto_baseline->SetLineColor(kAzure + 1);
    h_veto_baseline->SetLineWidth(2);
    h_veto_baseline->Draw("HIST");
    c_veto->cd(4); gPad->SetGrid();
    h_veto_t_peak->SetLineColor(kSpring + 3);
    h_veto_t_peak->SetLineWidth(2);
    h_veto_t_peak->Draw("HIST");
    c_veto->Write("Canvas_Veto_Diagnostics");

// ---- Salvataggio di tutti gli istogrammi ----
    // 4 pipeline comparabili + Cconst, via WritePipeline (Sezione 3c).
    for (int p = 0; p < N_PIPES; p++) WritePipeline(pipes[p]);
    WritePipeline(pipe_cconst);
    // Diagnostica veto (invariata).
    h_veto_amp_all->Write();      h_veto_sig_all->Write();
    h_veto_baseline->Write();     h_veto_baseline_rms->Write();
    h_veto_t_peak->Write();

    // ======================================================================
    //  FIT DELLE OSSERVABILI FISICHE (beta, TOF, 1/v) — 3 metodi per pipeline
    // ======================================================================
    //  Per ciascuna pipeline (4 comparabili + Cconst separata) applichiamo i 3
    //  metodi di fit (gaussiano-picco, EMG, RooKeysPdf) a beta, TOF e 1/v.
    //  Riferimenti fisici marcati sui canvas:
    //    beta -> linea a 1 (causalita'); 1/v -> linea a 1/c; TOF -> nessuna.
    //  beta ha coda a SINISTRA (left_tail=true): EMG riflessa. TOF e 1/v hanno
    //  coda a DESTRA (left_tail=false): EMG diretta.
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  FIT DELLE OSSERVABILI (3 metodi per pipeline)" << std::endl;
    std::cout << "=============================================" << std::endl;

    const double INV_C = 1.0 / C_LIGHT;          // 1/c [ns/cm], riferimento 1/v
    std::vector<ObsFitResult> fit_table;         // raccolta per la tabella finale

    // Lambda che esegue i 3 fit sulle 3 osservabili di UNA pipeline e accumula
    // i risultati in fit_table. var_label include il tag per essere univoco.
    auto fitPipeline = [&](PipelineHistos& ph) {
        std::string tg = ph.tag;
        // beta: coda sinistra, riferimento 1.0
        fit_table.push_back( FitObservable(
            ph.h_beta, BETA_LO, BETA_HI,
            Form("beta_%s", tg.c_str()),
            ph.label.c_str(), "beta", "",
            /*left_tail=*/true,  /*ref=*/1.0, fout) );
        // TOF: coda destra, nessun riferimento (NAN)
        fit_table.push_back( FitObservable(
            ph.h_tof, TOF_LO, TOF_HI,
            Form("tof_%s", tg.c_str()),
            ph.label.c_str(), "TOF", "ns",
            /*left_tail=*/false, /*ref=*/NAN, fout) );
        // 1/v: coda destra, riferimento 1/c
        fit_table.push_back( FitObservable(
            ph.h_invv, INVV_LO, INVV_HI,
            Form("invv_%s", tg.c_str()),
            ph.label.c_str(), "1/v", "ns/cm",
            /*left_tail=*/false, /*ref=*/INV_C, fout) );
    };

    // Le 4 pipeline comparabili.
    for (int p = 0; p < N_PIPES; p++) fitPipeline(pipes[p]);
    // Cconst (accessoria): fittata e aggiunta in coda, separata nella tabella.
    int idx_cconst_start = (int)fit_table.size();
    fitPipeline(pipe_cconst);

    // ======================================================================
    //  TABELLA RIEPILOGATIVA DEI FIT
    // ======================================================================
    //  Stampa allineata: per ogni (pipeline, osservabile) i valori centrali dei
    //  3 metodi con errore, la sigma del picco (gauss/EMG) e la mediana KDE.
    //  Le 4 pipeline comparabili sono stampate per prime; Cconst e' separata da
    //  una riga di demarcazione (e' una diagnostica di sistematica, non una
    //  delle 4 pipeline ufficiali).
    auto printFitTable = [&]() {
        std::cout << "\n";
        std::cout << "==========================================================="
                  << "===========================" << std::endl;
        std::cout << "  TABELLA RIEPILOGATIVA DEI FIT  (valore centrale #pm errore)"
                  << std::endl;
        std::cout << "==========================================================="
                  << "===========================" << std::endl;
        // Intestazione colonne.
        std::cout << "  "
                  << std::left  << std::setw(13) << "Pipeline"
                  << std::left  << std::setw(7)  << "Oss."
                  << std::right << std::setw(20) << "Gauss-picco"
                  << std::right << std::setw(20) << "EMG (moda)"
                  << std::right << std::setw(20) << "RooKeys (moda)"
                  << std::right << std::setw(22) << "KDE mediana"
                  << std::endl;
        std::cout << "  "
                  << std::string(13 + 7 + 20 + 20 + 20 + 22, '-') << std::endl;

        auto printRow = [&](const ObsFitResult& r) {
// Formattatori compatti "val±err" a precisione adattiva (vedi FmtVE):
            // mostrano sempre ~2 cifre significative sull'errore, risolvendo il
            // caso 1/v (~0.035, err ~3e-5) che con "%.4f" dava "0.0000". Per la
            // tabella usiamo il segno "+-" ASCII (niente LaTeX) e niente spazi.
            auto fmtVE = [](bool ok, double v, double e) -> std::string {
                if (!ok) return std::string("---");
                int nd = 4;
                if (e > 0.0) {
                    int k = (int)std::floor(std::log10(e));
                    int need = -k + 1;
                    if (need > nd) nd = need;
                    if (nd > 8) nd = 8;
                }
                char b[96];
                snprintf(b, sizeof(b), "%.*f+-%.*f", nd, v, nd, e);
                return std::string(b);
            };
            auto fmtV = [](bool ok, double v, double e) -> std::string {
                if (!ok) return std::string("---");
                int nd = 4;
                if (e > 0.0) {
                    int k = (int)std::floor(std::log10(e));
                    int need = -k + 1;
                    if (need > nd) nd = need;
                    if (nd > 8) nd = 8;
                }
                char b[96];
                snprintf(b, sizeof(b), "%.*f+-%.*f", nd, v, nd, e);
                return std::string(b);
            };
            std::cout << "  "
                      << std::left  << std::setw(13) << r.pipeline.c_str()
                      << std::left  << std::setw(7)  << r.observable.c_str()
                      << std::right << std::setw(20)
                      << fmtVE(r.gpeak.fit_ok, r.gpeak.center, r.gpeak.center_err).c_str()
                      << std::right << std::setw(20)
                      << fmtVE(r.emg.fit_ok,   r.emg.center,   r.emg.center_err).c_str()
                      << std::right << std::setw(20)
                      << fmtVE(r.rkeys.fit_ok, r.rkeys.center, r.rkeys.center_err).c_str()
                    << std::right << std::setw(22)
                      << fmtV(r.rkeys.fit_ok, r.rkeys_median,
                              r.rkeys_median_err).c_str()
                      << std::endl;
        };

        // 4 pipeline comparabili.
        for (int i = 0; i < idx_cconst_start; i++) {
            printRow(fit_table[i]);
            // Riga vuota tra una pipeline e l'altra (ogni 3 osservabili).
            if ((i + 1) % 3 == 0 && (i + 1) < idx_cconst_start)
                std::cout << std::endl;
        }
        // Separatore + Cconst (diagnostica di sistematica).
        std::cout << "  " << std::string(13 + 7 + 20 + 20 + 20 + 22, '=')
                  << std::endl;
        std::cout << "  [diagnostica sistematica C(x): pipeline accessoria]"
                  << std::endl;
        for (int i = idx_cconst_start; i < (int)fit_table.size(); i++)
            printRow(fit_table[i]);

        std::cout << "==========================================================="
                  << "===========================" << std::endl;
        std::cout << "  Note: le sigma del picco (gauss/EMG) e i chi2/ndf sono nei"
                  << " canvas Canvas_fit_*." << std::endl;
        std::cout << "  beta: coda sinistra (EMG riflessa); TOF e 1/v: coda destra."
                  << std::endl;
        std::cout << "  Per 1/v il riferimento fisico e' 1/c = " << INV_C
                  << " ns/cm (v=c)." << std::endl;
        std::cout << "==========================================================="
                  << "===========================" << std::endl;
    };
    printFitTable();

    fout->Close();

    // ======================================================================
    //  RIEPILOGO FINALE
    // ======================================================================
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  ANALISI TOF v9 COMPLETATA                   " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  File processati:           " << file_list.size() << std::endl;
    std::cout << "  Eventi totali:             " << total_events << std::endl;
    std::cout << "  Eventi buoni (hybrid_tot): " << total_good << std::endl;
    std::cout << "  Eventi buoni (noclip):     " << total_good_noclip << std::endl;
    std::cout << "  Eventi buoni (Cconst):     " << total_good_Cconst << std::endl;
    std::cout << "  Eventi buoni (rise_corr):  " << total_good_risecorr << std::endl;
    std::cout << "  rise_corr nel range lineare [" << RISE_LINEAR_X_LO
              << "," << RISE_LINEAR_X_HI << "] cm: "
              << total_good_risecorr_lin << std::endl;
    std::cout << "  Parametri usati:" << std::endl;
    std::cout << "    h       = " << PAR_H << " cm" << std::endl;
    std::cout << "    x_PMT3  = " << PAR_X_PMT3 << " cm" << std::endl;
    std::cout << "    m (cal) = " << CAL_M << " ns/cm" << std::endl;
    std::cout << "    q (cal) = " << CAL_Q << " ns" << std::endl;
    std::cout << "    C(x):   "
              << (gC_poly_loaded  ? "QUADRATICO (principale)" :
                  gC_lin_loaded   ? "LINEARE (fallback)"      :
                                    "INTERPOLAZIONE (backup)") << std::endl;
    if (gC_poly_loaded)
        std::cout << "            C(x) = " << gC_poly_p0 << " + ("
                  << gC_poly_p1 << ")*x + (" << gC_poly_p2 << ")*x^2" << std::endl;
    if (gC_const_loaded)
        std::cout << "    <C> costante = " << gC_const_val << " +/- "
                  << gC_const_err << " ns" << std::endl;
    if (gRiseCorrCal_loaded)
        std::cout << "    m_corr/q_corr = " << CAL_M_CORR << " ns/cm, "
                  << CAL_Q_CORR << " ns" << std::endl;
    if (gC_corr_poly_loaded)
        std::cout << "    C_corr(x) = " << gC_corr_poly_p0 << " + ("
                  << gC_corr_poly_p1 << ")*x + (" << gC_corr_poly_p2
                  << ")*x^2" << std::endl;
std::cout << "  Output: " << outname << std::endl;
    std::cout << "=============================================" << std::endl;
}


// ==========================================================================
//  SEZIONE 8b: ANALISI TOF CON FILTRO DI DECADIMENTO FIFO (PIOMBI)
// ==========================================================================
//
//  Lead_analysis(): variante di TOF_Analysis() per il run con i PIOMBI.
//
//  Esegue ESATTAMENTE la stessa ricostruzione e gli stessi tagli di
//  TOF_Analysis() (4 pipeline: hybrid_tot, noclip, Cconst, rise_corr; veto
//  offline PMT4; cache XML->ROOT), ma in piu' applica il FILTRO DI DECADIMENTO
//  basato sul FIFO della DE10-Nano:
//    - per ogni file .txt costruisce le coppie di decadimento START->STOP;
//    - per ogni evento DRS verifica se il suo timestamp coincide (entro la
//      finestra di matching) con lo START di una coppia non ancora usata;
//    - marca l'evento con fifo_ok = 1 (decadimento confermato) o 0.
//
//  Ogni file .xml e' accoppiato POSIZIONALMENTE al rispettivo .txt: il primo
//  xml col primo txt, il secondo col secondo, ecc. Questo e' essenziale perche'
//  il matching temporale e' valido solo tra DRS e FIFO della STESSA acquisizione
//  (stessa origine temporale, stesso orologio FPGA).
//
//  L'output ROOT contiene, per ogni grandezza, l'istogramma su TUTTI gli eventi
//  buoni (riferimento, suffisso _all) e quello sui soli eventi di decadimento
//  (suffisso _lead). Il TTree salva tutti gli eventi con i branch diagnostici
//  del FIFO (fifo_ok, fifo_matched_dt, fifo_dt_decay_us, drs_abs_s) per il
//  tuning a posteriori della finestra di matching.
//
//  ARGOMENTI
//    xml_files : file .xml DRS separati da virgola (uno o piu')
//    fifo_files: file .txt FIFO separati da virgola, STESSO NUMERO e STESSO
//                ORDINE degli .xml (accoppiamento posizionale)
//    outname   : nome del file ROOT di output
//    cal_file  : (opzionale) file ROOT di calibrazione (retta, C(x), A(TOT))
//
//  UTILIZZO
//    root -l 'TOF_Analysis_v10.cpp'
//    Lead_analysis("run1.xml,run2.xml", "fifo1.txt,fifo2.txt",
//                  "tof_lead.root", "cal_v13.root")
// ==========================================================================
void Lead_analysis(const char* xml_files,
                   const char* fifo_files,
                   const char* outname  = "TOF_Lead_output.root",
                   const char* cal_file = "") {

    gROOT->SetBatch(kTRUE);

    std::cout << "=============================================" << std::endl;
    std::cout << "  TOF Analysis con FILTRO DECADIMENTO (PIOMBI)" << std::endl;
    std::cout << "=============================================" << std::endl;

    // --- Caricamento della calibrazione (retta, C(x), polinomio A(TOT)) ---
    if (strlen(cal_file) > 0) {
        std::cout << "[INFO] Caricamento calibrazione da: " << cal_file << std::endl;
        LoadCalibrationFromFile(cal_file);
    } else if (!gC_poly_loaded && !gC_lin_loaded) {
        std::cerr << "[WARNING] Nessuna calibrazione fornita/caricata: parametri "
                  << "hardcoded, recupero clippati DISABILITATO." << std::endl;
    }

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "  Geometria:    h = " << PAR_H << " cm, x_PMT3 = "
              << PAR_X_PMT3 << " cm" << std::endl;
    std::cout << "  Retta cal.:   m = " << CAL_M << " ns/cm, q = " << CAL_Q
              << " ns  (v_eff = " << 2.0 / fabs(CAL_M) << " cm/ns)" << std::endl;
    std::cout << "  Finestra match DRS<->FIFO: +-" << FIFO_MATCH_WINDOW_S * 1e3
              << " ms" << std::endl;
    std::cout << "  GATE decadimento: [" << FIFO_DECAY_DT_MIN << ", "
              << FIFO_DECAY_DT_MAX << "] clock = ["
              << FIFO_DECAY_DT_MIN * FIFO_CLK_PERIOD_S * 1e6 << ", "
              << FIFO_DECAY_DT_MAX * FIFO_CLK_PERIOD_S * 1e6 << "] us" << std::endl;
    std::cout << "  Offset fuso DRS-FPGA: " << FIFO_TZ_OFFSET_S << " s" << std::endl;
    std::cout << "=============================================" << std::endl;

    // ----------------------------------------------------------------------
    //  PARSING E ACCOPPIAMENTO POSIZIONALE DELLE LISTE .xml / .txt
    // ----------------------------------------------------------------------
    auto split_csv = [](const char* s) -> std::vector<std::string> {
        std::vector<std::string> out;
        std::string str(s);
        std::istringstream ss(str);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            size_t a = tok.find_first_not_of(" \t");
            size_t b = tok.find_last_not_of(" \t");
            if (a != std::string::npos) out.push_back(tok.substr(a, b - a + 1));
        }
        return out;
    };

    std::vector<std::string> xml_list  = split_csv(xml_files);
    std::vector<std::string> fifo_list = split_csv(fifo_files);

    if (xml_list.empty()) {
        std::cerr << "[ERRORE] Nessun file XML specificato." << std::endl;
        return;
    }
    // CONTROLLO RICHIESTO: stesso numero di .xml e .txt.
    if (xml_list.size() != fifo_list.size()) {
        std::cerr << "[ERRORE] Numero di file XML (" << xml_list.size()
                  << ") diverso dal numero di file FIFO (" << fifo_list.size()
                  << "). Ogni .xml deve avere il suo .txt (accoppiamento "
                  << "posizionale)." << std::endl;
        return;
    }
    std::cout << "[INFO] " << xml_list.size()
              << " coppie (XML, FIFO) da processare:" << std::endl;
    for (size_t i = 0; i < xml_list.size(); i++)
        std::cout << "   [" << i << "]  " << xml_list[i]
                  << "  <->  " << fifo_list[i] << std::endl;

    // --- Setup grafico (identico a TOF_Analysis) ---
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleSize(0.05, "t");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetTitleSize(0.045, "xy");

    TFile *fout = new TFile(outname, "RECREATE");
    if (!fout->IsOpen()) {
        std::cerr << "[ERRORE] Impossibile creare " << outname << std::endl;
        return;
    }

    // ======================================================================
    //  TTree DI OUTPUT "tof_data" — stesso layout di TOF_Analysis + branch FIFO
    // ======================================================================
    TTree *tree = new TTree("tof_data",
                            "Dati TOF evento per evento (run piombi, filtro FIFO)");

    // ---- Pipeline hybrid_tot (principale, senza suffisso) ----
    Float_t t_cfd[3], amp[3];
    tree->Branch("t_cfd", t_cfd, "t_cfd[3]/F");
    tree->Branch("amp",   amp,   "amp[3]/F");
    Float_t dt12, x_imp, T_meas, C_used, tof;
    Float_t path_len, theta, vel, beta, inv_vel;
    Int_t   good;
    tree->Branch("dt12",     &dt12,     "dt12/F");
    tree->Branch("x_imp",    &x_imp,    "x_imp/F");
    tree->Branch("T_meas",   &T_meas,   "T_meas/F");
    tree->Branch("C_used",   &C_used,   "C_used/F");
    tree->Branch("tof",      &tof,      "tof/F");
    tree->Branch("path_len", &path_len, "path_len/F");
    tree->Branch("theta",    &theta,    "theta/F");
    tree->Branch("vel",      &vel,      "vel/F");
    tree->Branch("beta",     &beta,     "beta/F");
    tree->Branch("inv_vel",  &inv_vel,  "inv_vel/F");
    tree->Branch("good",     &good,     "good/I");

    // ---- Pipeline noclip ----
    Float_t dt12_noclip, x_imp_noclip, T_meas_noclip, C_used_noclip, tof_noclip;
    Float_t path_len_noclip, theta_noclip, vel_noclip, beta_noclip, inv_vel_noclip;
    Int_t   good_noclip;
    tree->Branch("dt12_noclip",     &dt12_noclip,     "dt12_noclip/F");
    tree->Branch("x_imp_noclip",    &x_imp_noclip,    "x_imp_noclip/F");
    tree->Branch("T_meas_noclip",   &T_meas_noclip,   "T_meas_noclip/F");
    tree->Branch("C_used_noclip",   &C_used_noclip,   "C_used_noclip/F");
    tree->Branch("tof_noclip",      &tof_noclip,      "tof_noclip/F");
    tree->Branch("path_len_noclip", &path_len_noclip, "path_len_noclip/F");
    tree->Branch("theta_noclip",    &theta_noclip,    "theta_noclip/F");
    tree->Branch("vel_noclip",      &vel_noclip,      "vel_noclip/F");
    tree->Branch("beta_noclip",     &beta_noclip,     "beta_noclip/F");
    tree->Branch("inv_vel_noclip",  &inv_vel_noclip,  "inv_vel_noclip/F");
    tree->Branch("good_noclip",     &good_noclip,     "good_noclip/I");

    // ---- Pipeline Cconst ----
    Float_t tof_Cconst, C_used_Cconst, beta_Cconst, inv_vel_Cconst, vel_Cconst;
    Int_t   good_Cconst;
    tree->Branch("tof_Cconst",     &tof_Cconst,     "tof_Cconst/F");
    tree->Branch("C_used_Cconst",  &C_used_Cconst,  "C_used_Cconst/F");
    tree->Branch("beta_Cconst",    &beta_Cconst,    "beta_Cconst/F");
    tree->Branch("inv_vel_Cconst", &inv_vel_Cconst, "inv_vel_Cconst/F");
    tree->Branch("vel_Cconst",     &vel_Cconst,     "vel_Cconst/F");
    tree->Branch("good_Cconst",    &good_Cconst,    "good_Cconst/I");

    // ---- Pipeline rise_corr ----
    Float_t dt12_risecorr, x_imp_risecorr, T_meas_risecorr, C_used_risecorr;
    Float_t tof_risecorr, path_len_risecorr, theta_risecorr;
    Float_t vel_risecorr, beta_risecorr, inv_vel_risecorr;
    Int_t   good_risecorr, good_risecorr_linrange;
    tree->Branch("dt12_risecorr",     &dt12_risecorr,     "dt12_risecorr/F");
    tree->Branch("x_imp_risecorr",    &x_imp_risecorr,    "x_imp_risecorr/F");
    tree->Branch("T_meas_risecorr",   &T_meas_risecorr,   "T_meas_risecorr/F");
    tree->Branch("C_used_risecorr",   &C_used_risecorr,   "C_used_risecorr/F");
    tree->Branch("tof_risecorr",      &tof_risecorr,      "tof_risecorr/F");
    tree->Branch("path_len_risecorr", &path_len_risecorr, "path_len_risecorr/F");
    tree->Branch("theta_risecorr",    &theta_risecorr,    "theta_risecorr/F");
    tree->Branch("vel_risecorr",      &vel_risecorr,      "vel_risecorr/F");
    tree->Branch("beta_risecorr",     &beta_risecorr,     "beta_risecorr/F");
    tree->Branch("inv_vel_risecorr",  &inv_vel_risecorr,  "inv_vel_risecorr/F");
    tree->Branch("good_risecorr",     &good_risecorr,     "good_risecorr/I");
    tree->Branch("good_risecorr_linrange",
                 &good_risecorr_linrange, "good_risecorr_linrange/I");

    // ---- Diagnostici clipping / rise time ----
    Int_t clipped_chan[3], recov_tot[3], file_idx;
    Float_t rise_time[3], t10_frac[3], t30_frac[3];
    Int_t   rise_time_ok[3], t1030_ok[3];
    tree->Branch("clipped_chan", clipped_chan, "clipped_chan[3]/I");
    tree->Branch("recov_tot",    recov_tot,    "recov_tot[3]/I");
    tree->Branch("file_idx",     &file_idx,    "file_idx/I");
    tree->Branch("rise_time",    rise_time,    "rise_time[3]/F");
    tree->Branch("t10_frac",     t10_frac,     "t10_frac[3]/F");
    tree->Branch("t30_frac",     t30_frac,     "t30_frac[3]/F");
    tree->Branch("rise_time_ok", rise_time_ok, "rise_time_ok[3]/I");
    tree->Branch("t1030_ok",     t1030_ok,     "t1030_ok[3]/I");

    // ---- Diagnostici veto PMT4 ----
    Float_t veto_amp, veto_t_peak, veto_sig, veto_thr;
    Float_t veto_baseline, veto_baseline_rms;
    Int_t   veto_nsamp, veto_fired, veto_avail;
    tree->Branch("veto_amp",          &veto_amp,          "veto_amp/F");
    tree->Branch("veto_t_peak",       &veto_t_peak,       "veto_t_peak/F");
    tree->Branch("veto_sig",          &veto_sig,          "veto_sig/F");
    tree->Branch("veto_thr",          &veto_thr,          "veto_thr/F");
    tree->Branch("veto_baseline",     &veto_baseline,     "veto_baseline/F");
    tree->Branch("veto_baseline_rms", &veto_baseline_rms, "veto_baseline_rms/F");
    tree->Branch("veto_nsamp",        &veto_nsamp,        "veto_nsamp/I");
    tree->Branch("veto_fired",        &veto_fired,        "veto_fired/I");
    tree->Branch("veto_avail",        &veto_avail,        "veto_avail/I");

    // ---- Diagnostici FIFO (NUOVI) ----
    //   fifo_ok          = 1 se l'evento e' un decadimento confermato dal FIFO
    //   fifo_matched_dt  = |t_DRS - t_START| della coppia associata [ms] (-1 se no match)
    //   fifo_dt_decay_us = dt START->STOP della coppia associata [us] (-1 se no match)
    //   drs_abs_s        = istante DRS sull'asse FPGA [s] (per diagnostica/tuning)
    Int_t   fifo_ok;
    Float_t fifo_matched_dt, fifo_dt_decay_us;
    Double_t drs_abs_s_branch;
    tree->Branch("fifo_ok",          &fifo_ok,          "fifo_ok/I");
    tree->Branch("fifo_matched_dt",  &fifo_matched_dt,  "fifo_matched_dt/F");
    tree->Branch("fifo_dt_decay_us", &fifo_dt_decay_us, "fifo_dt_decay_us/F");
    tree->Branch("drs_abs_s",        &drs_abs_s_branch, "drs_abs_s/D");

    // ======================================================================
    //  ISTOGRAMMI — per ogni grandezza: "_all" (riferimento) e "_lead" (FIFO)
    // ======================================================================
    // Macro-aiuto per creare la coppia di TH1D con stesso binning.
    auto mkTH1 = [&](const char* base, const char* suf, const char* titax,
                     int nb, double lo, double hi) -> TH1D* {
        return new TH1D(Form("%s%s", base, suf), titax, nb, lo, hi);
    };

    // Pipeline hybrid_tot: TOF, beta, beta_norm, 1/v, theta, x_imp (all & lead)
    TH1D *h_tof_all   = mkTH1("h_tof",  "_all",  "TOF (tutti);TOF [ns];Conteggi",   NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_tof_lead  = mkTH1("h_tof",  "_lead", "TOF (decadimenti);TOF [ns];Conteggi", NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta_all  = mkTH1("h_beta", "_all",  "#beta (tutti);#beta;Conteggi",    NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_lead = mkTH1("h_beta", "_lead", "#beta (decadimenti);#beta;Conteggi", NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_norm_all  = mkTH1("h_beta_norm", "_all",  "PDF #beta (tutti);#beta;PDF",  NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_norm_lead = mkTH1("h_beta_norm", "_lead", "PDF #beta (decadimenti);#beta;PDF", NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_invv_all  = mkTH1("h_invv", "_all",  "1/v (tutti);1/v [ns/cm];Conteggi", NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_invv_lead = mkTH1("h_invv", "_lead", "1/v (decadimenti);1/v [ns/cm];Conteggi", NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_theta_all  = mkTH1("h_theta", "_all",  "#theta (tutti);#theta [rad];Conteggi", NBINS_THETA, THETA_LO, THETA_HI);
    TH1D *h_theta_lead = mkTH1("h_theta", "_lead", "#theta (decadimenti);#theta [rad];Conteggi", NBINS_THETA, THETA_LO, THETA_HI);
    TH1D *h_x_imp_all  = mkTH1("h_x_imp", "_all",  "x_{imp} (tutti);x_{imp} [cm];Conteggi", NBINS_X, X_LO, X_HI);
    TH1D *h_x_imp_lead = mkTH1("h_x_imp", "_lead", "x_{imp} (decadimenti);x_{imp} [cm];Conteggi", NBINS_X, X_LO, X_HI);
    TH1D *h_path_all  = mkTH1("h_path", "_all",  "Distanza (tutti);l [cm];Conteggi", 100, 90, 250);
    TH1D *h_path_lead = mkTH1("h_path", "_lead", "Distanza (decadimenti);l [cm];Conteggi", 100, 90, 250);
    TH1D *h_T_meas_all  = mkTH1("h_T_meas", "_all",  "T_{meas} (tutti);T_{meas} [ns];Conteggi", 150, -5, 25);
    TH1D *h_T_meas_lead = mkTH1("h_T_meas", "_lead", "T_{meas} (decadimenti);T_{meas} [ns];Conteggi", 150, -5, 25);
    TH2D *h2_tof_x_all  = new TH2D("h2_tof_x_all",  "TOF vs x (tutti);x_{imp} [cm];TOF [ns]", 70, X_LO, X_HI, 75, TOF_LO, TOF_HI);
    TH2D *h2_tof_x_lead = new TH2D("h2_tof_x_lead", "TOF vs x (decadimenti);x_{imp} [cm];TOF [ns]", 70, X_LO, X_HI, 75, TOF_LO, TOF_HI);
    TH2D *h2_beta_x_all  = new TH2D("h2_beta_x_all",  "#beta vs x (tutti);x_{imp} [cm];#beta", 70, X_LO, X_HI, 100, BETA_LO, BETA_HI);
    TH2D *h2_beta_x_lead = new TH2D("h2_beta_x_lead", "#beta vs x (decadimenti);x_{imp} [cm];#beta", 70, X_LO, X_HI, 100, BETA_LO, BETA_HI);

    // Pipeline noclip (solo all & lead per le grandezze chiave)
    TH1D *h_tof_nc_all   = mkTH1("h_tof_noclip",  "_all",  "TOF noclip (tutti);TOF [ns];Conteggi", NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_tof_nc_lead  = mkTH1("h_tof_noclip",  "_lead", "TOF noclip (decadimenti);TOF [ns];Conteggi", NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta_nc_all  = mkTH1("h_beta_noclip", "_all",  "#beta noclip (tutti);#beta;Conteggi", NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_nc_lead = mkTH1("h_beta_noclip", "_lead", "#beta noclip (decadimenti);#beta;Conteggi", NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_invv_nc_all  = mkTH1("h_invv_noclip", "_all",  "1/v noclip (tutti);1/v [ns/cm];Conteggi", NBINS_INVV, INVV_LO, INVV_HI);
    TH1D *h_invv_nc_lead = mkTH1("h_invv_noclip", "_lead", "1/v noclip (decadimenti);1/v [ns/cm];Conteggi", NBINS_INVV, INVV_LO, INVV_HI);

    // Pipeline Cconst (all & lead, grandezze chiave)
    TH1D *h_tof_cc_all   = mkTH1("h_tof_Cconst",  "_all",  "TOF Cconst (tutti);TOF [ns];Conteggi", NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_tof_cc_lead  = mkTH1("h_tof_Cconst",  "_lead", "TOF Cconst (decadimenti);TOF [ns];Conteggi", NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta_cc_all  = mkTH1("h_beta_Cconst", "_all",  "#beta Cconst (tutti);#beta;Conteggi", NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_cc_lead = mkTH1("h_beta_Cconst", "_lead", "#beta Cconst (decadimenti);#beta;Conteggi", NBINS_BETA, BETA_LO, BETA_HI);

    // Pipeline rise_corr (all & lead, grandezze chiave)
    TH1D *h_tof_rc_all   = mkTH1("h_tof_risecorr",  "_all",  "TOF rise_corr (tutti);TOF [ns];Conteggi", NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_tof_rc_lead  = mkTH1("h_tof_risecorr",  "_lead", "TOF rise_corr (decadimenti);TOF [ns];Conteggi", NBINS_TOF, TOF_LO, TOF_HI);
    TH1D *h_beta_rc_all  = mkTH1("h_beta_risecorr", "_all",  "#beta rise_corr (tutti);#beta;Conteggi", NBINS_BETA, BETA_LO, BETA_HI);
    TH1D *h_beta_rc_lead = mkTH1("h_beta_risecorr", "_lead", "#beta rise_corr (decadimenti);#beta;Conteggi", NBINS_BETA, BETA_LO, BETA_HI);

    // Diagnostici veto PMT4 (come in TOF_Analysis)
    TH1D *h_veto_amp_all = new TH1D("h_veto_amp_all",
        "Ampiezza CH4 (finestra veto);A_{4} [mV];Conteggi", 300, 0, 60);
    TH1D *h_veto_sig_all = new TH1D("h_veto_sig_all",
        "Significativit#grave{a} CH4;A_{4}/#sigma_{bl};Conteggi", 300, 0, 60);

    // Diagnostici FIFO: distribuzioni per il tuning della finestra
    TH1D *h_fifo_matched_dt = new TH1D("h_fifo_matched_dt",
        "Scarto temporale DRS-FIFO degli eventi matchati;|#Deltat| [ms];Conteggi",
        200, 0, FIFO_MATCH_WINDOW_S * 1e3);
    TH1D *h_fifo_dt_decay = new TH1D("h_fifo_dt_decay",
        "Tempo di decadimento START#rightarrowSTOP (coppie matchate);#Deltat [#mus];Conteggi",
        100, 0, FIFO_DECAY_DT_MAX * FIFO_CLK_PERIOD_S * 1e6);

    const double tof_min_phys = PAR_H / C_LIGHT;
    std::cout << "  TOF minimo fisico (h/c) = " << tof_min_phys << " ns" << std::endl;
    std::cout << "=============================================" << std::endl;

    // ======================================================================
    //  CONTATORI
    // ======================================================================
    int total_events     = 0;
    int total_good       = 0, total_good_lead       = 0;   // hybrid_tot
    int total_good_nc    = 0, total_good_nc_lead    = 0;   // noclip
    int total_good_cc    = 0, total_good_cc_lead    = 0;   // Cconst
    int total_good_rc    = 0, total_good_rc_lead    = 0;   // rise_corr
    int rej_no_cfd = 0, rej_oscillating = 0, rej_clipped = 0;
    int rej_x_out = 0, rej_tof_unphys = 0, rej_veto = 0;
    int n_veto_available = 0;
    int n_fifo_pairs_tot = 0;     // coppie di decadimento totali su tutti i .txt
    int n_fifo_matched   = 0;     // eventi DRS con match FIFO

    // ======================================================================
    //  LOOP SULLE COPPIE (XML, FIFO)
    // ======================================================================
    for (size_t fi = 0; fi < xml_list.size(); fi++) {

        // --- 1) Carica e prepara il FIFO di questa acquisizione ---
        FifoData fd;
        if (!ParseFifoStartTime(fifo_list[fi], fd)) {
            std::cerr << "[WARNING] FIFO " << fifo_list[fi]
                      << " : orario di start non leggibile, coppia saltata."
                      << std::endl;
            continue;
        }
        int npairs = LoadFifoAndBuildPairs(fifo_list[fi], fd, /*discard_buffer0=*/true);
        if (npairs < 0) {
            std::cerr << "[WARNING] FIFO " << fifo_list[fi]
                      << " : caricamento fallito, coppia saltata." << std::endl;
            continue;
        }
        n_fifo_pairs_tot += npairs;

        // --- 2) Carica gli eventi DRS (con cache) ---
        std::vector<EventData> events;
        int n_parsed = LoadOrParseDataset(xml_list[fi].c_str(), events);
        if (n_parsed <= 0) {
            std::cerr << "[WARNING] Nessun evento in " << xml_list[fi] << std::endl;
            continue;
        }

        // ------------------------------------------------------------------
        //  LOOP SUGLI EVENTI DEL FILE
        // ------------------------------------------------------------------
        for (size_t ev = 0; ev < events.size(); ev++) {

            EventData &e = events[ev];
            total_events++;
            file_idx = (Int_t)fi;

            // ---- Default sentinella di TUTTI i branch ----
            for (int k = 0; k < 3; k++) {
                t_cfd[k] = -999.0f; amp[k] = 0.0f;
                clipped_chan[k] = 0; recov_tot[k] = 0;
                rise_time[k] = -999.0f; t10_frac[k] = -999.0f; t30_frac[k] = -999.0f;
                rise_time_ok[k] = 0; t1030_ok[k] = 0;
            }
            dt12 = x_imp = T_meas = C_used = tof = -999.0f;
            path_len = theta = vel = beta = inv_vel = -999.0f; good = 0;
            dt12_noclip = x_imp_noclip = T_meas_noclip = -999.0f;
            C_used_noclip = tof_noclip = -999.0f;
            path_len_noclip = theta_noclip = -999.0f;
            vel_noclip = beta_noclip = inv_vel_noclip = -999.0f; good_noclip = 0;
            tof_Cconst = C_used_Cconst = -999.0f;
            beta_Cconst = inv_vel_Cconst = vel_Cconst = -999.0f; good_Cconst = 0;
            dt12_risecorr = x_imp_risecorr = T_meas_risecorr = -999.0f;
            C_used_risecorr = tof_risecorr = -999.0f;
            path_len_risecorr = theta_risecorr = -999.0f;
            vel_risecorr = beta_risecorr = inv_vel_risecorr = -999.0f;
            good_risecorr = 0; good_risecorr_linrange = 0;
            veto_amp = veto_t_peak = veto_sig = veto_thr = -999.0f;
            veto_baseline = veto_baseline_rms = -999.0f;
            veto_nsamp = 0; veto_fired = 0; veto_avail = 0;
            fifo_ok = 0; fifo_matched_dt = -1.0f; fifo_dt_decay_us = -1.0f;
            drs_abs_s_branch = -1.0;

            // ============================================================
            //  GATE FIFO — calcolato per OGNI evento (marca, non scarta qui)
            // ============================================================
            // Il timestamp DRS (campo <Time>) viene portato sull'asse FPGA della
            // STESSA acquisizione (origine = inizio giorno di avvio del FIFO fd).
            // Poi si cerca una coppia di decadimento il cui START coincida entro
            // la finestra di matching. fifo_ok = 1 se trovata.
            double drs_abs = -1.0;
            if (ParseDRSTimestampAbs(e.timestamp, fd.start_day_index, drs_abs)) {
                drs_abs_s_branch = drs_abs;
                double matched_dt = -1.0;
                int idx = MatchDRStoFifo(fd, drs_abs, matched_dt);
                if (idx >= 0) {
                    fifo_ok          = 1;
                    fifo_matched_dt  = (Float_t)(matched_dt * 1e3);   // ms
                    fifo_dt_decay_us = (Float_t)(fd.pairs[idx].dt_clk
                                                 * FIFO_CLK_PERIOD_S * 1e6); // us
                    n_fifo_matched++;
                    h_fifo_matched_dt->Fill(fifo_matched_dt);
                    h_fifo_dt_decay->Fill(fifo_dt_decay_us);
                }
            }
            // NOTA: NON facciamo "continue" qui. L'evento prosegue nella pipeline
            // come in TOF_Analysis; gli istogrammi _lead vengono riempiti solo se
            // fifo_ok==1. Cosi' il TTree contiene tutti gli eventi (con il flag)
            // e si possono confrontare direttamente "tutti" vs "decadimenti".

            // ============================================================
            //  TAGLIO 0: almeno 3 canali con impulso e CFD valido
            // ============================================================
            if (e.nchannels < 3) { tree->Fill(); rej_no_cfd++; continue; }
            bool all_cfd_ok = true;
            for (int k = 0; k < 3; k++) {
                t_cfd[k] = (Float_t)e.ch[k].t_cfd;
                amp[k]   = (Float_t)e.ch[k].amplitude;
                rise_time[k]    = e.ch[k].rise_time_ok ? (Float_t)e.ch[k].rise_time : -999.0f;
                t10_frac[k]     = e.ch[k].t10_ok       ? (Float_t)e.ch[k].t_10      : -999.0f;
                t30_frac[k]     = e.ch[k].t30_ok       ? (Float_t)e.ch[k].t_30      : -999.0f;
                rise_time_ok[k] = e.ch[k].rise_time_ok ? 1 : 0;
                t1030_ok[k]     = (e.ch[k].t10_ok && e.ch[k].t30_ok) ? 1 : 0;
                if (!e.ch[k].has_pulse || !e.ch[k].cfd_ok) all_cfd_ok = false;
            }
            if (!all_cfd_ok) { tree->Fill(); rej_no_cfd++; continue; }

            // ============================================================
            //  VETO OFFLINE PMT4
            // ============================================================
            if (e.nchannels >= 4) {
                veto_avail = 1; n_veto_available++;
                double t_mu = 0.5 * (e.ch[0].t_cfd + e.ch[1].t_cfd);
                VetoResult vr = AnalyzeVeto(e.ch[VETO_CH_INDEX], t_mu);
                veto_amp = (Float_t)vr.amplitude;  veto_t_peak = (Float_t)vr.t_peak;
                veto_sig = (Float_t)vr.significance; veto_thr = (Float_t)vr.v_thr_used;
                veto_baseline = (Float_t)vr.baseline;
                veto_baseline_rms = (Float_t)vr.baseline_rms;
                veto_nsamp = (Int_t)vr.n_samples_below;
                veto_fired = vr.has_signal ? 1 : 0;
                if (vr.window_ok) {
                    h_veto_amp_all->Fill(vr.amplitude);
                    h_veto_sig_all->Fill(vr.significance);
                }
                if (vr.has_signal) { tree->Fill(); rej_veto++; continue; }
            }

            // ============================================================
            //  TAGLIO OSCILLAZIONE
            // ============================================================
            bool any_osc = false;
            for (int k = 0; k < 3; k++)
                if (e.ch[k].is_oscillating) { any_osc = true; break; }
            if (any_osc) { tree->Fill(); rej_oscillating++; continue; }

            // ============================================================
            //  RECUPERO CLIPPATI VIA TOT POLINOMIALE
            // ============================================================
            for (int k = 0; k < 3; k++)
                if (e.ch[k].is_clipped) RecoverClippedCFD_TOT(e.ch[k], k);
            for (int k = 0; k < 3; k++) {
                clipped_chan[k] = e.ch[k].is_clipped        ? 1 : 0;
                recov_tot[k]    = e.ch[k].cfd_recovered_tot ? 1 : 0;
            }

            // ============================================================
            //  LAMBDA RunPipeline (identica a TOF_Analysis)
            // ============================================================
            auto RunPipeline = [&](const Float_t t_use[3], const bool ok_use[3],
                                   double cal_m, double cal_q, bool use_corr_C,
                                   Float_t &dt12_o,  Float_t &x_imp_o,
                                   Float_t &T_meas_o, Float_t &C_used_o,
                                   Float_t &tof_o,   Float_t &path_o,
                                   Float_t &theta_o, Float_t &vel_o,
                                   Float_t &beta_o,  Float_t &invv_o) -> int {
                dt12_o = x_imp_o = T_meas_o = C_used_o = tof_o = -999.0f;
                path_o = theta_o = vel_o = beta_o = invv_o = -999.0f;
                if (!ok_use[0] || !ok_use[1] || !ok_use[2]) return -1;
                dt12_o  = t_use[0] - t_use[1];
                x_imp_o = (Float_t)(((double)dt12_o - cal_q) / cal_m);
                if (x_imp_o < X_CUT_LO || x_imp_o > X_CUT_HI) return -2;
                T_meas_o = t_use[2] - (t_use[0] + t_use[1]) / 2.0f;
                C_used_o = (Float_t)(use_corr_C ? GetC_corr((double)x_imp_o)
                                                : GetC((double)x_imp_o));
                tof_o    = T_meas_o - C_used_o;
                if (tof_o < tof_min_phys) return -3;
                double dxp = (double)x_imp_o - PAR_X_PMT3;
                path_o  = (Float_t)sqrt(dxp * dxp + PAR_H * PAR_H);
                theta_o = (Float_t)atan2(dxp, PAR_H);
                if (tof_o > 0.1f) {
                    vel_o  = (Float_t)((double)path_o / (double)tof_o);
                    beta_o = (Float_t)((double)vel_o / C_LIGHT);
                    invv_o = (Float_t)((double)tof_o / (double)path_o);
                }
                return 0;
            };

            // ---- Pipeline hybrid_tot ----
            Float_t t_hy[3]; bool ok_hy[3];
            for (int k = 0; k < 3; k++) {
                if (!e.ch[k].is_clipped) {
                    t_hy[k] = (Float_t)e.ch[k].t_cfd; ok_hy[k] = true;
                    amp[k] = (Float_t)e.ch[k].amplitude;
                } else if (e.ch[k].cfd_recovered_tot) {
                    t_hy[k] = (Float_t)e.ch[k].t_cfd_rec_tot; ok_hy[k] = true;
                    amp[k] = (Float_t)e.ch[k].amplitude_tot;
                } else {
                    t_hy[k] = -999.0f; ok_hy[k] = false;
                    amp[k] = (Float_t)e.ch[k].amplitude;
                }
                t_cfd[k] = t_hy[k];
            }
            int st_hy = RunPipeline(t_hy, ok_hy, CAL_M, CAL_Q, false,
                                    dt12, x_imp, T_meas, C_used, tof,
                                    path_len, theta, vel, beta, inv_vel);
            switch (st_hy) {
                case  0: good = 1; total_good++;     break;
                case -1: good = 0; rej_clipped++;    break;
                case -2: good = 0; rej_x_out++;      break;
                case -3: good = 0; rej_tof_unphys++; break;
                default: good = 0;                   break;
            }

            // ---- Pipeline noclip ----
            Float_t t_nc[3]; bool ok_nc[3];
            for (int k = 0; k < 3; k++) {
                if (!e.ch[k].is_clipped) { t_nc[k] = (Float_t)e.ch[k].t_cfd; ok_nc[k] = true; }
                else                     { t_nc[k] = -999.0f; ok_nc[k] = false; }
            }
            int st_nc = RunPipeline(t_nc, ok_nc, CAL_M, CAL_Q, false,
                                    dt12_noclip, x_imp_noclip, T_meas_noclip,
                                    C_used_noclip, tof_noclip, path_len_noclip,
                                    theta_noclip, vel_noclip, beta_noclip,
                                    inv_vel_noclip);
            if (st_nc == 0) { good_noclip = 1; total_good_nc++; }

            // ---- Pipeline Cconst ----
            if (good == 1 && gC_const_loaded) {
                C_used_Cconst = (Float_t)GetC_const((double)x_imp);
                tof_Cconst    = T_meas - C_used_Cconst;
                if (tof_Cconst >= tof_min_phys && tof_Cconst > 0.1f) {
                    double dxp = (double)x_imp - PAR_X_PMT3;
                    double l   = sqrt(dxp * dxp + PAR_H * PAR_H);
                    vel_Cconst     = (Float_t)(l / (double)tof_Cconst);
                    beta_Cconst    = (Float_t)((double)vel_Cconst / C_LIGHT);
                    inv_vel_Cconst = (Float_t)((double)tof_Cconst / l);
                    good_Cconst = 1; total_good_cc++;
                }
            }

            // ---- Pipeline rise_corr ----
            {
                bool no_clip_any = !e.ch[0].is_clipped && !e.ch[1].is_clipped
                                && !e.ch[2].is_clipped;
                bool t1030_all = e.ch[0].t10_ok && e.ch[0].t30_ok
                              && e.ch[1].t10_ok && e.ch[1].t30_ok
                              && e.ch[2].t10_ok && e.ch[2].t30_ok;
                if (good == 1 && no_clip_any && t1030_all) {
                    Float_t t_rc[3]; bool ok_rc[3] = { true, true, true };
                    for (int k = 0; k < 3; k++) {
                        double dT = e.ch[k].t_30 - e.ch[k].t_10;
                        t_rc[k] = (Float_t)(e.ch[k].t_cfd - dT);
                    }
                    double m_use = gRiseCorrCal_loaded ? CAL_M_CORR : CAL_M;
                    double q_use = gRiseCorrCal_loaded ? CAL_Q_CORR : CAL_Q;
                    int st_rc = RunPipeline(t_rc, ok_rc, m_use, q_use, true,
                                            dt12_risecorr, x_imp_risecorr,
                                            T_meas_risecorr, C_used_risecorr,
                                            tof_risecorr, path_len_risecorr,
                                            theta_risecorr, vel_risecorr,
                                            beta_risecorr, inv_vel_risecorr);
                    if (st_rc == 0) {
                        good_risecorr = 1; total_good_rc++;
                        if (x_imp_risecorr >= RISE_LINEAR_X_LO &&
                            x_imp_risecorr <= RISE_LINEAR_X_HI)
                            good_risecorr_linrange = 1;
                    }
                }
            }

            // ============================================================
            //  RIEMPIMENTO ISTOGRAMMI — sempre "_all", e "_lead" se fifo_ok
            // ============================================================
            // Helper: riempie la coppia (all, lead) di una grandezza scalare.
            auto fillPair = [&](TH1D* h_all, TH1D* h_lead, double val) {
                h_all->Fill(val);
                if (fifo_ok == 1) h_lead->Fill(val);
            };

            if (good == 1) {
                fillPair(h_T_meas_all, h_T_meas_lead, T_meas);
                fillPair(h_x_imp_all,  h_x_imp_lead,  x_imp);
                if (tof > 0.1f) {
                    fillPair(h_tof_all,  h_tof_lead,  tof);
                    fillPair(h_path_all, h_path_lead, path_len);
                    fillPair(h_theta_all,h_theta_lead,theta);
                    h2_tof_x_all->Fill(x_imp, tof);
                    if (fifo_ok == 1) h2_tof_x_lead->Fill(x_imp, tof);
                    if (beta > 0.0f && beta < 3.0f) {
                        fillPair(h_beta_all,      h_beta_lead,      beta);
                        fillPair(h_beta_norm_all, h_beta_norm_lead, beta);
                        h2_beta_x_all->Fill(x_imp, beta);
                        if (fifo_ok == 1) h2_beta_x_lead->Fill(x_imp, beta);
                    }
                    if (inv_vel > 0.0f && inv_vel < 0.2f)
                        fillPair(h_invv_all, h_invv_lead, inv_vel);
                }
                if (fifo_ok == 1) total_good_lead++;
            }
            if (good_noclip == 1 && tof_noclip > 0.1f) {
                fillPair(h_tof_nc_all, h_tof_nc_lead, tof_noclip);
                if (beta_noclip > 0.0f && beta_noclip < 3.0f)
                    fillPair(h_beta_nc_all, h_beta_nc_lead, beta_noclip);
                if (inv_vel_noclip > 0.0f && inv_vel_noclip < 0.2f)
                    fillPair(h_invv_nc_all, h_invv_nc_lead, inv_vel_noclip);
                if (fifo_ok == 1) total_good_nc_lead++;
            }
            if (good_Cconst == 1 && tof_Cconst > 0.1f) {
                fillPair(h_tof_cc_all, h_tof_cc_lead, tof_Cconst);
                if (beta_Cconst > 0.0f && beta_Cconst < 3.0f)
                    fillPair(h_beta_cc_all, h_beta_cc_lead, beta_Cconst);
                if (fifo_ok == 1) total_good_cc_lead++;
            }
            if (good_risecorr == 1 && tof_risecorr > 0.1f) {
                fillPair(h_tof_rc_all, h_tof_rc_lead, tof_risecorr);
                if (beta_risecorr > 0.0f && beta_risecorr < 3.0f)
                    fillPair(h_beta_rc_all, h_beta_rc_lead, beta_risecorr);
                if (fifo_ok == 1) total_good_rc_lead++;
            }

            tree->Fill();
        }   // fine loop eventi
    }       // fine loop coppie (XML,FIFO)

    // ======================================================================
    //  NORMALIZZAZIONE PDF beta (all & lead)
    // ======================================================================
    auto normPDF = [&](TH1D* h) {
        if (h->Integral() > 0.0) {
            h->Scale(1.0 / h->Integral("width"));
            h->GetYaxis()->SetTitle("PDF [1/unit. #beta]");
        }
    };
    normPDF(h_beta_norm_all);
    normPDF(h_beta_norm_lead);

    // ======================================================================
    //  CANVAS DI CONFRONTO all vs lead (prodotto fisico principale)
    // ======================================================================
    // Sovrappone tutti gli eventi buoni e i soli decadimenti confermati dal
    // FIFO, per le grandezze chiave. Il confronto evidenzia l'effetto del
    // filtro di decadimento sulla distribuzione di beta (i muoni che si fermano
    // e decadono sono a bassa quantita' di moto -> beta in media piu' bassa).
    {
        TCanvas *c_cmp = new TCanvas("c_compare_all_lead",
                                     "Confronto tutti vs decadimenti", 1200, 800);
        c_cmp->Divide(2, 2);
        auto drawPair = [&](int pad, TH1D* h_all, TH1D* h_lead, const char* tit) {
            c_cmp->cd(pad); gPad->SetGrid();
            h_all->SetLineColor(kGray + 2);  h_all->SetLineWidth(2);
            h_lead->SetLineColor(kRed + 1);  h_lead->SetLineWidth(2);
            // Per confrontare la FORMA, disegno entrambi normalizzati ad area 1
            // su cloni temporanei (non tocco gli istogrammi salvati).
            TH1D* a = (TH1D*)h_all->Clone(Form("%s_cln", h_all->GetName()));
            TH1D* l = (TH1D*)h_lead->Clone(Form("%s_cln", h_lead->GetName()));
            if (a->Integral() > 0) a->Scale(1.0 / a->Integral());
            if (l->Integral() > 0) l->Scale(1.0 / l->Integral());
            double ymax = std::max(a->GetMaximum(), l->GetMaximum());
            a->SetTitle(tit); a->SetMaximum(1.18 * ymax);
            a->Draw("HIST"); l->Draw("HIST same");
            TLegend *lg = new TLegend(0.58, 0.74, 0.89, 0.89);
            lg->SetTextFont(42); lg->SetTextSize(0.034);
            lg->AddEntry(a, "tutti (norm.)", "l");
            lg->AddEntry(l, "decadimenti (norm.)", "l");
            lg->Draw();
        };
        drawPair(1, h_tof_all,  h_tof_lead,  "TOF: tutti vs decadimenti;TOF [ns];frazione");
        drawPair(2, h_beta_all, h_beta_lead, "#beta: tutti vs decadimenti;#beta;frazione");
        drawPair(3, h_invv_all, h_invv_lead, "1/v: tutti vs decadimenti;1/v [ns/cm];frazione");
        drawPair(4, h_x_imp_all,h_x_imp_lead,"x_{imp}: tutti vs decadimenti;x_{imp} [cm];frazione");
        c_cmp->Write("Canvas_Compare_all_lead");
    }

    // Canvas diagnostico FIFO: scarto di matching e tempo di decadimento.
    {
        TCanvas *c_fifo = new TCanvas("c_fifo_diag",
                                      "Diagnostica FIFO", 1100, 450);
        c_fifo->Divide(2, 1);
        c_fifo->cd(1); gPad->SetGrid();
        h_fifo_matched_dt->SetLineColor(kAzure + 1);
        h_fifo_matched_dt->SetLineWidth(2); h_fifo_matched_dt->Draw("HIST");
        c_fifo->cd(2); gPad->SetGrid();
        h_fifo_dt_decay->SetLineColor(kGreen + 2);
        h_fifo_dt_decay->SetLineWidth(2); h_fifo_dt_decay->Draw("HIST");
        c_fifo->Write("Canvas_FIFO_Diagnostics");
    }

    // ======================================================================
    //  SCRITTURA ISTOGRAMMI
    // ======================================================================
    h_tof_all->Write();   h_tof_lead->Write();
    h_beta_all->Write();  h_beta_lead->Write();
    h_beta_norm_all->Write(); h_beta_norm_lead->Write();
    h_invv_all->Write();  h_invv_lead->Write();
    h_theta_all->Write(); h_theta_lead->Write();
    h_x_imp_all->Write(); h_x_imp_lead->Write();
    h_path_all->Write();  h_path_lead->Write();
    h_T_meas_all->Write();h_T_meas_lead->Write();
    h2_tof_x_all->Write(); h2_tof_x_lead->Write();
    h2_beta_x_all->Write();h2_beta_x_lead->Write();
    h_tof_nc_all->Write();  h_tof_nc_lead->Write();
    h_beta_nc_all->Write(); h_beta_nc_lead->Write();
    h_invv_nc_all->Write(); h_invv_nc_lead->Write();
    h_tof_cc_all->Write();  h_tof_cc_lead->Write();
    h_beta_cc_all->Write(); h_beta_cc_lead->Write();
    h_tof_rc_all->Write();  h_tof_rc_lead->Write();
    h_beta_rc_all->Write(); h_beta_rc_lead->Write();
    h_veto_amp_all->Write(); h_veto_sig_all->Write();
    h_fifo_matched_dt->Write(); h_fifo_dt_decay->Write();

    fout->Close();

    // ======================================================================
    //  RIEPILOGO FINALE
    // ======================================================================
    auto pct = [&](int n) -> double {
        return (total_events > 0) ? 100.0 * n / total_events : 0.0;
    };
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  ANALISI TOF PIOMBI (filtro FIFO) COMPLETATA " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  Coppie (XML,FIFO) processate: " << xml_list.size() << std::endl;
    std::cout << "  Coppie di decadimento FIFO:   " << n_fifo_pairs_tot << std::endl;
    std::cout << "  Eventi DRS totali:            " << total_events << std::endl;
    std::cout << "  Eventi con match FIFO:        " << n_fifo_matched
              << " (" << pct(n_fifo_matched) << "%)" << std::endl;
    std::cout << "  ---- Tagli di qualita' (comuni) ----" << std::endl;
    std::cout << "  Scartati (no CFD):       " << rej_no_cfd << std::endl;
    std::cout << "  Scartati (oscillazione): " << rej_oscillating << std::endl;
    std::cout << "  Scartati (VETO PMT4):    " << rej_veto
              << "  [su " << n_veto_available << " con CH4]" << std::endl;
    std::cout << "  Scartati (clipped):      " << rej_clipped << std::endl;
    std::cout << "  Scartati (x fuori barra):" << rej_x_out << std::endl;
    std::cout << "  Scartati (TOF<h/c):      " << rej_tof_unphys << std::endl;
    std::cout << "  ---- Eventi buoni: TUTTI / DECADIMENTI ----" << std::endl;
    std::cout << "  hybrid_tot: " << total_good << " / " << total_good_lead << std::endl;
    std::cout << "  noclip:     " << total_good_nc << " / " << total_good_nc_lead << std::endl;
    std::cout << "  Cconst:     " << total_good_cc << " / " << total_good_cc_lead << std::endl;
    std::cout << "  rise_corr:  " << total_good_rc << " / " << total_good_rc_lead << std::endl;
    std::cout << "  Output: " << outname << std::endl;
    std::cout << "=============================================" << std::endl;
}


// ==========================================================================
//  SEZIONE 9: STUDIO DELLE DISTRIBUZIONI vs FRAZIONE CFD
// ==========================================================================

/// CFD_study(): ricalcola la ricostruzione TOF al variare della frazione CFD.
///
/// Le frazioni studiate sono {10, 15, 20, 30, 40, 50} %. La funzione usa
/// LoadOrParseDataset(), quindi beneficia automaticamente della cache XML->ROOT.
/// Per mantenere il confronto pulito fra frazioni usa lo stesso campione di
/// eventi non clippati e non oscillanti; sui clippati non viene applicato il
/// recupero TOT perche' il polinomio A(TOT) e' calibrato alla frazione nominale.
void CFD_study(const char* xml_files,
               const char* outname  = "TOF_CFD_study.root",
               const char* cal_file = "") {

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    if (strlen(cal_file) > 0) {
        std::cout << "[INFO] CFD_study: carico calibrazione da " << cal_file
                  << std::endl;
        LoadCalibrationFromFile(cal_file);
    }

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
        std::cerr << "[ERRORE] CFD_study: nessun file XML specificato."
                  << std::endl;
        return;
    }

    const int    NF = 6;
    const double fracs[NF] = { 0.10, 0.15, 0.20, 0.30, 0.40, 0.50 };
    const int    pct  [NF] = {   10,   15,   20,   30,   40,   50 };
    const int    cols [NF] = { kRed + 1, kOrange + 1, kGreen + 2,
                               kCyan + 2, kBlue + 1,  kViolet + 1 };
    const int    styles[NF] = { 1, 1, 1, 2, 2, 2 };

    TFile* fout = new TFile(outname, "RECREATE");
    if (!fout || !fout->IsOpen()) {
        std::cerr << "[ERRORE] CFD_study: impossibile creare " << outname
                  << std::endl;
        delete fout;
        return;
    }

    TH1D* h_beta_cfd[NF];
    TH1D* h_tof_cfd[NF];
    TH1D* h_invv_cfd[NF];
    for (int f = 0; f < NF; f++) {
        h_beta_cfd[f] = new TH1D(Form("h_beta_cfd%02d", pct[f]),
            Form("#beta (CFD %d%%);#beta;Conteggi", pct[f]),
            NBINS_BETA, BETA_LO, BETA_HI);
        h_tof_cfd[f] = new TH1D(Form("h_tof_cfd%02d", pct[f]),
            Form("TOF (CFD %d%%);TOF [ns];Conteggi", pct[f]),
            NBINS_TOF, TOF_LO, TOF_HI);
        h_invv_cfd[f] = new TH1D(Form("h_invv_cfd%02d", pct[f]),
            Form("1/v (CFD %d%%);1/v [ns/cm];Conteggi", pct[f]),
            NBINS_INVV, INVV_LO, INVV_HI);
    }

    const double tof_min_phys = PAR_H / C_LIGHT;
    int total_events = 0;
    int clean_events = 0;
    int accepted_by_frac[NF] = {0, 0, 0, 0, 0, 0};

    for (size_t fi = 0; fi < file_list.size(); fi++) {
        std::vector<EventData> events;
        int n = LoadOrParseDataset(file_list[fi].c_str(), events);
        if (n <= 0) continue;

        for (const auto& e : events) {
            total_events++;
            if (e.nchannels < 3) continue;

            const ChannelData& c1 = e.ch[0];
            const ChannelData& c2 = e.ch[1];
            const ChannelData& c3 = e.ch[2];

            if (!c1.has_pulse || !c2.has_pulse || !c3.has_pulse) continue;
            if (c1.is_clipped || c2.is_clipped || c3.is_clipped) continue;
            if (c1.is_oscillating || c2.is_oscillating || c3.is_oscillating) continue;

            // Veto offline identico alla macro principale, se il CH4 e' presente.
            if (e.nchannels >= 4 && c1.cfd_ok && c2.cfd_ok) {
                double t_mu = 0.5 * (c1.t_cfd + c2.t_cfd);
                VetoResult vr = AnalyzeVeto(e.ch[VETO_CH_INDEX], t_mu);
                if (vr.has_signal) continue;
            }

            clean_events++;

            for (int f = 0; f < NF; f++) {
                bool ok1, ok2, ok3;
                double t1 = ComputeCFDTimeAtFraction(c1, fracs[f], ok1);
                double t2 = ComputeCFDTimeAtFraction(c2, fracs[f], ok2);
                double t3 = ComputeCFDTimeAtFraction(c3, fracs[f], ok3);
                if (!ok1 || !ok2 || !ok3) continue;

                double dt12 = t1 - t2;
                double x    = (dt12 - CAL_Q) / CAL_M;
                if (x < X_CUT_LO || x > X_CUT_HI) continue;

                double Tm   = t3 - 0.5 * (t1 + t2);
                double tof  = Tm - GetC(x);
                if (tof < tof_min_phys || tof <= 0.1) continue;

                double dxp  = x - PAR_X_PMT3;
                double path = sqrt(dxp * dxp + PAR_H * PAR_H);
                double beta = (path / tof) / C_LIGHT;
                double invv = tof / path;

                h_tof_cfd[f]->Fill(tof);
                if (beta > 0.0 && beta < 3.0) h_beta_cfd[f]->Fill(beta);
                if (invv > 0.0 && invv < 0.2) h_invv_cfd[f]->Fill(invv);
                accepted_by_frac[f]++;
            }
        }
    }

    fout->cd();

    TCanvas* c_all = new TCanvas("c_cfd_study_distributions",
        "Distribuzioni TOF/#beta/1v vs frazione CFD", 1350, 450);
    c_all->Divide(3, 1);

    auto drawCFDOverlay = [&](int pad, TH1D* hs[NF], const char* title) {
        c_all->cd(pad);
        gPad->SetGrid();
        double ymax = 0.0;
        for (int f = 0; f < NF; f++) ymax = std::max(ymax, hs[f]->GetMaximum());

        TLegend* leg = new TLegend(0.58, 0.62, 0.89, 0.89);
        leg->SetTextFont(42);
        leg->SetTextSize(0.032);

        for (int f = 0; f < NF; f++) {
            hs[f]->SetLineColor(cols[f]);
            hs[f]->SetLineWidth(2);
            hs[f]->SetLineStyle(styles[f]);
            hs[f]->SetTitle(title);
            hs[f]->SetMaximum((ymax > 0.0) ? 1.20 * ymax : 1.0);
            hs[f]->Draw(f == 0 ? "HIST" : "HIST same");
            leg->AddEntry(hs[f], Form("CFD %d%%", pct[f]), "l");
        }
        leg->Draw();
    };

    drawCFDOverlay(1, h_beta_cfd, "#beta vs frazione CFD;#beta;Conteggi");
    drawCFDOverlay(2, h_tof_cfd,  "TOF vs frazione CFD;TOF [ns];Conteggi");
    drawCFDOverlay(3, h_invv_cfd, "1/v vs frazione CFD;1/v [ns/cm];Conteggi");
    c_all->Write("Canvas_CFD_study_distributions");

    for (int f = 0; f < NF; f++) {
        h_beta_cfd[f]->Write();
        h_tof_cfd[f]->Write();
        h_invv_cfd[f]->Write();
    }

    fout->Close();

    std::cout << "\n=============================================" << std::endl;
    std::cout << "  CFD_study completata" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "  File processati:       " << file_list.size() << std::endl;
    std::cout << "  Eventi totali letti:    " << total_events << std::endl;
    std::cout << "  Eventi puliti comuni:   " << clean_events << std::endl;
    for (int f = 0; f < NF; f++) {
        std::cout << "  CFD " << pct[f] << "%: eventi accettati = "
                  << accepted_by_frac[f] << std::endl;
    }
    std::cout << "  Output: " << outname << std::endl;
    std::cout << "=============================================" << std::endl;
}
