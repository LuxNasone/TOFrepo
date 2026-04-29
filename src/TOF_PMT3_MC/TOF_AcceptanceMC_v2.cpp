/*******************************************************************************
 *
 *  TOF_AcceptanceMC.cpp
 *  =====================
 *
 *  SCOPO GENERALE
 *  ===============
 *  Simulazione Monte Carlo dell'accettanza geometrica dello scintillatore
 *  piccolo accoppiato a PMT3, nell'esperienza di Time-of-Flight (TOF)
 *  del corso di Laboratorio di Interazioni Fondamentali, Universita' di Pisa.
 *
 *  L'accettanza geometrica e' definita come la probabilita' che un raggio
 *  cosmico che attraversa la barra di scintillatore (causando una doppia
 *  coincidenza PMT1 & PMT2) attraversi ANCHE lo scintillatore di PMT3
 *  (causando una tripla coincidenza PMT1 & PMT2 & PMT3):
 *
 *     A_geom = N_triple / N_doppie = P(colpisce PMT3 | colpisce barra)
 *
 *
 *  GEOMETRIA DEL SETUP
 *  ====================
 *  CONFIGURAZIONE A (calibrazione):
 *    PMT3 SOPRA la barra (guida A, ~3.5 cm sopra).
 *    Il cosmico attraversa prima PMT3, poi la barra.
 *
 *  CONFIGURAZIONE B (misura TOF):
 *    PMT3 SOTTO la barra (guida B, ~156 cm sotto).
 *    Il cosmico attraversa prima la barra, poi PMT3.
 *
 *
 *  SISTEMA DI COORDINATE
 *  =====================
 *    x : lungo l'asse della barra  (0 = PMT1, L = PMT2)
 *    y : orizzontale perpendicolare (0 = centro sezione)
 *    z : verticale, positivo ALTO   (0 = faccia top barra)
 *
 *  FORMULA DI PROIEZIONE (unificata sopra/sotto):
 *    x(z_t) = x0  -  z_t * tan(theta) * cos(phi)
 *    y(z_t) = y0  -  z_t * tan(theta) * sin(phi)
 *
 *
 *  COMPILAZIONE E USO
 *  ==================
 *    .L TOF_AcceptanceMC.cpp++
 *
 *    // --- Simulazione principale ---
 *    TOF_AcceptanceMC(1000000, 2.0, 140.0, -160.6);   // guida B (sotto)
 *    TOF_AcceptanceMC(1000000, 2.0, 140.0, +4.1);     // guida A (sopra)
 *
 *    // --- Debug: verifica distribuzione angolare e generazione ---
 *    DebugDistributions(500000, 2.0, 140.0, -160.6);
 *
 *    // --- Scan in posizione / altezza ---
 *    ScanPosition(500000, 2.0, -160.6);
 *    ScanHeight(500000, 2.0, 140.0);
 *
 *    // --- Confronto esponenti ---
 *    CompareExponents(300000, 140.0, -160.6);
 *
 *
 *  NOTA SULL'ORIENTAZIONE DI PMT3
 *  ===============================
 *  Dimensioni: 1.2 x 5.8 x 20 cm.
 *  Orientazione:  5.8 cm lungo x (asse barra),  20 cm lungo y.
 *
 *
 *******************************************************************************/


// ============================================================================
//  INCLUSIONI
// ============================================================================

#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>        // std::max

#include <TRandom3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <fstream>

// ============================================================================
//  COSTANTI GEOMETRICHE DEL SETUP  (tutte in cm)
// ============================================================================

// Barra di scintillatore BICRON BC408
const double L_bar   = 280.0;    // lunghezza lungo x
const double W_bar   = 4.0;      // larghezza lungo y
const double T_bar   = 4.0;      // spessore  lungo z

// Scintillatore accoppiato a PMT3 (ORIENTAZIONE DA VERIFICARE!)
const double scint3_dx    = 5.8;    // lungo x (asse barra)
const double scint3_dy    = 20.0;   // lungo y (perpendicolare)
const double scint3_thick = 1.2;    // spessore z

// Quote z delle guide (centro scintillatore PMT3)
const double z_guidaA =    4.1;     // sopra la barra (calibrazione)
const double z_guidaB = -160.6;     // sotto la barra (TOF)

// Velocita' della luce
const double c_light_cm_ns = 29.9792;  // cm/ns


// ============================================================================
//  CLASSE CosmicRay
// ============================================================================
//
//  Genera raggi cosmici sulla faccia superiore della barra (z = 0).
//  Distribuzione angolare:  dN/dOmega ~ cos^alpha(theta).
//  Generazione:  cos(theta) = u^{1/(alpha+1)},  u ~ Uniform(0,1).
//
class CosmicRay {
private:
    double m_cosTheta, m_phi;
    double m_x0, m_y0;
    double m_exponent;
    double m_Lx, m_Ly;

public:
    CosmicRay(double exponent, double Lx, double Ly)
        : m_exponent(exponent), m_Lx(Lx), m_Ly(Ly) { Throw(); }

    void Throw() {
        m_x0 = gRandom->Uniform(0.0, m_Lx);
        m_y0 = gRandom->Uniform(-m_Ly / 2.0, m_Ly / 2.0);
        m_phi = gRandom->Uniform(0.0, 2.0 * M_PI);
        m_cosTheta = pow(gRandom->Uniform(), 1.0 / (m_exponent + 1.0));
    }

    double CosTheta() const { return m_cosTheta; }
    double Theta()    const { return acos(m_cosTheta); }
    double Phi()      const { return m_phi; }
    double X0()       const { return m_x0; }
    double Y0()       const { return m_y0; }

    // Proiezione unificata (z_t > 0 = sopra, z_t < 0 = sotto)
    double X_at_z(double z_t) const {
        double tanTh = sqrt(1.0 - m_cosTheta * m_cosTheta) / m_cosTheta;
        return m_x0 - z_t * tanTh * cos(m_phi);
    }
    double Y_at_z(double z_t) const {
        double tanTh = sqrt(1.0 - m_cosTheta * m_cosTheta) / m_cosTheta;
        return m_y0 - z_t * tanTh * sin(m_phi);
    }

    // Lunghezza percorso tra due quote
    double PathLength(double z_from, double z_to) const {
        return fabs(z_from - z_to) / m_cosTheta;
    }
};


// ============================================================================
//  CLASSE Scintillator
// ============================================================================
//
//  Rettangolo piano a quota z fissa. Verifica intersezione con un raggio.
//
class Scintillator {
private:
    double m_xmin, m_xmax, m_ymin, m_ymax, m_z;

public:
    Scintillator(double cx, double cy, double z, double dx, double dy)
        : m_xmin(cx - dx/2.0), m_xmax(cx + dx/2.0),
          m_ymin(cy - dy/2.0), m_ymax(cy + dy/2.0), m_z(z) {}

    bool CheckIntersect(const CosmicRay &ray) const {
        double x = ray.X_at_z(m_z);
        double y = ray.Y_at_z(m_z);
        return (x > m_xmin && x < m_xmax && y > m_ymin && y < m_ymax);
    }

    double Z()    const { return m_z; }
    double Xmin() const { return m_xmin; }
    double Xmax() const { return m_xmax; }
    double Ymin() const { return m_ymin; }
    double Ymax() const { return m_ymax; }
};


// ============================================================================
//  DICHIARAZIONI ANTICIPATE
// ============================================================================

void TOF_AcceptanceMC(int n = 1000000, double exponent = 2.0,
                       double x_pmt3 = 140.0, double z_pmt3 = -160.6);

void DebugDistributions(int n = 500000, double exponent = 2.0,
                         double x_pmt3 = 140.0, double z_pmt3 = -160.6);

void ScanPosition(int n_per_point = 500000, double exponent = 2.0,
                  double z_pmt3 = -160.6, int n_points = 28);

void ScanHeight(int n_per_point = 500000, double exponent = 2.0,
                double x_pmt3 = 140.0, int n_points = 30);

void CompareExponents(int n_per_point = 300000, double x_pmt3 = 140.0,
                      double z_pmt3 = -160.6);


// ============================================================================
//  TOF_AcceptanceMC  -  FUNZIONE PRINCIPALE
// ============================================================================
//
//  Produce 5 canvas indipendenti (uno per grafico):
//    1) Posizione impatto x sulla barra (accettati)
//    2) Footprint 2D sulla barra (accettati)
//    3) Mappa 2D intercetta sul piano PMT3
//    4) Distribuzione lunghezze di percorso
//    5) Distribuzione TOF ideale (beta = 1)
//
//  Per i grafici di debug/verifica della generazione: DebugDistributions().
//
void TOF_AcceptanceMC(int n, double exponent, double x_pmt3, double z_pmt3)
{
    gRandom = new TRandom3(0);

    bool pmt3_sopra = (z_pmt3 > 0);
    std::string config_name = pmt3_sopra
        ? "GUIDA A (PMT3 sopra, calibrazione)"
        : "GUIDA B (PMT3 sotto, TOF)";

    Scintillator pmt3(x_pmt3, 0.0, z_pmt3, scint3_dx, scint3_dy);

    // --- Stampa configurazione ---
    std::cout << "\n============================================================\n";
    std::cout << "  TOF Acceptance Monte Carlo\n";
    std::cout << "  Configurazione: " << config_name << "\n";
    std::cout << "============================================================\n";
    std::cout << "  Raggi generati:      " << n << "\n";
    std::cout << "  Esponente alpha:     " << exponent << "\n";
    std::cout << "  Barra:               " << L_bar << " x " << W_bar
              << " x " << T_bar << " cm^3\n";
    std::cout << "  Scint. PMT3:         " << scint3_dx << " x " << scint3_dy
              << " x " << scint3_thick << " cm^3\n";
    std::cout << "  Posiz. PMT3 (x):     " << x_pmt3 << " cm\n";
    std::cout << "  Quota PMT3 (z):      " << z_pmt3 << " cm  ("
              << (pmt3_sopra ? "SOPRA" : "SOTTO") << ")\n";
    std::cout << "============================================================\n\n";


    // --- ISTOGRAMMI ---

    // 1) x_imp sulla barra (solo accettati)
    TH1F *h_ximp_acc = new TH1F("h_ximp_acc",
        "Posizione impatto sulla barra;x [cm];Conteggi",
        140, 0.0, L_bar);

    // 2) Footprint 2D sulla barra (accettati)
    TH2F *h_xy_bar_acc = new TH2F("h_xy_bar_acc",
        "Heatmap di intercetta sulla barra;x [cm];y [cm]",
        280, 0.0, L_bar, 40, -W_bar/2.0, W_bar/2.0);

    // 3) Intercetta sul piano PMT3 (accettati)
    double x_range = std::max(100.0, fabs(z_pmt3) * 1.0);
    TH2F *h_xy_pmt3 = new TH2F("h_xy_pmt3",
        "Heatmap sul piano PMT3 ;x [cm];y [cm]",
        200, x_pmt3 - x_range/2, x_pmt3 + x_range/2,
        200, -x_range/3, x_range/3);

    // 4) Path length (accettati)
    double path_max = fabs(z_pmt3) * 2.5 + 50.0;
    TH1F *h_pathlen = new TH1F("h_pathlen",
        "Lunghezza percorso barra #leftrightarrow PMT3;L [cm];Conteggi",
        200, 0.0, path_max);

    // 5) TOF ideale (beta = 1)
    // double tof_max = path_max / c_light_cm_ns;
    // TH1F *h_tof_ideal = new TH1F("h_tof_ideal",
    //     "TOF ideale (#beta = 1);TOF [ns];Conteggi",
    //     200, 0.0, tof_max);


    // --- CICLO MC ---
    CosmicRay ray(exponent, L_bar, W_bar);
    int n_accepted = 0;

    for (int i = 0; i < n; i++) {
        ray.Throw();

        if (pmt3.CheckIntersect(ray)) {
            n_accepted++;

            h_ximp_acc->Fill(ray.X0());
            h_xy_bar_acc->Fill(ray.X0(), ray.Y0());
            h_xy_pmt3->Fill(ray.X_at_z(z_pmt3), ray.Y_at_z(z_pmt3));

            double z_face = pmt3_sopra ? 0.0 : -T_bar;
            double path = ray.PathLength(z_face, z_pmt3);
            h_pathlen->Fill(path);
            // h_tof_ideal->Fill(path / c_light_cm_ns);
        }
    }


    // --- RISULTATI ---
    double acceptance = static_cast<double>(n_accepted) / static_cast<double>(n);
    double sigma_acc  = sqrt(acceptance * (1.0 - acceptance) / static_cast<double>(n));

    std::cout << "  RISULTATI\n";
    std::cout << "  ---------------------------------------------------------\n";
    std::cout << "  Eventi generati (doppie):    " << n << "\n";
    std::cout << "  Eventi accettati (triple):   " << n_accepted << "\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Accettanza geometrica:       " << acceptance
              << " +/- " << sigma_acc << "\n";
    std::cout << "  (= " << std::setprecision(3) << acceptance * 100 << " %)\n";
    std::cout << "  ---------------------------------------------------------\n\n";

    double rate_doppie = 20.0;  // Hz (indicativo, da misurare!)
    double rate_triple = rate_doppie * acceptance;
    std::cout << "  Stima rate triple (rate doppie ~ " << rate_doppie << " Hz):\n";
    std::cout << "    Rate triple ~ " << std::setprecision(4)
              << rate_triple << " Hz\n";
    if (rate_triple > 0)
        std::cout << "    Tempo per 1000 eventi ~ "
                  << std::setprecision(0) << 1000.0/rate_triple
                  << " s = " << std::setprecision(1)
                  << 1000.0/rate_triple/60.0 << " min\n";
    std::cout << "\n";


    // --- VISUALIZZAZIONE (un canvas per grafico) ---
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPalette(kBird);

    // Canvas 1: posizione impatto x
    TCanvas *c1 = new TCanvas("c_ximp", "x impatto (accettati)", 800, 600);
    h_ximp_acc->SetLineColor(kRed);
    h_ximp_acc->SetLineWidth(2);
    h_ximp_acc->SetFillColorAlpha(kRed, 0.15);
    h_ximp_acc->Draw();
    TLine *lp = new TLine(x_pmt3, 0, x_pmt3, h_ximp_acc->GetMaximum());
    lp->SetLineColor(kBlue); lp->SetLineStyle(2); lp->SetLineWidth(2);
    lp->Draw("same");
    TLatex *tp = new TLatex(x_pmt3+5, h_ximp_acc->GetMaximum()*0.9,
                             Form("x_{PMT3} = %.0f cm", x_pmt3));
    tp->SetTextSize(0.035); tp->SetTextColor(kBlue); tp->Draw();
    c1->Update();

    // Canvas 2: footprint 2D barra
    TCanvas *c2 = new TCanvas("c_footprint", "Footprint barra (accettati)", 900, 400);
    h_xy_bar_acc->Draw("COLZ");
    c2->Update();

    // Canvas 3: mappa piano PMT3
    TCanvas *c3 = new TCanvas("c_pmt3_map", "Intercetta piano PMT3", 800, 600);
    h_xy_pmt3->Draw("COLZ");
    // Rettangolo bordi scintillatore
    TLine *box[4];
    box[0] = new TLine(pmt3.Xmin(), pmt3.Ymin(), pmt3.Xmax(), pmt3.Ymin());
    box[1] = new TLine(pmt3.Xmax(), pmt3.Ymin(), pmt3.Xmax(), pmt3.Ymax());
    box[2] = new TLine(pmt3.Xmax(), pmt3.Ymax(), pmt3.Xmin(), pmt3.Ymax());
    box[3] = new TLine(pmt3.Xmin(), pmt3.Ymax(), pmt3.Xmin(), pmt3.Ymin());
    for (int i = 0; i < 4; i++) {
        box[i]->SetLineColor(kRed); box[i]->SetLineWidth(2);
        box[i]->Draw("same");
    }
    c3->Update();

    // // Canvas 4: path length
    // TCanvas *c4 = new TCanvas("c_pathlen", "Path length", 800, 600);
    // h_pathlen->SetLineColor(kMagenta+1);
    // h_pathlen->SetLineWidth(2);
    // h_pathlen->SetFillColorAlpha(kMagenta+1, 0.15);
    // h_pathlen->Draw();
    // Linea per il percorso minimo (verticale)
    // double z_face = pmt3_sopra ? 0.0 : -T_bar;
    // double L_vert = fabs(z_face - z_pmt3);
    // TLine *lv = new TLine(L_vert, 0, L_vert, h_pathlen->GetMaximum()*0.7);
    // lv->SetLineColor(kBlue); lv->SetLineStyle(2); lv->SetLineWidth(2);
    // lv->Draw("same");
    // TLatex *tv = new TLatex(L_vert+2, h_pathlen->GetMaximum()*0.72,
    //                          Form("L_{min} = %.1f cm", L_vert));
    // tv->SetTextSize(0.03); tv->SetTextColor(kBlue); tv->Draw();
    // c4->Update();

    // // Canvas 5: TOF ideale
    // TCanvas *c5 = new TCanvas("c_tof_ideal", "TOF ideale (beta=1)", 800, 600);
    // h_tof_ideal->SetLineColor(kBlue+1);
    // h_tof_ideal->SetLineWidth(2);
    // h_tof_ideal->SetFillColorAlpha(kBlue+1, 0.15);
    // h_tof_ideal->Draw();
    // c5->Update();

    // std::cout << "  Simulazione completata.\n";
    // std::cout << "  Per grafici di debug: DebugDistributions(...)\n\n";
}


// ============================================================================
//  DebugDistributions  -  GRAFICI DI CONTROLLO (funzione separata)
// ============================================================================
//
//  Verifica la correttezza del generatore MC. Da chiamare solo quando serve.
//
//  Produce 4 canvas indipendenti:
//    DBG 1:  cos(theta) generato  +  curva teorica sovrapposta
//    DBG 2:  cos(theta) accettati (effetto selezione geometrica)
//    DBG 3:  theta [rad] accettati
//    DBG 4:  x_imp generato (verifica uniformita')
//
//  USO:
//    DebugDistributions(500000, 2.0, 140.0, -160.6);
//
void DebugDistributions(int n, double exponent, double x_pmt3, double z_pmt3)
{
    gRandom = new TRandom3(0);

    Scintillator pmt3(x_pmt3, 0.0, z_pmt3, scint3_dx, scint3_dy);

    std::cout << "\n============================================================\n";
    std::cout << "  DEBUG DISTRIBUZIONI\n";
    std::cout << "  n = " << n << ",  alpha = " << exponent
              << ",  x = " << x_pmt3 << " cm,  z = " << z_pmt3 << " cm\n";
    std::cout << "============================================================\n\n";

    // --- Istogrammi di debug ---

    // cos(theta) generato: deve seguire cos^alpha(theta)
    TH1F *h_ct_gen = new TH1F("h_dbg_ct_gen",
        "cos(#theta) generato (tutti);cos(#theta);Normalizzato",
        100, 0.0, 1.0);

    // cos(theta) accettati: mostra la selezione geometrica
    TH1F *h_ct_acc = new TH1F("h_dbg_ct_acc",
        "cos(#theta) accettati;cos(#theta);Normalizzato",
        100, 0.0, 1.0);

    // theta accettati
    TH1F *h_th_acc = new TH1F("h_dbg_th_acc",
        "#theta accettati;#theta [rad];Normalizzato",
        100, 0.0, M_PI / 2.0);

    // x_imp generato: deve essere uniforme
    TH1F *h_x_gen = new TH1F("h_dbg_x_gen",
        "x_{imp} generato (verifica uniformita');x [cm];Conteggi",
        140, 0.0, L_bar);


    // --- Ciclo di generazione ---
    CosmicRay ray(exponent, L_bar, W_bar);

    for (int i = 0; i < n; i++) {
        ray.Throw();
        h_ct_gen->Fill(ray.CosTheta());
        h_x_gen->Fill(ray.X0());

        if (pmt3.CheckIntersect(ray)) {
            h_ct_acc->Fill(ray.CosTheta());
            h_th_acc->Fill(ray.Theta());
        }
    }


    // --- Stile ---
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);


    // ===== DBG 1: cos(theta) generato + curva teorica =====
    //
    //  DrawNormalized() di TH1 divide ogni bin per N_entries.
    //  Il contenuto normalizzato di ogni bin e':  p(c) * w
    //  dove p(c) = (alpha+1)*c^alpha  e  w = bin_width.
    //
    //  Per sovrapporre una TF1 alla stessa scala, impostiamo:
    //    f(c) = [0] * c^[1]  con  [0] = (alpha+1) * w
    //
    TCanvas *cd1 = new TCanvas("c_dbg_ct_gen",
        "DEBUG: cos(theta) generato", 800, 600);
    h_ct_gen->SetLineColor(kBlue);
    h_ct_gen->SetLineWidth(2);
    h_ct_gen->DrawNormalized();
    double bw = h_ct_gen->GetBinWidth(1);
    TF1 *f_cos = new TF1("f_cos_dbg", "[0]*pow(x,[1])", 0, 1);
    f_cos->SetParameters((exponent + 1.0) * bw, exponent);
    f_cos->SetLineColor(kRed);
    f_cos->SetLineStyle(2);
    f_cos->SetLineWidth(2);
    f_cos->Draw("same");
    TLegend *leg = new TLegend(0.15, 0.72, 0.58, 0.88);
    leg->AddEntry(h_ct_gen, "MC generato", "l");
    leg->AddEntry(f_cos, Form("PDF: (#alpha+1) cos^{#alpha}(#theta) #times w"), "l");
    leg->Draw();
    cd1->Update();


    // ===== DBG 2: cos(theta) accettati =====
    TCanvas *cd2 = new TCanvas("c_dbg_ct_acc",
        "DEBUG: cos(theta) accettati", 800, 600);
    h_ct_acc->SetLineColor(kRed);
    h_ct_acc->SetLineWidth(2);
    h_ct_acc->SetFillColorAlpha(kRed, 0.2);
    h_ct_acc->DrawNormalized();
    cd2->Update();


    // ===== DBG 3: theta accettati =====
    TCanvas *cd3 = new TCanvas("c_dbg_th_acc",
        "DEBUG: theta accettati", 800, 600);
    h_th_acc->SetLineColor(kGreen + 2);
    h_th_acc->SetLineWidth(2);
    h_th_acc->SetFillColorAlpha(kGreen + 2, 0.2);
    h_th_acc->DrawNormalized();
    cd3->Update();


    // ===== DBG 4: x_imp generato (uniformita') =====
    TCanvas *cd4 = new TCanvas("c_dbg_x_gen",
        "DEBUG: x generato (uniformita')", 800, 600);
    h_x_gen->SetLineColor(kBlue);
    h_x_gen->SetLineWidth(2);
    h_x_gen->Draw();
    // Linea al valore atteso per bin (N / n_bins)
    double expected = static_cast<double>(n) / 140.0;
    TLine *lf = new TLine(0, expected, L_bar, expected);
    lf->SetLineColor(kRed); lf->SetLineStyle(2); lf->SetLineWidth(2);
    lf->Draw("same");
    TLatex *tf = new TLatex(10, expected * 1.06,
                             Form("Atteso: %.0f / bin", expected));
    tf->SetTextColor(kRed); tf->SetTextSize(0.03); tf->Draw();
    cd4->Update();

    std::cout << "  Debug completato.\n\n";
}


// ============================================================================
//  ScanPosition:  accettanza vs posizione x
// ============================================================================
void ScanPosition(int n_per_point, double exponent,
                  double z_pmt3, int n_points)
{
    gRandom = new TRandom3(0);
    std::string config = (z_pmt3 > 0) ? "sopra" : "sotto";

    std::vector<double> vx, vy, vex, vey;

    std::cout << "\n============================================================\n";
    std::cout << "  SCAN IN POSIZIONE  (z = " << z_pmt3 << " cm, "
              << config << ", alpha = " << exponent << ")\n";
    std::cout << "============================================================\n";
    std::cout << std::setw(12) << "x [cm]"
              << std::setw(15) << "Accettanza"
              << std::setw(15) << "Sigma\n";
    std::cout << "--------------------------------------------\n";

    double x_start = 10.0, x_end = 270.0;
    double dx = (x_end - x_start) / (n_points - 1);
    CosmicRay ray(exponent, L_bar, W_bar);

    for (int ip = 0; ip < n_points; ip++) {
        double x_pos = x_start + ip * dx;
        Scintillator pmt3(x_pos, 0.0, z_pmt3, scint3_dx, scint3_dy);

        int n_acc = 0;
        for (int i = 0; i < n_per_point; i++) {
            ray.Throw();
            if (pmt3.CheckIntersect(ray)) n_acc++;
        }
        double acc = static_cast<double>(n_acc) / static_cast<double>(n_per_point);
        double sig = sqrt(acc * (1.0-acc) / static_cast<double>(n_per_point));

        vx.push_back(x_pos); vy.push_back(acc);
        vex.push_back(0.0);  vey.push_back(sig);

        std::cout << std::fixed << std::setprecision(1) << std::setw(12) << x_pos
                  << std::setprecision(6)
                  << std::setw(15) << acc << std::setw(15) << sig << "\n";
    }

    TGraphErrors *gr = new TGraphErrors(vx.size(), vx.data(), vy.data(),
                                         vex.data(), vey.data());
    gr->SetTitle(Form("Accettanza PMT3 vs posizione (z=%.1f cm, #alpha=%.1f)"
                      ";x_{PMT3} [cm];Accettanza", z_pmt3, exponent));
    gr->SetMarkerStyle(20); gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(kBlue+1); gr->SetLineColor(kBlue+1);

    TCanvas *cs = new TCanvas("c_scan_pos",
        Form("Scan posizione (z=%.0f)", z_pmt3), 900, 600);
    gStyle->SetOptStat(0);
    gr->Draw("APE");

    // Media regione centrale
    double acc_mean = 0.0; int nc = 0;
    for (size_t i = 0; i < vx.size(); i++)
        if (vx[i] > 50 && vx[i] < 230) { acc_mean += vy[i]; nc++; }
    if (nc > 0) acc_mean /= nc;

    TLine *lm = new TLine(0, acc_mean, L_bar, acc_mean);
    lm->SetLineColor(kRed); lm->SetLineStyle(2); lm->SetLineWidth(2);
    lm->Draw("same");
    TLatex *tm = new TLatex(15, acc_mean*1.08,
                             Form("Media (50-230 cm): %.4f", acc_mean));
    tm->SetTextColor(kRed); tm->SetTextSize(0.035); tm->Draw();
    cs->Update();

    std::cout << "\n  Media (50-230 cm): " << std::setprecision(5)
              << acc_mean << "\n";
    std::cout << "============================================================\n\n";
}


// ============================================================================
//  ScanHeight:  accettanza vs distanza verticale (sotto la barra)
// ============================================================================
void ScanHeight(int n_per_point, double exponent,
                double x_pmt3, int n_points)
{
    gRandom = new TRandom3(0);
    std::vector<double> vh, va, veh, vea;

    std::cout << "\n============================================================\n";
    std::cout << "  SCAN IN ALTEZZA  (x = " << x_pmt3 << " cm, alpha = "
              << exponent << ")\n";
    std::cout << "============================================================\n";
    std::cout << std::setw(15) << "h_gap [cm]" << std::setw(15) << "z [cm]"
              << std::setw(15) << "Accettanza" << std::setw(15) << "Sigma\n";
    std::cout << "------------------------------------------------------------\n";

    double h0 = 0.0, h1 = 200.0;
    double dh = (h1 - h0) / (n_points - 1);
    CosmicRay ray(exponent, L_bar, W_bar);

    for (int ip = 0; ip < n_points; ip++) {
        double h = h0 + ip * dh;
        double z = -(T_bar + h + scint3_thick / 2.0);
        Scintillator pmt3(x_pmt3, 0.0, z, scint3_dx, scint3_dy);

        int n_acc = 0;
        for (int i = 0; i < n_per_point; i++) {
            ray.Throw();
            if (pmt3.CheckIntersect(ray)) n_acc++;
        }
        double acc = static_cast<double>(n_acc) / static_cast<double>(n_per_point);
        double sig = sqrt(acc*(1.0-acc) / static_cast<double>(n_per_point));

        vh.push_back(h); va.push_back(acc);
        veh.push_back(0.0); vea.push_back(sig);

        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(15) << h << std::setw(15) << z
                  << std::setprecision(6)
                  << std::setw(15) << acc << std::setw(15) << sig << "\n";
    }

    TGraphErrors *gr = new TGraphErrors(vh.size(), vh.data(), va.data(),
                                         veh.data(), vea.data());
    gr->SetTitle(Form("Accettanza vs h_{gap} (x=%.0f cm, #alpha=%.1f)"
                      ";h_{gap} [cm];Accettanza", x_pmt3, exponent));
    gr->SetMarkerStyle(20); gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(kRed+1); gr->SetLineColor(kRed+1);

    TCanvas *ch = new TCanvas("c_scan_height", "Scan altezza", 900, 600);
    gStyle->SetOptStat(0);
    gr->Draw("APE");

    // Linee guide
    struct G { double h; const char* n; int c; };
    G gs[] = {{0.0,"Contatto",kGreen+2},{105.0,"Guida int.",kOrange+1},
              {156.0,"Guida inf.",kBlue}};
    for (int ig = 0; ig < 3; ig++) {
        TLine *lg = new TLine(gs[ig].h, gr->GetYaxis()->GetXmin(),
                              gs[ig].h, gr->GetYaxis()->GetXmax());
        lg->SetLineColor(gs[ig].c); lg->SetLineStyle(2); lg->SetLineWidth(2);
        lg->Draw("same");
        TLatex *lt = new TLatex(gs[ig].h+3,
                                gr->GetYaxis()->GetXmax()*(0.92-ig*0.08),
                                gs[ig].n);
        lt->SetTextColor(gs[ig].c); lt->SetTextSize(0.028); lt->Draw();
    }
    ch->Update();
    std::cout << "============================================================\n\n";
}


// ============================================================================
//  CompareExponents:  confronto per diversi alpha
// ============================================================================
void CompareExponents(int n_per_point, double x_pmt3, double z_pmt3)
{
    gRandom = new TRandom3(0);
    double alphas[] = {1.5, 2.0, 2.5, 3.0, 3.5};
    int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kOrange+1};
    int na = 5, np = 20;
    double x0 = 10.0, x1 = 270.0, dx = (x1-x0)/(np-1);

    TCanvas *cc = new TCanvas("c_compare", "Confronto esponenti", 900, 600);
    TLegend *leg = new TLegend(0.65, 0.62, 0.88, 0.88);
    bool first = true;

    for (int ia = 0; ia < na; ia++) {
        CosmicRay ray(alphas[ia], L_bar, W_bar);
        std::vector<double> vx, vy;
        for (int ip = 0; ip < np; ip++) {
            double xp = x0 + ip * dx;
            Scintillator pmt3(xp, 0.0, z_pmt3, scint3_dx, scint3_dy);
            int n_acc = 0;
            for (int i = 0; i < n_per_point; i++) {
                ray.Throw();
                if (pmt3.CheckIntersect(ray)) n_acc++;
            }
            vx.push_back(xp);
            vy.push_back(static_cast<double>(n_acc)/static_cast<double>(n_per_point));
        }
        TGraph *gr = new TGraph(vx.size(), vx.data(), vy.data());
        gr->SetLineColor(colors[ia]); gr->SetLineWidth(2);
        gr->SetMarkerColor(colors[ia]); gr->SetMarkerStyle(20); gr->SetMarkerSize(0.6);
        if (first) {
            gr->SetTitle(Form("Accettanza vs posizione per diversi #alpha "
                              "(z=%.0f cm);x_{PMT3} [cm];Accettanza", z_pmt3));
            gr->Draw("ALP"); first = false;
        } else gr->Draw("LP same");
        leg->AddEntry(gr, Form("#alpha = %.1f", alphas[ia]), "lp");
    }
    leg->Draw(); cc->Update();
    std::cout << "\n  Confronto esponenti completato.\n\n";
}
// ============================================================================
//  ============================================================================
//      SEZIONE PARALLASSE — STIMA DELL'ERRORE σ_x DA MC
//  ============================================================================
//
//  Funzioni per stimare l'errore stocastico σ_x sulla posizione di passaggio
//  del muone NELLA BARRA, dato che il PMT3 è centrato in x_PMT3 (la posizione
//  letta col metro).
//
//  L'errore "di parallasse" σ_x include:
//    - dimensione finita dello scintillatore di PMT3 (~ ℓ/√12)
//    - parallasse geometrica dovuta alla distribuzione angolare dei cosmici
//    - smearing verticale dovuto allo spessore della barra (opzionale)
//
//  Il MC genera questi tre contributi insieme nel modo fisicamente corretto.
//
//  CONVENZIONE COORDINATE: queste funzioni usano coordinate CENTRATE
//  (x ∈ [-L_bar/2, +L_bar/2], cioè x=0 al centro della barra), consistenti
//  con il codice di calibrazione TOF_Calibration_v3.cpp. La conversione
//  alle coordinate native del MC (0..L_bar) è fatta internamente.
//
// ============================================================================


// ----------------------------------------------------------------------------
//  Conversione coordinate centrate ↔ native MC
// ----------------------------------------------------------------------------
//  Coordinate centrate (codice calibrazione): x ∈ [-140, +140] cm
//  Coordinate native MC:                       x ∈ [0, 280] cm
//  Relazione: x_native = x_centered + L_bar/2

inline double CenteredToNative(double x_centered) {
    return x_centered + L_bar / 2.0;
}
inline double NativeToCentered(double x_native) {
    return x_native - L_bar / 2.0;
}


// ----------------------------------------------------------------------------
//  ParallaxScan
// ----------------------------------------------------------------------------
//
//  Fa lo scan delle posizioni del PMT3 e per ciascuna calcola:
//    bias_x  = <x_barra - x_PMT3>           (valor medio, bias sistematico)
//    sigma_x = sqrt(Var[x_barra - x_PMT3])  (dev. std. = errore stocastico)
//    mse_x   = sqrt(<(x_barra - x_PMT3)²>)  (sqrt MSE = bias² + sigma²)
//
//  L'errore stocastico σ_x è quello che ENTRA come errore su x nel fit
//  Δt vs x_PMT3 (insieme ad altri termini come la lettura del metro).
//  Il bias è una correzione SISTEMATICA che andrebbe sottratta a x_PMT3
//  prima del fit, oppure trattata come incertezza separata se piccola.
//
//  PARAMETRI:
//    n_per_point  — numero di muoni generati per ogni posizione
//    exponent     — esponente α della distribuzione cos^α(θ)
//    z_pmt3       — quota verticale del PMT3 [cm]:
//                     +4.1   → Guida A (sopra, calibrazione)
//                     -160.6 → Guida B (sotto, TOF)
//    x_min, x_max, dx_step — griglia di posizioni in COORDINATE CENTRATE [cm]
//    use_bar_thickness     — se true, simula z_scint uniforme in [-T_bar, 0]
//                            (cattura lo smearing verticale dello spessore barra)
//                            Per Guida A irrilevante; per Guida B può contribuire.
//    output_file  — nome del file ROOT di output. Salva:
//                     TTree "parallax" con bias, sigma, mse, accettanza
//                     TGraphErrors "g_sigma_x" per uso diretto
//
//  USO TIPICO:
//    // Per la calibrazione (Guida A):
//    ParallaxScan(500000, 2.0, +4.1, -135.0, +135.0, 5.0,
//                 false, "parallax_GuidaA.root");
//
//    // Per la misura TOF (Guida B):
//    ParallaxScan(500000, 2.0, -160.6, -135.0, +135.0, 5.0,
//                 true, "parallax_GuidaB.root");

void ParallaxScan(int n_per_point = 500000,
                  double exponent = 2.0,
                  double z_pmt3 = +4.1,
                  double x_min = -135.0,
                  double x_max = +135.0,
                  double dx_step = 5.0,
                  bool use_bar_thickness = false,
                  const char* output_file = "parallax_scan.root")
{
    gRandom = new TRandom3(0);

    bool pmt3_sopra = (z_pmt3 > 0);
    std::string config_name = pmt3_sopra
        ? "GUIDA A (PMT3 sopra, calibrazione)"
        : "GUIDA B (PMT3 sotto, TOF)";

    // --- Costruzione griglia posizioni ---
    int n_points = (int)round((x_max - x_min) / dx_step) + 1;

    std::cout << "\n============================================================\n";
    std::cout << "  PARALLAX SCAN — " << config_name << "\n";
    std::cout << "============================================================\n";
    std::cout << "  Eventi per punto:    " << n_per_point << "\n";
    std::cout << "  Esponente α:         " << exponent << "\n";
    std::cout << "  z_PMT3:              " << z_pmt3 << " cm\n";
    std::cout << "  Griglia x:           [" << x_min << ", " << x_max
              << "] cm, passo " << dx_step << " cm  ("
              << n_points << " punti)\n";
    std::cout << "  Smearing barra:      "
              << (use_bar_thickness ? "ON" : "OFF") << "\n";
    std::cout << "  Output:              " << output_file << "\n";
    std::cout << "============================================================\n\n";


    // --- File ROOT di output ---
    TFile *fout = new TFile(output_file, "RECREATE");

    // TTree: una entry per posizione di scan
    TTree *tree = new TTree("parallax",
        "Statistiche di parallasse per posizione del PMT3");
    Double_t t_x_pmt3, t_bias_x, t_sigma_x, t_mse_x, t_accept;
    Int_t    t_n_acc;
    tree->Branch("x_pmt3",  &t_x_pmt3, "x_pmt3/D");   // posizione centrata [cm]
    tree->Branch("bias_x",  &t_bias_x, "bias_x/D");   // <x_barra - x_PMT3> [cm]
    tree->Branch("sigma_x", &t_sigma_x,"sigma_x/D");  // dev std stocastica [cm]
    tree->Branch("mse_x",   &t_mse_x,  "mse_x/D");    // sqrt(MSE) [cm]
    tree->Branch("accept",  &t_accept, "accept/D");   // accettanza
    tree->Branch("n_acc",   &t_n_acc,  "n_acc/I");    // n eventi accettati


    // --- Vettori per il TGraphErrors ---
    std::vector<double> v_x, v_sigma, v_ex, v_esigma;
    std::vector<double> v_bias;


    // --- Stampa header tabella ---
    std::cout << std::setw(10) << "x [cm]"
              << std::setw(13) << "bias [cm]"
              << std::setw(13) << "sigma [cm]"
              << std::setw(13) << "MSE [cm]"
              << std::setw(13) << "Acc [%]"
              << std::setw(10) << "N_acc\n";
    std::cout << "----------------------------------------------------------------\n";


    // --- Loop sulle posizioni ---
    CosmicRay ray(exponent, L_bar, W_bar);

    for (int ip = 0; ip < n_points; ip++) {

        // Posizione del PMT3 in coordinate centrate
        double x_centered = x_min + ip * dx_step;

        // Conversione a coordinate native del MC (0..L_bar) per Scintillator
        double x_native   = CenteredToNative(x_centered);

        Scintillator pmt3(x_native, 0.0, z_pmt3, scint3_dx, scint3_dy);

        // Accumulatori per le statistiche di (x_barra - x_PMT3)
        double sum_dx  = 0.0;   // somma differenze
        double sum_dx2 = 0.0;   // somma quadrati
        int    n_acc   = 0;

        for (int i = 0; i < n_per_point; i++) {

            ray.Throw();

            if (pmt3.CheckIntersect(ray)) {
                // x_barra = posizione di scintillazione nella barra (in coord native)
                // Nel codice originale è semplicemente ray.X0() (sulla superficie z=0).
                //
                // Se use_bar_thickness=true, simuliamo z_scint uniforme dentro
                // lo spessore della barra: z_scint ∈ [-T_bar, 0]. La correzione
                // a x_barra è: -z_scint * tan(θ) * cos(φ).
                double x_barra_native = ray.X0();

                if (use_bar_thickness) {
                    double z_scint = gRandom->Uniform(-T_bar, 0.0);
                    double cosTh   = ray.CosTheta();
                    double tanTh   = sqrt(1.0 - cosTh*cosTh) / cosTh;
                    x_barra_native = ray.X0() - z_scint * tanTh * cos(ray.Phi());
                }

                // Differenza in coordinate centrate (è invariante per traslazione)
                double dx_event = x_barra_native - x_native;

                sum_dx  += dx_event;
                sum_dx2 += dx_event * dx_event;
                n_acc++;
            }
        }

        // --- Calcolo statistiche ---
        double accept = (double)n_acc / (double)n_per_point;
        double bias_x = 0.0, sigma_x = 0.0, mse_x = 0.0;

        if (n_acc > 1) {
            bias_x = sum_dx / n_acc;
            // Varianza campionaria (con correzione di Bessel non strettamente
            // necessaria per n_acc grande, ma ben venga essere precisi)
            double var_x = (sum_dx2 - n_acc * bias_x * bias_x) / (n_acc - 1);
            if (var_x < 0) var_x = 0;  // protezione errori numerici
            sigma_x = sqrt(var_x);
            mse_x   = sqrt(sum_dx2 / n_acc);  // sqrt(MSE) = sqrt(bias² + var)
        }

        // --- Salva nel TTree ---
        t_x_pmt3 = x_centered;
        t_bias_x = bias_x;
        t_sigma_x = sigma_x;
        t_mse_x  = mse_x;
        t_accept = accept;
        t_n_acc  = n_acc;
        tree->Fill();

        // --- Salva nei vettori per il TGraph ---
        v_x.push_back(x_centered);
        v_sigma.push_back(sigma_x);
        v_bias.push_back(bias_x);
        v_ex.push_back(0.0);
        // Errore su sigma stesso (errore della deviazione standard di un
        // campione gaussiano): σ(σ) ≈ σ / √(2(n-1))
        double sig_err = (n_acc > 1) ? sigma_x / sqrt(2.0 * (n_acc - 1)) : 0;
        v_esigma.push_back(sig_err);

        // --- Stampa riga ---
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << x_centered
                  << std::setw(13) << bias_x
                  << std::setw(13) << sigma_x
                  << std::setw(13) << mse_x
                  << std::setw(13) << std::setprecision(3) << accept * 100.0
                  << std::setw(10) << n_acc << "\n";
    }


    // --- Crea e salva TGraphErrors di sigma_x ---
    TGraphErrors *g_sigma = new TGraphErrors(v_x.size(),
                                              v_x.data(), v_sigma.data(),
                                              v_ex.data(), v_esigma.data());
    g_sigma->SetName("g_sigma_x");
    g_sigma->SetTitle(Form("#sigma_{x} di parallasse vs x_{PMT3} (%s);"
                            "x_{PMT3} [cm];#sigma_{x} [cm]",
                            config_name.c_str()));
    g_sigma->SetMarkerStyle(20);
    g_sigma->SetMarkerColor(pmt3_sopra ? kBlue+1 : kRed+1);
    g_sigma->SetLineColor(pmt3_sopra ? kBlue+1 : kRed+1);

    // --- TGraph per il bias (errori su bias trascurabili a queste statistiche) ---
    TGraph *g_bias = new TGraph(v_x.size(), v_x.data(), v_bias.data());
    g_bias->SetName("g_bias_x");
    g_bias->SetTitle(Form("Bias di parallasse vs x_{PMT3} (%s);"
                           "x_{PMT3} [cm];#LT x_{barra} - x_{PMT3} #GT [cm]",
                           config_name.c_str()));
    g_bias->SetMarkerStyle(21);
    g_bias->SetMarkerColor(kGreen+2);
    g_bias->SetLineColor(kGreen+2);

    // --- Salva tutto su file ---
    fout->cd();
    tree->Write();
    g_sigma->Write();
    g_bias->Write();
    fout->Close();

    std::cout << "\n  Risultati salvati in: " << output_file << "\n";
    std::cout << "============================================================\n\n";
}


// ----------------------------------------------------------------------------
//  PlotParallax
// ----------------------------------------------------------------------------
//
//  Carica uno o due file di output di ParallaxScan e disegna i grafici di
//  σ_x e bias in funzione di x_PMT3. Utile per confronto Guida A vs B.
//
//  USO:
//    PlotParallax("parallax_GuidaA.root");
//    PlotParallax("parallax_GuidaA.root", "parallax_GuidaB.root");

void PlotParallax(const char* file_a,
                  const char* file_b = nullptr)
{
    TFile *fa = TFile::Open(file_a);
    if (!fa || fa->IsZombie()) {
        std::cerr << "[ERRORE] Impossibile aprire " << file_a << "\n";
        return;
    }
    TGraphErrors *g_sa = (TGraphErrors*)fa->Get("g_sigma_x");
    TGraph       *g_ba = (TGraph*)fa->Get("g_bias_x");

    TFile *fb = nullptr;
    TGraphErrors *g_sb = nullptr;
    TGraph       *g_bb = nullptr;
    if (file_b) {
        fb = TFile::Open(file_b);
        if (fb && !fb->IsZombie()) {
            g_sb = (TGraphErrors*)fb->Get("g_sigma_x");
            g_bb = (TGraph*)fb->Get("g_bias_x");
        }
    }

    gStyle->SetOptStat(0);

    // Canvas σ_x
    TCanvas *cs = new TCanvas("c_parallax_sigma",
                               "Sigma_x parallasse vs posizione", 900, 600);
    cs->SetGrid();
    g_sa->SetTitle(";x_{PMT3} [cm];#sigma_{x} parallasse [cm]");
    g_sa->Draw("APE");
    if (g_sb) {
        g_sb->Draw("PE same");
        TLegend *leg = new TLegend(0.62, 0.74, 0.88, 0.88);
        leg->AddEntry(g_sa, "Guida A (calibrazione)", "lp");
        leg->AddEntry(g_sb, "Guida B (TOF)",         "lp");
        leg->Draw();
    }
    cs->Update();

    // Canvas bias
    TCanvas *cb = new TCanvas("c_parallax_bias",
                               "Bias x parallasse vs posizione", 900, 600);
    cb->SetGrid();
    g_ba->SetTitle(";x_{PMT3} [cm];#LT x_{barra} - x_{PMT3} #GT [cm]");
    g_ba->Draw("APL");
    if (g_bb) {
        g_bb->Draw("PL same");
        TLegend *leg = new TLegend(0.62, 0.74, 0.88, 0.88);
        leg->AddEntry(g_ba, "Guida A (calibrazione)", "lp");
        leg->AddEntry(g_bb, "Guida B (TOF)",         "lp");
        leg->Draw();
    }
    cb->Update();
}


// ----------------------------------------------------------------------------
//  GetParallaxSigma
// ----------------------------------------------------------------------------
//
//  Helper che il codice di calibrazione può chiamare per ottenere σ_x
//  parallasse a una posizione specifica, leggendola dal file di scan.
//  Se la posizione richiesta non è esattamente nella griglia, fa
//  interpolazione lineare tra i due punti più vicini.
//
//  Restituisce -1.0 in caso di errore (file non trovato, posizione fuori range).
//
//  USO (dal codice di calibrazione):
//    double sx = GetParallaxSigma("parallax_GuidaA.root", -56.0);

double GetParallaxSigma(const char* file, double x_pmt3_centered)
{
    TFile *f = TFile::Open(file);
    if (!f || f->IsZombie()) {
        std::cerr << "[ERRORE] GetParallaxSigma: file non leggibile: "
                  << file << "\n";
        return -1.0;
    }

    TGraph *g = (TGraph*)f->Get("g_sigma_x");
    if (!g) {
        std::cerr << "[ERRORE] GetParallaxSigma: grafico g_sigma_x assente in "
                  << file << "\n";
        f->Close();
        return -1.0;
    }

    // TGraph::Eval fa interpolazione lineare automaticamente
    double sigma = g->Eval(x_pmt3_centered);
    f->Close();
    return sigma;
}


// ----------------------------------------------------------------------------
//  GetParallaxBias
// ----------------------------------------------------------------------------
//
//  Analogo di GetParallaxSigma per il bias di parallasse.
//  Utile se vuoi correggere x_PMT3 prima del fit:
//    x_corretto = x_PMT3 + bias_parallax(x_PMT3)
//
//  USO:
//    double bx = GetParallaxBias("parallax_GuidaA.root", -56.0);

double GetParallaxBias(const char* file, double x_pmt3_centered)
{
    TFile *f = TFile::Open(file);
    if (!f || f->IsZombie()) return 0.0;

    TGraph *g = (TGraph*)f->Get("g_bias_x");
    if (!g) { f->Close(); return 0.0; }

    double bias = g->Eval(x_pmt3_centered);
    f->Close();
    return bias;
}