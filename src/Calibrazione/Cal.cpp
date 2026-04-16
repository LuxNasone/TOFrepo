#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include <TSystem.h> 
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TLegend.h>

std::vector<std::string> point_names = {"N130", "N112", "N84", "N56", "N28", "X0", "X28", "X56", "X84", "X112", "X130"};

const int n_h = 3;

struct Points{

    std::vector<double> mean;
    std::vector<double> mean_err;
    std::vector<double> std;
    std::vector<double> std_err;

};

std::vector<double> H_Fit(const char* treename, const char* path, int kw){

    ROOT::RDataFrame df(treename, path);

    auto new_df = df.Define("dt_12", [](ROOT::RVecF a){ return a[0] - a[1]; }, {"DropTime"})
                    .Define("dt_13", [](ROOT::RVecF a){ return a[0] - a[2]; }, {"DropTime"})
                    .Define("dt_23", [](ROOT::RVecF a){ return a[1] - a[2]; }, {"DropTime"});

    auto h_12 = new_df.Histo1D({"dt12", "Differenze temporali 01-02", 100, -20, 20}, "dt_12");
    auto h_13 = new_df.Histo1D({"dt13", "Differenze temporali 01-03", 100, -20, 20}, "dt_13");
    auto h_23 = new_df.Histo1D({"dt23", "Differenze temporali 02-03", 100, -20, 20}, "dt_23");

    auto h12 = h_12.GetValue();
    auto h13 = h_13.GetValue();
    auto h23 = h_23.GetValue();

    TCanvas* C = new TCanvas(treename, treename, 800, 600);

    h12.Draw();
    h13.Draw("SAME");
    h23.Draw("SAME");

    C->SaveAs(Form("%s.png", treename));

    std::vector<TH1D*> h_v = {h_12.GetPtr(), h_13.GetPtr(), h_23.GetPtr()};

    TF1 *f = new TF1("f", "breitwigner", -20, 20);

    f->SetParameters(h_v[kw]->GetMaximum(), h_v[kw]->GetBinCenter(h_v[kw]->GetMaximumBin()), h_v[kw]->GetStdDev());

    h_v[kw]->Fit(f, "ILSQ");

    std::vector<double> result = {f->GetParameter(1), f->GetParError(1), f->GetParameter(2), f->GetParError(2)};

    return result;

}

void Lin_Fit(const char* path, const char* out){

    gROOT->SetBatch(kTRUE);

    Points Lines[n_h];

    for(size_t i = 0; i < point_names.size(); i++){

        for(int kw = 0; kw < n_h; kw++){

            std::vector<double> result = H_Fit(point_names[i].c_str(), path, kw);

            Points &current_line = Lines[kw];

            current_line.mean.push_back(result[0]);

            current_line.mean_err.push_back(result[1]);

            current_line.std.push_back(result[2]);

            current_line.std_err.push_back(result[3]);

        }
    }

    std::vector<double> x;
    std::vector<double> dx;

    for (size_t j = 0; j < point_names.size(); j++){

        if(point_names[j][0]  == 'X'){x.push_back(std::stod(point_names[j].substr(1)));}

        if(point_names[j][0]  == 'N'){x.push_back( - std::stod(point_names[j].substr(1)));}
        
        dx.push_back(2.45);

    }

    TFile* output = new TFile(out, "RECREATE");

    TCanvas* c_mean12 = new TCanvas("Differenze_temporali(1-2)", "Delta_t(1-2)", 800, 600);

    TCanvas* c_std12 = new TCanvas("Risoluzione_temporale(1-2)", "Sigma_t(1-2)", 800, 600);

    Points &line_to_fit12 = Lines[0];

    c_mean12->Divide(1,2);

    c_mean12->cd(1);

    TGraphErrors* g_mean12 = new TGraphErrors(x.size(), line_to_fit12.mean.data(), x.data(), line_to_fit12.mean_err.data(), dx.data());

    g_mean12->SetMarkerColor(1);
    g_mean12->SetLineColor(1);
    g_mean12->SetMarkerStyle(20);

    g_mean12->SetTitle("Ritardo (1-2) vs punto sulla sbarra;#Delta t_{12}[ns];x[cm]");

    g_mean12->Fit("pol1");
    TF1* f12 = g_mean12->GetFunction("pol1");

    g_mean12->Draw("AP");

    c_mean12->Clear();
    c_mean12->Divide(1,2);

    c_mean12->cd(1);
    gPad->SetPad(0, 0.3, 1, 1);
    gPad->SetBottomMargin(0.02);

    g_mean12->Draw("AP");

    TPaveText *pt = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.8);
    pt->SetLineColor(kBlack);
    pt->SetTextFont(42);
    pt->SetTextSize(0.06);
    pt->SetTextAlign(12); 

    pt->AddText(Form("Parametri di fit:"));
    pt->AddText(Form("c : %.3f +/- %.3f cm/ns", f12->GetParameter(1), f12->GetParError(1)));
    pt->AddText(Form("q : %.3f +/- %.3f cm", f12->GetParameter(0), f12->GetParError(0)));

    pt->Draw();

    c_mean12->cd(2);
    gPad->SetPad(0, 0.0, 1, 0.3);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.3);

    // Calcolo residui
    std::vector<double> residui(x.size());
    std::vector<double> residui_err(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        double fit_val = f12->Eval(line_to_fit12.mean[i]);
        residui[i] = x[i] - fit_val;   // y - fit(x)
        residui_err[i] = dx[i];        // errore su y
    }

    // Grafico residui
    TGraphErrors* g_res = new TGraphErrors(
        x.size(),
        line_to_fit12.mean.data(),
        residui.data(),
        line_to_fit12.mean_err.data(),
        residui_err.data()
    );

    g_res->SetMarkerStyle(20);
    g_res->SetTitle(";#Delta t_{12}[ns];Residui [cm]");

    g_res->Draw("AP");

    // Linea a zero
    TF1* zero = new TF1("zero", "0", 
        *std::min_element(line_to_fit12.mean.begin(), line_to_fit12.mean.end()),
        *std::max_element(line_to_fit12.mean.begin(), line_to_fit12.mean.end())
    );
    zero->SetLineColor(kRed);
    zero->Draw("same");

    c_mean12->Update();

    c_mean12->Write("Mean_graphs(1-2)");

    c_std12->cd();

    TGraphErrors* g_std12 = new TGraphErrors(x.size(), x.data(), line_to_fit12.std.data(), dx.data(), line_to_fit12.std_err.data());

    g_std12->SetMarkerColor(1);
    g_std12->SetLineColor(1);
    g_std12->SetMarkerStyle(20);

    g_std12->SetTitle("Risoluzione temporale (1-2) vs piunto sulla sbarra;x[cm];#sigma_{t_{12}}[ns]");

    g_std12->Draw("AP");

    c_std12->Update();

    c_std12->Write("STD_graphs(1-2)");

    TCanvas* c_mean = new TCanvas("Differenze_temporali", "Delta_t", 800, 600);

    TCanvas* c_std = new TCanvas("Risoluzione_temporale", "Sigma_t", 800, 600);

    TLegend* legend  = new TLegend(0.6, 0.7, 0.9, 0.9);

    std::vector<std::string> label = {"PMT01 - PMT03", "PMT02 - PMT03"};

    for (int i = 1; i < n_h; i++){

        Points &line = Lines[i];

        c_mean->cd();

        TGraphErrors* g_mean = new TGraphErrors(x.size(), line.mean.data(), x.data(), line.mean_err.data(), dx.data());

        g_mean->SetMarkerColor(i+1);
        g_mean->SetLineColor(i+1);
        g_mean->SetMarkerStyle(20);

        g_mean->SetTitle("Punto vs Intervallo di tempo;#Delta t[ns];x[cm]");

        g_mean->Fit("pol1");

        if (i == 1){g_mean->Draw("AP");}
        else{g_mean->Draw("P SAME");}

        legend->AddEntry(g_mean, label[i - 1].c_str(), "lep");

        legend->Draw();

        c_mean->Update();

        c_std->cd();

        TGraphErrors* g_std = new TGraphErrors(x.size(), x.data(), line.std.data(), dx.data(), line.std_err.data());

        g_std->SetMarkerColor(i+1);
        g_std->SetLineColor(i+1);
        g_std->SetMarkerStyle(20);

        g_std->SetTitle("Punto vs Risoluzione temporale;x[cm];#sigma_{t}[ns]");

        if (i == 1){g_std->Draw("AP");}
        else{g_std->Draw("P SAME");}

        legend->Draw();

        c_std->Update();

    }

    c_mean->Write("Mean_graphs");
    c_std->Write("STD_graphs");

    output->Close();

    gROOT->SetBatch(kFALSE);

}