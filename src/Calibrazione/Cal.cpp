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
#include <TGraph.h>
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

    auto new_df = df.Define("dt_12", [](ROOT::RVecF a){ return a[0] - a[1]; }, {"CFT"})
                    .Define("dt_13", [](ROOT::RVecF a){ return a[0] - a[2]; }, {"CFT"})
                    .Define("dt_23", [](ROOT::RVecF a){ return a[1] - a[2]; }, {"CFT"});

    auto h_12 = new_df.Histo1D({"dt12", "Differenze temporali 01-02", 100, -30, 30}, "dt_12");
    auto h_13 = new_df.Histo1D({"dt13", "Differenze temporali 01-03", 100, -30, 30}, "dt_13");
    auto h_23 = new_df.Histo1D({"dt23", "Differenze temporali 02-03", 100, -30, 30}, "dt_23");

    auto h12 = h_12.GetValue();
    auto h13 = h_13.GetValue();
    auto h23 = h_23.GetValue();

    std::vector<TH1D*> h_v = {h_12.GetPtr(), h_13.GetPtr(), h_23.GetPtr()};

    //TCanvas* c = new TCanvas(Form("%d", kw), Form("%d", kw), 800, 600);

    TF1 *f = new TF1("f", "gaus", -30, 30);

    f->SetParameters(h_v[kw]->GetMaximum(), h_v[kw]->GetBinCenter(h_v[kw]->GetMaximumBin()), h_v[kw]->GetStdDev());

    h_v[kw]->Fit(f, "ILSQ");

    //h_v[0]->Draw();

    //c->Update();

    //c->SaveAs(Form("%s.png", treename));

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

        int k = j - 5;

        if(point_names[j][0]  == 'X'){x.push_back(std::stod(point_names[j].substr(1)));}

        if(point_names[j][0]  == 'N'){x.push_back( - std::stod(point_names[j].substr(1)));}
        
        dx.push_back(2.45 + abs(k) * 0.1);

    }

    TFile* output = new TFile(out, "RECREATE");

    TCanvas* c_mean12 = new TCanvas("Differenze_temporali(1-2)", "Delta_t(1-2)", 800, 600);

    TCanvas* c_std12 = new TCanvas("Risoluzione_temporale(1-2)", "Sigma_t(1-2)", 800, 600);

    Points &line_to_fit12 = Lines[0];

    c_mean12->Divide(1,2,0,0);

    c_mean12->cd(1);

    TGraphErrors* g_mean12 = new TGraphErrors(x.size(), line_to_fit12.mean.data(), x.data(), line_to_fit12.mean_err.data(), dx.data());

    g_mean12->SetMarkerColor(1);
    g_mean12->SetLineColor(1);
    g_mean12->SetMarkerStyle(20);

    g_mean12->SetTitle("Ritardo (1-2) vs punto sulla sbarra;#Delta t_{12}[ns];x[cm]");

    g_mean12->Fit("pol1");
    TF1* f12 = g_mean12->GetFunction("pol1");

    g_mean12->GetXaxis()->SetTitle("");
    g_mean12->GetXaxis()->SetLabelSize(0);

    g_mean12->Draw("AP");

    c_mean12->Clear();
    c_mean12->Divide(1,2);

    c_mean12->cd(1);
    gPad->SetPad(0, 0.3, 1, 1);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.02);

    g_mean12->Draw("AP");

    TPaveText *pt = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.8);
    pt->SetLineColor(kBlack);
    pt->SetTextFont(42);
    pt->SetTextSize(0.06);
    pt->SetTextAlign(12); 

    pt->AddText(Form("Parametri di fit:"));
    pt->AddText(Form("m : %.2f +/- %.2f cm/ns", f12->GetParameter(1), f12->GetParError(1)));
    pt->AddText(Form("q : %.1f +/- %.1f cm", f12->GetParameter(0), f12->GetParError(0)));

    pt->Draw();

    c_mean12->cd(2);
    gPad->SetPad(0, 0.0, 1, 0.3);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.3);

    std::vector<double> residui(x.size());
    std::vector<double> residui_err(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        double fit_val = f12->Eval(line_to_fit12.mean[i]);
        residui[i] = x[i] - fit_val; 
        residui_err[i] = dx[i];       
    }

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

        TGraphErrors* g_mean = new TGraphErrors(x.size(), x.data(), line.mean.data(), dx.data(), line.mean_err.data());

        g_mean->SetMarkerColor(i+1);
        g_mean->SetLineColor(i+1);
        g_mean->SetMarkerStyle(20);

        g_mean->SetTitle("Punto vs Intervallo di tempo;#Delta t[ns];x[cm]");

        TF1* l = new TF1("line", "[0] + [1] * x", 0, 22);

        g_mean->Fit("line");

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

void W_bias(const char* path = "/home/lux_n/TOFrepo/OutFiles/scanW"){

    std::vector<double> res;
    std::vector<double> res_err;
    std::vector<double> w_v;
    std::vector<double> w_err;

    for (int i = 0; i < 10; i++){

        float w = 0.1 * (1 + i);

        w = std::round(w * 10) / 10.0f;

        w_v.push_back(w);
        w_err.push_back(1e-4);

        std::vector<double> v = H_Fit(Form("%.1f", w), Form("%s/%.1f.root", path, w), 0);

        res.push_back(v[2]);
        res_err.push_back(v[3]);

    }

    TCanvas* c_w = new TCanvas("W vs Sigma", "W vs Sigma", 800, 600);

    TGraphErrors* g_w = new TGraphErrors(w_v.size(), w_v.data(), res.data(), w_err.data(), res_err.data());

    g_w->GetXaxis()->SetLimits(0, 1.2);
    g_w->GetYaxis()->SetRangeUser(0, 5);
    g_w->SetMarkerStyle(20);
    g_w->SetMarkerSize(1.2);

    g_w->SetTitle("W vs Sigma;W;#sigma_{t} [ns]");

    g_w->Draw("AP");

    c_w->Update();

    c_w->SaveAs("W_graph.png");

    delete g_w;
    delete c_w;

}

void IntrinsicRes(const char* path, const char* out){

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

        int k = j - 5;

        if(point_names[j][0]  == 'X'){x.push_back(std::stod(point_names[j].substr(1)));}

        if(point_names[j][0]  == 'N'){x.push_back( - std::stod(point_names[j].substr(1)));}
        
        dx.push_back(2.8 + abs(k) * 0.1);

    }

    std::vector<double> sigma_1;
    std::vector<double> sigma_2;
    std::vector<double> sigma_3;

    for (size_t j = 0; j < point_names.size(); j++){

        double s_1 = sqrt(0.5 * (pow(Lines[0].std[j], 2) + pow(Lines[2].std[j],2) - pow(Lines[1].std[j], 2))); 

        std::cout << s_1 <<std::endl;

        sigma_1.push_back(s_1);

        double s_2 = sqrt(0.5 * (pow(Lines[1].std[j], 2) + pow(Lines[0].std[j], 2) - pow(Lines[2].std[j], 2)));

        std::cout << s_2 <<std::endl; 

        sigma_2.push_back(s_2);

        double s_3 = sqrt(0.5 * (abs(pow(Lines[1].std[j], 2) + pow(Lines[2].std[j], 2) - pow(Lines[0].std[j], 2))));

        std::cout << s_3 <<std::endl;

        sigma_3.push_back(s_3);

    }

    std::vector<std::vector<double>> sigma = {sigma_1, sigma_2, sigma_3};

    TFile* output = new TFile(out, "RECREATE");

    for (size_t j = 0; j < sigma.size(); j++){

        TCanvas* c = new TCanvas(Form("PMT0%d", (int)j+1), Form("PMT0%d", (int)j+1), 800, 600);

        TGraph* g = new TGraph(x.size(), x.data(), sigma[j].data());

        g->SetMarkerColor(1);
        g->SetLineColor(1);
        g->SetMarkerStyle(20);

        g->SetTitle(Form("Punto vs Incertezza intrinseca PMT0%d;x[cm]; #sigma_t", (int)j + 1));

        g->Draw("AP");

        c->Update();

        c->Write();

        delete c;

        delete g;

    }

    output->Close();

    gROOT->SetBatch(kFALSE);

}