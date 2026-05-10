#include <TSystem.h> 
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

#include <cmath>

void TOF(const char* path, const char* tree){

    TFile* out = new TFile("/home/lux_n/TOFrepo/TOF_MC/out.root", "RECREATE");

    ROOT::RDataFrame df(tree, path);

    auto tof_df = df.Filter("t_PMT03 > 0 && t_PMT02 > 0 && t_PMT01 > 0")
                    .Filter("std::abs(x_Bar) < 139.75 && std::abs(y_Bar) < 2")
                    .Filter("std::abs(x_Surf) < 2.95 && std::abs(y_Bar - 10.3) < 10.3")
                    .Define("t1", "t_PMT01 - std::abs(139.75 + x_Bar)/18.9873")
                    .Define("t2", "t_PMT02 - std::abs(139.75 - x_Bar)/18.9873")
                    .Define("t3", "t_PMT03 - std::abs(10.3 - y_Surf)/18.9873")
                    .Define("Time_of_Flight", "t3 - 0.5 * (t1 + t2)")
                    .Filter("Time_of_Flight > 5.3")
                    .Define("Distance", "std::pow(std::pow(x_Bar - x_Surf, 2) + std::pow(y_Bar - y_Surf, 2) + std::pow(z_Bar - z_Surf, 2), 0.5)")
                    .Define("Beta", "Distance/(30 * Time_of_Flight)");
    
    auto h_tof = tof_df.Histo1D({"Time_of_Flight", "Time_of_Flight", 100, 5, 25}, "Time_of_Flight");

    auto h_D = tof_df.Histo1D({"Distance", "Distance", 100, 150, 220}, "Distance");

    auto h_beta = tof_df.Histo1D({"Beta", "Beta", 100, 0.2, 1.5}, "Beta");

    h_tof->Write();

    h_D->Write();
    
    h_beta->Write();

    out->Close();

}