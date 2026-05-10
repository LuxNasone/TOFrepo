#include <TSystem.h> 
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

void TOF(const char* path, const char* tree){

    TFile* out = new TFile("out.root", "RECREATE");

    ROOT::RDataFrame df(tree, path);

    auto tof_df = df.Define("TOF", "t3 - 0.5 * (t1 + t2)");
    
    auto h_tof = tof_df.Histo1D({"TOF", "Time of Flight", 50, 0, 20}, "TOF");

    h_tof->Write();

    out->Close();

}