/*

Codice per la letturra dei file del DRS:
    1)Salva la forma d'onda in due vettori tempo-voltaggio
    2)Implementa metodi per calcolo dell'area, istante di discesa e altre grandezze fisiche di interesse

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream> 
#include <algorithm>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include "WF.h"

std::vector<std::string> point_names = {"N130", "N112", "N84", "N56", "N28", "X0", "X28", "X56", "X84", "X112", "X130"};

//Reading function, it converts a DRS file in a TTree

void DRSread(const char* fname, 
            const std::string &outname,
             const char* outfile = "OutFiles/Points.root",
             const char* opt = "UPDATE",
             int n = 100,
             double w = 0.5,
             double k = 0.6,
             double s = 1.25){

    //Import del file con check di apertura

    std::ifstream file(fname);

    if (!file.is_open()){ std::cout << "Errore di apertura" << std::endl; }

    //Loop principale

    std::string line; 

    //Gruppo di eventi

    std::vector<Event> all_events;
    Event current_event;
    WF* current_wf = nullptr;

    while (std::getline(file, line)) {

        if (line.find("<Event>") != std::string::npos) {

            current_event = Event(); 

            current_wf = nullptr;

        }

        if (line.find("<Time>") != std::string::npos){

            size_t start = line.find(">") + 1;
            size_t end = line.find("</");
            std::string stamp = line.substr(start, end-start);

            if (current_wf) current_wf->SetTimeStamp(stamp);
            
        }

        if (line.find("<CHN") != std::string::npos){

            size_t start = line.find("N") + 1;
            size_t end = line.find(">", start);
            std::string ch_str = line.substr(start, end-start);

            int ch = std::stoi(ch_str);

            current_wf = &current_event.waveforms[ch-1]; 
            current_wf->SetChannel(ch);
        
        }

        if (line.find("<Data>") != std::string::npos){

            size_t start = line.find(">") + 1;
            size_t end = line.find("</");
            std::string content = line.substr(start, end - start);

            std::stringstream ss(content);
            std::string t_str, V_str;

            getline(ss, t_str, ',');
            getline(ss, V_str, ',');

            double t_val = std::stod(t_str);
            double V_val = std::stod(V_str);

            if (current_wf) current_wf->PushBack(t_val, V_val);

        }
        
        else if (line.find("</Event>") != std::string::npos) { all_events.push_back(current_event); }
    }

    int N = all_events.size();

    std::vector<double> BL[3];
    
    std::vector<double> A[3];

    std::vector<double> CFT[3];

    std::vector<double> ARC[3];

    std::vector<double> MDer[3];

    for(int i = 0; i < N; i++){

        Event& e = all_events[i];

        if(e.IsBad(n, 500, -500)){
        
            for (int ev_ind = 0; ev_ind < 3; ev_ind++){

                const WF& wf = e.GetChannel(ev_ind);

                BL[ev_ind].push_back(wf.Baseline(n));

                A[ev_ind].push_back(wf.Amp(n));

                CFT[ev_ind].push_back(wf.CFT(w, n));

                double rise_time = wf.CFT(0.9, n) - wf.CFT(0.1, n);

                ARC[ev_ind].push_back(wf.ARC(k, s * rise_time));

                MDer[ev_ind].push_back(wf.maxDer());

            }
        }
    }

    TFile* f = new TFile(outfile, opt);

    TTree* t =  new TTree(outname.c_str(), outname.c_str());

    Float_t bl[3], amp[3], cft[3], arc[3], mder[3];

    t->Branch("BL", bl, "BL[3]/F");
    t->Branch("Amp", amp, "Amp[3]/F");
    t->Branch("CFT", cft, "CFT[3]/F");
    t->Branch("ARC", arc, "ARC[3]/F");
    t->Branch("MDer", mder, "MDer[3]/F");

    for (int i = 0; i < N; i++){

        for (int ev_ind = 0; ev_ind < 3; ev_ind++) {


            bl[ev_ind] = BL[ev_ind][i];

            amp[ev_ind] = A[ev_ind][i];

            cft[ev_ind] = CFT[ev_ind][i];

            arc[ev_ind] = ARC[ev_ind][i];

            mder[ev_ind] = MDer[ev_ind][i];

        }

        t->Fill();

    }

    t->Write();

    f->Close();

}

//It converts all file with name in point_names in TTree. Just run ReadFolder() to have an analysis, but it is unoptimized 

void ReadFolder(const char* outname = "/home/lux_n/TOFrepo/OutFiles/Points.root", int n = 10, double w = 0.5, double k = 0.6, double s = 0.75){

    for(size_t i = 0; i < point_names.size(); i++){

        const char* opt;

        if (i == 0){opt = "RECREATE";}
        else {opt = "UPDATE";}

        DRSread(Form("/home/lux_n/TOF_DRS/%s.xml", point_names[i].c_str()), point_names[i].c_str(), outname, opt, n, w, k, s);
    }

}

