/*

Codice per la letturra dei file del DRS:
    1)Salva la forma d'onda in due vettori tempo-voltaggio
    2)Implementa metodi per calcolo dell'area, istante di discesa e altre grandezze fisiche di interesse
    3)implementazione successiva : creare una funzione che ne faccia derivata e integrale

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
#include "WF.h"

void DRSread(const char* fname, const std::string &outname, const char* opt = "UPDATE"){

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

    std::vector<double> Dtime[3];

    std::vector<double> Rtime[3];

    for(int i = 0; i < N; i++){

        Event e = all_events[i];
        
        for (int ev_ind = 0; ev_ind < 3; ev_ind++){

            WF wf = e.GetChannel(ev_ind);

            BL[ev_ind].push_back(wf.Baseline(10));

            A[ev_ind].push_back(wf.Amp(10));

            Dtime[ev_ind].push_back(wf.DropTime(0.5, 10));

            Rtime[ev_ind].push_back(wf.RiseTime(10, 10));

        }

    }

    TFile* f = new TFile("OutFiles/Points.root", opt);

    TTree* t =  new TTree(outname.c_str(), outname.c_str());

    Float_t bl[3], amp[3], dt[3], rt[3];

    t->Branch("BL", bl, "BL[3]/F");
    t->Branch("Amp", amp, "Amp[3]/F");
    t->Branch("DropTime", dt, "DropTime[3]/F");
    t->Branch("RiseTime", rt, "RiseTime[3]/F");


    for (int i = 0; i < N; i++){

        for (int ev_ind = 0; ev_ind < 3; ev_ind++) {


            bl[ev_ind] = BL[ev_ind][i];

            amp[ev_ind] = A[ev_ind][i];

            dt[ev_ind] = Dtime[ev_ind][i];

            rt[ev_ind] = Rtime[ev_ind][i];

        }

        t->Fill();

    }

    t->Write();

    f->Close();

}