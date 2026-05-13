#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"

RunAction::RunAction(){

    auto analysis = G4AnalysisManager::Instance();

    analysis->SetVerboseLevel(1);

}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*){

    auto analysis = G4AnalysisManager::Instance();

    analysis->OpenFile("/home/lux_n/TOFrepo/TOF_MC/TOF.root");

    analysis->CreateNtuple("TOF", "TOF data");

    analysis->CreateNtupleDColumn("t_PMT01");

    analysis->CreateNtupleDColumn("t_PMT02");

    analysis->CreateNtupleDColumn("t_PMT03");

    analysis->CreateNtupleDColumn("x_Bar");

    analysis->CreateNtupleDColumn("y_Bar");

    analysis->CreateNtupleDColumn("z_Bar");

    analysis->CreateNtupleDColumn("x_Surf");

    analysis->CreateNtupleDColumn("y_Surf");

    analysis->CreateNtupleDColumn("z_Surf");

    analysis->CreateNtupleDColumn("t_PMT04");

    analysis->CreateNtupleDColumn("E_Pb");

    analysis->CreateNtupleDColumn("TrackLengthPb");

    analysis->FinishNtuple();
    
}

void RunAction::EndOfRunAction(const G4Run*){

    auto analysis = G4AnalysisManager::Instance();

    analysis->Write();

    analysis->CloseFile();

}