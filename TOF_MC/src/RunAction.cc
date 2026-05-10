#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"

RunAction::RunAction()
{
    auto analysis = G4AnalysisManager::Instance();

    analysis->SetVerboseLevel(1);

}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*)
{
    auto analysis = G4AnalysisManager::Instance();

    analysis->OpenFile("/home/lux_n/TOFrepo/TOF_MC/tof.root");

    analysis->CreateNtuple("tof", "TOF data");

    analysis->CreateNtupleDColumn("t1");
    analysis->CreateNtupleDColumn("t2");
    analysis->CreateNtupleDColumn("t3");

    analysis->FinishNtuple();
}

void RunAction::EndOfRunAction(const G4Run*)
{
    auto analysis = G4AnalysisManager::Instance();

    analysis->Write();
    analysis->CloseFile();
}