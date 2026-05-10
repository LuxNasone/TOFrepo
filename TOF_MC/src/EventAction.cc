#include "EventAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

EventAction::EventAction()
: G4UserEventAction(),
  fTimeBar1(-1.),
  fTimeBar2(-1.),
  fTimeBarSmall(-1.)
{}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*){
    
    fTimeBar1 = -1.;
    fTimeBar2 = -1.;
    fTimeBarSmall = -1.;
}

void EventAction::EndOfEventAction(const G4Event*)
{

    if (fTimeBar1 < 0 || fTimeBar2 < 0 || fTimeBarSmall < 0) return;
    
    auto analysis = G4AnalysisManager::Instance();

    analysis->FillNtupleDColumn(0, fTimeBar1 / ns);

    analysis->FillNtupleDColumn(1, fTimeBar2 / ns);

    analysis->FillNtupleDColumn(2, fTimeBarSmall / ns);
    
    analysis->AddNtupleRow();
}