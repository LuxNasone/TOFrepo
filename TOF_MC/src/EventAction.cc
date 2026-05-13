#include "EventAction.hh"
#include "PbSD.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"

EventAction::EventAction()
: G4UserEventAction(), fTimeBar1(-1.), fTimeBar2(-1.), fTimeBarSmall(-1.) {}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*){
    
    fTimeBar1 = -1.;

    fTimeBar2 = -1.;

    fTimeBarSmall = -1.;

    fTimeVeto = -1.;

    fEdepPb = 0.; 

    fTrackLenPb  = 0.;

    fHitBarPos  = G4ThreeVector(-999,-999,-999);

    fHitSurfPos = G4ThreeVector(-999,-999,-999);
}

void EventAction::EndOfEventAction(const G4Event*){

    if ((fTimeBar1 < 0 || fTimeBar2 < 0 || fTimeBarSmall < 0) && (fTimeBarSmall < fTimeBar1) && (fTimeBarSmall < fTimeBar2)) return;
    if (fHitBarPos.x() == -999 || fHitSurfPos.x() == -999) return;
    if (fHitBarPos.y() == -999 || fHitSurfPos.y() == -999) return;
    if (fHitBarPos.z() == -999 || fHitSurfPos.z() == -999) return;

    auto sdMan = G4SDManager::GetSDMpointer();

    auto pbSD  = static_cast<PbSD*>(sdMan->FindSensitiveDetector("PbSD"));

    fEdepPb = 0.;

    fTrackLenPb  = 0.;

    if (pbSD) {

        fEdepPb = pbSD->GetEdep();

        fTrackLenPb = pbSD->GetTrackLength();

        pbSD->Reset(); 

    }
    
    auto analysis = G4AnalysisManager::Instance();

    analysis->FillNtupleDColumn(0, fTimeBar1 / ns);

    analysis->FillNtupleDColumn(1, fTimeBar2 / ns);

    analysis->FillNtupleDColumn(2, fTimeBarSmall / ns);

    analysis->FillNtupleDColumn(3, fHitBarPos.x()  / cm);

    analysis->FillNtupleDColumn(4, fHitBarPos.y()  / cm);

    analysis->FillNtupleDColumn(5, fHitBarPos.z()  / cm);

    analysis->FillNtupleDColumn(6, fHitSurfPos.x() / cm);

    analysis->FillNtupleDColumn(7, fHitSurfPos.y() / cm);

    analysis->FillNtupleDColumn(8, fHitSurfPos.z() / cm);

    analysis->FillNtupleDColumn(9, fTimeVeto / ns);

    analysis->FillNtupleDColumn(10, fEdepPb / MeV);

    analysis->FillNtupleDColumn(11, fTrackLenPb / cm); 

    analysis->AddNtupleRow();

}