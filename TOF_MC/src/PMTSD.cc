#include "PMTSD.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

PMTSD::PMTSD(const G4String& name, EventAction* eventAction)
: G4VSensitiveDetector(name), fEventAction(eventAction) {}

PMTSD::~PMTSD() {}
 
G4bool PMTSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
 
    if (!fEventAction) return false;

    auto track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::Definition())
        return false;

    G4double time = step->GetPreStepPoint()->GetGlobalTime();

    auto touch = step->GetPreStepPoint()->GetTouchableHandle();
    G4int copyNo = touch->GetCopyNumber(1);

    if (copyNo == 0) fEventAction->SetTimeBar1(time);
    else if (copyNo == 1) fEventAction->SetTimeBar2(time);
    else if (copyNo == 2) fEventAction->SetTimeBarSmall(time);

    track->SetTrackStatus(fStopAndKill);
    return true;
}

