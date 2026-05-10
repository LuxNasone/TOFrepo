#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4OpticalPhoton.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction) {}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step){
    
    auto track = step->GetTrack();

    if (track->GetDefinition() == G4OpticalPhoton::Definition()) return;

    auto preVol = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!preVol) return;

    G4String vol = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

    G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();

    if (vol == "Bar")  {fEventAction->SetHitBarPos(pos);}

    if (vol == "Surf") {fEventAction->SetHitSurfPos(pos);}
    
}