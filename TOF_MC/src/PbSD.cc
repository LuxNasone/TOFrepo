#include "PbSD.hh"
#include "G4SystemOfUnits.hh"

PbSD::PbSD(const G4String& name)
: G4VSensitiveDetector(name) {}

G4bool PbSD::ProcessHits(G4Step* step, G4TouchableHistory*) {

    if (step->GetTrack()->GetTrackID() != 1) return false;

    G4double edep = step->GetTotalEnergyDeposit();

    if (edep <= 0.) return false;

    fEdepPb += edep;

    fTrackLength += step->GetStepLength();
    
    return true;
}