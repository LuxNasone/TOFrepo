#include "SensitiveDetector.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhoton.hh"

SensitiveDetector::SensitiveDetector(const G4String& name)
: G4VSensitiveDetector(name),
  fTotalEdep(0.)
{}

SensitiveDetector::~SensitiveDetector() {}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{

    G4double edep = step->GetTotalEnergyDeposit();
    if(edep <= 0.) return false;

    fTotalEdep += edep;
    return true;
}