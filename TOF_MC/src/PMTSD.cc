#include "PMTSD.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

PMTSD::PMTSD(const G4String& name)
: G4VSensitiveDetector(name) {}

PMTSD::~PMTSD() {}

G4bool PMTSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4Track* track = step->GetTrack();

    // only optical photons
    if(track->GetDefinition() != G4OpticalPhoton::Definition())
        return false;

    // only entry into sensitive volume
    if(step->GetPreStepPoint()->GetStepStatus() != fGeomBoundary)
        return false;

    G4double time = step->GetPreStepPoint()->GetGlobalTime();

    // PMT id
    G4int copyNo =
        step->GetPreStepPoint()
        ->GetTouchableHandle()
        ->GetCopyNumber();

    // ---- FIRST PHOTON LOGIC (IMPORTANT FIX) ----
    static std::map<G4int, G4double> firstHitTime;

    if(firstHitTime.find(copyNo) == firstHitTime.end())
    {
        firstHitTime[copyNo] = time;
    }
    else
    {
        return false; // ignore later photons
    }

    // ---- SAVE TO ROOT ----
    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillNtupleIColumn(0, copyNo);
    analysisManager->FillNtupleDColumn(1, time/ns);
    analysisManager->AddNtupleRow();

    track->SetTrackStatus(fStopAndKill);

    return true;
}