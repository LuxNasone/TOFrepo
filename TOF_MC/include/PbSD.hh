#ifndef PbSD_hh
#define PbSD_hh

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class PbSD : public G4VSensitiveDetector {

    public:

        PbSD(const G4String& name);

        ~PbSD() override = default;

        G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

        G4double GetEdep() const { return fEdepPb; }

        G4double GetTrackLength() const { return fTrackLength; }

        void Reset() { fEdepPb = 0.; fTrackLength = 0.; }

    private:

        G4double fEdepPb = 0.;
        
        G4double fTrackLength = 0.;
};

#endif