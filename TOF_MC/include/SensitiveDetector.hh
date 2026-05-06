#ifndef SensitiveDetector_h
#define SensitiveDetector_h

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4TouchableHistory;

class SensitiveDetector : public G4VSensitiveDetector{

    public:

        SensitiveDetector(const G4String& name);

        virtual ~SensitiveDetector();

        virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

    private:
    
        G4double fTotalEdep;
};

#endif