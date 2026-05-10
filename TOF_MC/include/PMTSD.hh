#ifndef PMTSD_HH
#define PMTSD_HH

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;
class G4TouchableHistory;
class EventAction;

class PMTSD : public G4VSensitiveDetector{
    
    public:

        PMTSD(const G4String& name, EventAction* eventAction);

        virtual ~PMTSD();

        virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    
    private:

        EventAction* fEventAction;
        
};

#endif