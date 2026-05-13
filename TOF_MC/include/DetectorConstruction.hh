#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"

class EventAction;

class DetectorConstruction : public G4VUserDetectorConstruction{

    public:

        DetectorConstruction(EventAction* eventAction);

        virtual ~DetectorConstruction();

        virtual G4VPhysicalVolume* Construct() override;

        virtual void ConstructSDandField() override;

    private:
    
        G4LogicalVolume* fScoringVolume;

        G4LogicalVolume* fLogicCathode;

        G4LogicalVolume* fLogicPb = nullptr;

        EventAction* fEventAction;

};

#endif