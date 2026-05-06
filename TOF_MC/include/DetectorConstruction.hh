#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction{

    public:

        DetectorConstruction();

        virtual ~DetectorConstruction();

        virtual G4VPhysicalVolume* Construct() override;

        virtual void ConstructSDandField() override;

    private:
    
        G4LogicalVolume* fScoringVolume = nullptr;

        G4LogicalVolume* fLogicPMT = nullptr;

};

#endif