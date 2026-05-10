#include "PhysicsList.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"

PhysicsList::PhysicsList() : FTFP_BERT()
{
    auto opticalPhysics = new G4OpticalPhysics();
    RegisterPhysics(opticalPhysics);

    auto optParams = G4OpticalParameters::Instance();
    optParams->SetProcessActivation("Scintillation", true);
    optParams->SetProcessActivation("Cerenkov", false);
    
}