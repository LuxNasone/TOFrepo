#include "PhysicsList.hh"
#include "G4OpticalPhysics.hh"

PhysicsList::PhysicsList()
: FTFP_BERT()
{

    auto opticalPhysics = new G4OpticalPhysics();

    RegisterPhysics(opticalPhysics);
}