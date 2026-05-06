#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
    
    public:

        PrimaryGeneratorAction();

        virtual ~PrimaryGeneratorAction();

        virtual void GeneratePrimaries(G4Event* anEvent);

    private:

        G4ParticleGun* fParticleGun;

        G4double SampleMuonEnergy();
};

#endif