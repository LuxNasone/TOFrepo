#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "CLHEP/Units/PhysicalConstants.h"

PrimaryGeneratorAction::PrimaryGeneratorAction(){

    //Generating beam

    fParticleGun = new G4ParticleGun(1);

    //Beam energy

    fParticleGun->SetParticleEnergy(3.0*GeV);

    //Beam origin

    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 2.5*m));

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

G4double PrimaryGeneratorAction::SampleMuonEnergy(){

    //Minimum energy around muon mass

    G4double Emin = 110.0 * MeV;

    //Maximum energy

    G4double Emax = 100.0 * GeV;

    //Energy spectrum index

    G4double gamma = 2.7;

    //Sampling energies

    G4double u = G4UniformRand();

    G4double termMin = std::pow(Emin, 1.0 - gamma);

    G4double termMax = std::pow(Emax, 1.0 - gamma);

    return std::pow(termMin + u*(termMax - termMin), 1.0/(1.0 - gamma));

}

G4ThreeVector PrimaryGeneratorAction::SampleMuonDir(){

    G4double txBar = (G4UniformRand() - 0.5) * 2.795*m;

    G4double tyBar = (G4UniformRand() - 0.5) * 4*cm;

    G4ThreeVector startPos(txBar, tyBar, 165*cm);

    fParticleGun->SetParticlePosition(startPos);

    // Slab centrata in (0, sy/2, 0) con dimensioni sx x sy x sz

    G4double sx = 5.9*cm;

    G4double sy = 20.6*cm;
    
    // Target: punto casuale sulla faccia superiore della slab

    G4double txSlab = (G4UniformRand() - 0.5) * sx;

    G4double tySlab = (G4UniformRand()) * sy; // centro + offset
    
    G4double tzSlab = 1.1*cm / 2; // faccia superiore della slab

    G4ThreeVector targetPoint(txSlab, tySlab, tzSlab);

    G4ThreeVector dir = (targetPoint - startPos).unit();

    return dir;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

    //Negative muon generator

    if (G4UniformRand() > 0.44) {fParticleGun->SetParticleDefinition(G4MuonMinus::MuonMinusDefinition());}

    //Positive muon generator

    else {fParticleGun->SetParticleDefinition(G4MuonPlus::MuonPlusDefinition());}

    //Energy sampling

    fParticleGun->SetParticleEnergy(SampleMuonEnergy());

    //Direction sampling

    fParticleGun->SetParticleMomentumDirection(SampleMuonDir());

    //Generating event

    fParticleGun->GeneratePrimaryVertex(anEvent);
}