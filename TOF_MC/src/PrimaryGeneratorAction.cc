#include "PrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4MuonMinus.hh"
#include "CLHEP/Units/PhysicalConstants.h"

PrimaryGeneratorAction::PrimaryGeneratorAction(){

    fParticleGun = new G4ParticleGun(1);

    // default particle (ONLY ONCE)
    fParticleGun->SetParticleDefinition(G4MuonMinus::MuonMinusDefinition());

    fParticleGun->SetParticleEnergy(3.0*GeV);

    // start above detector (safe position)
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 1.2*m));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {delete fParticleGun;}

G4double PrimaryGeneratorAction::SampleMuonEnergy(){
    G4double Emin = 0.5*GeV;
    G4double Emax = 100.0*GeV;
    G4double gamma = 2.7;

    G4double u = G4UniformRand();

    G4double termMin = std::pow(Emin, 1.0 - gamma);
    G4double termMax = std::pow(Emax, 1.0 - gamma);

    return std::pow(termMin + u*(termMax - termMin), 1.0/(1.0 - gamma));
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

    // energy sampling
    fParticleGun->SetParticleEnergy(SampleMuonEnergy());

    // isotropic-ish angular distribution (cos²-like)
    G4double u = G4UniformRand();
    G4double theta = std::acos(std::sqrt(1.0 - u));
    G4double phi = 2.0 * CLHEP::pi * G4UniformRand();

    G4ThreeVector dir(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));

    fParticleGun->SetParticleMomentumDirection(dir);

    // uniform over detector area
    G4double x = (G4UniformRand() - 0.5) * 3.0 * m;
    G4double y = (G4UniformRand() - 0.5) * 1.0 * m;

    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, 1.2*m));

    fParticleGun->GeneratePrimaryVertex(anEvent);
}