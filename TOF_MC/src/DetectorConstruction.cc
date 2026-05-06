#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"
#include "PMTSD.hh"
#include "EventAction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SDManager.hh"

DetectorConstruction::DetectorConstruction() {}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct(){

    auto nist = G4NistManager::Instance();

    // WORLD
    G4double worldSize = 4*m;

    auto worldMat = nist->FindOrBuildMaterial("G4_AIR");

    auto solidWorld = new G4Box("World", worldSize/2, worldSize/2, worldSize/2);
    auto logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    auto physWorld  = new G4PVPlacement(nullptr, {}, logicWorld, "World", nullptr, false, 0);

    // ENVELOPE
    G4double envSize = 3*m;
    auto solidEnv = new G4Box("Env", envSize/2, envSize/2, envSize/2);
    auto logicEnv = new G4LogicalVolume(solidEnv, worldMat, "Env");

    new G4PVPlacement(nullptr, {}, logicEnv, "Env", logicWorld, false, 0);

    // SCINTILLATOR MATERIAL (BC408)
    auto BC408 = new G4Material("BC408", 1.032*g/cm3, 2);
    auto C = nist->FindOrBuildElement("C");
    auto H = nist->FindOrBuildElement("H");

    BC408->AddElement(C, 9);
    BC408->AddElement(H, 10);

    const G4int n = 2;

    G4double e[n] = {2.0*eV, 3.5*eV};

    G4double rindex[n] = {1.58, 1.58};
    G4double abslen[n] = {2.1*m, 2.1*m};

    auto mpt = new G4MaterialPropertiesTable();

    mpt->AddProperty("RINDEX", e, rindex, n);
    mpt->AddProperty("ABSLENGTH", e, abslen, n);

    mpt->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1*ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);

    BC408->SetMaterialPropertiesTable(mpt);

    // BAR
    G4double barX = 2.795*m;
    G4double barYZ = 4*cm;

    auto solidBar = new G4Box("Bar", barX/2, barYZ/2, barYZ/2);
    auto logicBar = new G4LogicalVolume(solidBar, BC408, "Bar");

    new G4PVPlacement(nullptr, {0,0,73.5*cm}, logicBar, "Bar", logicEnv, false, 0);

    fScoringVolume = logicBar;

    // SMALL SCINTILLATOR
    G4double sx = 5.9*cm;
    G4double sy = 20.6*cm;
    G4double sz = 1.1*cm;

    auto solidSurf = new G4Box("Surf", sx/2, sy/2, sz/2);
    auto logicSurf = new G4LogicalVolume(solidSurf, BC408, "Surf");

    new G4PVPlacement(nullptr, {0,30*cm,0}, logicSurf, "Surf", logicEnv, false, 0);

    // OPTICAL SURFACE (BETTER APPROACH)
    auto opSurface = new G4OpticalSurface("ScintSurface");
    opSurface->SetType(dielectric_dielectric);
    opSurface->SetModel(unified);
    opSurface->SetFinish(ground);

    new G4LogicalSkinSurface("BarSkin", logicBar, opSurface);
    new G4LogicalSkinSurface("SurfSkin", logicSurf, opSurface);

    // PMT MATERIAL
    auto glass = nist->FindOrBuildMaterial("G4_Pyrex_Glass");

    G4double r = 2.25*cm;
    G4double h = 12*cm;

    auto solidPMT = new G4Tubs("PMT", 0, r, h/2, 0, 360*deg);
    fLogicPMT = new G4LogicalVolume(solidPMT, glass, "PMT");

    G4double xEnd = barX/2 + h/2;

    new G4PVPlacement(nullptr, {-xEnd,0,73.5*cm}, fLogicPMT, "PMT1", logicEnv, false, 0);
    new G4PVPlacement(nullptr, {+xEnd,0,73.5*cm}, fLogicPMT, "PMT2", logicEnv, false, 1);
    new G4PVPlacement(nullptr, {0,30*cm+sy/2+h/2,0}, fLogicPMT, "PMT3", logicEnv, false, 2);

    return physWorld;

}

void DetectorConstruction::ConstructSDandField()
{
    auto sdMan = G4SDManager::GetSDMpointer();

    // scintillator SD
    auto scintSD = new SensitiveDetector("ScintSD");
    sdMan->AddNewDetector(scintSD);
    fScoringVolume->SetSensitiveDetector(scintSD);

    // PMT SD (CORRECT: eventAction passed safely)
    auto pmtSD = new PMTSD("PMTSD");
    sdMan->AddNewDetector(pmtSD);
    fLogicPMT->SetSensitiveDetector(pmtSD);
}