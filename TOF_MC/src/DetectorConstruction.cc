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
#include "G4RotationMatrix.hh"

DetectorConstruction::DetectorConstruction(EventAction* eventAction)
: fEventAction(eventAction) {}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct(){

    auto nist = G4NistManager::Instance();

    const G4int n = 2;
    G4double e[n] = {2.0*eV, 3.5*eV};

    // WORLD
    auto worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4double rindexAir[n] = {1.0, 1.0};
    auto airMPT = new G4MaterialPropertiesTable();
    airMPT->AddProperty("RINDEX", e, rindexAir, n);
    worldMat->SetMaterialPropertiesTable(airMPT);

    auto solidWorld = new G4Box("World", 5*m/2, 5*m/2, 5*m/2);
    auto logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
    auto physWorld  = new G4PVPlacement(nullptr, {}, logicWorld, "World", nullptr, false, 0);

    // ENVELOPE
    auto solidEnv = new G4Box("Env", 4*m/2, 4*m/2, 4*m/2);
    auto logicEnv = new G4LogicalVolume(solidEnv, worldMat, "Env");
    new G4PVPlacement(nullptr, {}, logicEnv, "Env", logicWorld, false, 0);

    // SCINTILLATOR MATERIAL BC408
    auto BC408 = new G4Material("BC408", 1.032*g/cm3, 2);
    BC408->AddElement(nist->FindOrBuildElement("C"), 9);
    BC408->AddElement(nist->FindOrBuildElement("H"), 10);

    G4double rindex[n] = {1.58, 1.58};
    G4double abslen[n] = {2.1*m, 2.1*m};
    auto mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", e, rindex, n);
    mpt->AddProperty("ABSLENGTH", e, abslen, n);
    mpt->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1*ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    BC408->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
    BC408->SetMaterialPropertiesTable(mpt);

    // GLASS per PMT
    auto glass = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
    G4double rindexGlass[n] = {1.47, 1.47};
    auto glassMPT = new G4MaterialPropertiesTable();
    glassMPT->AddProperty("RINDEX", e, rindexGlass, n);
    glass->SetMaterialPropertiesTable(glassMPT);

    // BAR (lungo X, centrata in z=73.5cm)
    G4double barX = 2.795*m;
    G4double barYZ = 4*cm;
    auto solidBar = new G4Box("Bar", barX/2, barYZ/2, barYZ/2);
    auto logicBar = new G4LogicalVolume(solidBar, BC408, "Bar");
    auto physBar = new G4PVPlacement(nullptr, {0,0,160*cm}, logicBar, "Bar", logicEnv, false, 0);
    fScoringVolume = logicBar;

    // PICCOLO SCINTILLATORE (y=30cm, z=0)
    G4double sx = 5.9*cm, sy = 20.6*cm, sz = 1.1*cm;
    auto solidSurf = new G4Box("Surf", sx/2, sy/2, sz/2);
    auto logicSurf = new G4LogicalVolume(solidSurf, BC408, "Surf");
    auto physSurf = new G4PVPlacement(nullptr, {0,24*cm + sy/2,0}, logicSurf, "Surf", logicEnv, false, 0);

    // PMT geometry
    G4double r = 2.25*cm;
    G4double h = 12*cm;
    auto solidPMT = new G4Tubs("PMT", 0, r, h/2, 0, 360*deg);
    fLogicPMT = new G4LogicalVolume(solidPMT, glass, "PMT");

    // Catodo (faccia frontale del PMT, in coordinate locali del PMT)
    // Il PMT ha asse Z locale. La faccia frontale è a -h/2 in Z locale.

    // Rotazioni
    G4RotationMatrix* roty = new G4RotationMatrix();
    roty->rotateY(90.*deg);
    G4RotationMatrix* rotx = new G4RotationMatrix();
    rotx->rotateX(90.*deg);

    auto solidCathode = new G4Tubs("Cathode", 0, r, 0.5*mm, 0, 360*deg);

    // Catodo per PMT1 e PMT2 (asse lungo X, faccia verso -Z locale)
    fLogicCathode = new G4LogicalVolume(solidCathode, glass, "Cathode");
    new G4PVPlacement(nullptr, {0, 0, -(h/2 - 0.5*mm)}, fLogicCathode, "Cathode", fLogicPMT, false, 0);

    G4double xEnd = barX/2 + h/2;
    G4double ySmall = 24*cm + sy/2;
    G4double yEnd = ySmall + sy/2 + h/2; 

    new G4PVPlacement(roty, {-xEnd, 0., 160*cm}, fLogicPMT, "PMT1", logicEnv, false, 0);
    new G4PVPlacement(roty, {+xEnd, 0., 160*cm}, fLogicPMT, "PMT2", logicEnv, false, 1);
    new G4PVPlacement(rotx, {0, yEnd, 0}, fLogicPMT, "PMT3", logicEnv, false, 2);

    // SUPERFICIE OTTICA scintillatori (riflettente sui lati)
    auto opSurface = new G4OpticalSurface("ScintSurface");
    opSurface->SetType(dielectric_dielectric);
    opSurface->SetModel(unified);
    opSurface->SetFinish(polished);
    new G4LogicalSkinSurface("BarSkin", logicBar, opSurface);
    new G4LogicalSkinSurface("SurfSkin", logicSurf, opSurface);

    return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
    auto sdMan = G4SDManager::GetSDMpointer();

    // SD scintillatore
    auto scintSD = new SensitiveDetector("ScintSD");
    sdMan->AddNewDetector(scintSD);
    fScoringVolume->SetSensitiveDetector(scintSD);

    // SD catodo - un solo logical, copyNo del PMT parent distingue quale
    auto pmtSD = new PMTSD("PMTSD", fEventAction);
    sdMan->AddNewDetector(pmtSD);
    fLogicCathode->SetSensitiveDetector(pmtSD);
}