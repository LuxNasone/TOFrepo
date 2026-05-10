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

    //Manager for creating/selecting materials

    auto nist = G4NistManager::Instance();

    //Number of points for optical properties

    const G4int n = 2;

    //Scintillating photon energy range, to red to UV

    G4double e[n] = {2.0*eV, 3.5*eV};

    // World material : air

    auto worldMat = nist->FindOrBuildMaterial("G4_AIR");

    //Air refractive index, needed for photon transport to an interface to another
    //Also creation of optical properties table

    G4double rindexAir[n] = {1.0, 1.0};

    auto airMPT = new G4MaterialPropertiesTable();

    airMPT->AddProperty("RINDEX", e, rindexAir, n);

    worldMat->SetMaterialPropertiesTable(airMPT);

    //World creation

    auto solidWorld = new G4Box("World", 7*m/2, 7*m/2, 7*m/2);

    auto logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");

    auto physWorld  = new G4PVPlacement(nullptr, {}, logicWorld, "World", nullptr, false, 0);

    // Envelope creation

    auto solidEnv = new G4Box("Env", 6*m/2, 6*m/2, 6*m/2);

    auto logicEnv = new G4LogicalVolume(solidEnv, worldMat, "Env");

    new G4PVPlacement(nullptr, {}, logicEnv, "Env", logicWorld, false, 0);

    // Scintillator material : BC408 
    // Density

    auto BC408 = new G4Material("BC408", 1.032*g/cm3, 2);

    //Molecular composition

    BC408->AddElement(nist->FindOrBuildElement("C"), 9);

    BC408->AddElement(nist->FindOrBuildElement("H"), 10);

    //Refractive index

    G4double rindex[n] = {1.58, 1.58};

    //Absorption length

    G4double abslen[n] = {2.1*m, 2.1*m};

    //Optical properties table

    auto mpt = new G4MaterialPropertiesTable();

    mpt->AddProperty("RINDEX", e, rindex, n);

    mpt->AddProperty("ABSLENGTH", e, abslen, n);

    G4double scintSpectrum[n] = {1.0, 1.0};
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", e, scintSpectrum, n);

    //Scintillation yield

    mpt->AddConstProperty("SCINTILLATIONYIELD", 108./MeV);
    

    //Number of decay component

    mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);

    mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);

    //Time constant (decay)

    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1*ns);

    //Birks Costant and saving material stuff

    BC408->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

    BC408->SetMaterialPropertiesTable(mpt);

    // GLASS per PMT

    auto glass = nist->FindOrBuildMaterial("G4_Pyrex_Glass");

    G4double rindexGlass[n] = {1.47, 1.47};

    auto glassMPT = new G4MaterialPropertiesTable();

    glassMPT->AddProperty("RINDEX", e, rindexGlass, n);

    glass->SetMaterialPropertiesTable(glassMPT);

    // Bar, x dimension

    G4double barX = 2.795*m;

    //Bar, transverse dimension

    G4double barYZ = 4*cm;

    //Bar creation

    auto solidBar = new G4Box("Bar", barX/2, barYZ/2, barYZ/2);

    auto logicBar = new G4LogicalVolume(solidBar, BC408, "Bar");

    auto physBar = new G4PVPlacement(nullptr, {0, 0, 160 * cm}, logicBar, "Bar", logicEnv, false, 0);

    //Setting scoring volume

    fScoringVolume = logicBar;

    // Little scintillating slab dimensions

    G4double sx = 5.9 * cm, sy = 20.6 * cm, sz = 1.1 * cm;

    //Creating scintillating slab

    auto solidSurf = new G4Box("Surf", sx/2, sy/2, sz/2);

    auto logicSurf = new G4LogicalVolume(solidSurf, BC408, "Surf");

    auto physSurf = new G4PVPlacement(nullptr, {0, sy/2, 0}, logicSurf, "Surf", logicEnv, false, 0);

    //Cathode (PMT) geometry (cylindrical)

    G4double r = 2.25 * cm;

    G4double h = 0.5 * mm;

    //Rotation for PMT (one for bar readout other for slab)

    G4RotationMatrix* roty = new G4RotationMatrix();

    roty->rotateY(90.*deg);

    G4RotationMatrix* rotx = new G4RotationMatrix();

    rotx->rotateX(90.*deg);

    //Cathode creation

    auto solidCathode = new G4Tubs("Cathode", 0, r, 0.5*mm, 0, 360*deg);

    fLogicCathode = new G4LogicalVolume(solidCathode, glass, "Cathode");

    new G4PVPlacement(roty, {-barX/2 - h/2, 0, 1.6 * m}, fLogicCathode, "Cathode1", logicEnv, false, 0);

    new G4PVPlacement(roty, {barX/2 +  h/2, 0, 1.6 * m}, fLogicCathode, "Cathode2", logicEnv, false, 1);

    new G4PVPlacement(rotx, {0, sy + h/2, 0}, fLogicCathode, "Cathode3", logicEnv, false, 2);

    //Some optical properties (reflectivity and quantum eff)

    G4double reflectivity[n] = {0.0, 0.0};

    G4double qe[n] = {0.28, 0.28};

    auto cathodeMPT = new G4MaterialPropertiesTable();

    cathodeMPT->AddProperty("EFFICIENCY", e, qe, n);

    cathodeMPT->AddProperty("REFLECTIVITY", e, reflectivity, n);

    //Optical surface for cathode
    
    auto cathodeSurface = new G4OpticalSurface("CathodeSurface");

    cathodeSurface->SetType(dielectric_metal);

    cathodeSurface->SetModel(glisur);

    cathodeSurface->SetFinish(polished);

    cathodeSurface->SetMaterialPropertiesTable(cathodeMPT);

    new G4LogicalSkinSurface("CathodeSkin", fLogicCathode, cathodeSurface);

    // Scintillators optical surface

    G4double reflectivity_s[n] = {0.9, 0.9};

    auto opMPT = new G4MaterialPropertiesTable();

    opMPT->AddProperty("REFLECTIVITY", e, reflectivity_s, n);

    auto opSurface = new G4OpticalSurface("ScintSurface");

    opSurface->SetType(dielectric_metal);

    opSurface->SetModel(glisur);

    opSurface->SetFinish(polished);

    opSurface->SetMaterialPropertiesTable(opMPT);

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