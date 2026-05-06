#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SensitiveDetector.hh"
#include "PhysicsList.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc, char** argv)
{
    auto runManager = new G4RunManager();

    auto* physicsList = new PhysicsList();

    runManager->SetUserInitialization(physicsList);

    auto* eventAction = new EventAction();

    runManager->SetUserAction(eventAction);

    runManager->SetUserInitialization(new DetectorConstruction());

    runManager->SetUserAction(new PrimaryGeneratorAction());

    auto* runAction = new RunAction();

    runManager->SetUserAction(runAction);

    runManager->Initialize();

    G4VisManager* visManager = new G4VisExecutive();

    visManager->Initialize();

    G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    UImanager->ApplyCommand("/control/execute init_vis.mac");

    UImanager->ApplyCommand("/run/verbose 0");

    UImanager->ApplyCommand("/event/verbose 0");

    UImanager->ApplyCommand("/tracking/verbose 0");

    ui->SessionStart();

    delete ui;

    delete visManager;

    delete runManager;

    return 0;

}