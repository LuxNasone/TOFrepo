cmake_minimum_required(VERSION 3.16)
project(TOF_MC)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

include_directories(include)

file(GLOB sources src/*.cc)
file(GLOB headers include/*.hh)

add_executable(TOF_MC
    main.cc
    src/DetectorConstruction.cc
    src/PrimaryGeneratorAction.cc
    src/RunAction.cc
    src/EventAction.cc 
    src/SensitiveDetector.cc
    src/PhysicsList.cc
    src/PMTSD.cc
)

target_link_libraries(TOF_MC ${Geant4_LIBRARIES})