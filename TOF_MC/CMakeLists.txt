cmake_minimum_required(VERSION 3.16)
project(TOF_MC)

# =========================
# Geant4 setup
# =========================
find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# =========================
# Sources
# =========================
file(GLOB sources
    ${PROJECT_SOURCE_DIR}/src/*.cc
)

file(GLOB headers
    ${PROJECT_SOURCE_DIR}/include/*.hh
)

# =========================
# Executable
# =========================
add_executable(TOF_MC main.cc ${sources} ${headers})

# =========================
# Linking
# =========================
target_link_libraries(TOF_MC ${Geant4_LIBRARIES})

# =========================
# ROOT / analysis (IMPORTANT)
# =========================
add_definitions(-DG4ANALYSIS_USE)

# =========================
# Optional: compile flags
# =========================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# =========================
# Copy macro files (optional)
# =========================
file(GLOB MACRO_FILES
    "${PROJECT_SOURCE_DIR}/*.mac"
)

foreach(_file ${MACRO_FILES})
    configure_file(
        ${_file}
        ${PROJECT_BINARY_DIR}
        COPYONLY
    )
endforeach()