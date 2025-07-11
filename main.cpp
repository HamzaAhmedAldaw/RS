// main.cpp - Black Oil Reservoir Simulator Main Program
#include "simulator.h"
#include "eclipse_input_reader.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

// Function to detect deck format
enum class DeckFormat {
    UNKNOWN,
    ECLIPSE,
    SIMPLE
};

DeckFormat detectDeckFormat(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return DeckFormat::UNKNOWN;
    }
    
    // Read first few lines to detect format
    std::string line;
    int lineCount = 0;
    bool hasEclipseKeywords = false;
    
    while (std::getline(file, line) && lineCount < 50) {
        lineCount++;
        
        // Convert to uppercase for comparison
        std::transform(line.begin(), line.end(), line.begin(), ::toupper);
        
        // Check for Eclipse section keywords
        if (line.find("RUNSPEC") != std::string::npos ||
            line.find("GRID") != std::string::npos ||
            line.find("PROPS") != std::string::npos ||
            line.find("SOLUTION") != std::string::npos ||
            line.find("SCHEDULE") != std::string::npos) {
            hasEclipseKeywords = true;
            break;
        }
    }
    
    file.close();
    
    return hasEclipseKeywords ? DeckFormat::ECLIPSE : DeckFormat::SIMPLE;
}

// Function to create a simple test deck
void createSimpleTestDeck(const std::string& filename) {
    std::ofstream test(filename);
    test << "-- Simple test deck for Black Oil Simulator\n";
    test << "RUNSPEC\n";
    test << "DIMENS\n";
    test << "  3 3 2 /\n";
    test << "FIELD\n";
    test << "OIL\n";
    test << "WATER\n";
    test << "GAS\n";
    test << "\n";
    test << "GRID\n";
    test << "DX\n";
    test << "  18*100 /\n";
    test << "DY\n";
    test << "  18*100 /\n";
    test << "DZ\n";
    test << "  18*20 /\n";
    test << "PORO\n";
    test << "  18*0.25 /\n";
    test << "PERMX\n";
    test << "  18*200 /\n";
    test << "\n";
    test << "PROPS\n";
    test << "PVTW\n";
    test << "  3600 1.0 3e-6 0.5 0 /\n";
    test << "ROCK\n";
    test << "  3600 4e-6 /\n";
    test << "DENSITY\n";
    test << "  45 63 0.07 /\n";
    test << "\n";
    test << "SOLUTION\n";
    test << "EQUIL\n";
    test << "  8000 3600 9000 0 7000 0 /\n";
    test << "\n";
    test << "SCHEDULE\n";
    test << "WELSPECS\n";
    test << "  'P1' 'G' 2 2 8000 'OIL' /\n";
    test << "  'I1' 'G' 3 3 8000 'WATER' /\n";
    test << "/\n";
    test << "COMPDAT\n";
    test << "  'P1' 2 2 1 2 'OPEN' 1* 1* 0.5 /\n";
    test << "  'I1' 3 3 1 2 'OPEN' 1* 1* 0.5 /\n";
    test << "/\n";
    test << "WCONPROD\n";
    test << "  'P1' 'OPEN' 'BHP' 1* 1* 1* 1* 1* 1000 /\n";
    test << "/\n";
    test << "WCONINJE\n";
    test << "  'I1' 'WATER' 'OPEN' 'RATE' 100 1* 4000 /\n";
    test << "/\n";
    test << "TSTEP\n";
    test << "  10*1 /\n";
    test << "END\n";
    test.close();
}

// Function to print usage
void printUsage(const char* programName) {
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Usage: %s <deck_file> [options]\n", programName);
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Options:\n");
    PetscPrintf(PETSC_COMM_WORLD, "  --eclipse     Force Eclipse format parser\n");
    PetscPrintf(PETSC_COMM_WORLD, "  --simple      Force simple format parser\n");
    PetscPrintf(PETSC_COMM_WORLD, "  --help        Show this help message\n");
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Examples:\n");
    PetscPrintf(PETSC_COMM_WORLD, "  %s SPE1.DATA\n", programName);
    PetscPrintf(PETSC_COMM_WORLD, "  %s SPE9.DATA --eclipse\n", programName);
    PetscPrintf(PETSC_COMM_WORLD, "  %s test_simple.DATA --simple\n", programName);
    PetscPrintf(PETSC_COMM_WORLD, "\n");
}

int main(int argc, char **argv) {
    // Initialize PETSc
    PetscInitialize(&argc, &argv, NULL, NULL);
    
    // Parse command line arguments
    std::string deckFile = "";
    bool forceEclipse = false;
    bool forceSimple = false;
    bool showHelp = false;
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--eclipse") {
            forceEclipse = true;
        } else if (arg == "--simple") {
            forceSimple = true;
        } else if (arg == "--help" || arg == "-h") {
            showHelp = true;
        } else if (arg[0] != '-' && deckFile.empty()) {
            deckFile = arg;
        }
    }
    
    // Show help if requested or no arguments
    if (showHelp || argc < 2 || deckFile.empty()) {
        printUsage(argv[0]);
        
        if (argc < 2) {
            // Create a simple test deck for demonstration
            PetscPrintf(PETSC_COMM_WORLD, "No deck file specified. Creating test_simple.DATA for demonstration...\n");
            createSimpleTestDeck("test_simple.DATA");
            PetscPrintf(PETSC_COMM_WORLD, "Created test_simple.DATA\n");
            PetscPrintf(PETSC_COMM_WORLD, "Run with: %s test_simple.DATA\n", argv[0]);
        }
        
        PetscFinalize();
        return showHelp ? 0 : 1;
    }
    
    // Print header
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "================================================\n");
    PetscPrintf(PETSC_COMM_WORLD, "   Black Oil Reservoir Simulator                \n");
    PetscPrintf(PETSC_COMM_WORLD, "   SPE1/SPE9 Compatible Version                 \n");
    PetscPrintf(PETSC_COMM_WORLD, "================================================\n");
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    
    // Detect deck format
    DeckFormat format = DeckFormat::UNKNOWN;
    if (forceEclipse) {
        format = DeckFormat::ECLIPSE;
        PetscPrintf(PETSC_COMM_WORLD, "Forcing Eclipse format parser\n");
    } else if (forceSimple) {
        format = DeckFormat::SIMPLE;
        PetscPrintf(PETSC_COMM_WORLD, "Forcing simple format parser\n");
    } else {
        format = detectDeckFormat(deckFile);
        PetscPrintf(PETSC_COMM_WORLD, "Auto-detected deck format: %s\n", 
                   format == DeckFormat::ECLIPSE ? "Eclipse" : "Simple");
    }
    
    // Create simulator
    BlackOil::BlackOilSimulator simulator;
    
    // Load deck
    PetscPrintf(PETSC_COMM_WORLD, "Loading deck file: %s\n", deckFile.c_str());
    
    bool loadSuccess = false;
    if (format == DeckFormat::ECLIPSE) {
        // Use Eclipse parser
        loadSuccess = simulator.loadEclipseDeck(deckFile);
    } else {
        // Use simple built-in parser
        loadSuccess = simulator.loadDeck(deckFile);
    }
    
    if (!loadSuccess) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Failed to load deck file\n");
        
        // If Eclipse format failed, try simple format as fallback
        if (format == DeckFormat::ECLIPSE && !forceEclipse) {
            PetscPrintf(PETSC_COMM_WORLD, "Trying simple format parser as fallback...\n");
            loadSuccess = simulator.loadDeck(deckFile);
        }
        
        if (!loadSuccess) {
            PetscFinalize();
            return 1;
        }
    }
    
    // Initialize simulator
    PetscPrintf(PETSC_COMM_WORLD, "\nInitializing simulator...\n");
    simulator.initialize();
    
    // Run simulation
    PetscPrintf(PETSC_COMM_WORLD, "Starting simulation...\n");
    simulator.run();
    
    // Finalize PETSc
    PetscFinalize();
    
    PetscPrintf(PETSC_COMM_WORLD, "\nSimulation completed successfully.\n");
    return 0;
}