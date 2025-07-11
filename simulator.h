// simulator.h - Black Oil Reservoir Simulator Header
#ifndef BLACK_OIL_SIMULATOR_H
#define BLACK_OIL_SIMULATOR_H

#include <petsc.h>
#include <petscsnes.h>
#include <petscdmda.h>
#include <string>
#include <vector>
#include <map>
#include <memory>

// Forward declarations
namespace BlackOil {
    class DeckParser;
    class BlackOilSimulator;
}

// ==================== DATA STRUCTURES ====================

// Field structure for solution vector (p, Sw, Sg)
typedef struct {
    PetscScalar p;   // Pressure
    PetscScalar Sw;  // Water saturation
    PetscScalar Sg;  // Gas saturation
} Field;

// Grid structure
struct Grid3D {
    PetscInt Nx, Ny, Nz;          // Number of cells
    PetscReal dx, dy, dz;         // Cell dimensions (if uniform)
    std::vector<PetscReal> DX;     // Cell dimensions in X (per cell)
    std::vector<PetscReal> DY;     // Cell dimensions in Y (per cell)
    std::vector<PetscReal> DZ;     // Cell dimensions in Z (per cell)
    std::vector<PetscReal> TOPS;   // Top depths
    std::vector<PetscReal> porosity;
    std::vector<PetscReal> permX;
    std::vector<PetscReal> permY;
    std::vector<PetscReal> permZ;
    std::vector<PetscReal> swat;   // Initial water saturation
    std::vector<PetscReal> sgas;   // Initial gas saturation
};

// Well structure
struct Well {
    std::string name;
    PetscInt i, j, k1, k2;        // Cell location (k1 to k2 for multilayer)
    PetscReal refDepth;           // Reference depth for BHP
    PetscReal BHP;                // Bottom hole pressure
    PetscReal rate;               // Flow rate
    PetscInt type;                // 0: producer, 1: injector
    PetscInt control;             // 0: BHP control, 1: rate control
    std::string phase;            // OIL, WATER, GAS
    PetscReal WI;                 // Well index
    PetscReal rw;                 // Wellbore radius
    bool isOpen;                  // Well status
};

// PVT Table structures
struct PVTTable {
    std::vector<PetscReal> x;      // Independent variable
    std::vector<PetscReal> y;      // Dependent variable
};

struct PVTOTable {
    std::vector<PetscReal> Rs;     // Dissolved GOR
    std::vector<PetscReal> pbub;   // Bubble point pressure
    std::vector<PetscReal> Bo;     // Oil FVF
    std::vector<PetscReal> muo;    // Oil viscosity
    std::vector<std::vector<PetscReal>> p_usat;   // Undersaturated pressures
    std::vector<std::vector<PetscReal>> Bo_usat;  // Undersaturated Bo
    std::vector<std::vector<PetscReal>> muo_usat; // Undersaturated viscosity
};

// Fluid properties structure
struct FluidProps {
    // Water properties
    PetscReal pwref, Bwref, cw, muw, muw_p;
    
    // Oil properties
    PVTOTable pvto;
    
    // Gas properties
    std::vector<PetscReal> pvdg_p;    // Pressure
    std::vector<PetscReal> pvdg_Bg;   // Gas FVF
    std::vector<PetscReal> pvdg_mug;  // Gas viscosity
    
    // Relative permeability tables
    std::vector<PetscReal> swof_Sw, swof_krw, swof_krow, swof_pcow;
    std::vector<PetscReal> sgof_Sg, sgof_krg, sgof_krog, sgof_pcog;
    
    // Densities at surface conditions
    PetscReal rhoO_surf, rhoW_surf, rhoG_surf;
};

// Rock properties
struct RockProps {
    PetscReal crock;              // Rock compressibility
    PetscReal pref;               // Reference pressure
};

// Equilibration data
struct EquilData {
    PetscReal datum;              // Datum depth
    PetscReal pres_datum;         // Pressure at datum
    PetscReal zwoc;               // Water-oil contact
    PetscReal pcow_woc;           // Pc at WOC
    PetscReal zgoc;               // Gas-oil contact
    PetscReal pcog_goc;           // Pc at GOC
    std::vector<PetscReal> rsvd_depth;  // Depth for RSVD
    std::vector<PetscReal> rsvd_rs;     // Rs values for RSVD
};

// Application context
typedef struct {
    DM da;                         // Distributed array
    Grid3D grid;                   // Grid structure
    FluidProps fluid;              // Fluid properties
    RockProps rock;                // Rock properties
    std::vector<Well> wells;       // Well list
    EquilData equil;              // Equilibration data
    
    // Simulation parameters
    PetscReal dt;                  // Time step
    PetscReal t_current;           // Current time
    PetscReal t_end;               // End time
    PetscInt step;                 // Current step
    std::vector<PetscReal> timesteps; // Scheduled timesteps
    PetscInt timestep_idx;         // Current timestep index
    
    // Unit system
    bool useField;                 // true: field units, false: metric
    
    // PETSc objects
    Vec solution;                  // Current solution
    Vec solution_old;              // Previous timestep solution
    Mat J;                         // Jacobian matrix
    SNES snes;                     // Nonlinear solver
} AppCtx;

// ==================== FUNCTION DECLARATIONS ====================

// PETSc callback functions
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx);
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat J, Mat P, void *ctx);

// Utility functions
std::string toUpper(const std::string& str);
std::string trim(const std::string& str);
std::string trim(const std::string& str, const std::string& chars);
bool isComment(const std::string& line);
bool getDataLine(std::ifstream& file, std::string& line);
std::vector<PetscReal> parseNumbers(const std::string& line);
PetscReal interpolate(PetscReal x, const std::vector<PetscReal>& xData, 
                     const std::vector<PetscReal>& yData);

// PVT functions
void getOilProperties(const FluidProps& fluid, PetscReal p, PetscReal Rs,
                     PetscReal& Bo, PetscReal& muo);
void getWaterProperties(const FluidProps& fluid, PetscReal p,
                       PetscReal& Bw, PetscReal& muw);
void getGasProperties(const FluidProps& fluid, PetscReal p,
                     PetscReal& Bg, PetscReal& mug);
void getRelPerm(const FluidProps& fluid, PetscReal Sw, PetscReal Sg,
               PetscReal& krw, PetscReal& kro, PetscReal& krg);

// Initialization
void initializeEquilibration(AppCtx* ctx);

// Well model
PetscReal calculateWellIndex(const Well& well, const Grid3D& grid, 
                           PetscInt i, PetscInt j, PetscInt k);

// ==================== DECK PARSER CLASS ====================

namespace BlackOil {

class DeckParser {
private:
    AppCtx* ctx;
    std::ifstream* file;
    std::string currentSection;
    std::string deckFilename;
    
    // Unit conversion factors
    PetscReal ft_to_m = 0.3048;
    PetscReal psi_to_Pa = 6894.757;
    PetscReal rb_to_m3 = 0.1589873;
    PetscReal stb_to_m3 = 0.1589873;
    PetscReal mD_to_m2 = 9.869233e-16;
    PetscReal cp_to_Pas = 1e-3;
    PetscReal lbft3_to_kgm3 = 16.01846;
    PetscReal Mscf_to_sm3 = 28.31685;
    PetscReal day_to_sec = 86400.0;
    
    // Parsing methods
    bool parseRunspec();
    bool parseGrid();
    bool parseProps();
    bool parseSolution();
    bool parseSchedule();
    bool parseDimens();
    bool parsePVTW();
    bool parseROCK();
    bool parseDENSITY();
    bool parsePVTO();
    bool parsePVDG();
    bool parseSGOF();
    bool parseSWOF();
    bool parseEQUIL();
    bool parseRSVD();
    bool parseWELSPECS();
    bool parseCOMPDAT();
    bool parseWCONINJE();
    bool parseWCONPROD();
    bool parseTSTEP();
    
    std::string getKeyword(const std::string& line);
    std::vector<PetscReal> readArrayData(PetscInt expectedSize);
    bool matchesWildcard(const std::string& pattern, const std::string& str);
    
public:
    DeckParser(AppCtx* context);
    bool parseDeck(const std::string& filename);
};

// ==================== SIMULATOR CLASS ====================

class BlackOilSimulator {
private:
    AppCtx ctx;
    bool initialized;
    
    void setupPETScObjects();
    void outputResults();
    
public:
    BlackOilSimulator();
    ~BlackOilSimulator();
    
    bool loadDeck(const std::string& filename);
    bool loadEclipseDeck(const std::string& filename);  // New method for Eclipse format
    void initialize();
    void run();
    
    // Getters for Eclipse integration
    AppCtx& getContext() { return ctx; }
    const AppCtx& getContext() const { return ctx; }
};

} // namespace BlackOil

#endif // BLACK_OIL_SIMULATOR_H