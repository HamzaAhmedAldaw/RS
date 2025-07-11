#ifndef ECLIPSE_INPUT_READER_H
#define ECLIPSE_INPUT_READER_H

#include "simulator.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

// ============================================================================
// Data Structures
// ============================================================================

struct RunspecData {
    std::string title;
    int nx = 0, ny = 0, nz = 0;
    bool oil_phase = false;
    bool water_phase = false;
    bool gas_phase = false;
    bool dissolved_gas = false;
    bool vaporized_oil = false;
    std::string unit_system = "METRIC";
    std::string start_date;
};

struct GridData {
    std::vector<double> dx, dy, dz;
    std::vector<double> tops;
    std::vector<double> poro;
    std::vector<double> permx, permy, permz;
    std::vector<double> ntg;
    std::vector<int> actnum;
};

struct PropsData {
    struct PVTData {
        std::vector<double> pressure;
        std::vector<double> fvf;
        std::vector<double> viscosity;
        std::vector<double> rs;
        double compressibility = 0.0;
        double viscosibility = 0.0;
    };
    
    struct RelPermData {
        std::vector<double> sw, sg, so;
        std::vector<double> krw, krg, krow, krog;
        std::vector<double> pcow, pcog;
    };
    
    struct RockData {
        double pref;
        double compressibility;
    };
    
    std::vector<PVTData> pvto;
    std::vector<PVTData> pvtw;
    std::vector<PVTData> pvtg;
    std::vector<RelPermData> swof;
    std::vector<RelPermData> sgof;
    std::vector<RockData> rock;
    double oil_density = 0.0;
    double water_density = 0.0;
    double gas_density = 0.0;
};

struct SolutionData {
    struct EquilData {
        double datum_depth = 0.0;
        double datum_pressure = 0.0;
        double owc_depth = 0.0;
        double pc_owc = 0.0;
        double goc_depth = 0.0;
        double pc_goc = 0.0;
        int rsvd_table = 0;
        int rvvd_table = 0;
        int n = 0;
    };
    
    std::vector<double> pressure;
    std::vector<double> swat;
    std::vector<double> sgas;
    std::vector<double> rs;
    std::vector<EquilData> equil;
};

struct WellData {
    struct Completion {
        int i = 0, j = 0, k1 = 0, k2 = 0;
        std::string status = "OPEN";
        int sat_table = 0;
        double conn_factor = 0.0;
        double wellbore_diam = 0.5;
        double kh = 0.0;
        double skin = 0.0;
        double dfact = 0.0;
    };
    
    std::string name;
    std::string group;
    int i = 0, j = 0, k1 = 0, k2 = 0;
    double ref_depth = 0.0;
    std::string phase;
    double radius = 0.5;
    
    // Control data
    std::string control_mode;
    std::string status = "OPEN";
    double oil_rate = 0.0;
    double water_rate = 0.0;
    double gas_rate = 0.0;
    double liquid_rate = 0.0;
    double resv_rate = 0.0;
    double bhp_limit = 0.0;
    double thp_limit = 0.0;
    
    // Well type
    bool is_producer = false;
    std::string injection_type = "WATER";
    
    // Completions
    std::vector<Completion> completions;
};

struct ScheduleData {
    struct TimeStep {
        double days = 0.0;
        std::vector<WellData> wells;
    };
    
    std::vector<TimeStep> timesteps;
};

// ============================================================================
// Eclipse Input Reader Class
// ============================================================================

class EclipseInputReader {
private:
    // File handling
    std::string filename;
    std::ifstream file;
    int line_number = 0;
    
    // Data storage
    RunspecData runspec;
    GridData grid_data;
    PropsData props_data;
    SolutionData solution_data;
    ScheduleData schedule_data;
    
    // Current section
    enum class Section {
        NONE, RUNSPEC, GRID, EDIT, PROPS, REGIONS, SOLUTION, SUMMARY, SCHEDULE
    };
    Section current_section = Section::NONE;
    
    // Helper methods
    bool getline_with_continuation(std::string& line);
    bool is_section_keyword(const std::string& line);
    void handle_section_keyword(const std::string& line);
    void process_keyword(const std::string& line);
    
    // Section processors
    void process_runspec_keyword(const std::string& keyword);
    void process_grid_keyword(const std::string& keyword);
    void process_props_keyword(const std::string& keyword);
    void process_solution_keyword(const std::string& keyword);
    void process_schedule_keyword(const std::string& keyword);
    
    // Data readers
    std::vector<std::string> read_tokens_until_slash();
    std::vector<double> read_numeric_array(int expected_count);
    void process_copy_keyword();
    
    // Table readers
    void read_pvto_table();
    void read_pvtw_table();
    void read_pvtg_table();
    void read_swof_table();
    void read_sgof_table();
    void read_equil_data();
    
    // Schedule readers
    void read_welspecs();
    void read_compdat();
    void read_wconprod();
    void read_wconinje();
    void read_tstep();
    
    // Unit conversion
    void apply_unit_conversions();
    
    // Validation
    bool validate_data();

    // INCLUDE file support
    void process_include(const std::string& filename);
    
public:
    explicit EclipseInputReader(const std::string& fname);
    
    bool parse();
    
    // Getters
    const RunspecData& get_runspec() const;
    const GridData& get_grid_data() const;
    const PropsData& get_props_data() const;
    const SolutionData& get_solution_data() const;
    const ScheduleData& get_schedule_data() const;
};

// Global parse function
bool parse_eclipse_deck(const std::string& filename, 
                       RunspecData& runspec, 
                       GridData& grid,
                       PropsData& props, 
                       SolutionData& solution, 
                       ScheduleData& schedule);

#endif // ECLIPSE_INPUT_READER_H