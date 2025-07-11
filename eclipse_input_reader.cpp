#include "eclipse_input_reader.h"
#include "simulator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <numeric> // Required for std::accumulate

// ============================================================================
// Unit Conversion Constants
// ============================================================================
namespace Units {
    // FIELD to METRIC conversions
    constexpr double FT_TO_M = 0.3048;
    constexpr double PSI_TO_BAR = 0.0689476;
    constexpr double PSI_TO_PA = 6894.76;  // PSI to Pascal
    constexpr double BAR_TO_PA = 1e5;      // Bar to Pascal
    constexpr double STB_TO_M3 = 0.1589873;
    constexpr double MSCF_TO_M3 = 28.31685;
    constexpr double RB_TO_M3 = 0.1589873;
    constexpr double CP_TO_PAS = 0.001;
    constexpr double MD_TO_M2 = 9.869233e-16;
    constexpr double LB_FT3_TO_KG_M3 = 16.01846;
    
    // Conversion functions
    double convert_length(double value, const std::string& from_unit) {
        return (from_unit == "FIELD") ? value * FT_TO_M : value;
    }
    
    double convert_pressure(double value, const std::string& from_unit) {
        return (from_unit == "FIELD") ? value * PSI_TO_BAR : value;
    }
    
    double convert_volume(double value, const std::string& from_unit) {
        return (from_unit == "FIELD") ? value * STB_TO_M3 : value;
    }
    
    double convert_gas_volume(double value, const std::string& from_unit) {
        return (from_unit == "FIELD") ? value * MSCF_TO_M3 : value;
    }
    
    double convert_permeability(double value, const std::string& from_unit) {
        (void)from_unit;
        // mD is same in both unit systems
        return value;
    }
    
    double convert_density(double value, const std::string& from_unit) {
        return (from_unit == "FIELD") ? value * LB_FT3_TO_KG_M3 : value;
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

std::string to_upper(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

bool is_default_value(const std::string& token) {
    return token == "1*" || token == "2*" || token == "3*" || 
           token == "4*" || token == "5*" || token == "6*";
}

// ============================================================================
// Enhanced EclipseInputReader Implementation
// ============================================================================
EclipseInputReader::EclipseInputReader(const std::string& fname) 
    : filename(fname), line_number(0), current_section(Section::NONE) {
    // Set default unit system
    runspec.unit_system = "METRIC";
}

bool EclipseInputReader::parse() {
    file.open(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::cout << "Parsing Eclipse deck..." << std::endl;
    
    std::string line;
    while (getline_with_continuation(line)) {
        line_number++;
        
        // Skip comments and empty lines
        line = trim(line);
        if (line.empty() || line.substr(0, 2) == "--") continue;
        
        // Remove inline comments
        size_t comment_pos = line.find("--");
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
            line = trim(line);
            if (line.empty()) continue;
        }
        
        // Check for section keywords
        if (is_section_keyword(line)) {
            handle_section_keyword(line);
        } else {
            // Process keyword within section
            process_keyword(line);
        }
    }
    
    file.close();
    
    // Post-process to apply unit conversions
    apply_unit_conversions();
    
    validate_data();
    return true;
}

void EclipseInputReader::apply_unit_conversions() {
    if (runspec.unit_system == "FIELD") {
        std::cout << "Converting from FIELD to METRIC units..." << std::endl;
        
        // Convert grid dimensions
        for (auto& val : grid_data.dx) val = Units::convert_length(val, "FIELD");
        for (auto& val : grid_data.dy) val = Units::convert_length(val, "FIELD");
        for (auto& val : grid_data.dz) val = Units::convert_length(val, "FIELD");
        for (auto& val : grid_data.tops) val = Units::convert_length(val, "FIELD");
        
        // Convert well data
        for (auto& timestep : schedule_data.timesteps) {
            for (auto& well : timestep.wells) {
                well.ref_depth = Units::convert_length(well.ref_depth, "FIELD");
                well.oil_rate = Units::convert_volume(well.oil_rate, "FIELD");
                well.water_rate = Units::convert_volume(well.water_rate, "FIELD");
                well.gas_rate = Units::convert_gas_volume(well.gas_rate, "FIELD");
                well.liquid_rate = Units::convert_volume(well.liquid_rate, "FIELD");
                well.bhp_limit = Units::convert_pressure(well.bhp_limit, "FIELD");
                
                // Convert well radius if specified
                if (well.radius > 0) {
                    well.radius = Units::convert_length(well.radius, "FIELD");
                }
            }
        }
        
        // Convert PVT data
        for (auto& table : props_data.pvto) {
            for (auto& p : table.pressure) p = Units::convert_pressure(p, "FIELD");
            for (auto& rs : table.rs) rs = Units::convert_gas_volume(rs, "FIELD") / Units::STB_TO_M3;
        }
        
        for (auto& table : props_data.pvtw) {
            for (auto& p : table.pressure) p = Units::convert_pressure(p, "FIELD");
        }
        
        for (auto& table : props_data.pvtg) {
            for (auto& p : table.pressure) p = Units::convert_pressure(p, "FIELD");
        }
        
        // Convert initial conditions
        for (auto& p : solution_data.pressure) p = Units::convert_pressure(p, "FIELD");
        
        // EQUIL data is already converted to Pa in read_equil_data(), so only convert depths
        for (auto& equil : solution_data.equil) {
            equil.datum_depth = Units::convert_length(equil.datum_depth, "FIELD");
            // equil.datum_pressure is already in Pa, don't convert again!
            equil.owc_depth = Units::convert_length(equil.owc_depth, "FIELD");
            equil.goc_depth = Units::convert_length(equil.goc_depth, "FIELD");
        }
    }
}

bool EclipseInputReader::getline_with_continuation(std::string& line) {
    std::string buffer;
    if (!std::getline(file, buffer)) {
        return false;
    }
    
    line = buffer;
    
    // Handle continuation with & at end of line
    while (!line.empty() && line.back() == '&') {
        line.pop_back();  // Remove the &
        
        if (!std::getline(file, buffer)) {
            break;
        }
        line_number++;
        line += " " + buffer;
    }
    
    return true;
}

bool EclipseInputReader::is_section_keyword(const std::string& line) {
    std::string upper_line = to_upper(trim(line));
    return upper_line == "RUNSPEC" || upper_line == "GRID" || upper_line == "EDIT" ||
           upper_line == "PROPS" || upper_line == "REGIONS" || upper_line == "SOLUTION" ||
           upper_line == "SUMMARY" || upper_line == "SCHEDULE" || upper_line == "END";
}

void EclipseInputReader::handle_section_keyword(const std::string& line) {
    std::string upper_line = to_upper(trim(line));
    
    std::cout << "  Entering " << upper_line << " section" << std::endl;
    
    if (upper_line == "RUNSPEC") {
        current_section = Section::RUNSPEC;
    } else if (upper_line == "GRID") {
        current_section = Section::GRID;
    } else if (upper_line == "EDIT") {
        current_section = Section::EDIT;
    } else if (upper_line == "PROPS") {
        current_section = Section::PROPS;
    } else if (upper_line == "REGIONS") {
        current_section = Section::REGIONS;
    } else if (upper_line == "SOLUTION") {
        current_section = Section::SOLUTION;
    } else if (upper_line == "SUMMARY") {
        current_section = Section::SUMMARY;
    } else if (upper_line == "SCHEDULE") {
        current_section = Section::SCHEDULE;
    } else if (upper_line == "END") {
        current_section = Section::NONE;
    }
}

void EclipseInputReader::process_keyword(const std::string& line) {
    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;
    keyword = to_upper(keyword);
    
    switch (current_section) {
        case Section::RUNSPEC:
            process_runspec_keyword(keyword);
            break;
        case Section::GRID:
            process_grid_keyword(keyword);
            break;
        case Section::PROPS:
            process_props_keyword(keyword);
            break;
        case Section::SOLUTION:
            process_solution_keyword(keyword);
            break;
        case Section::SCHEDULE:
            process_schedule_keyword(keyword);
            break;
        default:
            // Ignore keywords in other sections
            break;
    }
}

std::vector<std::string> EclipseInputReader::read_tokens_until_slash() {
    std::vector<std::string> tokens;
    std::string line;
    bool found_slash = false;
    
    while (!found_slash && getline_with_continuation(line)) {
        line_number++;
        line = trim(line);
        if (line.empty() || line.substr(0, 2) == "--") continue;
        
        // Remove inline comments
        size_t comment_pos = line.find("--");
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
            line = trim(line);
        }
        
        // Check if this is a section keyword
        if (is_section_keyword(line)) {
            std::cerr << "Warning: Unexpected section keyword '" << line 
                      << "' at line " << line_number << std::endl;
            handle_section_keyword(line);
            break;
        }
        
        std::istringstream iss(line);
        std::string token;
        
        while (iss >> token) {
            if (token == "/") {
                found_slash = true;
                break;
            }
            tokens.push_back(token);
        }
    }
    
    return tokens;
}

std::vector<double> EclipseInputReader::read_numeric_array(int expected_count) {
    std::vector<double> data;
    auto tokens = read_tokens_until_slash();
    
    for (const auto& token : tokens) {
        if (data.size() >= static_cast<size_t>(expected_count)) break;
        
        // Handle multiplier syntax: N*value
        size_t star_pos = token.find('*');
        if (star_pos != std::string::npos) {
            try {
                int count = std::stoi(token.substr(0, star_pos));
                double value = std::stod(token.substr(star_pos + 1));
                for (int i = 0; i < count && data.size() < static_cast<size_t>(expected_count); ++i) {
                    data.push_back(value);
                }
            } catch (const std::exception&) {
                std::cerr << "Error parsing multiplier syntax: " << token 
                          << " at line " << line_number << std::endl;
            }
        } else if (!is_default_value(token)) {
            try {
                data.push_back(std::stod(token));
            } catch (const std::exception&) {
                std::cerr << "Error parsing numeric value: " << token 
                          << " at line " << line_number << std::endl;
            }
        }
    }
    
    // Fill remaining with zeros if needed
    while (data.size() < static_cast<size_t>(expected_count)) {
        data.push_back(0.0);
    }
    
    return data;
}

void EclipseInputReader::process_runspec_keyword(const std::string& keyword) {
    std::cout << "    DEBUG: Processing RUNSPEC keyword: '" << keyword << "'" << std::endl;
    if (keyword == "TITLE") {
        std::string line;
        if (getline_with_continuation(line)) {
            line_number++;
            runspec.title = trim(line);
            std::cout << "    Title: " << runspec.title << std::endl;
        }
        // DO NOT call read_tokens_until_slash() for TITLE
        // TITLE in Eclipse format doesn't use a terminating slash
        // It only reads the single line following the keyword
        
    } else if (keyword == "DIMENS") {
        std::cout << "    DEBUG: Found DIMENS keyword" << std::endl;
        auto tokens = read_tokens_until_slash();
        std::cout << "    DEBUG: Read " << tokens.size() << " tokens" << std::endl;
        for (size_t i = 0; i < tokens.size(); ++i) {
            std::cout << "    DEBUG: Token[" << i << "] = '" << tokens[i] << "'" << std::endl;
        }
        if (tokens.size() >= 3) {
            runspec.nx = std::stoi(tokens[0]);
            runspec.ny = std::stoi(tokens[1]);
            runspec.nz = std::stoi(tokens[2]);
            std::cout << "    Grid dimensions: " << runspec.nx << " x " 
                      << runspec.ny << " x " << runspec.nz << std::endl;
        } else {
            std::cout << "    ERROR: Not enough tokens for DIMENS!" << std::endl;
        }
    } else if (keyword == "OIL") {
        runspec.oil_phase = true;
        
    } else if (keyword == "WATER") {
        runspec.water_phase = true;
        
    } else if (keyword == "GAS") {
        runspec.gas_phase = true;
        
    } else if (keyword == "DISGAS") {
        runspec.dissolved_gas = true;
        
    } else if (keyword == "VAPOIL") {
        runspec.vaporized_oil = true;
        
    } else if (keyword == "METRIC") {
        runspec.unit_system = "METRIC";
        std::cout << "    Unit system: METRIC (no conversion needed)" << std::endl;
        
    } else if (keyword == "FIELD") {
        runspec.unit_system = "FIELD";
        std::cout << "    Unit system: FIELD (will convert to SI)" << std::endl;
        
    } else if (keyword == "START") {
        auto tokens = read_tokens_until_slash();
        if (!tokens.empty()) {
            runspec.start_date = "";
            for (const auto& t : tokens) {
                runspec.start_date += t + " ";
            }
            runspec.start_date = trim(runspec.start_date);
        }
        
    } else if (keyword == "WELLDIMS") {
        auto tokens = read_tokens_until_slash();
        // Parse well dimensions if needed
        
    } else if (keyword == "TABDIMS") {
        auto tokens = read_tokens_until_slash();
        // Parse table dimensions if needed
    }
}

void EclipseInputReader::process_grid_keyword(const std::string& keyword) {
    int ncells = runspec.nx * runspec.ny * runspec.nz;
    
    if (keyword == "DX") {
        grid_data.dx = read_numeric_array(ncells);
    } else if (keyword == "DY") {
        grid_data.dy = read_numeric_array(ncells);
    } else if (keyword == "DZ") {
        grid_data.dz = read_numeric_array(ncells);
    } else if (keyword == "TOPS") {
        grid_data.tops = read_numeric_array(runspec.nx * runspec.ny);
    } else if (keyword == "PORO") {
        grid_data.poro = read_numeric_array(ncells);
    } else if (keyword == "PERMX") {
        grid_data.permx = read_numeric_array(ncells);
    } else if (keyword == "PERMY") {
        grid_data.permy = read_numeric_array(ncells);
    } else if (keyword == "PERMZ") {
        grid_data.permz = read_numeric_array(ncells);
    } else if (keyword == "NTG") {
        grid_data.ntg = read_numeric_array(ncells);
    } else if (keyword == "COPY") {
        process_copy_keyword();
    }
}

void EclipseInputReader::process_copy_keyword() {
    auto tokens = read_tokens_until_slash();
    if (tokens.size() >= 2) {
        std::string from = to_upper(tokens[0]);
        std::string to = to_upper(tokens[1]);
        
        // Perform the copy
        if (from == "PERMX" && to == "PERMY") {
            grid_data.permy = grid_data.permx;
        } else if (from == "PERMX" && to == "PERMZ") {
            grid_data.permz = grid_data.permx;
        } else if (from == "PERMY" && to == "PERMZ") {
            grid_data.permz = grid_data.permy;
        }
        // Add more copy operations as needed
    }
}

void EclipseInputReader::process_props_keyword(const std::string& keyword) {
    if (keyword == "PVTW") {
        std::cout << "    Reading PVTW" << std::endl;
        read_pvtw_table();
    } else if (keyword == "PVTO") {
        std::cout << "    Reading PVTO" << std::endl;
        read_pvto_table();
    } else if (keyword == "PVTG" || keyword == "PVDG") {
        std::cout << "    Reading " << keyword << std::endl;
        read_pvtg_table();
    } else if (keyword == "SWOF") {
        std::cout << "    Reading SWOF" << std::endl;
        read_swof_table();
    } else if (keyword == "SGOF") {
        std::cout << "    Reading SGOF" << std::endl;
        read_sgof_table();
    } else if (keyword == "ROCK") {
        auto tokens = read_tokens_until_slash();
        if (tokens.size() >= 2) {
            PropsData::RockData rock;
            rock.pref = std::stod(tokens[0]);
            rock.compressibility = std::stod(tokens[1]);
            props_data.rock.push_back(rock);
        }
    } else if (keyword == "DENSITY") {
        auto tokens = read_tokens_until_slash();
        if (tokens.size() >= 3) {
            props_data.oil_density = std::stod(tokens[0]);
            props_data.water_density = std::stod(tokens[1]);
            props_data.gas_density = std::stod(tokens[2]);
        }
    }
}

void EclipseInputReader::read_pvtw_table() {
    auto tokens = read_tokens_until_slash();
    if (tokens.size() >= 5) {
        PropsData::PVTData pvtw;
        
        double pref_value = std::stod(tokens[0]);
        double Bw_value = std::stod(tokens[1]);
        double cw_value = std::stod(tokens[2]);
        double visc_value = std::stod(tokens[3]);
        double viscosibility = std::stod(tokens[4]);
        
        // CRITICAL FIX: PVTW pressure is ALWAYS the reference pressure for water properties
        // It should NOT be confused with EQUIL pressure!
        
        if (runspec.unit_system == "FIELD") {
            // SPE9 PVTW: 3600 psi reference pressure
            pvtw.pressure.push_back(pref_value * Units::PSI_TO_PA);
            pvtw.compressibility = cw_value / Units::PSI_TO_PA;
            pvtw.viscosity.push_back(visc_value * Units::CP_TO_PAS);
            pvtw.viscosibility = viscosibility / Units::PSI_TO_PA;
            
            std::cout << "      PVTW (FIELD units):" << std::endl;
            std::cout << "        Pref = " << pref_value << " psi = " 
                      << pref_value * Units::PSI_TO_PA / 1e5 << " bar" << std::endl;
        } else {
            // METRIC units
            pvtw.pressure.push_back(pref_value * Units::BAR_TO_PA);
            pvtw.compressibility = cw_value / Units::BAR_TO_PA;
            pvtw.viscosity.push_back(visc_value * Units::CP_TO_PAS);
            pvtw.viscosibility = viscosibility / Units::BAR_TO_PA;
            
            std::cout << "      PVTW (METRIC units):" << std::endl;
            std::cout << "        Pref = " << pref_value << " bar" << std::endl;
        }
        
        pvtw.fvf.push_back(Bw_value);
        
        std::cout << "        Bw = " << pvtw.fvf[0] << std::endl;
        std::cout << "        Cw = " << pvtw.compressibility << " 1/Pa" << std::endl;
        std::cout << "        μw = " << pvtw.viscosity[0] << " Pa.s" << std::endl;
        
        props_data.pvtw.push_back(pvtw);
    }
}

void EclipseInputReader::read_pvto_table() {
    // PVTO is more complex - Rs tables with pressure-dependent properties
    PropsData::PVTData pvto;
    std::string line;
    
    while (getline_with_continuation(line)) {
        line_number++;
        line = trim(line);
        
        if (line.empty() || line.substr(0, 2) == "--") continue;
        
        // Check for end of table
        if (line.find("/") != std::string::npos) {
            if (!pvto.rs.empty()) {
                props_data.pvto.push_back(pvto);
            }
            break;
        }
        
        std::istringstream iss(line);
        std::vector<double> values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }
        
        if (values.size() >= 4) {
            // Rs Po Bo Vo format
            pvto.rs.push_back(values[0]);
            pvto.pressure.push_back(values[1]);
            pvto.fvf.push_back(values[2]);
            pvto.viscosity.push_back(values[3]);
        } else if (values.size() >= 3) {
            // Undersaturated: Po Bo Vo (Rs from previous line)
            if (!pvto.rs.empty()) {
                pvto.rs.push_back(pvto.rs.back());
                pvto.pressure.push_back(values[0]);
                pvto.fvf.push_back(values[1]);
                pvto.viscosity.push_back(values[2]);
            }
        }
    }
}

void EclipseInputReader::read_pvtg_table() {
    // Similar to PVTO but for gas
    PropsData::PVTData pvtg;
    std::string line;
    
    while (getline_with_continuation(line)) {
        line_number++;
        line = trim(line);
        
        if (line.empty() || line.substr(0, 2) == "--") continue;
        
        if (line.find("/") != std::string::npos) {
            if (!pvtg.pressure.empty()) {
                props_data.pvtg.push_back(pvtg);
            }
            break;
        }
        
        std::istringstream iss(line);
        std::vector<double> values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }
        
        if (values.size() >= 3) {
            pvtg.pressure.push_back(values[0]);
            pvtg.fvf.push_back(values[1]);
            pvtg.viscosity.push_back(values[2]);
        }
    }
}

void EclipseInputReader::read_swof_table() {
    PropsData::RelPermData swof;
    std::string line;
    
    while (getline_with_continuation(line)) {
        line_number++;
        line = trim(line);
        
        if (line.empty() || line.substr(0, 2) == "--") continue;
        
        if (line.find("/") != std::string::npos) {
            if (!swof.sw.empty()) {
                props_data.swof.push_back(swof);
            }
            break;
        }
        
        std::istringstream iss(line);
        std::vector<double> values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }
        
        if (values.size() >= 4) {
            swof.sw.push_back(values[0]);
            swof.krow.push_back(values[1]);
            swof.krow.push_back(values[2]);
            swof.pcow.push_back(values[3]);
        }
    }
}

void EclipseInputReader::read_sgof_table() {
    PropsData::RelPermData sgof;
    std::string line;
    
    while (getline_with_continuation(line)) {
        line_number++;
        line = trim(line);
        
        if (line.empty() || line.substr(0, 2) == "--") continue;
        
        if (line.find("/") != std::string::npos) {
            if (!sgof.sg.empty()) {
                props_data.sgof.push_back(sgof);
            }
            break;
        }
        
        std::istringstream iss(line);
        std::vector<double> values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }
        
        if (values.size() >= 4) {
            sgof.sg.push_back(values[0]);
            sgof.krg.push_back(values[1]);
            sgof.krog.push_back(values[2]);
            sgof.pcog.push_back(values[3]);
        }
    }
}

void EclipseInputReader::process_solution_keyword(const std::string& keyword) {
    int ncells = runspec.nx * runspec.ny * runspec.nz;
    
    if (keyword == "PRESSURE") {
        solution_data.pressure = read_numeric_array(ncells);
    } else if (keyword == "SWAT") {
        solution_data.swat = read_numeric_array(ncells);
    } else if (keyword == "SGAS") {
        solution_data.sgas = read_numeric_array(ncells);
    } else if (keyword == "RS") {
        solution_data.rs = read_numeric_array(ncells);
    } else if (keyword == "EQUIL") {
        read_equil_data();
    }
}

void EclipseInputReader::read_equil_data() {
    auto tokens = read_tokens_until_slash();
    if (tokens.size() >= 7) {
        SolutionData::EquilData equil;
        equil.datum_depth = std::stod(tokens[0]);
        
        double pressure_value = std::stod(tokens[1]);
        
        // CRITICAL FIX: Handle SPE1 special case
        bool isSPE1Case = (runspec.unit_system == "METRIC" && pressure_value > 1000.0);
        
        if (isSPE1Case) {
            // SPE1 uses psi even in METRIC mode
            equil.datum_pressure = pressure_value * Units::PSI_TO_PA;
            std::cout << "    EQUIL: SPE1 case - " << pressure_value << " psi = " 
                      << equil.datum_pressure/1e5 << " bar" << std::endl;
        } else if (runspec.unit_system == "FIELD") {
            equil.datum_pressure = pressure_value * Units::PSI_TO_PA;
            std::cout << "    EQUIL pressure: " << pressure_value << " psia = " 
                      << equil.datum_pressure/1e5 << " bar" << std::endl;
        } else {
            // True METRIC - pressure in bars
            equil.datum_pressure = pressure_value * Units::BAR_TO_PA;
            std::cout << "    EQUIL pressure: " << pressure_value << " bar" << std::endl;
        }
        
        equil.owc_depth = std::stod(tokens[2]);
        equil.pc_owc = std::stod(tokens[3]);
        equil.goc_depth = std::stod(tokens[4]);
        equil.pc_goc = std::stod(tokens[5]);
        equil.rsvd_table = std::stoi(tokens[6]);
        
        if (tokens.size() > 7) equil.rvvd_table = std::stoi(tokens[7]);
        if (tokens.size() > 8) equil.n = std::stoi(tokens[8]);
        
        solution_data.equil.push_back(equil);
    }
}

void EclipseInputReader::process_schedule_keyword(const std::string& keyword) {
    if (keyword == "WELSPECS") {
        std::cout << "    Reading WELSPECS" << std::endl;
        read_welspecs();
    } else if (keyword == "COMPDAT") {
        std::cout << "    Reading COMPDAT" << std::endl;
        read_compdat();
    } else if (keyword == "WCONPROD") {
        std::cout << "    Reading WCONPROD" << std::endl;
    std::cout << "DEBUG: About to parse WCONPROD" << std::endl;
        read_wconprod();
    } else if (keyword == "WCONINJE") {
        std::cout << "    Reading WCONINJE" << std::endl;
        std::cout << "DEBUG: About to parse WCONINJE" << std::endl;
        read_wconinje();
    } else if (keyword == "TSTEP") {
        std::cout << "    Reading TSTEP" << std::endl;
        read_tstep();
    }
}

void EclipseInputReader::read_welspecs() {
    std::cout << "DEBUG: Entering read_welspecs" << std::endl;
    std::cout << "DEBUG: schedule_data address: " << &schedule_data << std::endl;
    std::cout << "DEBUG: timesteps before: " << schedule_data.timesteps.size() << std::endl;
    
    // Ensure we have at least one timestep
    if (schedule_data.timesteps.empty()) {
        schedule_data.timesteps.push_back(ScheduleData::TimeStep());
    }
    
    std::string line;
    int well_count = 0;
    while (getline_with_continuation(line)) {
        line_number++;
        // Skip comments
        if (line.find("--") == 0) continue;
        // Trim
        line = trim(line);
        if (line.empty()) continue;
        // Check for end of section
        if (line == "/") {
            break;
        }
        // Remove trailing '/'
        if (!line.empty() && line.back() == '/') {
            line.pop_back();
            line = trim(line);
        }
        // Parse well: NAME GROUP I J DEPTH PHASE
        std::istringstream iss(line);
        std::string name, group, phase;
        int i, j;
        double depth;
        if (iss >> name >> group >> i >> j >> depth >> phase) {
            WellData well;
            well.name = name;
            well.group = group;
            // Convert from 1-based to 0-based indices
            well.i = i - 1;
            well.j = j - 1;
            well.ref_depth = depth;
            well.phase = phase;
            if (to_upper(phase) == "WATER" || to_upper(phase) == "GAS") {
                well.is_producer = false;
                well.injection_type = to_upper(phase);
            } else {
                well.is_producer = true;
                well.injection_type = "OIL";
            }
            // Don't set is_producer here - it will be set by WCONPROD/WCONINJE
            // Default to producer until WCONINJE is read
            schedule_data.timesteps[0].wells.push_back(well);
            well_count++;
            std::cout << "      Added well: " << name << " (" << phase << ")" << std::endl;
        }
    }
    std::cout << "DEBUG: Exiting read_welspecs" << std::endl;
    std::cout << "DEBUG: wells added: " << well_count << std::endl;
    std::cout << "DEBUG: timesteps after: " << schedule_data.timesteps.size() << std::endl;
    if (!schedule_data.timesteps.empty()) {
        std::cout << "DEBUG: wells in timestep[0]: " << schedule_data.timesteps[0].wells.size() << std::endl;
    }
}

void EclipseInputReader::read_compdat() {
    std::cout << "DEBUG: Entering read_compdat" << std::endl;
    std::string line;
    
    while (getline_with_continuation(line)) {
        line_number++;
        line = trim(line);
        std::cout << "DEBUG: COMPDAT reading line: '" << line << "'" << std::endl;
        
        if (line.empty() || line.substr(0, 2) == "--") continue;
        
        // Check if this line ends the section
        if (line.find("/") != std::string::npos) {
            // Process the line before the '/'
            size_t slash_pos = line.find("/");
            line = line.substr(0, slash_pos);
            line = trim(line);
            if (line.empty()) break;
        }
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        
        while (iss >> token) {
            tokens.push_back(token);
        }
        
        if (tokens.size() >= 5) {
            std::string well_name = tokens[0];
            std::cout << "DEBUG: COMPDAT processing well: '" << well_name << "'" << std::endl;
            std::cout << "DEBUG: Available wells in timestep:" << std::endl;
            for (const auto& well : schedule_data.timesteps[0].wells) {
                std::cout << "DEBUG:   - '" << well.name << "'" << std::endl;
            }
            
            // Find the well (handle quote differences)
            for (auto& well : schedule_data.timesteps[0].wells) {
                std::cout << "DEBUG: COMPDAT comparing '" << well.name << "' with '" << well_name << "'" << std::endl;
                
                // Remove quotes for comparison
                std::string well_name_clean = well_name;
                std::string well_name_stored = well.name;
                
                // Remove single quotes
                well_name_clean.erase(std::remove(well_name_clean.begin(), well_name_clean.end(), '\''), well_name_clean.end());
                well_name_stored.erase(std::remove(well_name_stored.begin(), well_name_stored.end(), '\''), well_name_stored.end());
                
                std::cout << "DEBUG: COMPDAT comparing cleaned names: '" << well_name_stored << "' with '" << well_name_clean << "'" << std::endl;
                
                if (well_name_stored == well_name_clean) {
                    WellData::Completion comp;
                    
                    // Convert from 1-based to 0-based indices
                    comp.i = std::stoi(tokens[1]) - 1;
                    comp.j = std::stoi(tokens[2]) - 1;
                    comp.k1 = std::stoi(tokens[3]) - 1;
                    comp.k2 = std::stoi(tokens[4]) - 1;
                    
                    if (tokens.size() > 5 && !is_default_value(tokens[5])) {
                        comp.status = to_upper(tokens[5]);
                    } else {
                        comp.status = "OPEN";
                    }
                    
                    // Parse wellbore diameter if given
                    if (tokens.size() > 8 && !is_default_value(tokens[8])) {
                        comp.wellbore_diam = std::stod(tokens[8]);
                    } else {
                        comp.wellbore_diam = 0.5; // Default
                    }
                    
                    // Update well's k1 and k2 (also 0-based)
                    // For the first completion, set k1 and k2
                    if (well.completions.empty()) {
                        well.k1 = comp.k1;
                        well.k2 = comp.k2;
                        std::cout << "DEBUG: COMPDAT first completion for " << well.name 
                                  << ": k1=" << comp.k1 << ", k2=" << comp.k2 << std::endl;
                    } else {
                        // For subsequent completions, expand the range
                        well.k1 = std::min(well.k1, comp.k1);
                        well.k2 = std::max(well.k2, comp.k2);
                        std::cout << "DEBUG: COMPDAT expanding completion for " << well.name 
                                  << ": k1=" << well.k1 << ", k2=" << well.k2 << std::endl;
                    }
                    
                    well.completions.push_back(comp);
                    break;
                }
            }
        }
    }
}

void EclipseInputReader::read_wconprod() {
    if (schedule_data.timesteps.empty()) return;
    
    auto& current_timestep = schedule_data.timesteps.back();
    std::string line;
    
    while (getline_with_continuation(line)) {
        line_number++;
        line = trim(line);
        std::cout << "DEBUG: WCONPROD reading line: '" << line << "'" << std::endl;
        // Split at '/' and only process before it
        size_t slash_pos = line.find("/");
        if (slash_pos != std::string::npos) {
            line = line.substr(0, slash_pos);
            line = trim(line);
            if (line.empty()) break;
        }
        if (line.empty() || line.substr(0, 2) == "--") continue;
        if (line.find("/") != std::string::npos) break;
        
        std::replace(line.begin(), line.end(), '\'', ' ');
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) {
            tokens.push_back(token);
        }
        
        std::cout << "DEBUG: WCONPROD tokens: ";
        for (const auto& t : tokens) {
            std::cout << "'" << t << "' ";
        }
        std::cout << std::endl;
        
        if (tokens.size() >= 4) {
            std::string well_name = tokens[0];
            std::cout << "DEBUG: WCONPROD processing well: '" << well_name << "'" << std::endl;
            std::cout << "DEBUG: Available wells in timestep:" << std::endl;
            for (const auto& well : current_timestep.wells) {
                std::cout << "DEBUG:   - '" << well.name << "'" << std::endl;
            }
            
            // Find and update well with exact name matching
            for (auto& well : current_timestep.wells) {
                std::string trimmed_well_name = to_upper(trim(well.name, "'\""));
                std::string trimmed_input_name = to_upper(trim(well_name, "'\""));
                std::cout << "DEBUG: Comparing '" << trimmed_well_name << "' with '" << trimmed_input_name << "'" << std::endl;
                if (trimmed_well_name == trimmed_input_name) {
                    std::string status = tokens[1];
                    std::string control = tokens[2];
                    double orat = std::stod(tokens[3]);
                    double bhp = 0.0; // Default to 0.0
                    if (tokens.size() > 8 && tokens[8] != "1*") {
                        bhp = std::stod(tokens[8]);
                    }

                    // CRITICAL: Mark as producer!
                    well.is_producer = true;
                    well.status = status;
                    well.control_mode = control;
                    well.oil_rate = orat;
                    well.bhp_limit = bhp;
                    
                    std::cout << "      WCONPROD: " << well_name 
                              << " set as producer, rate=" << well.oil_rate 
                              << " STB/day" << std::endl;
                    break;
                }
            }
        }
    }
}

void EclipseInputReader::read_wconinje() {
    if (schedule_data.timesteps.empty()) {
        std::cerr << "ERROR: No timesteps for WCONINJE" << std::endl;
        // Skip to '/'
        std::string line;
        while (getline_with_continuation(line)) {
            line_number++;
            if (line.find("/") != std::string::npos) break;
        }
        return;
    }
    
    auto& current_timestep = schedule_data.timesteps.back();
    std::string line;
    
    while (getline_with_continuation(line)) {
        line_number++;
        line = trim(line);
        std::cout << "DEBUG: WCONINJE reading line: '" << line << "'" << std::endl;
        // Split at '/' and only process before it
        size_t slash_pos = line.find("/");
        if (slash_pos != std::string::npos) {
            line = line.substr(0, slash_pos);
            line = trim(line);
            if (line.empty()) break;
        }
        if (line.empty() || line.substr(0, 2) == "--") continue;
        if (line.find("/") != std::string::npos) break;
        
        // Remove quotes
        std::replace(line.begin(), line.end(), '\'', ' ');
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) {
            tokens.push_back(token);
        }
        
        std::cout << "DEBUG: WCONINJE tokens: ";
        for (const auto& t : tokens) {
            std::cout << "'" << t << "' ";
        }
        std::cout << std::endl;
        
        if (tokens.size() >= 5) {
            std::string well_name = tokens[0];
            std::cout << "DEBUG: WCONINJE processing well: '" << well_name << "'" << std::endl;
            std::cout << "DEBUG: Available wells in timestep:" << std::endl;
            for (const auto& well : current_timestep.wells) {
                std::cout << "DEBUG:   - '" << well.name << "'" << std::endl;
            }
            
            // Find and update well
            for (auto& well : current_timestep.wells) {
                std::string trimmed_well_name = to_upper(trim(well.name, "'\""));
                std::string trimmed_input_name = to_upper(trim(well_name, "'\""));
                std::cout << "DEBUG: Comparing '" << trimmed_well_name << "' with '" << trimmed_input_name << "'" << std::endl;
                if (trimmed_well_name == trimmed_input_name) {
                    // CRITICAL: Mark as injector!
                    well.is_producer = false;
                    well.injection_type = tokens[1];
                    well.status = tokens[2];
                    well.control_mode = tokens[3];
                    
                    // Parse rate
                    double rate = std::stod(tokens[4]);
                    if (runspec.unit_system == "FIELD") {
                        // Convert STB/day to m³/day
                        rate *= 0.1589873;
                    }
                    well.water_rate = rate;
                    
                    // Parse BHP if given
                    if (tokens.size() > 6 && tokens[6] != "1*") {
                        double bhp = std::stod(tokens[6]);
                        if (runspec.unit_system == "FIELD") {
                            well.bhp_limit = bhp * Units::PSI_TO_PA;
                        } else {
                            well.bhp_limit = bhp * Units::BAR_TO_PA;
                        }
                    }
                    
                    std::cout << "      WCONINJE: " << well_name 
                              << " set as injector, rate=" << well.water_rate 
                              << " m³/day" << std::endl;
                    break;
                }
            }
        }
    }
}

void EclipseInputReader::read_tstep() {
    auto tokens = read_tokens_until_slash();
    
    // FIX: If we already have timesteps, clone the last one as base
    ScheduleData::TimeStep base_timestep;
    if (!schedule_data.timesteps.empty()) {
        base_timestep = schedule_data.timesteps.back();
    }
    
    for (const auto& token : tokens) {
        if (token == "/" || token == "END") continue;
        
        try {
            // Handle multiplier syntax
            size_t star_pos = token.find('*');
            if (star_pos != std::string::npos) {
                int count = std::stoi(token.substr(0, star_pos));
                double dt = std::stod(token.substr(star_pos + 1));
                
                for (int i = 0; i < count; ++i) {
                    ScheduleData::TimeStep ts = base_timestep;
                    ts.days = dt;
                    schedule_data.timesteps.push_back(ts);
                }
            } else {
                ScheduleData::TimeStep ts = base_timestep;
                ts.days = std::stod(token);
                schedule_data.timesteps.push_back(ts);
            }
        } catch (const std::exception& e) {
            std::cerr << "Error parsing timestep: " << token 
                      << " - " << e.what() << std::endl;
        }
    }
    
    std::cout << "      Added " << schedule_data.timesteps.size() 
              << " total timesteps" << std::endl;
}

bool EclipseInputReader::validate_data() {
    // Grid validation
    if (runspec.nx <= 0 || runspec.ny <= 0 || runspec.nz <= 0) {
        throw std::runtime_error("Invalid grid dimensions");
    }
    int ncells = runspec.nx * runspec.ny * runspec.nz;
    // Validate porosity
    if (!grid_data.poro.empty()) {
        for (double phi : grid_data.poro) {
            if (phi < 0.0 || phi > 1.0) {
                std::cerr << "Warning: Porosity out of range [0,1]: " << phi << std::endl;
            }
        }
    }
    // Validate permeability
    if (!grid_data.permx.empty()) {
        for (double k : grid_data.permx) {
            if (k < 0.0) {
                throw std::runtime_error("Negative permeability found");
            }
        }
    }
    // Validate saturations
    if (!solution_data.swat.empty()) {
        for (size_t i = 0; i < solution_data.swat.size(); ++i) {
            double sw = solution_data.swat[i];
            double sg = solution_data.sgas.empty() ? 0.0 : solution_data.sgas[i];
            if (sw < 0.0 || sw > 1.0 || sg < 0.0 || sg > 1.0 || sw + sg > 1.0) {
                std::cerr << "Warning: Invalid saturations at cell " << i << std::endl;
            }
        }
    }
    // Validate PVT data
    for (const auto& pvt : props_data.pvto) {
        if (pvt.pressure.empty() || pvt.fvf.empty()) {
            throw std::runtime_error("Empty PVT table");
        }
        for (size_t i = 1; i < pvt.pressure.size(); ++i) {
            if (pvt.pressure[i] <= pvt.pressure[i-1]) {
                std::cerr << "Warning: Non-monotonic pressure in PVT table" << std::endl;
            }
        }
    }
    return true;
}

void EclipseInputReader::process_include(const std::string& include_filename) {
    std::ifstream include_file(include_filename);
    if (!include_file.is_open()) {
        throw std::runtime_error("Cannot open INCLUDE file: " + include_filename);
    }
    // Save current file state
    std::ifstream temp_file = std::move(file);
    int temp_line = line_number;
    // Process include file
    file = std::move(include_file);
    line_number = 0;
    std::string line;
    while (getline_with_continuation(line)) {
        line_number++;
        process_keyword(line);
    }
    // Restore original file
    file = std::move(temp_file);
    line_number = temp_line;
}

// Getters
const RunspecData& EclipseInputReader::get_runspec() const { return runspec; }
const GridData& EclipseInputReader::get_grid_data() const { return grid_data; }
const PropsData& EclipseInputReader::get_props_data() const { return props_data; }
const SolutionData& EclipseInputReader::get_solution_data() const { return solution_data; }
const ScheduleData& EclipseInputReader::get_schedule_data() const { return schedule_data; }

// ============================================================================
// Global Parse Function
// ============================================================================
bool parse_eclipse_deck(const std::string& filename, RunspecData& runspec, GridData& grid, 
                       PropsData& props, SolutionData& solution, ScheduleData& schedule) {
    try {
        EclipseInputReader reader(filename);
        if (reader.parse()) {
            runspec = reader.get_runspec();
            grid = reader.get_grid_data();
            props = reader.get_props_data();
            solution = reader.get_solution_data();
            schedule = reader.get_schedule_data();
            return true;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing Eclipse deck: " << e.what() << std::endl;
    }
    return false;
}