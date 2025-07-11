// simulator.cpp - Black Oil Reservoir Simulator Implementation
#include "simulator.h"
#include "eclipse_input_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cctype>
#include <map>

// ==================== UTILITY FUNCTIONS ====================

// Convert string to uppercase
std::string toUpper(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

// Trim whitespace from string
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

// Trim specific characters from string
std::string trim(const std::string& str, const std::string& chars) {
    size_t first = str.find_first_not_of(chars);
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(chars);
    return str.substr(first, (last - first + 1));
}

// Check if line is a comment
bool isComment(const std::string& line) {
    std::string trimmed = trim(line);
    return trimmed.empty() || trimmed.substr(0, 2) == "--";
}

// Read next data line (skip comments and empty lines)
bool getDataLine(std::ifstream& file, std::string& line) {
    std::string originalLine;
    while (std::getline(file, originalLine)) {
        line = originalLine;
        
        // Skip empty lines
        if (line.empty()) continue;
        
        // Trim leading whitespace
        size_t firstNonSpace = line.find_first_not_of(" \t");
        if (firstNonSpace == std::string::npos) continue;
        
        // Skip comment lines
        if (line.substr(firstNonSpace, 2) == "--") continue;
        
        // Remove inline comments
        size_t commentPos = line.find("--");
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }
        
        // Trim whitespace
        line = trim(line);
        if (!line.empty()) {
            return true;
        }
    }
    return false;
}

// Parse a line of numbers
std::vector<PetscReal> parseNumbers(const std::string& line) {
    std::vector<PetscReal> numbers;
    std::istringstream iss(line);
    std::string token;
    
    while (iss >> token) {
        // Skip slashes and terminators
        if (token == "/" || token == "1*") continue;
        
        // Check for repeat syntax (e.g., "600*0.087")
        size_t starPos = token.find('*');
        if (starPos != std::string::npos && starPos > 0) {
            try {
                PetscInt count = std::stoi(token.substr(0, starPos));
                PetscReal val = std::stod(token.substr(starPos + 1));
                for (PetscInt i = 0; i < count; i++) {
                    numbers.push_back(val);
                }
            } catch (...) {
                // Skip invalid tokens
            }
        } else {
            try {
                numbers.push_back(std::stod(token));
            } catch (...) {
                // Skip non-numeric tokens
            }
        }
    }
    return numbers;
}

// Linear interpolation
PetscReal interpolate(PetscReal x, const std::vector<PetscReal>& xData, 
                     const std::vector<PetscReal>& yData) {
    if (xData.size() != yData.size() || xData.empty()) return 0.0;
    
    // Handle extrapolation
    if (x <= xData[0]) return yData[0];
    if (x >= xData.back()) return yData.back();
    
    // Find interval
    for (size_t i = 1; i < xData.size(); i++) {
        if (x <= xData[i]) {
            PetscReal t = (x - xData[i-1]) / (xData[i] - xData[i-1]);
            return yData[i-1] + t * (yData[i] - yData[i-1]);
        }
    }
    return yData.back();
}

// ==================== PVT FUNCTIONS ====================

// Get oil properties at given pressure and Rs
void getOilProperties(const FluidProps& fluid, PetscReal p, PetscReal Rs,
                     PetscReal& Bo, PetscReal& muo) {
    // Default values if no data
    if (fluid.pvto.Rs.empty()) {
        Bo = 1.2;
        muo = 2e-3;
        return;
    }
    
    // Find correct Rs interval
    size_t idx = 0;
    for (size_t i = 0; i < fluid.pvto.Rs.size() - 1; i++) {
        if (Rs <= fluid.pvto.Rs[i+1]) {
            idx = i;
            break;
        }
    }
    if (Rs > fluid.pvto.Rs.back()) idx = fluid.pvto.Rs.size() - 1;
    
    // Check if saturated or undersaturated
    PetscReal pbub = interpolate(Rs, fluid.pvto.Rs, fluid.pvto.pbub);
    
    if (p <= pbub) {
        // Saturated
        Bo = interpolate(Rs, fluid.pvto.Rs, fluid.pvto.Bo);
        muo = interpolate(Rs, fluid.pvto.Rs, fluid.pvto.muo);
    } else {
        // Undersaturated - use the table for current Rs
        if (idx < fluid.pvto.p_usat.size() && 
            idx < fluid.pvto.Bo_usat.size() && 
            idx < fluid.pvto.muo_usat.size()) {
            Bo = interpolate(p, fluid.pvto.p_usat[idx], fluid.pvto.Bo_usat[idx]);
            muo = interpolate(p, fluid.pvto.p_usat[idx], fluid.pvto.muo_usat[idx]);
        } else {
            // Fallback to saturated values if undersaturated data not available
            Bo = fluid.pvto.Bo[idx];
            muo = fluid.pvto.muo[idx];
        }
    }
}

// Get water properties
void getWaterProperties(const FluidProps& fluid, PetscReal p,
                       PetscReal& Bw, PetscReal& muw) {
    Bw = fluid.Bwref * (1.0 - fluid.cw * (p - fluid.pwref));
    muw = fluid.muw * (1.0 + fluid.muw_p * (p - fluid.pwref));
    
    // Ensure reasonable values
    if (Bw <= 0) Bw = 1.0;
    if (muw <= 0) muw = 5e-4;
}

// Get gas properties
void getGasProperties(const FluidProps& fluid, PetscReal p,
                     PetscReal& Bg, PetscReal& mug) {
    if (fluid.pvdg_p.empty()) {
        // Default gas properties
        PetscReal pref = 1e5;  // 1 bar
        Bg = pref / p;
        mug = 2e-5;
        return;
    }
    
    Bg = interpolate(p, fluid.pvdg_p, fluid.pvdg_Bg);
    mug = interpolate(p, fluid.pvdg_p, fluid.pvdg_mug);
}

// Get relative permeabilities
void getRelPerm(const FluidProps& fluid, PetscReal Sw, PetscReal Sg,
               PetscReal& krw, PetscReal& kro, PetscReal& krg) {
    PetscReal So = 1.0 - Sw - Sg;
    
    // Default relative perms if no tables
    if (fluid.swof_Sw.empty()) {
        krw = Sw * Sw;
        kro = So * So;
        krg = Sg * Sg;
        return;
    }
    
    // Water relative permeability
    krw = interpolate(Sw, fluid.swof_Sw, fluid.swof_krw);
    
    // Gas relative permeability
    if (!fluid.sgof_Sg.empty()) {
        krg = interpolate(Sg, fluid.sgof_Sg, fluid.sgof_krg);
    } else {
        krg = Sg * Sg;
    }
    
    // Oil relative permeability (Stone's model II)
    PetscReal krow = interpolate(Sw, fluid.swof_Sw, fluid.swof_krow);
    PetscReal krog = 1.0;
    if (!fluid.sgof_Sg.empty()) {
        krog = interpolate(Sg, fluid.sgof_Sg, fluid.sgof_krog);
    }
    
    if (krow > 0.0 && krog > 0.0) {
        PetscReal kro_max = fluid.swof_krow.empty() ? 1.0 : fluid.swof_krow[0];
        kro = krow * krog / kro_max;
    } else {
        kro = 0.0;
    }
}

// ==================== INITIALIZATION ====================

// Initialize with equilibration
void initializeEquilibration(AppCtx* ctx) {
    Vec local;
    Field ***x;
    
    DMGetLocalVector(ctx->da, &local);
    DMDAVecGetArray(ctx->da, local, &x);
    
    DMDALocalInfo info;
    DMDAGetLocalInfo(ctx->da, &info);
    
    PetscReal g = 9.81;  // Gravity
    
    // Check if essential arrays are initialized
    PetscInt ncells = ctx->grid.Nx * ctx->grid.Ny * ctx->grid.Nz;
    
    // Ensure arrays are populated
    if (ctx->grid.DZ.empty()) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: DZ not specified, using default 10m\n");
        ctx->grid.DZ.resize(ncells, 10.0);
    }
    if (ctx->grid.permX.empty()) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: PERMX not specified, using default 100mD\n");
        ctx->grid.permX.resize(ncells, 100.0 * 9.869233e-16);
        ctx->grid.permY = ctx->grid.permX;
        ctx->grid.permZ.resize(ncells, 10.0 * 9.869233e-16);
    }
    if (ctx->grid.porosity.empty()) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: PORO not specified, using default 0.2\n");
        ctx->grid.porosity.resize(ncells, 0.2);
    }
    
    // Set default equilibration data if not specified
    if (ctx->equil.datum == 0.0) {
        ctx->equil.datum = 2000.0;  // 2000m depth
        ctx->equil.pres_datum = 2e7;  // 200 bar
        ctx->equil.zwoc = 2100.0;  // 100m below datum
        ctx->equil.zgoc = 1900.0;  // 100m above datum
        PetscPrintf(PETSC_COMM_WORLD, "Using default equilibration data\n");
    }
    
    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                PetscInt idx = k * ctx->grid.Ny * ctx->grid.Nx + j * ctx->grid.Nx + i;
                
                // Calculate cell center depth
                PetscReal depth = 0.0;
                
                // Handle TOPS array
                if (!ctx->grid.TOPS.empty()) {
                    size_t top_idx = j * ctx->grid.Nx + i;
                    if (top_idx < ctx->grid.TOPS.size()) {
                        depth = ctx->grid.TOPS[top_idx];
                    }
                } else {
                    // If no TOPS, assume top of reservoir at datum depth minus some offset
                    depth = ctx->equil.datum - 100.0;  // Start 100m above datum
                }
                
                // Add depths of layers above
                for (PetscInt kk = 0; kk <= k; kk++) {
                    size_t kidx = kk * ctx->grid.Ny * ctx->grid.Nx + j * ctx->grid.Nx + i;
                    if (kidx < ctx->grid.DZ.size()) {
                        if (kk < k) {
                            depth += ctx->grid.DZ[kidx];
                        } else {
                            depth += 0.5 * ctx->grid.DZ[kidx];
                        }
                    }
                }
                
                // Calculate pressure using OWC/GOC
                PetscReal p = ctx->equil.pres_datum;
                
                if (depth < ctx->equil.datum) {
                    // Above datum - integrate upward
                    PetscReal rho_avg = ctx->fluid.rhoO_surf;  // Approximate
                    p -= rho_avg * g * (ctx->equil.datum - depth);
                } else {
                    // Below datum - integrate downward
                    PetscReal rho_avg = ctx->fluid.rhoO_surf;
                    p += rho_avg * g * (depth - ctx->equil.datum);
                }
                
                x[k][j][i].p = p;
                
                // Set saturations from Eclipse deck data if available, otherwise use depth-based defaults
                if (idx < static_cast<PetscInt>(ctx->grid.swat.size()) && 
                    idx < static_cast<PetscInt>(ctx->grid.sgas.size())) {
                    // Use SWAT and SGAS from Eclipse deck
                    x[k][j][i].Sw = ctx->grid.swat[idx];
                    x[k][j][i].Sg = ctx->grid.sgas[idx];
                } else {
                    // Fallback to depth-based saturations
                    if (depth < ctx->equil.zgoc) {
                        // Gas cap
                        x[k][j][i].Sg = 0.8;
                        x[k][j][i].Sw = 0.2;
                    } else if (depth > ctx->equil.zwoc) {
                        // Water zone
                        x[k][j][i].Sw = 0.8;
                        x[k][j][i].Sg = 0.0;
                    } else {
                        // Oil zone
                        x[k][j][i].Sw = 0.2;  // Connate water
                        x[k][j][i].Sg = 0.0;
                    }
                }
            }
        }
    }
    
    DMDAVecRestoreArray(ctx->da, local, &x);
    DMLocalToGlobalBegin(ctx->da, local, INSERT_VALUES, ctx->solution);
    DMLocalToGlobalEnd(ctx->da, local, INSERT_VALUES, ctx->solution);
    DMRestoreLocalVector(ctx->da, &local);
    
    VecCopy(ctx->solution, ctx->solution_old);
}

// ==================== WELL MODEL ====================

// Calculate well index
PetscReal calculateWellIndex(const Well& well, const Grid3D& grid, 
                           PetscInt i, PetscInt j, PetscInt k) {
    PetscInt idx = k * grid.Ny * grid.Nx + j * grid.Nx + i;
    
    // Get permeabilities with bounds checking
    PetscReal kx = (idx < static_cast<PetscInt>(grid.permX.size())) ? grid.permX[idx] : 100.0 * 9.869233e-16;
    PetscReal ky = (idx < static_cast<PetscInt>(grid.permY.size())) ? grid.permY[idx] : kx;
    PetscReal kz = (idx < static_cast<PetscInt>(grid.permZ.size())) ? grid.permZ[idx] : kx * 0.1;
    
    // Get cell dimensions with defaults
    PetscReal dx = (!grid.DX.empty() && idx < static_cast<PetscInt>(grid.DX.size())) ? grid.DX[idx] : 50.0;
    PetscReal dy = (!grid.DY.empty() && idx < static_cast<PetscInt>(grid.DY.size())) ? grid.DY[idx] : 50.0;
    PetscReal dz = (!grid.DZ.empty() && idx < static_cast<PetscInt>(grid.DZ.size())) ? grid.DZ[idx] : 10.0;
    
    // Peaceman's formula for vertical well
    PetscReal kh = std::sqrt(kx * ky);
    PetscReal req = 0.28 * std::sqrt(std::sqrt(ky/kx) * dx * dx + 
                                     std::sqrt(kx/ky) * dy * dy) /
                   (std::pow(ky/kx, 0.25) + std::pow(kx/ky, 0.25));
    
    PetscReal rw = well.rw;
    if (rw <= 0.0) rw = 0.1;  // Default 0.1m
    
    // Ensure req > rw to avoid log(0) or negative values
    if (req <= rw) req = rw * 1.1;
    
    PetscReal WI = 2.0 * PETSC_PI * kh * dz / std::log(req / rw);
    
    // Ensure reasonable well index
    if (WI <= 0.0 || std::isnan(WI) || std::isinf(WI)) {
        WI = 1e-12; // Small but finite value
    }
    
    return WI;
}

// ==================== RESIDUAL FUNCTION ====================

PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx_ptr) {
    AppCtx *ctx = (AppCtx*)ctx_ptr;
    DM da = ctx->da;
    DMDALocalInfo info;
    Field ***xlocal, ***flocal, ***xold;
    
    DMDAGetLocalInfo(da, &info);
    
    // Get local vectors
    Vec x_local, f_local, xold_local;
    DMGetLocalVector(da, &x_local);
    DMGetLocalVector(da, &f_local);
    DMGetLocalVector(da, &xold_local);
    
    DMGlobalToLocalBegin(da, x, INSERT_VALUES, x_local);
    DMGlobalToLocalEnd(da, x, INSERT_VALUES, x_local);
    
    DMGlobalToLocalBegin(da, ctx->solution_old, INSERT_VALUES, xold_local);
    DMGlobalToLocalEnd(da, ctx->solution_old, INSERT_VALUES, xold_local);
    
    DMDAVecGetArrayRead(da, x_local, &xlocal);
    DMDAVecGetArray(da, f_local, &flocal);
    DMDAVecGetArrayRead(da, xold_local, &xold);
    
    // Clear residual
    VecSet(f_local, 0.0);
    
    // Compute residual for each cell
    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                PetscInt idx = k * ctx->grid.Ny * ctx->grid.Nx + j * ctx->grid.Nx + i;
                
                // Bounds check
                if (idx >= static_cast<PetscInt>(ctx->grid.porosity.size())) {
                    PetscPrintf(PETSC_COMM_WORLD, "Error: idx=%d out of bounds (size=%d)\n", 
                               idx, static_cast<PetscInt>(ctx->grid.porosity.size()));
                    continue;
                }
                
                // Cell properties
                PetscReal phi0 = ctx->grid.porosity[idx];
                PetscReal V = 1.0;  // Default volume
                
                // Calculate cell volume
                if (!ctx->grid.DX.empty() && idx < static_cast<PetscInt>(ctx->grid.DX.size())) {
                    V = ctx->grid.DX[idx] * ctx->grid.DY[idx] * ctx->grid.DZ[idx];
                } else {
                    // Use uniform grid spacing if not specified
                    V = 50.0 * 50.0 * 10.0;  // Default 50m x 50m x 10m cells
                }
                
                // Current state
                PetscReal p = xlocal[k][j][i].p;
                PetscReal Sw = xlocal[k][j][i].Sw;
                PetscReal Sg = xlocal[k][j][i].Sg;
                PetscReal So = 1.0 - Sw - Sg;
                
                // Saturation constraints
                if (Sw < -0.1 || Sw > 1.1 || Sg < -0.1 || Sg > 1.1 || So < -0.1) {
                    PetscReal penalty = 1e8;
                    flocal[k][j][i].p = penalty * (PetscMax(-Sw, 0.0) + PetscMax(Sw - 1.0, 0.0) +
                                                  PetscMax(-Sg, 0.0) + PetscMax(Sg - 1.0, 0.0) +
                                                  PetscMax(-So, 0.0));
                    flocal[k][j][i].Sw = penalty * (PetscMax(-Sw, 0.0) + PetscMax(Sw - 1.0, 0.0));
                    flocal[k][j][i].Sg = penalty * (PetscMax(-Sg, 0.0) + PetscMax(Sg - 1.0, 0.0));
                    continue;
                }
                
                // Clamp saturations
                Sw = PetscMax(0.0, PetscMin(1.0, Sw));
                Sg = PetscMax(0.0, PetscMin(1.0, Sg));
                So = 1.0 - Sw - Sg;
                
                // Old state
                PetscReal p_old = xold[k][j][i].p;
                PetscReal Sw_old = xold[k][j][i].Sw;
                PetscReal Sg_old = xold[k][j][i].Sg;
                PetscReal So_old = 1.0 - Sw_old - Sg_old;
                
                // Rock compressibility
                PetscReal phi = phi0 * (1.0 + ctx->rock.crock * (p - ctx->rock.pref));
                PetscReal phi_old = phi0 * (1.0 + ctx->rock.crock * (p_old - ctx->rock.pref));
                
                // Get Rs from RSVD
                PetscReal Rs = 0.0, Rs_old = 0.0;
                if (!ctx->equil.rsvd_depth.empty()) {
                    PetscReal depth = 0.0;  // Should calculate actual depth
                    Rs = Rs_old = interpolate(depth, ctx->equil.rsvd_depth, ctx->equil.rsvd_rs);
                }
                
                // Fluid properties
                PetscReal Bo, Bw, Bg, muo, muw, mug;
                PetscReal Bo_old, Bw_old, Bg_old, muo_old, muw_old, mug_old;
                
                getOilProperties(ctx->fluid, p, Rs, Bo, muo);
                getWaterProperties(ctx->fluid, p, Bw, muw);
                getGasProperties(ctx->fluid, p, Bg, mug);
                
                getOilProperties(ctx->fluid, p_old, Rs_old, Bo_old, muo_old);
                getWaterProperties(ctx->fluid, p_old, Bw_old, muw_old);
                getGasProperties(ctx->fluid, p_old, Bg_old, mug_old);
                
                // Accumulation terms
                PetscReal acc_o = V * phi * So / Bo;
                PetscReal acc_w = V * phi * Sw / Bw;
                PetscReal acc_g = V * phi * (Sg / Bg + Rs * So / Bo);
                
                PetscReal acc_o_old = V * phi_old * So_old / Bo_old;
                PetscReal acc_w_old = V * phi_old * Sw_old / Bw_old;
                PetscReal acc_g_old = V * phi_old * (Sg_old / Bg_old + Rs_old * So_old / Bo_old);
                
                // Start with accumulation terms
                flocal[k][j][i].p = (acc_o - acc_o_old + acc_w - acc_w_old + acc_g - acc_g_old) / ctx->dt;
                flocal[k][j][i].Sw = (acc_w - acc_w_old) / ctx->dt;
                flocal[k][j][i].Sg = (acc_g - acc_g_old) / ctx->dt;
                
                // Flux contributions
                const int offsets[6][3] = {
                    {+1,0,0}, {-1,0,0}, {0,+1,0}, {0,-1,0}, {0,0,+1}, {0,0,-1}
                };
                
                for (int d = 0; d < 6; d++) {
                    PetscInt ni = i + offsets[d][0];
                    PetscInt nj = j + offsets[d][1];
                    PetscInt nk = k + offsets[d][2];
                    
                    if (ni < 0 || ni >= info.mx || nj < 0 || nj >= info.my || 
                        nk < 0 || nk >= info.mz) continue;
                    
                    // Transmissibility
                    PetscReal T = 0.0;
                    PetscInt nidx = nk * ctx->grid.Ny * ctx->grid.Nx + nj * ctx->grid.Nx + ni;
                    
                    // Bounds check
                    if (nidx >= static_cast<PetscInt>(ctx->grid.permX.size())) continue;
                    
                    // Get cell dimensions
                    PetscReal dx = (!ctx->grid.DX.empty() && idx < static_cast<PetscInt>(ctx->grid.DX.size())) ? ctx->grid.DX[idx] : 50.0;
                    PetscReal dy = (!ctx->grid.DY.empty() && idx < static_cast<PetscInt>(ctx->grid.DY.size())) ? ctx->grid.DY[idx] : 50.0;
                    PetscReal dz = (!ctx->grid.DZ.empty() && idx < static_cast<PetscInt>(ctx->grid.DZ.size())) ? ctx->grid.DZ[idx] : 10.0;
                    
                    PetscReal dx_n = (!ctx->grid.DX.empty() && nidx < static_cast<PetscInt>(ctx->grid.DX.size())) ? ctx->grid.DX[nidx] : 50.0;
                    PetscReal dy_n = (!ctx->grid.DY.empty() && nidx < static_cast<PetscInt>(ctx->grid.DY.size())) ? ctx->grid.DY[nidx] : 50.0;
                    PetscReal dz_n = (!ctx->grid.DZ.empty() && nidx < static_cast<PetscInt>(ctx->grid.DZ.size())) ? ctx->grid.DZ[nidx] : 10.0;
                    
                    if (d < 2) { // X-direction
                        PetscReal area = dy * dz;
                        PetscReal dist = 0.5 * (dx + dx_n);
                        PetscReal perm = 2.0 * ctx->grid.permX[idx] * ctx->grid.permX[nidx] /
                                        (ctx->grid.permX[idx] + ctx->grid.permX[nidx]);
                        T = perm * area / dist;
                    } else if (d < 4) { // Y-direction
                        PetscReal area = dx * dz;
                        PetscReal dist = 0.5 * (dy + dy_n);
                        PetscReal perm = 2.0 * ctx->grid.permY[idx] * ctx->grid.permY[nidx] /
                                        (ctx->grid.permY[idx] + ctx->grid.permY[nidx]);
                        T = perm * area / dist;
                    } else { // Z-direction
                        PetscReal area = dx * dy;
                        PetscReal dist = 0.5 * (dz + dz_n);
                        PetscReal perm = 2.0 * ctx->grid.permZ[idx] * ctx->grid.permZ[nidx] /
                                        (ctx->grid.permZ[idx] + ctx->grid.permZ[nidx]);
                        T = perm * area / dist;
                    }
                    
                    // Neighbor state
                    PetscReal pn = xlocal[nk][nj][ni].p;
                    PetscReal Swn = PetscMax(0.0, PetscMin(1.0, xlocal[nk][nj][ni].Sw));
                    PetscReal Sgn = PetscMax(0.0, PetscMin(1.0, xlocal[nk][nj][ni].Sg));
                    PetscReal Son = 1.0 - Swn - Sgn;
                    
                    // Neighbor fluid properties
                    PetscReal Bon, Bwn, Bgn, muon, muwn, mugn;
                    PetscReal Rsn = Rs;  // Simplified
                    
                    getOilProperties(ctx->fluid, pn, Rsn, Bon, muon);
                    getWaterProperties(ctx->fluid, pn, Bwn, muwn);
                    getGasProperties(ctx->fluid, pn, Bgn, mugn);
                    
                    // Phase densities
                    PetscReal rho_o = ctx->fluid.rhoO_surf / Bo;
                    PetscReal rho_w = ctx->fluid.rhoW_surf / Bw;
                    PetscReal rho_g = ctx->fluid.rhoG_surf / Bg;
                    
                    // Gravity
                    PetscReal g = 9.81;
                    PetscReal dz_gravity = 0.0;
                    if (d == 4) dz_gravity = 0.5 * (dz + dz_n);  // Up
                    if (d == 5) dz_gravity = -0.5 * (dz + dz_n); // Down
                    
                    // Potential differences
                    PetscReal dpot_o = (p - pn) - rho_o * g * dz_gravity;
                    PetscReal dpot_w = (p - pn) - rho_w * g * dz_gravity;
                    PetscReal dpot_g = (p - pn) - rho_g * g * dz_gravity;
                    
                    // Upstream weighting
                    PetscReal krw, kro, krg;
                    PetscReal krwn, kron, krgn;
                    
                    getRelPerm(ctx->fluid, Sw, Sg, krw, kro, krg);
                    getRelPerm(ctx->fluid, Swn, Sgn, krwn, kron, krgn);
                    
                    PetscReal krw_up = (dpot_w > 0) ? krw : krwn;
                    PetscReal kro_up = (dpot_o > 0) ? kro : kron;
                    PetscReal krg_up = (dpot_g > 0) ? krg : krgn;
                    
                    PetscReal Bo_up = (dpot_o > 0) ? Bo : Bon;
                    PetscReal Bw_up = (dpot_w > 0) ? Bw : Bwn;
                    PetscReal Bg_up = (dpot_g > 0) ? Bg : Bgn;
                    
                    PetscReal Rs_up = (dpot_o > 0) ? Rs : Rsn;
                    
                    // Phase fluxes
                    PetscReal qo = T * kro_up / muo * dpot_o / Bo_up;
                    PetscReal qw = T * krw_up / muw * dpot_w / Bw_up;
                    PetscReal qg = T * krg_up / mug * dpot_g / Bg_up + Rs_up * qo;
                    
                    // Add to residual
                    flocal[k][j][i].p += qo + qw + qg;
                    flocal[k][j][i].Sw += qw;
                    flocal[k][j][i].Sg += qg;
                }
                
                // Well contributions
                for (const auto& well : ctx->wells) {
                    if (!well.isOpen) continue;
                    if (well.i == i && well.j == j && k >= well.k1 && k <= well.k2) {
                        PetscReal WI = calculateWellIndex(well, ctx->grid, i, j, k);
                        
                        // Get relative permeabilities
                        PetscReal krw, kro, krg;
                        getRelPerm(ctx->fluid, Sw, Sg, krw, kro, krg);
                        
                        if (well.type == 0) { // Producer
                            PetscReal qo, qw, qg;
                            
                            if (well.control == 1) { // Rate control
                                // Calculate total mobility
                                PetscReal lambda_t = kro / muo + krw / muw + krg / mug;
                                if (lambda_t > 1e-12) {
                                    // Distribute rate according to mobility
                                    PetscReal q_total = well.rate; // Negative for production
                                    qo = q_total * (kro / muo) / lambda_t;
                                    qw = q_total * (krw / muw) / lambda_t;
                                    qg = q_total * (krg / mug) / lambda_t;
                                } else {
                                    qo = qw = qg = 0.0;
                                }
                            } else { // BHP control
                                PetscReal dp = p - well.BHP;
                                qo = WI * kro / muo * dp / Bo;
                                qw = WI * krw / muw * dp / Bw;
                                qg = WI * krg / mug * dp / Bg;
                            }
                            
                            // Add dissolved gas
                            qg += Rs * qo;
                            
                            flocal[k][j][i].p += qo + qw + qg;
                            flocal[k][j][i].Sw += qw;
                            flocal[k][j][i].Sg += qg;
                            
                        } else { // Injector
                            PetscReal qi;
                            if (well.control == 1) { // Rate control
                                qi = well.rate; // Positive for injection
                            } else { // BHP control
                                PetscReal dp = well.BHP - p;
                                if (well.phase == "WATER") {
                                    qi = WI / muw * dp / Bw;
                                } else if (well.phase == "GAS") {
                                    qi = WI / mug * dp / Bg;
                                } else {
                                    qi = 0.0;
                                }
                            }
                            
                            if (well.phase == "WATER") {
                                flocal[k][j][i].p -= qi;
                                flocal[k][j][i].Sw -= qi;
                            } else if (well.phase == "GAS") {
                                flocal[k][j][i].p -= qi;
                                flocal[k][j][i].Sg -= qi;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Restore arrays
    DMDAVecRestoreArrayRead(da, x_local, &xlocal);
    DMDAVecRestoreArray(da, f_local, &flocal);
    DMDAVecRestoreArrayRead(da, xold_local, &xold);
    
    // Local to global
    DMLocalToGlobalBegin(da, f_local, INSERT_VALUES, f);
    DMLocalToGlobalEnd(da, f_local, INSERT_VALUES, f);
    
    DMRestoreLocalVector(da, &x_local);
    DMRestoreLocalVector(da, &f_local);
    DMRestoreLocalVector(da, &xold_local);
    
    return 0;
}

// ==================== JACOBIAN FUNCTION ====================

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat J, Mat P, void *ctx_ptr) {
    AppCtx *ctx = (AppCtx*)ctx_ptr;
    
    // Clear matrix
    MatZeroEntries(J);
    
    // For now, use finite differences
    SNESComputeJacobianDefault(snes, x, J, P, ctx);
    
    return 0;
}

// ==================== DECK PARSER IMPLEMENTATION ====================

namespace BlackOil {

// DeckParser constructor
DeckParser::DeckParser(AppCtx* context) : ctx(context), file(nullptr) {}

// DeckParser method implementations
bool DeckParser::parseDeck(const std::string& filename) {
    deckFilename = filename;  // Store the filename for INCLUDE file resolution
    std::ifstream deckFile(filename);
    if (!deckFile.is_open()) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Cannot open deck file %s\n", filename.c_str());
        return false;
    }
    
    file = &deckFile;
    std::string line;
    bool foundRunspec = false, foundGrid = false, foundProps = false;
    bool foundSolution = false, foundSchedule = false;
    
    while (getDataLine(*file, line)) {
        std::string keyword = getKeyword(line);
        
        if (keyword == "RUNSPEC") {
            currentSection = "RUNSPEC";
            foundRunspec = true;
            if (!parseRunspec()) return false;
        } else if (keyword == "GRID") {
            currentSection = "GRID";
            foundGrid = true;
            if (!parseGrid()) return false;
        } else if (keyword == "PROPS") {
            currentSection = "PROPS";
            foundProps = true;
            if (!parseProps()) return false;
        } else if (keyword == "SOLUTION") {
            currentSection = "SOLUTION";
            foundSolution = true;
            if (!parseSolution()) return false;
        } else if (keyword == "SCHEDULE") {
            currentSection = "SCHEDULE";
            foundSchedule = true;
            if (!parseSchedule()) return false;
        }
    }
    
    file = nullptr;
    deckFile.close();
    
    PetscPrintf(PETSC_COMM_WORLD, "\nSections found: RUNSPEC=%d, GRID=%d, PROPS=%d, SOLUTION=%d, SCHEDULE=%d\n",
               foundRunspec, foundGrid, foundProps, foundSolution, foundSchedule);
    
    // Validate data
    if (ctx->grid.Nx <= 0 || ctx->grid.Ny <= 0 || ctx->grid.Nz <= 0) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Invalid grid dimensions\n");
        return false;
    }
    
    // Set defaults for missing data
    PetscInt ncells = ctx->grid.Nx * ctx->grid.Ny * ctx->grid.Nz;
    
    PetscPrintf(PETSC_COMM_WORLD, "\nData Summary:\n");
    PetscPrintf(PETSC_COMM_WORLD, "  Grid: %d x %d x %d = %d cells\n", 
               ctx->grid.Nx, ctx->grid.Ny, ctx->grid.Nz, ncells);
    PetscPrintf(PETSC_COMM_WORLD, "  DX: %d values\n", static_cast<PetscInt>(ctx->grid.DX.size()));
    PetscPrintf(PETSC_COMM_WORLD, "  DY: %d values\n", static_cast<PetscInt>(ctx->grid.DY.size()));
    PetscPrintf(PETSC_COMM_WORLD, "  DZ: %d values\n", static_cast<PetscInt>(ctx->grid.DZ.size()));
    PetscPrintf(PETSC_COMM_WORLD, "  PORO: %d values\n", static_cast<PetscInt>(ctx->grid.porosity.size()));
    PetscPrintf(PETSC_COMM_WORLD, "  PERMX: %d values\n", static_cast<PetscInt>(ctx->grid.permX.size()));
    PetscPrintf(PETSC_COMM_WORLD, "  TOPS: %d values\n", static_cast<PetscInt>(ctx->grid.TOPS.size()));
    PetscPrintf(PETSC_COMM_WORLD, "  Wells: %d\n", static_cast<PetscInt>(ctx->wells.size()));
    
    // Default cell dimensions if not specified
    if (ctx->grid.DX.empty()) {
        PetscPrintf(PETSC_COMM_WORLD, "Info: Using default cell size 50m x 50m x 10m\n");
        ctx->grid.DX.resize(ncells, 50.0);
        ctx->grid.DY.resize(ncells, 50.0);
        ctx->grid.DZ.resize(ncells, 10.0);
    }
    
    // Default porosity if not specified
    if (ctx->grid.porosity.empty()) {
        PetscPrintf(PETSC_COMM_WORLD, "Info: Using default porosity 0.2\n");
        ctx->grid.porosity.resize(ncells, 0.2);
    }
    
    // Default permeability if not specified
    if (ctx->grid.permX.empty()) {
        PetscPrintf(PETSC_COMM_WORLD, "Info: Using default permeability 100mD\n");
        ctx->grid.permX.resize(ncells, 100.0 * 9.869233e-16);
        ctx->grid.permY = ctx->grid.permX;
        ctx->grid.permZ.resize(ncells, 10.0 * 9.869233e-16);
    }
    
    // Default fluid properties if not specified
    if (ctx->fluid.rhoO_surf <= 0) {
        ctx->fluid.rhoO_surf = 800.0;
        ctx->fluid.rhoW_surf = 1000.0;
        ctx->fluid.rhoG_surf = 1.0;
    }
    
    return true;
}

std::string DeckParser::getKeyword(const std::string& line) {
    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;
    std::string upper = toUpper(keyword);
    
    // Only log major section keywords
    if (upper == "RUNSPEC" || upper == "GRID" || upper == "PROPS" || 
        upper == "SOLUTION" || upper == "SCHEDULE" || upper == "END") {
        PetscPrintf(PETSC_COMM_WORLD, "Found section keyword: %s\n", upper.c_str());
    }
    
    return upper;
}

bool DeckParser::parseRunspec() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        std::string keyword = getKeyword(line);
        
        if (keyword == "TITLE") {
            getDataLine(*file, line);
            PetscPrintf(PETSC_COMM_WORLD, "Title: %s\n", line.c_str());
            
        } else if (keyword == "DIMENS") {
            if (!parseDimens()) return false;
            
        } else if (keyword == "FIELD") {
            ctx->useField = true;
            PetscPrintf(PETSC_COMM_WORLD, "  Using FIELD units\n");
            
        } else if (keyword == "METRIC") {
            ctx->useField = false;
            PetscPrintf(PETSC_COMM_WORLD, "  Using METRIC units\n");
            
        } else if (keyword == "OIL" || keyword == "WATER" || keyword == "GAS" || 
                  keyword == "DISGAS" || keyword == "VAPOIL") {
            // Phase keywords - noted but not used directly
            
        } else if (keyword == "WELLDIMS" || keyword == "TABDIMS" || 
                  keyword == "EQLDIMS" || keyword == "START") {
            // Skip these for now
            
        } else if (keyword == "GRID" || keyword == "PROPS" || 
                  keyword == "SOLUTION" || keyword == "SCHEDULE") {
            // End of RUNSPEC - rewind to allow main parser to see this keyword
            file->seekg(-static_cast<std::streamoff>(line.length())-1, std::ios_base::cur);
            break;
        }
    }
    
    return true;
}

bool DeckParser::parseDimens() {
    std::string line;
    if (!getDataLine(*file, line)) return false;
    
    std::vector<PetscReal> dims = parseNumbers(line);
    if (dims.size() >= 3) {
        ctx->grid.Nx = (PetscInt)dims[0];
        ctx->grid.Ny = (PetscInt)dims[1];
        ctx->grid.Nz = (PetscInt)dims[2];
        PetscPrintf(PETSC_COMM_WORLD, "Grid dimensions: %d x %d x %d\n", 
                   ctx->grid.Nx, ctx->grid.Ny, ctx->grid.Nz);
        return true;
    }
    return false;
}

bool DeckParser::parseGrid() {
    std::string line;
    PetscInt ncells = ctx->grid.Nx * ctx->grid.Ny * ctx->grid.Nz;
    
    PetscPrintf(PETSC_COMM_WORLD, "Parsing GRID section (expecting %d cells)...\n", ncells);
    
    // Track what we've found
    bool hasData = false;
    
    while (getDataLine(*file, line)) {
        std::string keyword = getKeyword(line);
        
        // Special keywords to skip
        if (keyword == "NOECHO" || keyword == "ECHO") {
            continue;
        }
        
        // Check for end of section - all these keywords mark the end of GRID
        if (keyword == "PROPS" || keyword == "SOLUTION" || 
            keyword == "SCHEDULE" || keyword == "RUNSPEC" || 
            keyword == "END" || keyword == "EDIT") {
            file->seekg(-static_cast<std::streamoff>(line.length())-1, std::ios_base::cur);
            break;
        }
        
        // Handle INCLUDE files
        if (keyword == "INCLUDE") {
            if (getDataLine(*file, line)) {
                // Extract filename from quotes
                std::string includeFile = trim(line, "'\"/ ");
                PetscPrintf(PETSC_COMM_WORLD, "  Found INCLUDE: %s\n", includeFile.c_str());
                
                // For SPE9, the grid data is often in an include file
                // Try to extract directory from the main deck filename
                std::string deckDir = "";
                size_t lastSlash = deckFilename.find_last_of("/\\");
                if (lastSlash != std::string::npos) {
                    deckDir = deckFilename.substr(0, lastSlash + 1);
                }
                
                // For now, skip INCLUDE files and use defaults
                PetscPrintf(PETSC_COMM_WORLD, "  Note: INCLUDE file parsing not implemented, using defaults\n");
                continue;
            }
        }
        
        // Handle data keywords
        if (keyword == "DX" || keyword == "DY" || keyword == "DZ") {
            hasData = true;
            std::vector<PetscReal> data = readArrayData(ncells);
            
            if (ctx->useField && !data.empty()) {
                // Convert feet to meters
                for (auto& val : data) val *= ft_to_m;
            }
            
            if (keyword == "DX") {
                ctx->grid.DX = data;
                PetscPrintf(PETSC_COMM_WORLD, "  Loaded DX: %d values\n", static_cast<PetscInt>(data.size()));
            } else if (keyword == "DY") {
                ctx->grid.DY = data;
                PetscPrintf(PETSC_COMM_WORLD, "  Loaded DY: %d values\n", static_cast<PetscInt>(data.size()));
            } else if (keyword == "DZ") {
                ctx->grid.DZ = data;
                PetscPrintf(PETSC_COMM_WORLD, "  Loaded DZ: %d values (first=%g)\n", 
                           static_cast<PetscInt>(data.size()), data.empty() ? 0.0 : data[0]);
            }
            
        } else if (keyword == "TOPS") {
            hasData = true;
            ctx->grid.TOPS = readArrayData(ctx->grid.Nx * ctx->grid.Ny);
            if (ctx->useField && !ctx->grid.TOPS.empty()) {
                for (auto& val : ctx->grid.TOPS) val *= ft_to_m;
            }
            PetscPrintf(PETSC_COMM_WORLD, "  Loaded TOPS: %d values\n", static_cast<PetscInt>(ctx->grid.TOPS.size()));
            
        } else if (keyword == "PORO") {
            hasData = true;
            ctx->grid.porosity = readArrayData(ncells);
            PetscPrintf(PETSC_COMM_WORLD, "  Loaded PORO: %d values\n", static_cast<PetscInt>(ctx->grid.porosity.size()));
            
        } else if (keyword == "PERMX") {
            hasData = true;
            ctx->grid.permX = readArrayData(ncells);
            if (!ctx->grid.permX.empty()) {
                for (auto& val : ctx->grid.permX) val *= mD_to_m2;
            }
            PetscPrintf(PETSC_COMM_WORLD, "  Loaded PERMX: %d values\n", static_cast<PetscInt>(ctx->grid.permX.size()));
            
        } else if (keyword == "PERMY") {
            hasData = true;
            ctx->grid.permY = readArrayData(ncells);
            if (!ctx->grid.permY.empty()) {
                for (auto& val : ctx->grid.permY) val *= mD_to_m2;
            }
            PetscPrintf(PETSC_COMM_WORLD, "  Loaded PERMY: %d values\n", static_cast<PetscInt>(ctx->grid.permY.size()));
            
        } else if (keyword == "PERMZ") {
            hasData = true;
            ctx->grid.permZ = readArrayData(ncells);
            if (!ctx->grid.permZ.empty()) {
                for (auto& val : ctx->grid.permZ) val *= mD_to_m2;
            }
            PetscPrintf(PETSC_COMM_WORLD, "  Loaded PERMZ: %d values\n", static_cast<PetscInt>(ctx->grid.permZ.size()));
        }
    }
    
    // Set defaults if not specified
    if (ctx->grid.permY.empty() && !ctx->grid.permX.empty()) {
        ctx->grid.permY = ctx->grid.permX;
    }
    if (ctx->grid.permZ.empty() && !ctx->grid.permX.empty()) {
        ctx->grid.permZ = ctx->grid.permX;
        for (auto& val : ctx->grid.permZ) val *= 0.1;  // Default Kv/Kh = 0.1
    }
    
    // For SPE9, if no DX/DY/DZ found, use standard SPE9 dimensions
    if (ctx->grid.DX.empty() && ctx->useField) {
        // SPE9 uses 300ft x 300ft cells
        PetscReal cellSizeXY = 300.0 * ft_to_m;
        ctx->grid.DX.resize(ncells, cellSizeXY);
        ctx->grid.DY.resize(ncells, cellSizeXY);
        PetscPrintf(PETSC_COMM_WORLD, "  Using SPE9 default DX/DY: 300 ft = %g m\n", cellSizeXY);
        
        // SPE9 has variable DZ - simplified here
        ctx->grid.DZ.resize(ncells, 25.0 * ft_to_m);
        PetscPrintf(PETSC_COMM_WORLD, "  Using default DZ: 25 ft = %g m\n", 25.0 * ft_to_m);
    }
    
    return true;
}

bool DeckParser::parseProps() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        std::string keyword = getKeyword(line);
        
        if (keyword == "PVTW") {
            if (!parsePVTW()) return false;
            
        } else if (keyword == "ROCK") {
            if (!parseROCK()) return false;
            
        } else if (keyword == "DENSITY") {
            if (!parseDENSITY()) return false;
            
        } else if (keyword == "PVTO") {
            if (!parsePVTO()) return false;
            
        } else if (keyword == "PVDG") {
            if (!parsePVDG()) return false;
            
        } else if (keyword == "SGOF") {
            if (!parseSGOF()) return false;
            
        } else if (keyword == "SWOF") {
            if (!parseSWOF()) return false;
            
        } else if (keyword == "SOLUTION" || keyword == "SCHEDULE" || 
                  keyword == "RUNSPEC" || keyword == "GRID" || 
                  keyword == "REGIONS" || keyword == "END") {
            file->seekg(-static_cast<std::streamoff>(line.length())-1, std::ios_base::cur);
            break;
        }
    }
    
    return true;
}

bool DeckParser::parsePVTW() {
    std::string line;
    if (!getDataLine(*file, line)) return false;
    
    std::vector<PetscReal> data = parseNumbers(line);
    if (data.size() >= 5) {
        ctx->fluid.pwref = data[0];
        ctx->fluid.Bwref = data[1];
        ctx->fluid.cw = data[2];
        ctx->fluid.muw = data[3];
        ctx->fluid.muw_p = data[4];
        
        // Convert units if field units
        if (ctx->useField) {
            ctx->fluid.pwref *= psi_to_Pa;
            ctx->fluid.cw /= psi_to_Pa;
            ctx->fluid.muw *= cp_to_Pas;
            ctx->fluid.muw_p /= psi_to_Pa;
        }
        
        PetscPrintf(PETSC_COMM_WORLD, "  Parsed PVTW\n");
        return true;
    }
    return false;
}

bool DeckParser::parseROCK() {
    std::string line;
    if (!getDataLine(*file, line)) return false;
    
    std::vector<PetscReal> data = parseNumbers(line);
    if (data.size() >= 2) {
        ctx->rock.pref = data[0];
        ctx->rock.crock = data[1];
        
        if (ctx->useField) {
            ctx->rock.pref *= psi_to_Pa;
            ctx->rock.crock /= psi_to_Pa;
        }
        
        PetscPrintf(PETSC_COMM_WORLD, "  Parsed ROCK\n");
        return true;
    }
    return false;
}

bool DeckParser::parseDENSITY() {
    std::string line;
    if (!getDataLine(*file, line)) return false;
    
    std::vector<PetscReal> data = parseNumbers(line);
    if (data.size() >= 3) {
        ctx->fluid.rhoO_surf = data[0];
        ctx->fluid.rhoW_surf = data[1];
        ctx->fluid.rhoG_surf = data[2];
        
        if (ctx->useField) {
            ctx->fluid.rhoO_surf *= lbft3_to_kgm3;
            ctx->fluid.rhoW_surf *= lbft3_to_kgm3;
            ctx->fluid.rhoG_surf *= lbft3_to_kgm3;
        }
        
        PetscPrintf(PETSC_COMM_WORLD, "  Parsed DENSITY\n");
        return true;
    }
    return false;
}

bool DeckParser::parsePVTO() {
    std::string line;
    ctx->fluid.pvto.Rs.clear();
    ctx->fluid.pvto.pbub.clear();
    ctx->fluid.pvto.Bo.clear();
    ctx->fluid.pvto.muo.clear();
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;  // End of table
        
        std::vector<PetscReal> data = parseNumbers(line);
        if (data.size() >= 4) {
            // First line of each Rs block
            PetscReal Rs = data[0];
            PetscReal pbub = data[1];
            PetscReal Bo = data[2];
            PetscReal muo = data[3];
            
            ctx->fluid.pvto.Rs.push_back(Rs);
            ctx->fluid.pvto.pbub.push_back(pbub);
            ctx->fluid.pvto.Bo.push_back(Bo);
            ctx->fluid.pvto.muo.push_back(muo);
            
            // Read undersaturated data
            std::vector<PetscReal> p_usat, Bo_usat, muo_usat;
            p_usat.push_back(pbub);
            Bo_usat.push_back(Bo);
            muo_usat.push_back(muo);
            
            while (getDataLine(*file, line)) {
                if (line.find("/") != std::string::npos) {
                    file->seekg(-static_cast<std::streamoff>(line.length())-1, std::ios_base::cur);
                    break;
                }
                
                data = parseNumbers(line);
                if (data.size() >= 3) {
                    p_usat.push_back(data[0]);
                    Bo_usat.push_back(data[1]);
                    muo_usat.push_back(data[2]);
                } else {
                    file->seekg(-static_cast<std::streamoff>(line.length())-1, std::ios_base::cur);
                    break;
                }
            }
            
            ctx->fluid.pvto.p_usat.push_back(p_usat);
            ctx->fluid.pvto.Bo_usat.push_back(Bo_usat);
            ctx->fluid.pvto.muo_usat.push_back(muo_usat);
        }
    }
    
    // Convert units
    if (ctx->useField) {
        for (auto& Rs : ctx->fluid.pvto.Rs) Rs *= Mscf_to_sm3 / stb_to_m3;
        for (auto& p : ctx->fluid.pvto.pbub) p *= psi_to_Pa;
        for (auto& Bo : ctx->fluid.pvto.Bo) Bo *= rb_to_m3 / stb_to_m3;
        for (auto& mu : ctx->fluid.pvto.muo) mu *= cp_to_Pas;
        
        for (auto& p_vec : ctx->fluid.pvto.p_usat) {
            for (auto& p : p_vec) p *= psi_to_Pa;
        }
        for (auto& Bo_vec : ctx->fluid.pvto.Bo_usat) {
            for (auto& Bo : Bo_vec) Bo *= rb_to_m3 / stb_to_m3;
        }
        for (auto& mu_vec : ctx->fluid.pvto.muo_usat) {
            for (auto& mu : mu_vec) mu *= cp_to_Pas;
        }
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "  Parsed PVTO (%d Rs values)\n", static_cast<PetscInt>(ctx->fluid.pvto.Rs.size()));
    return true;
}

bool DeckParser::parsePVDG() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;
        
        std::vector<PetscReal> data = parseNumbers(line);
        if (data.size() >= 3) {
            ctx->fluid.pvdg_p.push_back(data[0]);
            ctx->fluid.pvdg_Bg.push_back(data[1]);
            ctx->fluid.pvdg_mug.push_back(data[2]);
        }
    }
    
    // Convert units
    if (ctx->useField) {
        for (auto& p : ctx->fluid.pvdg_p) p *= psi_to_Pa;
        for (auto& Bg : ctx->fluid.pvdg_Bg) Bg *= rb_to_m3 / Mscf_to_sm3;
        for (auto& mu : ctx->fluid.pvdg_mug) mu *= cp_to_Pas;
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "  Parsed PVDG (%d entries)\n", static_cast<PetscInt>(ctx->fluid.pvdg_p.size()));
    return true;
}

bool DeckParser::parseSGOF() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;
        
        std::vector<PetscReal> data = parseNumbers(line);
        if (data.size() >= 4) {
            ctx->fluid.sgof_Sg.push_back(data[0]);
            ctx->fluid.sgof_krg.push_back(data[1]);
            ctx->fluid.sgof_krog.push_back(data[2]);
            ctx->fluid.sgof_pcog.push_back(data[3]);
        }
    }
    
    // Convert units
    if (ctx->useField) {
        for (auto& pc : ctx->fluid.sgof_pcog) pc *= psi_to_Pa;
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "  Parsed SGOF (%d entries)\n", static_cast<PetscInt>(ctx->fluid.sgof_Sg.size()));
    return true;
}

bool DeckParser::parseSWOF() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;
        
        std::vector<PetscReal> data = parseNumbers(line);
        if (data.size() >= 4) {
            ctx->fluid.swof_Sw.push_back(data[0]);
            ctx->fluid.swof_krw.push_back(data[1]);
            ctx->fluid.swof_krow.push_back(data[2]);
            ctx->fluid.swof_pcow.push_back(data[3]);
        }
    }
    
    // Convert units
    if (ctx->useField) {
        for (auto& pc : ctx->fluid.swof_pcow) pc *= psi_to_Pa;
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "  Parsed SWOF (%d entries)\n", static_cast<PetscInt>(ctx->fluid.swof_Sw.size()));
    return true;
}

bool DeckParser::parseSolution() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        std::string keyword = getKeyword(line);
        
        if (keyword == "EQUIL") {
            if (!parseEQUIL()) return false;
            
        } else if (keyword == "RSVD") {
            if (!parseRSVD()) return false;
            
        } else if (keyword == "SCHEDULE" || keyword == "RUNSPEC" || 
                  keyword == "GRID" || keyword == "PROPS" || 
                  keyword == "SUMMARY" || keyword == "END") {
            file->seekg(-static_cast<std::streamoff>(line.length())-1, std::ios_base::cur);
            break;
        }
    }
    
    return true;
}

bool DeckParser::parseEQUIL() {
    std::string line;
    if (!getDataLine(*file, line)) return false;
    
    std::vector<PetscReal> data = parseNumbers(line);
    if (data.size() >= 6) {
        ctx->equil.datum = data[0];
        ctx->equil.pres_datum = data[1];
        ctx->equil.zwoc = data[2];
        ctx->equil.pcow_woc = data[3];
        ctx->equil.zgoc = data[4];
        ctx->equil.pcog_goc = data[5];
        
        // Convert units
        if (ctx->useField) {
            ctx->equil.datum *= ft_to_m;
            ctx->equil.pres_datum *= psi_to_Pa;
            ctx->equil.zwoc *= ft_to_m;
            ctx->equil.pcow_woc *= psi_to_Pa;
            ctx->equil.zgoc *= ft_to_m;
            ctx->equil.pcog_goc *= psi_to_Pa;
        }
        
        PetscPrintf(PETSC_COMM_WORLD, "  Parsed EQUIL: datum=%g m, p=%g Pa\n", 
                   ctx->equil.datum, ctx->equil.pres_datum);
        return true;
    }
    return false;
}

bool DeckParser::parseRSVD() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;
        
        std::vector<PetscReal> data = parseNumbers(line);
        if (data.size() >= 2) {
            ctx->equil.rsvd_depth.push_back(data[0]);
            ctx->equil.rsvd_rs.push_back(data[1]);
        }
    }
    
    // Convert units
    if (ctx->useField) {
        for (auto& d : ctx->equil.rsvd_depth) d *= ft_to_m;
        for (auto& rs : ctx->equil.rsvd_rs) rs *= Mscf_to_sm3 / stb_to_m3;
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "  Parsed RSVD (%d entries)\n", static_cast<PetscInt>(ctx->equil.rsvd_depth.size()));
    return true;
}

bool DeckParser::parseSchedule() {
    std::string line;
    
    PetscPrintf(PETSC_COMM_WORLD, "Parsing SCHEDULE section...\n");
    
    while (getDataLine(*file, line)) {
        std::string keyword = getKeyword(line);
        
        // Skip report keywords
        if (keyword == "RPTRST" || keyword == "RPTSCHED") {
            // Skip until we find a slash
            while (getDataLine(*file, line)) {
                if (line.find("/") != std::string::npos) break;
            }
            continue;
        }
        
        if (keyword == "WELSPECS") {
            PetscPrintf(PETSC_COMM_WORLD, "  Parsing WELSPECS...\n");
            if (!parseWELSPECS()) return false;
            
        } else if (keyword == "COMPDAT") {
            PetscPrintf(PETSC_COMM_WORLD, "  Parsing COMPDAT...\n");
            if (!parseCOMPDAT()) return false;
            
        } else if (keyword == "WCONINJE") {
            PetscPrintf(PETSC_COMM_WORLD, "  Parsing WCONINJE...\n");
            if (!parseWCONINJE()) return false;
            
        } else if (keyword == "WCONPROD") {
            PetscPrintf(PETSC_COMM_WORLD, "  Parsing WCONPROD...\n");
            if (!parseWCONPROD()) return false;
            
        } else if (keyword == "TSTEP") {
            PetscPrintf(PETSC_COMM_WORLD, "  Parsing TSTEP...\n");
            if (!parseTSTEP()) return false;
            
        } else if (keyword == "END") {
            break;
        }
    }
    
    return true;
}

bool DeckParser::parseWELSPECS() {
    std::string line;
    int wellCount = 0;
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;
        
        std::istringstream iss(line);
        Well well;
        std::string wellName, group, preferredPhase;
        PetscReal drainageRadius = 0.5;  // Default
        
        // Read well name (might be quoted)
        if (!(iss >> wellName)) continue;
        well.name = trim(wellName, "'\"");
        
        // Read group (might be quoted)
        if (!(iss >> group)) continue;
        group = trim(group, "'\"");
        
        // Read I, J coordinates
        if (!(iss >> well.i >> well.j)) continue;
        well.i--;  // Convert to 0-based
        well.j--;
        
        // Read reference depth
        if (!(iss >> well.refDepth)) continue;
        
        // Read preferred phase (might be quoted) - this field is optional
        if (iss >> preferredPhase) {
            preferredPhase = trim(toUpper(preferredPhase), "'\"");
        } else {
            preferredPhase = "OIL";  // Default
        }
        
        // Try to read drainage radius
        iss >> drainageRadius;
        
        // Convert units
        if (ctx->useField) {
            well.refDepth *= ft_to_m;
            drainageRadius *= ft_to_m;
        }
        
        well.rw = 0.1;  // Default wellbore radius
        
        // Determine well type from name (common convention)
        if (well.name.find("INJ") != std::string::npos || 
            well.name.find("I") == 0) {
            well.type = 1;  // Injector
            well.phase = "WATER";
        } else {
            well.type = 0;  // Producer
            well.phase = "OIL";
        }
        
        well.isOpen = false;  // Will be set by COMPDAT
        well.control = 0;     // Default BHP control
        well.BHP = 1e7;       // Default 10 MPa
        well.rate = 0.0;
        
        ctx->wells.push_back(well);
        wellCount++;
        
        PetscPrintf(PETSC_COMM_WORLD, "    Added well %s at (%d,%d)\n", 
                   well.name.c_str(), well.i+1, well.j+1);
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "    Total wells defined: %d\n", wellCount);
    return true;
}

bool DeckParser::parseCOMPDAT() {
    std::string line;
    int completionCount = 0;
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;
        
        std::istringstream iss(line);
        std::string wellName, status;
        PetscInt i, j, k1, k2;
        PetscReal diameter = 0.5;
        
        // Read well name (might be quoted)
        if (!(iss >> wellName)) continue;
        wellName = trim(wellName, "'\"");
        
        // Read I, J, K1, K2
        if (!(iss >> i >> j >> k1 >> k2)) continue;
        
        // Read status (might be quoted)
        std::string token;
        if (iss >> token) {
            status = trim(toUpper(token), "'\"");
        } else {
            status = "OPEN";
        }
        
        // Skip next tokens and try to read diameter
        for (int skip = 0; skip < 2; skip++) {
            if (!(iss >> token)) break;
        }
        iss >> diameter;
        
        // Find matching wells
        for (auto& well : ctx->wells) {
            if (well.name == wellName || matchesWildcard(wellName, well.name)) {
                // Update location if specified (non-zero)
                if (i > 0) well.i = i - 1;  // Convert to 0-based
                if (j > 0) well.j = j - 1;
                if (k1 > 0) well.k1 = k1 - 1;  // Convert to 0-based
                if (k2 > 0) well.k2 = k2 - 1;  // Convert to 0-based
                
                // If k values are defaulted (0), complete in all layers
                if (k1 == 0 && k2 == 0) {
                    well.k1 = 0;
                    well.k2 = ctx->grid.Nz - 1;
                }
                
                // Debug output to verify layer assignment
                PetscPrintf(PETSC_COMM_WORLD, "    DEBUG: COMPDAT %s: i=%d, j=%d, k1=%d, k2=%d (0-based)\n",
                           well.name.c_str(), well.i, well.j, well.k1, well.k2);
                
                if (status == "OPEN") {
                    well.isOpen = true;
                }
                
                if (diameter > 0 && diameter < 10) {  // Reasonable diameter
                    well.rw = diameter / 2.0;
                    if (ctx->useField) well.rw *= ft_to_m;
                }
                
                completionCount++;
                PetscPrintf(PETSC_COMM_WORLD, "    Completed well %s at (%d,%d) layers %d-%d\n",
                           well.name.c_str(), well.i+1, well.j+1, well.k1+1, well.k2+1);
                break;
            }
        }
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "    Total completions: %d\n", completionCount);
    return true;
}

bool DeckParser::parseWCONINJE() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;
        
        std::istringstream iss(line);
        std::string wellName, phase, status, control;
        PetscReal rate = 0.0, bhp = 4000.0;  // Default BHP
        
        if (iss >> wellName >> phase >> status >> control) {
            wellName = trim(wellName, "'\"");
            phase = trim(toUpper(phase), "'\"");
            status = trim(toUpper(status), "'\"");
            control = trim(toUpper(control), "'\"");
            
            // Read rate
            iss >> rate;
            
            // Skip some fields and read BHP
            std::string token;
            for (int i = 0; i < 1; i++) {
                if (!(iss >> token)) break;
            }
            iss >> bhp;
            
            // Find matching wells
            for (auto& well : ctx->wells) {
                if (well.name == wellName || matchesWildcard(wellName, well.name)) {
                    well.type = 1;  // Injector
                    well.phase = phase;
                    
                    if (status == "SHUT") {
                        well.isOpen = false;
                    } else {
                        well.isOpen = true;
                    }
                    
                    if (control == "RATE") {
                        well.control = 1;  // Rate control
                        well.rate = rate;
                        if (ctx->useField) {
                            well.rate *= stb_to_m3 / day_to_sec;
                        }
                    } else {
                        well.control = 0;  // BHP control
                    }
                    
                    well.BHP = bhp;
                    if (ctx->useField) well.BHP *= psi_to_Pa;
                    
                    PetscPrintf(PETSC_COMM_WORLD, "    Set injector %s: %s injection, rate=%g, BHP=%g\n",
                               well.name.c_str(), phase.c_str(), rate, bhp);
                }
            }
        }
    }
    
    return true;
}

bool DeckParser::parseWCONPROD() {
    std::string line;
    
    while (getDataLine(*file, line)) {
        if (line.find("/") != std::string::npos) break;
        
        std::istringstream iss(line);
        std::string wellName, status, control;
        PetscReal orat = 0.0, wrat = 0.0, grat = 0.0, lrat = 0.0, resv = 0.0, bhp = 1000.0;
        
        if (iss >> wellName >> status >> control) {
            wellName = trim(wellName, "'\"");
            status = trim(toUpper(status), "'\"");
            control = trim(toUpper(control), "'\"");
            
            // Read rates
            iss >> orat >> wrat >> grat >> lrat >> resv >> bhp;
            
            // Find matching wells
            for (auto& well : ctx->wells) {
                if (well.name == wellName || matchesWildcard(wellName, well.name)) {
                    well.type = 0;  // Producer
                    
                    if (status == "SHUT") {
                        well.isOpen = false;
                    } else {
                        well.isOpen = true;
                    }
                    
                    if (control == "ORAT" && orat > 0) {
                        well.control = 1;  // Rate control
                        well.rate = -orat;  // Negative for production
                        if (ctx->useField) {
                            well.rate *= stb_to_m3 / day_to_sec;
                        }
                    } else {
                        well.control = 0;  // BHP control
                    }
                    
                    if (bhp > 0) {
                        well.BHP = bhp;
                        if (ctx->useField) well.BHP *= psi_to_Pa;
                    }
                    
                    PetscPrintf(PETSC_COMM_WORLD, "    Set producer %s: orat=%g, BHP=%g\n",
                               well.name.c_str(), orat, bhp);
                }
            }
        }
    }
    
    return true;
}

bool DeckParser::parseTSTEP() {
    std::string line;
    if (!getDataLine(*file, line)) return false;
    
    std::vector<PetscReal> steps = parseNumbers(line);
    for (auto dt : steps) {
        if (ctx->useField) dt *= day_to_sec;
        ctx->timesteps.push_back(dt);
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "    Added %d timesteps\n", static_cast<PetscInt>(steps.size()));
    return true;
}

std::vector<PetscReal> DeckParser::readArrayData(PetscInt expectedSize) {
    std::vector<PetscReal> data;
    std::string line;
    bool foundSlash = false;
    
    while (data.size() < expectedSize && getDataLine(*file, line)) {
        // Check if line contains a terminating slash
        if (line.find("/") != std::string::npos) {
            foundSlash = true;
            // Remove everything after the slash
            size_t slashPos = line.find("/");
            if (slashPos > 0) {
                line = line.substr(0, slashPos);
            } else {
                break;  // Line starts with slash, no more data
            }
        }
        
        // Parse numbers from this line
        std::vector<PetscReal> lineData = parseNumbers(line);
        data.insert(data.end(), lineData.begin(), lineData.end());
        
        // If we found a slash, stop reading
        if (foundSlash) {
            break;
        }
    }
    
    // Log what we read
    if (data.size() != expectedSize && expectedSize > 0) {
        if (data.empty()) {
            PetscPrintf(PETSC_COMM_WORLD, "    Warning: No data found (expected %d values)\n", 
                       expectedSize);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, "    Info: Read %d values (expected %d)\n", 
                       static_cast<PetscInt>(data.size()), expectedSize);
        }
    }
    
    // Resize to expected size
    if (data.size() > expectedSize) {
        data.resize(expectedSize);
    } else if (data.size() < expectedSize && !data.empty()) {
        // Fill remaining with last value
        PetscReal lastVal = data.back();
        PetscInt remaining = expectedSize - data.size();
        PetscPrintf(PETSC_COMM_WORLD, "    Info: Filling %d remaining values with %g\n", 
                   remaining, lastVal);
        while (data.size() < expectedSize) {
            data.push_back(lastVal);
        }
    }
    
    return data;
}

bool DeckParser::matchesWildcard(const std::string& pattern, const std::string& str) {
    if (pattern.find('*') == std::string::npos) {
        return pattern == str;
    }
    
    // Simple wildcard matching
    std::string prefix = pattern.substr(0, pattern.find('*'));
    return str.substr(0, prefix.length()) == prefix;
}

// ==================== SIMULATOR CLASS ====================

BlackOilSimulator::BlackOilSimulator() : initialized(false) {
    PetscMemzero(&ctx, sizeof(AppCtx));
    
    // Initialize default values
    ctx.fluid.Bwref = 1.0;
    ctx.fluid.cw = 3e-10;
    ctx.fluid.muw = 5e-4;
    ctx.fluid.muw_p = 0.0;
    ctx.rock.crock = 4e-10;
    ctx.rock.pref = 1e7;
}

BlackOilSimulator::~BlackOilSimulator() {
    if (initialized) {
        VecDestroy(&ctx.solution);
        VecDestroy(&ctx.solution_old);
        MatDestroy(&ctx.J);
        SNESDestroy(&ctx.snes);
        DMDestroy(&ctx.da);
    }
}

bool BlackOilSimulator::loadDeck(const std::string& filename) {
    DeckParser parser(&ctx);
    return parser.parseDeck(filename);
}

bool BlackOilSimulator::loadEclipseDeck(const std::string& filename) {
    // This method uses the Eclipse input reader instead of the built-in parser
    EclipseInputReader reader(filename);
    
    if (!reader.parse()) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Failed to parse Eclipse deck\n");
        return false;
    }
    
    // Convert Eclipse data to our internal format
    const auto& runspec = reader.get_runspec();
    const auto& grid_data = reader.get_grid_data();
    const auto& props_data = reader.get_props_data();
    const auto& solution_data = reader.get_solution_data();
    const auto& schedule_data = reader.get_schedule_data();
    
    // Set grid dimensions
    ctx.grid.Nx = runspec.nx;
    ctx.grid.Ny = runspec.ny;
    ctx.grid.Nz = runspec.nz;
    ctx.useField = (runspec.unit_system == "FIELD");
    
    // Set grid data
    ctx.grid.DX = grid_data.dx;
    ctx.grid.DY = grid_data.dy;
    ctx.grid.DZ = grid_data.dz;
    ctx.grid.TOPS = grid_data.tops;
    ctx.grid.porosity = grid_data.poro;
    ctx.grid.permX = grid_data.permx;
    ctx.grid.permY = grid_data.permy;
    ctx.grid.permZ = grid_data.permz;
    
    // Set initial saturations from Eclipse deck
    ctx.grid.swat = solution_data.swat;
    ctx.grid.sgas = solution_data.sgas;
    
    // Convert permeability from mD to m²
    for (auto& k : ctx.grid.permX) k *= 9.869233e-16;
    for (auto& k : ctx.grid.permY) k *= 9.869233e-16;
    for (auto& k : ctx.grid.permZ) k *= 9.869233e-16;
    
    // Set fluid properties
    if (!props_data.pvtw.empty()) {
        const auto& pvtw = props_data.pvtw[0];
        ctx.fluid.pwref = pvtw.pressure[0];
        ctx.fluid.Bwref = pvtw.fvf[0];
        ctx.fluid.cw = pvtw.compressibility;
        ctx.fluid.muw = pvtw.viscosity[0];
        ctx.fluid.muw_p = pvtw.viscosibility;
    }
    
    if (!props_data.rock.empty()) {
        ctx.rock.pref = props_data.rock[0].pref;
        ctx.rock.crock = props_data.rock[0].compressibility;
    }
    
    ctx.fluid.rhoO_surf = props_data.oil_density;
    ctx.fluid.rhoW_surf = props_data.water_density;
    ctx.fluid.rhoG_surf = props_data.gas_density;
    
    // Set equilibration data
    if (!solution_data.equil.empty()) {
        const auto& equil = solution_data.equil[0];
        ctx.equil.datum = equil.datum_depth;
        ctx.equil.pres_datum = equil.datum_pressure;
        ctx.equil.zwoc = equil.owc_depth;
        ctx.equil.pcow_woc = equil.pc_owc;
        ctx.equil.zgoc = equil.goc_depth;
        ctx.equil.pcog_goc = equil.pc_goc;
    }
    
    // Set wells
    if (!schedule_data.timesteps.empty() && !schedule_data.timesteps[0].wells.empty()) {
        for (const auto& well_data : schedule_data.timesteps[0].wells) {
            Well well;
            well.name = well_data.name;
            well.i = well_data.i;
            well.j = well_data.j;
            well.k1 = well_data.k1;
            well.k2 = well_data.k2;
            
            // Debug output for well conversion
            PetscPrintf(PETSC_COMM_WORLD, "DEBUG: Converting well %s: k1=%d, k2=%d (from Eclipse data)\n",
                       well.name.c_str(), well.k1, well.k2);
            well.refDepth = well_data.ref_depth;
            well.rw = well_data.radius;
            well.isOpen = (well_data.status == "OPEN");
            well.type = well_data.is_producer ? 0 : 1;
            well.phase = well_data.is_producer ? "OIL" : well_data.injection_type;
            
            if (well_data.control_mode == "RATE") {
                well.control = 1;
                if (well_data.is_producer) {
                    // For producers, use oil rate if available, otherwise use liquid rate
                    well.rate = -std::abs(well_data.oil_rate > 0 ? well_data.oil_rate : well_data.liquid_rate);
                } else {
                    // For injectors, use water rate
                    well.rate = well_data.water_rate;
                }
                // Convert from STB/day to m³/s (STB = 0.1589873 m³)
                well.rate = well.rate * 0.1589873 / 86400.0;
            } else {
                well.control = 0;
                well.BHP = well_data.bhp_limit;
            }
            
            ctx.wells.push_back(well);
        }
    }
    
    // Debug: Print wells after loading
    PetscPrintf(PETSC_COMM_WORLD, "DEBUG: Loaded %d wells\n", static_cast<PetscInt>(ctx.wells.size()));
    for (const auto& well : ctx.wells) {
        PetscPrintf(PETSC_COMM_WORLD, "DEBUG: Well %s: type=%d, control=%d, rate=%.2e, BHP=%.2e\n",
                   well.name.c_str(), well.type, well.control, well.rate, well.BHP);
    }
    
    // Set timesteps - each timestep.days is the individual timestep size
    ctx.timesteps.clear();
    for (const auto& ts : schedule_data.timesteps) {
        PetscReal dt = ts.days * 86400.0;  // Convert days to seconds
        if (dt > 0.0) {
            ctx.timesteps.push_back(dt);
        }
    }
    
    // If no timesteps were set, use a default
    if (ctx.timesteps.empty()) {
        ctx.timesteps.push_back(86400.0);  // 1 day default
    }
    
    // Limit the first few timesteps to improve convergence
    for (size_t i = 0; i < ctx.timesteps.size() && i < 5; i++) {
        if (ctx.timesteps[i] > 86400.0) {  // If greater than 1 day
            ctx.timesteps[i] = 86400.0;     // Limit to 1 day
        }
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "Timesteps: ");
    for (size_t i = 0; i < std::min(ctx.timesteps.size(), size_t(10)); i++) {
        PetscPrintf(PETSC_COMM_WORLD, "%.1f ", ctx.timesteps[i] / 86400.0);
    }
    if (ctx.timesteps.size() > 10) {
        PetscPrintf(PETSC_COMM_WORLD, "... (total %d steps)", static_cast<int>(ctx.timesteps.size()));
    }
    PetscPrintf(PETSC_COMM_WORLD, " days\n");
    
    return true;
}

void BlackOilSimulator::setupPETScObjects() {
    // Create DMDA
    DMDACreate3d(PETSC_COMM_WORLD,
                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_STAR,
                 ctx.grid.Nx, ctx.grid.Ny, ctx.grid.Nz,
                 PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                 3, 1, NULL, NULL, NULL, &ctx.da);
    DMSetUp(ctx.da);
    
    // Create vectors
    DMCreateGlobalVector(ctx.da, &ctx.solution);
    VecDuplicate(ctx.solution, &ctx.solution_old);
    
    // Create matrix
    DMCreateMatrix(ctx.da, &ctx.J);
    
    // Create SNES
    SNESCreate(PETSC_COMM_WORLD, &ctx.snes);
    SNESSetDM(ctx.snes, ctx.da);
    SNESSetFunction(ctx.snes, NULL, FormFunction, &ctx);
    SNESSetJacobian(ctx.snes, ctx.J, ctx.J, FormJacobian, &ctx);
    
    // Configure solver with more robust settings
    SNESSetType(ctx.snes, SNESNEWTONLS);
    SNESSetTolerances(ctx.snes, 1e-5, 1e-7, 1e-7, 100, 2000);  // More relaxed tolerances
    
    // Line search
    SNESLineSearch ls;
    SNESGetLineSearch(ctx.snes, &ls);
    SNESLineSearchSetType(ls, SNESLINESEARCHBT);
    SNESLineSearchSetOrder(ls, SNES_LINESEARCH_ORDER_CUBIC);
    
    // KSP with more robust settings
    KSP ksp;
    PC pc;
    SNESGetKSP(ctx.snes, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetTolerances(ksp, 1e-8, 1e-10, 1e6, 500);  // More relaxed KSP tolerances
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCILU);
    PCFactorSetLevels(pc, 2);  // Increase ILU levels
    
    SNESSetFromOptions(ctx.snes);
}

void BlackOilSimulator::outputResults() {
    PetscReal Qo_total = 0.0, Qw_total = 0.0, Qg_total = 0.0;
    PetscReal Qw_inj = 0.0;
    
    Vec local;
    Field ***x;
    DMGetLocalVector(ctx.da, &local);
    DMGlobalToLocalBegin(ctx.da, ctx.solution, INSERT_VALUES, local);
    DMGlobalToLocalEnd(ctx.da, ctx.solution, INSERT_VALUES, local);
    DMDAVecGetArrayRead(ctx.da, local, &x);
    
    DMDALocalInfo info;
    DMDAGetLocalInfo(ctx.da, &info);
    
    // Calculate well rates
    PetscPrintf(PETSC_COMM_WORLD, "  Processing %d wells...\n", static_cast<PetscInt>(ctx.wells.size()));
    for (const auto& well : ctx.wells) {
        if (!well.isOpen) {
            PetscPrintf(PETSC_COMM_WORLD, "    %s: SHUT\n", well.name.c_str());
            continue;
        }
        
        for (PetscInt k = well.k1; k <= well.k2; k++) {
            if (well.i >= info.xs && well.i < info.xs + info.xm &&
                well.j >= info.ys && well.j < info.ys + info.ym &&
                k >= info.zs && k < info.zs + info.zm) {
                
                PetscPrintf(PETSC_COMM_WORLD, "    %s: Processing cell (%d,%d,%d)\n", 
                           well.name.c_str(), well.i, well.j, k);
                
                // Bounds check
                if (well.i < 0 || well.i >= ctx.grid.Nx ||
                    well.j < 0 || well.j >= ctx.grid.Ny ||
                    k < 0 || k >= ctx.grid.Nz) {
                    continue;
                }
                
                PetscReal p = x[k][well.j][well.i].p;
                PetscReal Sw = x[k][well.j][well.i].Sw;
                PetscReal Sg = x[k][well.j][well.i].Sg;
                PetscReal So = 1.0 - Sw - Sg;
                
                // Calculate well index
                PetscReal WI = calculateWellIndex(well, ctx.grid, well.i, well.j, k);
                
                // Get fluid properties
                PetscReal Bo, Bw, Bg, muo, muw, mug;
                PetscReal Rs = 0.0;  // Simplified - should use RSVD
                if (!ctx.equil.rsvd_depth.empty()) {
                    PetscReal depth = 0.0;  // Calculate actual depth
                    Rs = interpolate(depth, ctx.equil.rsvd_depth, ctx.equil.rsvd_rs);
                }
                
                getOilProperties(ctx.fluid, p, Rs, Bo, muo);
                getWaterProperties(ctx.fluid, p, Bw, muw);
                getGasProperties(ctx.fluid, p, Bg, mug);
                
                // Get relative permeabilities
                PetscReal krw, kro, krg;
                getRelPerm(ctx.fluid, Sw, Sg, krw, kro, krg);
                
                PetscReal dp = p - well.BHP;
                
                if (well.type == 0) {  // Producer
                    PetscReal qo, qw, qg;
                    
                    if (well.control == 1) { // Rate control
                        // Calculate total mobility
                        PetscReal lambda_t = kro / muo + krw / muw + krg / mug;
                        if (lambda_t > 1e-12) {
                            // Distribute rate according to mobility
                            PetscReal q_total = well.rate; // Negative for production
                            qo = q_total * (kro / muo) / lambda_t;
                            qw = q_total * (krw / muw) / lambda_t;
                            qg = q_total * (krg / mug) / lambda_t;
                        } else {
                            qo = qw = qg = 0.0;
                        }
                        
                        // Debug output for rate control
                        if (well.name == "PROD1") {
                            PetscPrintf(PETSC_COMM_WORLD, "    DEBUG %s RATE: target=%.2e, lambda_t=%.2e, qo=%.2e, qw=%.2e, qg=%.2e\n",
                                       well.name.c_str(), well.rate, lambda_t, qo, qw, qg);
                        }
                    } else { // BHP control
                        PetscReal dp = p - well.BHP;
                        qo = WI * kro / muo * dp / Bo;
                        qw = WI * krw / muw * dp / Bw;
                        qg = WI * krg / mug * dp / Bg;
                        
                        // Debug output for BHP control
                        if (well.name == "PROD1") {
                            PetscPrintf(PETSC_COMM_WORLD, "    DEBUG %s BHP: p=%.1f, BHP=%.1f, dp=%.1f, WI=%.2e, kro=%.3f, qo=%.2e\n",
                                       well.name.c_str(), p/1e5, well.BHP/1e5, dp/1e5, WI, kro, qo);
                        }
                    }
                    
                    // Add dissolved gas
                    qg += Rs * qo;
                    
                    Qo_total += qo;
                    Qw_total += qw;
                    Qg_total += qg;
                } else {  // Injector
                    PetscReal qi;
                    if (well.control == 1) { // Rate control
                        qi = well.rate; // Positive for injection
                        
                        // Debug output for rate control injector
                        if (well.name == "INJ1") {
                            PetscPrintf(PETSC_COMM_WORLD, "    DEBUG %s RATE: target=%.2e, qi=%.2e\n",
                                       well.name.c_str(), well.rate, qi);
                        }
                    } else { // BHP control
                        PetscReal dp = well.BHP - p;
                        if (well.phase == "WATER") {
                            qi = WI / muw * dp / Bw;
                        } else if (well.phase == "GAS") {
                            qi = WI / mug * dp / Bg;
                        } else {
                            qi = 0.0;
                        }
                        
                        // Debug output for BHP control injector
                        if (well.name == "INJ1") {
                            PetscPrintf(PETSC_COMM_WORLD, "    DEBUG %s BHP: p=%.1f, BHP=%.1f, dp=%.1f, WI=%.2e, qi=%.2e\n",
                                       well.name.c_str(), p/1e5, well.BHP/1e5, dp/1e5, WI, qi);
                        }
                    }
                    
                    Qw_inj -= qi;
                }
            }
        }
    }
    
    DMDAVecRestoreArrayRead(ctx.da, local, &x);
    DMRestoreLocalVector(ctx.da, &local);
    
    // Print reservoir state changes for well cells
    DMGetLocalVector(ctx.da, &local);
    DMGlobalToLocalBegin(ctx.da, ctx.solution, INSERT_VALUES, local);
    DMGlobalToLocalEnd(ctx.da, ctx.solution, INSERT_VALUES, local);
    DMDAVecGetArrayRead(ctx.da, local, &x);
    
    PetscPrintf(PETSC_COMM_WORLD, "  Reservoir State at Wells:\n");
    for (const auto& well : ctx.wells) {
        if (well.isOpen) {
            for (PetscInt k = well.k1; k <= well.k2; k++) {
                if (well.i >= info.xs && well.i < info.xs + info.xm &&
                    well.j >= info.ys && well.j < info.ys + info.ym &&
                    k >= info.zs && k < info.zs + info.zm) {
                    
                    PetscReal p = x[k][well.j][well.i].p;
                    PetscReal Sw = x[k][well.j][well.i].Sw;
                    PetscReal Sg = x[k][well.j][well.i].Sg;
                    PetscReal So = 1.0 - Sw - Sg;
                    
                    // Get fluid properties for this cell
                    PetscReal Bo, Bw, Bg, muo, muw, mug;
                    PetscReal Rs = 0.0;
                    getOilProperties(ctx.fluid, p, Rs, Bo, muo);
                    getWaterProperties(ctx.fluid, p, Bw, muw);
                    getGasProperties(ctx.fluid, p, Bg, mug);
                    
                    // Get relative permeabilities
                    PetscReal krw, kro, krg;
                    getRelPerm(ctx.fluid, Sw, Sg, krw, kro, krg);
                    
                    PetscPrintf(PETSC_COMM_WORLD, "    %s(%d,%d,%d): p=%.1f bar, Sw=%.3f, So=%.3f, Sg=%.3f, krw=%.3f, kro=%.3f, krg=%.3f\n",
                               well.name.c_str(), well.i, well.j, k, p/1e5, Sw, So, Sg, krw, kro, krg);
                }
            }
        }
    }
    
    DMDAVecRestoreArrayRead(ctx.da, local, &x);
    DMRestoreLocalVector(ctx.da, &local);
    
    // Convert to field units for output
    PetscReal conv = 1.0;
    if (ctx.useField) {
        conv = 86400.0 / 0.1589873;  // m3/s to stb/day
    }
    
    PetscPrintf(PETSC_COMM_WORLD, 
               "  Production: Oil=%.1f, Water=%.1f, Gas=%.1f | Injection: Water=%.1f",
               Qo_total * conv, Qw_total * conv, 
               Qg_total * conv * 0.1589873 / 28.31685 * 1000.0,  // to Mscf/day
               Qw_inj * conv);
    
    if (Qo_total > 1e-10) {  // Check for meaningful oil production
        PetscReal GOR = Qg_total / Qo_total;
        if (ctx.useField) {
            GOR *= 0.1589873 / 28.31685 * 1000.0;  // Mscf/stb
        }
        PetscPrintf(PETSC_COMM_WORLD, " | GOR=%.1f", GOR);
    } else if (Qg_total > 1e-10) {
        PetscPrintf(PETSC_COMM_WORLD, " | GOR=INF (no oil production)");
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "\n");
}

void BlackOilSimulator::initialize() {
    setupPETScObjects();
    initializeEquilibration(&ctx);
    
    // Set initial timestep
    if (!ctx.timesteps.empty()) {
        ctx.dt = ctx.timesteps[0];
    } else {
        ctx.dt = 86400.0;  // 1 day default
    }
    
    // Calculate total time
    ctx.t_end = 0.0;
    for (auto dt : ctx.timesteps) {
        ctx.t_end += dt;
    }
    
    initialized = true;
}

void BlackOilSimulator::run() {
    ctx.t_current = 0.0;
    ctx.step = 0;
    ctx.timestep_idx = 0;

    // --- CSV Output Setup ---
    std::ofstream well_csv("well_output.csv");
    well_csv << "Step,Time_days,Well,Type,Layer,i,j,OilRate,WaterRate,GasRate,BHP,Pressure,Sw,So,Sg\n";
    // For summary
    struct WellTotals {
        double oil = 0, water = 0, gas = 0;
        double last_bhp = 0, last_p = 0;
    };
    std::map<std::string, WellTotals> well_totals;

    PetscPrintf(PETSC_COMM_WORLD, "\n=== Black Oil Reservoir Simulator ===\n");
    PetscPrintf(PETSC_COMM_WORLD, "Grid: %d x %d x %d = %d cells\n", 
                ctx.grid.Nx, ctx.grid.Ny, ctx.grid.Nz, 
                ctx.grid.Nx * ctx.grid.Ny * ctx.grid.Nz);
    PetscPrintf(PETSC_COMM_WORLD, "Wells: %d\n", static_cast<PetscInt>(ctx.wells.size()));
    
    // Print well information
    PetscPrintf(PETSC_COMM_WORLD, "Well Information:\n");
    for (const auto& well : ctx.wells) {
        PetscPrintf(PETSC_COMM_WORLD, "  %s: %s at (%d,%d) layers %d-%d, %s, %s=%.2f %s\n",
                   well.name.c_str(),
                   well.type == 0 ? "PRODUCER" : "INJECTOR",
                   well.i, well.j, well.k1, well.k2,
                   well.isOpen ? "OPEN" : "SHUT",
                   well.control == 1 ? "RATE" : "BHP",
                   well.control == 1 ? well.rate * 86400.0 / 0.1589873 : well.BHP / 1e5,
                   well.control == 1 ? "STB/day" : "bar");
    }
    
    // Print initial reservoir state summary
    PetscPrintf(PETSC_COMM_WORLD, "Initial Reservoir State Summary:\n");
    Vec local;
    Field ***x;
    DMGetLocalVector(ctx.da, &local);
    DMGlobalToLocalBegin(ctx.da, ctx.solution, INSERT_VALUES, local);
    DMGlobalToLocalEnd(ctx.da, ctx.solution, INSERT_VALUES, local);
    DMDAVecGetArrayRead(ctx.da, local, &x);
    
    PetscReal min_p = 1e20, max_p = -1e20;
    PetscReal min_Sw = 1.0, max_Sw = 0.0;
    PetscReal min_So = 1.0, max_So = 0.0;
    PetscReal min_Sg = 1.0, max_Sg = 0.0;
    
    for (PetscInt k = 0; k < ctx.grid.Nz; k++) {
        for (PetscInt j = 0; j < ctx.grid.Ny; j++) {
            for (PetscInt i = 0; i < ctx.grid.Nx; i++) {
                PetscReal p = x[k][j][i].p;
                PetscReal Sw = x[k][j][i].Sw;
                PetscReal Sg = x[k][j][i].Sg;
                PetscReal So = 1.0 - Sw - Sg;
                
                min_p = std::min(min_p, p);
                max_p = std::max(max_p, p);
                min_Sw = std::min(min_Sw, Sw);
                max_Sw = std::max(max_Sw, Sw);
                min_So = std::min(min_So, So);
                max_So = std::max(max_So, So);
                min_Sg = std::min(min_Sg, Sg);
                max_Sg = std::max(max_Sg, Sg);
            }
        }
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "  Pressure range: %.1f - %.1f bar\n", min_p/1e5, max_p/1e5);
    PetscPrintf(PETSC_COMM_WORLD, "  Water saturation range: %.3f - %.3f\n", min_Sw, max_Sw);
    PetscPrintf(PETSC_COMM_WORLD, "  Oil saturation range: %.3f - %.3f\n", min_So, max_So);
    PetscPrintf(PETSC_COMM_WORLD, "  Gas saturation range: %.3f - %.3f\n", min_Sg, max_Sg);
    
    DMDAVecRestoreArrayRead(ctx.da, local, &x);
    DMRestoreLocalVector(ctx.da, &local);
    
    PetscPrintf(PETSC_COMM_WORLD, "Simulation time: %.1f days\n\n", 
                ctx.t_end / 86400.0);
    
    // Initial output
    PetscPrintf(PETSC_COMM_WORLD, "Initial state:\n");
    
    // Print initial conditions for well cells
    DMGetLocalVector(ctx.da, &local);
    DMGlobalToLocalBegin(ctx.da, ctx.solution, INSERT_VALUES, local);
    DMGlobalToLocalEnd(ctx.da, ctx.solution, INSERT_VALUES, local);
    DMDAVecGetArrayRead(ctx.da, local, &x);
    
    for (const auto& well : ctx.wells) {
        if (well.isOpen) {
            PetscReal p = x[well.k1][well.j][well.i].p;
            PetscReal Sw = x[well.k1][well.j][well.i].Sw;
            PetscReal Sg = x[well.k1][well.j][well.i].Sg;
            PetscReal So = 1.0 - Sw - Sg;
            PetscPrintf(PETSC_COMM_WORLD, "  %s cell (%d,%d,%d): p=%.1f bar, Sw=%.3f, So=%.3f, Sg=%.3f\n",
                       well.name.c_str(), well.i, well.j, well.k1, p/1e5, Sw, So, Sg);
        }
    }
    
    DMDAVecRestoreArrayRead(ctx.da, local, &x);
    DMRestoreLocalVector(ctx.da, &local);
    
    outputResults();
    // Output initial well state to CSV
    {
        Vec local;
        Field ***x;
        DMGetLocalVector(ctx.da, &local);
        DMGlobalToLocalBegin(ctx.da, ctx.solution, INSERT_VALUES, local);
        DMGlobalToLocalEnd(ctx.da, ctx.solution, INSERT_VALUES, local);
        DMDAVecGetArrayRead(ctx.da, local, &x);
        for (const auto& well : ctx.wells) {
            if (well.isOpen) {
                for (PetscInt k = well.k1; k <= well.k2; k++) {
                    double p = x[k][well.j][well.i].p;
                    double Sw = x[k][well.j][well.i].Sw;
                    double Sg = x[k][well.j][well.i].Sg;
                    double So = 1.0 - Sw - Sg;
                    well_csv << ctx.step << "," << ctx.t_current/86400.0 << "," << well.name << ","
                             << (well.type==0?"PRODUCER":"INJECTOR") << "," << k << "," << well.i << "," << well.j
                             << ",0,0,0," << well.BHP/1e5 << "," << p/1e5 << "," << Sw << "," << So << "," << Sg << "\n";
                }
            }
        }
        DMDAVecRestoreArrayRead(ctx.da, local, &x);
        DMRestoreLocalVector(ctx.da, &local);
    }
    
    // Time stepping loop
    while (ctx.timestep_idx < ctx.timesteps.size()) {
        ctx.dt = ctx.timesteps[ctx.timestep_idx];
        ctx.t_current += ctx.dt;
        ctx.step++;
        
        PetscPrintf(PETSC_COMM_WORLD, "\nStep %d, Time = %.2f days (dt = %.2f days):\n", 
                   ctx.step, ctx.t_current / 86400.0, ctx.dt / 86400.0);
        
        // Copy current to old
        VecCopy(ctx.solution, ctx.solution_old);
        
        // Solve with progress monitoring
        PetscPrintf(PETSC_COMM_WORLD, "  Solving nonlinear system...\n");
        SNESSolve(ctx.snes, NULL, ctx.solution);
        
        // Check convergence
        SNESConvergedReason reason;
        SNESGetConvergedReason(ctx.snes, &reason);
        
        if (reason > 0) {
            PetscInt its;
            SNESGetIterationNumber(ctx.snes, &its);
            PetscPrintf(PETSC_COMM_WORLD, "  Converged in %d iterations\n", its);
            
            // Output results
            outputResults();
            // --- Output to CSV ---
            Vec local;
            Field ***x;
            DMGetLocalVector(ctx.da, &local);
            DMGlobalToLocalBegin(ctx.da, ctx.solution, INSERT_VALUES, local);
            DMGlobalToLocalEnd(ctx.da, ctx.solution, INSERT_VALUES, local);
            DMDAVecGetArrayRead(ctx.da, local, &x);
            for (const auto& well : ctx.wells) {
                if (!well.isOpen) continue;
                for (PetscInt k = well.k1; k <= well.k2; k++) {
                    double p = x[k][well.j][well.i].p;
                    double Sw = x[k][well.j][well.i].Sw;
                    double Sg = x[k][well.j][well.i].Sg;
                    double So = 1.0 - Sw - Sg;
                    // Calculate rates (reuse logic from outputResults)
                    double Bo, Bw, Bg, muo, muw, mug;
                    double Rs = 0.0;
                    getOilProperties(ctx.fluid, p, Rs, Bo, muo);
                    getWaterProperties(ctx.fluid, p, Bw, muw);
                    getGasProperties(ctx.fluid, p, Bg, mug);
                    double krw, kro, krg;
                    getRelPerm(ctx.fluid, Sw, Sg, krw, kro, krg);
                    double WI = calculateWellIndex(well, ctx.grid, well.i, well.j, k);
                    double qo=0, qw=0, qg=0;
                    if (well.type == 0) {
                        if (well.control == 1) {
                            double lambda_t = kro / muo + krw / muw + krg / mug;
                            if (lambda_t > 1e-12) {
                                double q_total = well.rate;
                                qo = q_total * (kro / muo) / lambda_t;
                                qw = q_total * (krw / muw) / lambda_t;
                                qg = q_total * (krg / mug) / lambda_t;
                            }
                        } else {
                            double dp = p - well.BHP;
                            qo = WI * kro / muo * dp / Bo;
                            qw = WI * krw / muw * dp / Bw;
                            qg = WI * krg / mug * dp / Bg;
                        }
                        qg += Rs * qo;
                        // Accumulate for summary
                        well_totals[well.name].oil += qo * ctx.dt;
                        well_totals[well.name].water += qw * ctx.dt;
                        well_totals[well.name].gas += qg * ctx.dt;
                        well_totals[well.name].last_bhp = well.BHP/1e5;
                        well_totals[well.name].last_p = p/1e5;
                    } else {
                        double qi = 0;
                        if (well.control == 1) {
                            qi = well.rate;
                        } else {
                            double dp = well.BHP - p;
                            if (well.phase == "WATER") qi = WI / muw * dp / Bw;
                            else if (well.phase == "GAS") qi = WI / mug * dp / Bg;
                        }
                        qw = qi;
                        well_totals[well.name].water += qi * ctx.dt;
                        well_totals[well.name].last_bhp = well.BHP/1e5;
                        well_totals[well.name].last_p = p/1e5;
                    }
                    // Convert to field units if needed
                    double conv = ctx.useField ? 86400.0 / 0.1589873 : 1.0;
                    well_csv << ctx.step << "," << ctx.t_current/86400.0 << "," << well.name << ","
                             << (well.type==0?"PRODUCER":"INJECTOR") << "," << k << "," << well.i << "," << well.j
                             << "," << qo*conv << "," << qw*conv << "," << qg*conv << "," << well.BHP/1e5 << "," << p/1e5 << "," << Sw << "," << So << "," << Sg << "\n";
                }
            }
            DMDAVecRestoreArrayRead(ctx.da, local, &x);
            DMRestoreLocalVector(ctx.da, &local);
            ctx.timestep_idx++;
        } else {
            PetscPrintf(PETSC_COMM_WORLD, "  Failed to converge (reason=%d)\n", reason);
            
            // Try to get more information about the failure
            PetscInt its;
            SNESGetIterationNumber(ctx.snes, &its);
            PetscPrintf(PETSC_COMM_WORLD, "  Iterations completed: %d\n", its);
            
            // Get residual norm
            Vec f;
            VecDuplicate(ctx.solution, &f);
            SNESComputeFunction(ctx.snes, ctx.solution, f);
            PetscReal norm;
            VecNorm(f, NORM_2, &norm);
            VecDestroy(&f);
            PetscPrintf(PETSC_COMM_WORLD, "  Final residual norm: %e\n", norm);
            
            // Try to continue with a smaller timestep if possible
            if (ctx.dt > 3600.0) {  // If timestep > 1 hour
                ctx.dt *= 0.5;  // Reduce timestep
                PetscPrintf(PETSC_COMM_WORLD, "  Reducing timestep to %.2f days and retrying...\n", ctx.dt / 86400.0);
                continue;  // Retry with smaller timestep
            } else {
                PetscPrintf(PETSC_COMM_WORLD, "  Timestep too small, stopping simulation\n");
                break;
            }
        }
    }
    // --- Write summary CSV ---
    std::ofstream summary_csv("well_summary.csv");
    summary_csv << "Well,Type,TotalOilProduced,TotalWaterProduced,TotalGasProduced,FinalBHP,FinalPressure\n";
    for (const auto& well : ctx.wells) {
        const auto& totals = well_totals[well.name];
        double conv = ctx.useField ? 86400.0 / 0.1589873 : 1.0;
        summary_csv << well.name << "," << (well.type==0?"PRODUCER":"INJECTOR") << ","
                    << totals.oil*conv << "," << totals.water*conv << "," << totals.gas*conv << ","
                    << totals.last_bhp << "," << totals.last_p << "\n";
    }
    summary_csv.close();
    well_csv.close();
    PetscPrintf(PETSC_COMM_WORLD, "\nSimulation completed. Results written to well_output.csv and well_summary.csv\n");
}

} // namespace BlackOil