#ifndef MD_REPORTER_H
#define MD_REPORTER_H

#include "md_system.h"
#include "config.h"
#include <string>
#include <string_view>
#include <fstream>

class MDReporter {
private:
    bool energy_init = false;
    bool traj_init = false;
    bool bonds_init = false;
    
    void ensure_directory(std::string_view path) const;

public:
    void export_topology(const OptimizedMolecularDynamics& engine, const std::string& filename) const;
    void export_xyz(const OptimizedMolecularDynamics& engine, const std::string& filename, int step) const;
    void export_energy(const OptimizedMolecularDynamics& engine, const std::string& filename, int step);
    void export_trajectory(const OptimizedMolecularDynamics& engine, const std::string& filename, int step);
    void export_bonds_angles(const OptimizedMolecularDynamics& engine, const std::string& filename, int step);
};

#endif
