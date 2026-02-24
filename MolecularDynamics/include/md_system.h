#ifndef MD_SYSTEM_H
#define MD_SYSTEM_H

#include "particle.h"
#include "integrator.h"
#include "molecule_registry.h"
#include "constants.h"
#include <vector>
#include <string>

class OptimizedMolecularDynamics {
private:
    std::vector<Particle> particles;
    std::vector<Molecule> molecules;
    OptimizedVerletIntegrator integrator;
    MoleculeRegistry registry;

    double box_size;
    double dt;
    double target_temperature;
    double langevin_friction;

    // File initialization flags
    bool energy_file_initialized;
    bool trajectory_file_initialized;
    bool bonds_angles_file_initialized;
    bool lj_forces_file_initialized;

public:
    OptimizedMolecularDynamics(double box);

    void update();
    void set_langevin_friction(double friction);

    double calculate_bond_length(const Particle& p1, const Particle& p2) const;
    double calculate_angle(const Particle& p1, const Particle& p2, const Particle& p3) const;

    void add_molecule(const std::string& type, double center_x, double center_y, double center_z);
    double calculate_kinetic_energy() const;
    void add_thermal_velocities(double temperature);

    double calculate_temperature() const;
    void print_status(int step) const;
    void analyze_stability() const;

    // Export functions
    void export_molecule_topology(const std::string& filename) const;
    void export_xyz(const std::string& filename, int step) const;
    void export_energy_data(const std::string& filename, int step);
    void export_bonds_angles(const std::string& filename, int step);
    void export_csv_trajectory(const std::string& filename, int step);
    void export_lj_forces(const std::string& filename, int step);
};

#endif