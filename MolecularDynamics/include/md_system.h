#ifndef MD_SYSTEM_H
#define MD_SYSTEM_H

#include "particle.h"
#include "integrator.h"
#include "molecule_registry.h"
#include "constants.h"
#include <vector>
#include <string>
#include <string_view>

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

public:
    explicit OptimizedMolecularDynamics(double box);

    void update();
    void set_langevin_friction(double friction);
    void set_target_temperature(double temp) { target_temperature = temp; }
    void set_dt(double delta_t) { dt = delta_t; }

    const std::vector<Particle>& get_particles() const { return particles; }
    const std::vector<Molecule>& get_molecules() const { return molecules; }
    double get_box_size() const { return box_size; }
    double get_dt() const { return dt; }

    double calculate_bond_length(const Particle& p1, const Particle& p2) const;
    double calculate_angle(const Particle& p1, const Particle& p2, const Particle& p3) const;

    void add_molecule(const std::string& type, double center_x, double center_y, double center_z);
    double calculate_kinetic_energy() const;
    double calculate_potential_energy() const;
    void add_thermal_velocities(double temperature);
    double calculate_temperature() const;
    void analyze_stability() const;
};

#endif