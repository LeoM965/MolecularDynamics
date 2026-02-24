#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "particle.h"
#include "lennard_jones.h"
#include "constants.h"
#include <vector>
#include <random>

class OptimizedVerletIntegrator {
private:
    std::mt19937 rng;
    std::normal_distribution<double> normal_dist;

public:
    OptimizedVerletIntegrator();

    LennardJones lj_calculator;

    Vec3 apply_pbc_vector(const Vec3& v, double box_size) const;
    double apply_pbc(double x, double box_size) const;

    void calculate_bond_forces(std::vector<Particle>& particles, double box_size,
        const std::vector<Molecule>& molecules);
    void calculate_angle_forces(std::vector<Particle>& particles, double box_size,
        const std::vector<Molecule>& molecules);
    void limit_forces(std::vector<Particle>& particles);

    void update_positions(std::vector<Particle>& particles, double box_size, double dt,
        const std::vector<Molecule>& molecules);

    void apply_langevin_thermostat(std::vector<Particle>& particles,
        double target_temp, double friction_coeff, double dt);

    double calculate_temperature(const std::vector<Particle>& particles) const;
    double compute_total_energy(const std::vector<Particle>& particles,
        double box_size, const std::vector<Molecule>& molecules) const;
};

#endif