#include "../include/md_system.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <cmath>
#include <algorithm>

OptimizedMolecularDynamics::OptimizedMolecularDynamics(double box)
    : box_size(box * Constants::ANGSTROM_TO_METERS),
      dt(0.02e-15),
      target_temperature(300.0),
      langevin_friction(1.0e10) {
    if (box <= 0) throw std::runtime_error("Invalid box size");
}

void OptimizedMolecularDynamics::set_langevin_friction(double friction) {
    if (friction <= 0) throw std::invalid_argument("Friction must be positive");
    langevin_friction = friction;
}

double OptimizedMolecularDynamics::calculate_bond_length(const Particle& p1, const Particle& p2) const {
    Vec3 dr = Vec3::apply_pbc_vector(p1.position - p2.position, box_size);
    return dr.magnitude() / Constants::ANGSTROM_TO_METERS;
}

double OptimizedMolecularDynamics::calculate_angle(const Particle& p1, const Particle& p2, const Particle& p3) const {
    Vec3 a = Vec3::apply_pbc_vector(p1.position - p2.position, box_size);
    Vec3 b = Vec3::apply_pbc_vector(p3.position - p2.position, box_size);
    double ra = a.magnitude();
    double rb = b.magnitude();
    if (ra < 1e-12 || rb < 1e-12) return 0.0;
    double cos_theta = std::max(-1.0, std::min(1.0, a.dot(b) / (ra * rb)));
    return std::acos(cos_theta) * 180.0 / Constants::PI;
}

void OptimizedMolecularDynamics::add_molecule(const std::string& type, double center_x, double center_y, double center_z) {
    const auto& def = registry.get_definition(std::string(type));
    int base = static_cast<int>(particles.size());

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-0.01, 0.01);

    Molecule mol;
    mol.type = type;
    for (size_t i = 0; i < def.atoms.size(); ++i) {
        particles.emplace_back(
            center_x + def.positions[i].x + dist(gen),
            center_y + def.positions[i].y + dist(gen),
            center_z + def.positions[i].z + dist(gen),
            def.atoms[i].second, def.atoms[i].first
        );
        mol.atom_indices.push_back(base + (int)i);
    }

    for (const auto& b : def.bonds) mol.bonds.emplace_back(base + b.first, base + b.second);
    for (const auto& a : def.angles) mol.angles.push_back({base + a[0], base + a[1], base + a[2]});
    
    mol.bond_lengths = def.bond_lengths;
    mol.angle_degrees = def.angle_degrees;
    molecules.push_back(std::move(mol));
}

double OptimizedMolecularDynamics::calculate_kinetic_energy() const {
    double kin = 0.0;
    double factor = Constants::ATOMIC_MASS_UNIT / Constants::KCAL_PER_MOL_TO_JOULES;
    for (const auto& p : particles) kin += 0.5 * p.mass * factor * p.velocity.magnitude_squared();
    return kin;
}

void OptimizedMolecularDynamics::add_thermal_velocities(double temp) {
    target_temperature = temp;
    std::mt19937 gen(std::random_device{}());
    for (auto& p : particles) {
        double sigma = std::sqrt(Constants::BOLTZMANN * temp / (p.mass * Constants::ATOMIC_MASS_UNIT));
        std::normal_distribution<double> d(0, sigma);
        p.velocity = Vec3(d(gen), d(gen), d(gen));
    }
}

void OptimizedMolecularDynamics::update() {
    integrator.update_positions(particles, box_size, dt, molecules);
    integrator.apply_langevin_thermostat(particles, target_temperature, langevin_friction, dt);
}

double OptimizedMolecularDynamics::calculate_temperature() const {
    return integrator.calculate_temperature(particles);
}

void OptimizedMolecularDynamics::analyze_stability() const {
    // Basic implementation for now
}