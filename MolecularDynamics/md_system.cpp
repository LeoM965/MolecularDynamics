#include "md_system.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <cmath>
#include <algorithm>

using namespace std;

OptimizedMolecularDynamics::OptimizedMolecularDynamics(double box)
    : box_size(box* Constants::ANGSTROM_TO_METERS),
    dt(0.02e-15),
    target_temperature(300.0),
    langevin_friction(1.0e10),
    energy_file_initialized(false),
    trajectory_file_initialized(false),
    bonds_angles_file_initialized(false),
    lj_forces_file_initialized(false) {

    if (box <= 0) throw runtime_error("Invalid box size");

    cout << fixed << setprecision(6);
    cout << "Box size: " << box << " A\n";
    cout << "Time step: " << dt * 1e15 << " fs\n";
    cout << "Langevin friction coefficient: " << langevin_friction << " s^-1\n";
}

void OptimizedMolecularDynamics::set_langevin_friction(double friction) {
    if (friction <= 0) throw invalid_argument("Friction coefficient must be positive");
    langevin_friction = friction;
    cout << "Langevin friction coefficient set to: " << friction << " s^-1\n";
}

inline double OptimizedMolecularDynamics::calculate_bond_length(const Particle& p1, const Particle& p2) const {
    Vec3 dr(
        integrator.apply_pbc(p1.position.x - p2.position.x, box_size),
        integrator.apply_pbc(p1.position.y - p2.position.y, box_size),
        integrator.apply_pbc(p1.position.z - p2.position.z, box_size)
    );
    return dr.magnitude() / Constants::ANGSTROM_TO_METERS;
}

double OptimizedMolecularDynamics::calculate_angle(const Particle& p1, const Particle& p2, const Particle& p3) const {
    Vec3 a(
        integrator.apply_pbc(p1.position.x - p2.position.x, box_size),
        integrator.apply_pbc(p1.position.y - p2.position.y, box_size),
        integrator.apply_pbc(p1.position.z - p2.position.z, box_size)
    );
    Vec3 b(
        integrator.apply_pbc(p3.position.x - p2.position.x, box_size),
        integrator.apply_pbc(p3.position.y - p2.position.y, box_size),
        integrator.apply_pbc(p3.position.z - p2.position.z, box_size)
    );

    const double ra = a.magnitude();
    const double rb = b.magnitude();

    if (ra < 1e-12 || rb < 1e-12) return 0.0;

    const double cos_theta = max(-1.0, min(1.0, a.dot(b) / (ra * rb)));
    return acos(cos_theta) * 180.0 / Constants::PI;
}

void OptimizedMolecularDynamics::add_molecule(const string& type, double center_x, double center_y, double center_z) {
    const MoleculeDefinition& def = registry.get_definition(type);
    const int base_index = static_cast<int>(particles.size());

    Molecule mol;
    mol.type = type;
    mol.atom_indices.reserve(def.atoms.size());
    mol.bonds.reserve(def.bonds.size());
    mol.angles.reserve(def.angles.size());

    // Add atoms with small random displacement to avoid perfect overlaps
    static mt19937 gen(random_device{}());
    static uniform_real_distribution<double> displacement(-0.01, 0.01); // 0.01 A displacement

    for (size_t i = 0; i < def.atoms.size(); ++i) {
        const double x = center_x + def.positions[i].x + displacement(gen);
        const double y = center_y + def.positions[i].y + displacement(gen);
        const double z = center_z + def.positions[i].z + displacement(gen);

        particles.emplace_back(x, y, z, def.atoms[i].second, def.atoms[i].first);
        mol.atom_indices.push_back(base_index + static_cast<int>(i));
    }

    // Add bonds
    for (const auto& bond : def.bonds) {
        mol.bonds.emplace_back(base_index + bond.first, base_index + bond.second);
    }

    // Add angles
    for (const auto& angle : def.angles) {
        mol.angles.push_back({ base_index + angle[0], base_index + angle[1], base_index + angle[2] });
    }

    mol.bond_lengths = def.bond_lengths;
    mol.angle_degrees = def.angle_degrees;
    molecules.push_back(move(mol));

    cout << "Molecule " << type << " added with " << def.bonds.size()
        << " bonds and " << def.angles.size() << " angles.\n";
}

double OptimizedMolecularDynamics::calculate_kinetic_energy() const {
    double kinetic_energy = 0.0;
    const double conversion_factor = Constants::ATOMIC_MASS_UNIT / Constants::KCAL_PER_MOL_TO_JOULES;

    for (const auto& p : particles) {
        kinetic_energy += 0.5 * p.mass * conversion_factor * p.velocity.magnitude_squared();
    }

    return kinetic_energy;
}

void OptimizedMolecularDynamics::add_thermal_velocities(double temperature) {
    target_temperature = temperature;
    mt19937 gen(random_device{}());

    // Add Maxwell-Boltzmann velocities
    for (auto& p : particles) {
        const double m = p.mass * Constants::ATOMIC_MASS_UNIT;
        const double sigma = sqrt(Constants::BOLTZMANN * temperature / m);
        normal_distribution<double> dist(0.0, sigma);

        p.velocity = Vec3(dist(gen), dist(gen), dist(gen));
    }

    // Remove center of mass velocity
    Vec3 cm_velocity(0.0, 0.0, 0.0);
    double total_mass = 0.0;
    for (const auto& p : particles) {
        cm_velocity += p.velocity * p.mass;
        total_mass += p.mass;
    }

    if (total_mass > 0) {
        cm_velocity *= (1.0 / total_mass);
        for (auto& p : particles) {
            p.velocity -= cm_velocity;
        }
    }

    // Scale to exact target temperature
    const double current_temp = calculate_temperature();
    if (current_temp > 0) {
        const double scale_factor = sqrt(temperature / current_temp);
        for (auto& p : particles) {
            p.velocity *= scale_factor;
        }
    }

    cout << "Thermal velocities added for T = " << temperature << " K\n";
    cout << "Initial kinetic energy: " << fixed << setprecision(4)
        << calculate_kinetic_energy() << " kcal/mol\n";
    cout << "Initial temperature: " << fixed << setprecision(2)
        << calculate_temperature() << " K\n";
}
 void OptimizedMolecularDynamics::update() {
    integrator.update_positions(particles, box_size, dt, molecules);
    integrator.apply_langevin_thermostat(particles, target_temperature, langevin_friction, dt);
}

 double OptimizedMolecularDynamics::calculate_temperature() const {
    return integrator.calculate_temperature(particles);
}

void OptimizedMolecularDynamics::print_status(int step) const {
    const double energy = integrator.compute_total_energy(particles, box_size, molecules);
    const double temp = calculate_temperature();

    cout << "\n=== Step " << step << " ===\n";
    cout << "Total energy: " << fixed << setprecision(4) << energy << " kcal/mol\n";
    cout << "Temperature: " << fixed << setprecision(2) << temp << " K\n";

    // Quick stability check
    if (temp > 2000.0 || abs(energy) > 10000.0) {
        cout << "WARNING: System may be unstable!\n";
    }

    for (size_t i = 0; i < molecules.size(); ++i) {
        const Molecule& mol = molecules[i];
        cout << "\nMolecule " << mol.type << " " << i << ":\n";

        // Print bonds
        cout << "  Bonds:\n";
        for (size_t j = 0; j < mol.bonds.size(); ++j) {
            const int idx1 = mol.bonds[j].first;
            const int idx2 = mol.bonds[j].second;
            const double length = calculate_bond_length(particles[idx1], particles[idx2]);
            const double deviation = abs(length - mol.bond_lengths[j]);

            cout << "    Bond " << j + 1 << " (" << idx1 << "-" << idx2
                << "): " << fixed << setprecision(4) << length
                << " A (equilibrium: " << mol.bond_lengths[j]
                << " A, deviation: " << deviation << " A)\n";
        }

        // Print angles
        cout << "  Angles:\n";
        for (size_t j = 0; j < mol.angles.size(); ++j) {
            const int idx1 = mol.angles[j][0];
            const int idx2 = mol.angles[j][1];
            const int idx3 = mol.angles[j][2];
            const double angle = calculate_angle(particles[idx1], particles[idx2], particles[idx3]);
            const double deviation = abs(angle - mol.angle_degrees[j]);

            cout << "    Angle " << j + 1 << " (" << idx1 << "-" << idx2 << "-" << idx3
                << "): " << fixed << setprecision(2) << angle
                << " degrees (equilibrium: " << mol.angle_degrees[j]
                << " degrees, deviation: " << deviation << " degrees)\n";
        }
    }
}

void OptimizedMolecularDynamics::analyze_stability() const {
    cout << "\n=== STABILITY ANALYSIS ===\n";

    double max_bond_deviation = 0.0;
    double max_angle_deviation = 0.0;
    double total_bond_deviation = 0.0;
    double total_angle_deviation = 0.0;
    size_t bond_count = 0;
    size_t angle_count = 0;

    // Calculate deviations efficiently
    for (const auto& mol : molecules) {
        for (size_t j = 0; j < mol.bonds.size(); ++j) {
            const int idx1 = mol.bonds[j].first;
            const int idx2 = mol.bonds[j].second;
            const double length = calculate_bond_length(particles[idx1], particles[idx2]);
            const double deviation = abs(length - mol.bond_lengths[j]);

            max_bond_deviation = max(max_bond_deviation, deviation);
            total_bond_deviation += deviation;
            ++bond_count;
        }

        for (size_t j = 0; j < mol.angles.size(); ++j) {
            const int idx1 = mol.angles[j][0];
            const int idx2 = mol.angles[j][1];
            const int idx3 = mol.angles[j][2];
            const double angle = calculate_angle(particles[idx1], particles[idx2], particles[idx3]);
            const double deviation = abs(angle - mol.angle_degrees[j]);

            max_angle_deviation = max(max_angle_deviation, deviation);
            total_angle_deviation += deviation;
            ++angle_count;
        }
    }

    const double avg_bond_deviation = bond_count > 0 ? total_bond_deviation / bond_count : 0.0;
    const double avg_angle_deviation = angle_count > 0 ? total_angle_deviation / angle_count : 0.0;

    cout << "Maximum bond deviation: " << fixed << setprecision(4) << max_bond_deviation << " A\n";
    cout << "Average bond deviation: " << fixed << setprecision(4) << avg_bond_deviation << " A\n";
    cout << "Maximum angle deviation: " << fixed << setprecision(2) << max_angle_deviation << " degrees\n";
    cout << "Average angle deviation: " << fixed << setprecision(2) << avg_angle_deviation << " degrees\n";

    // Stability assessment
    if (max_bond_deviation < 0.1 && max_angle_deviation < 15.0 &&
        avg_bond_deviation < 0.05 && avg_angle_deviation < 8.0) {
        cout << "Simulation is STABLE\n";
    }
    else if (max_bond_deviation < 0.5 && max_angle_deviation < 30.0) {
        cout << "Simulation is MODERATE\n";
    }
    else {
        cout << "Simulation is UNSTABLE\n";

        const double temp = calculate_temperature();
        const double energy = integrator.compute_total_energy(particles, box_size, molecules);

        cout << "Current temperature: " << temp << " K, Energy: " << energy << " kcal/mol\n";
        if (temp > 1000.0 || abs(energy) > 5000.0) {
            cout << "RECOMMENDATION: Reduce timestep to 0.01 fs or lower friction coefficient\n";
        }
    }
}