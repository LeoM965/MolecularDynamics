#include "integrator.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <iostream>

using namespace std;

OptimizedVerletIntegrator::OptimizedVerletIntegrator()
    : rng(random_device{}()), normal_dist(0.0, 1.0) {
}

Vec3 OptimizedVerletIntegrator::apply_pbc_vector(const Vec3& v, double box_size) const {
    if (box_size <= 0) throw invalid_argument("Box size must be positive");
    return Vec3(this->apply_pbc(v.x, box_size), this->apply_pbc(v.y, box_size), this->apply_pbc(v.z, box_size));
}

double OptimizedVerletIntegrator::apply_pbc(double x, double box_size) const {
    if (box_size <= 0) throw invalid_argument("Box size must be positive");
    return x - box_size * round(x / box_size);
}

void OptimizedVerletIntegrator::calculate_bond_forces(vector<Particle>& particles, double box_size,
    const vector<Molecule>& molecules) {
    if (box_size <= 0) throw invalid_argument("Box size must be positive");
    if (particles.empty()) throw invalid_argument("Particle list cannot be empty");

    Vec3 dr;
    for (size_t m = 0; m < molecules.size(); ++m) {
        const auto& mol = molecules[m];
        if (mol.bonds.size() != mol.bond_lengths.size()) throw invalid_argument("Mismatch in bonds and bond lengths");
        for (size_t i = 0; i < mol.bonds.size(); ++i) {
            int idx1 = mol.bonds[i].first;
            int idx2 = mol.bonds[i].second;
            if (idx1 >= static_cast<int>(particles.size()) || idx2 >= static_cast<int>(particles.size())) {
                throw out_of_range("Invalid bond indices");
            }
            double r0 = mol.bond_lengths[i] * Constants::ANGSTROM_TO_METERS;

            const auto& pos1 = particles[idx1].position;
            const auto& pos2 = particles[idx2].position;
            dr = this->apply_pbc_vector(pos1 - pos2, box_size);
            double r = dr.magnitude();

            if (r > 1e-12) {
                double f_magnitude = -Constants::BOND_FORCE_CONST * (r - r0) / r;
                Vec3 force = dr * f_magnitude;

                particles[idx1].force += force;
                particles[idx2].force -= force;
            }
        }
    }
}

void OptimizedVerletIntegrator::calculate_angle_forces(vector<Particle>& particles, double box_size,
    const vector<Molecule>& molecules) {
    if (box_size <= 0) throw invalid_argument("Box size must be positive");
    if (particles.empty()) throw invalid_argument("Particle list cannot be empty");

    Vec3 a, b;
    for (size_t m = 0; m < molecules.size(); ++m) {
        const auto& mol = molecules[m];
        if (mol.angles.size() != mol.angle_degrees.size()) throw invalid_argument("Mismatch in angles and angle degrees");
        for (size_t i = 0; i < mol.angles.size(); ++i) {
            int idx1 = mol.angles[i][0];
            int idx2 = mol.angles[i][1];
            int idx3 = mol.angles[i][2];
            if (idx1 >= static_cast<int>(particles.size()) || idx2 >= static_cast<int>(particles.size()) ||
                idx3 >= static_cast<int>(particles.size())) {
                throw out_of_range("Invalid angle indices");
            }
            double theta0 = mol.angle_degrees[i] * Constants::DEG_TO_RAD;

            const auto& pos1 = particles[idx1].position;
            const auto& pos2 = particles[idx2].position;
            const auto& pos3 = particles[idx3].position;
            a = this->apply_pbc_vector(pos1 - pos2, box_size);
            b = this->apply_pbc_vector(pos3 - pos2, box_size);

            double ra = a.magnitude();
            double rb = b.magnitude();

            if (ra > 1e-12 && rb > 1e-12) {
                double cos_theta = max(-1.0, min(1.0, a.dot(b) / (ra * rb)));
                double theta = acos(cos_theta);
                double dtheta = theta - theta0;
                double sin_theta = sin(theta);

                if (abs(sin_theta) > 1e-6) {
                    double f = -Constants::ANGLE_FORCE_CONST * dtheta;

                    // Scale angular deviation
                    double fa = f / (ra * sin_theta);
                    double fb = f / (rb * sin_theta);

                    // Direction of bonds
                    Vec3 a_norm = a * (1.0 / ra);
                    Vec3 b_norm = b * (1.0 / rb);

                    // Calculate forces
                    Vec3 force1 = (b_norm - a_norm * cos_theta) * fa;
                    Vec3 force3 = (a_norm - b_norm * cos_theta) * fb;

                    particles[idx1].force += force1;
                    particles[idx3].force += force3;
                    particles[idx2].force -= (force1 + force3);
                }
            }
        }
    }
}

void OptimizedVerletIntegrator::limit_forces(vector<Particle>& particles) {
    for (auto& p : particles) {
        double f_mag = p.force.magnitude();
        if (f_mag > Constants::MAX_FORCE) {
            p.force *= (Constants::MAX_FORCE / f_mag);
        }
    }
}

void OptimizedVerletIntegrator::update_positions(vector<Particle>& particles, double box_size, double dt,
    const vector<Molecule>& molecules) {
    if (box_size <= 0) throw invalid_argument("Box size must be positive");
    if (dt <= 0) throw invalid_argument("Time step must be positive");
    if (particles.empty()) throw invalid_argument("Particle list cannot be empty");

    // First force calculation
    for (auto& p : particles) {
        p.force = Vec3(0.0, 0.0, 0.0);
    }

    this->calculate_bond_forces(particles, box_size, molecules);
    this->calculate_angle_forces(particles, box_size, molecules);
    lj_calculator.calculate_lj_forces(particles, box_size);
    this->limit_forces(particles);

    // Update positions
    for (auto& p : particles) {
        double inv_mass = 1.0 / (p.mass * Constants::ATOMIC_MASS_UNIT);
        Vec3 acceleration = p.force * inv_mass;

        p.position += p.velocity * dt + acceleration * (0.5 * dt * dt);
        p.position.x = p.position.x - box_size * floor(p.position.x / box_size);
        p.position.y = p.position.y - box_size * floor(p.position.y / box_size);
        p.position.z = p.position.z - box_size * floor(p.position.z / box_size);
    }

    // Second force calculation
    for (auto& p : particles) {
        p.force = Vec3(0.0, 0.0, 0.0);
    }
    this->calculate_bond_forces(particles, box_size, molecules);
    this->calculate_angle_forces(particles, box_size, molecules);
    lj_calculator.calculate_lj_forces(particles, box_size);
    this->limit_forces(particles);

    // Update velocities
    for (auto& p : particles) {
        double inv_mass = 1.0 / (p.mass * Constants::ATOMIC_MASS_UNIT);
        p.velocity += p.force * (inv_mass * 0.5 * dt);
    }
}

void OptimizedVerletIntegrator::apply_langevin_thermostat(vector<Particle>& particles,
    double target_temp, double friction_coeff, double dt) {
    if (target_temp <= 0) throw invalid_argument("Target temperature must be positive");
    if (friction_coeff <= 0) throw invalid_argument("Friction coefficient must be positive");
    if (dt <= 0) throw invalid_argument("Time step must be positive");
    if (particles.empty()) throw invalid_argument("Particle list cannot be empty");

    double gamma = friction_coeff;
    for (auto& p : particles) {
        if (p.mass <= 0) throw invalid_argument("Particle mass must be positive");

        double mass_si = p.mass * Constants::ATOMIC_MASS_UNIT;
        double friction_factor = exp(-gamma * dt);
        double noise_variance = (Constants::BOLTZMANN * target_temp / mass_si) * (1.0 - friction_factor * friction_factor);
        double noise_amplitude = sqrt(noise_variance);

        double rand_x = normal_dist(rng);
        double rand_y = normal_dist(rng);
        double rand_z = normal_dist(rng);

        p.velocity.x = p.velocity.x * friction_factor + noise_amplitude * rand_x;
        p.velocity.y = p.velocity.y * friction_factor + noise_amplitude * rand_y;
        p.velocity.z = p.velocity.z * friction_factor + noise_amplitude * rand_z;
    }
}

double OptimizedVerletIntegrator::calculate_temperature(const vector<Particle>& particles) const {
    if (particles.empty()) throw invalid_argument("Particle list cannot be empty");

    double total_kinetic = 0.0;
    for (const auto& p : particles) {
        if (p.mass <= 0) throw invalid_argument("Particle mass must be positive");
        double m = p.mass * Constants::ATOMIC_MASS_UNIT;
        total_kinetic += 0.5 * m * p.velocity.magnitude_squared();
    }

    int dof = max(1, static_cast<int>(particles.size() * 3 - 6));
    return (2.0 * total_kinetic) / (dof * Constants::BOLTZMANN);
}

double OptimizedVerletIntegrator::compute_total_energy(const vector<Particle>& particles,
    double box_size, const vector<Molecule>& molecules) const {
    if (box_size <= 0) throw invalid_argument("Box size must be positive");
    if (particles.empty()) throw invalid_argument("Particle list cannot be empty");

    double potential_energy = 0.0;
    double kinetic_energy = 0.0;

    // Bond potential energy
    for (const auto& mol : molecules) {
        if (mol.bonds.size() != mol.bond_lengths.size()) throw invalid_argument("Mismatch in bonds and bond lengths");
        for (size_t i = 0; i < mol.bonds.size(); ++i) {
            int idx1 = mol.bonds[i].first;
            int idx2 = mol.bonds[i].second;
            if (idx1 >= static_cast<int>(particles.size()) || idx2 >= static_cast<int>(particles.size())) {
                throw out_of_range("Invalid bond indices");
            }
            double r0 = mol.bond_lengths[i] * Constants::ANGSTROM_TO_METERS;

            Vec3 dr = this->apply_pbc_vector(particles[idx1].position - particles[idx2].position, box_size);
            double r = dr.magnitude();
            double dr_val = r - r0;

            potential_energy += 0.5 * Constants::BOND_FORCE_CONST * dr_val * dr_val;
        }

        // Angle potential energy
        if (mol.angles.size() != mol.angle_degrees.size()) throw invalid_argument("Mismatch in angles and angle degrees");
        for (size_t i = 0; i < mol.angles.size(); ++i) {
            int idx1 = mol.angles[i][0];
            int idx2 = mol.angles[i][1];
            int idx3 = mol.angles[i][2];
            if (idx1 >= static_cast<int>(particles.size()) || idx2 >= static_cast<int>(particles.size()) ||
                idx3 >= static_cast<int>(particles.size())) {
                throw out_of_range("Invalid angle indices");
            }
            double theta0 = mol.angle_degrees[i] * Constants::DEG_TO_RAD;

            Vec3 a = this->apply_pbc_vector(particles[idx1].position - particles[idx2].position, box_size);
            Vec3 b = this->apply_pbc_vector(particles[idx3].position - particles[idx2].position, box_size);

            double ra = a.magnitude();
            double rb = b.magnitude();

            if (ra > 1e-12 && rb > 1e-12) {
                double cos_theta = max(-1.0, min(1.0, a.dot(b) / (ra * rb)));
                double theta = acos(cos_theta);
                double dtheta = theta - theta0;

                potential_energy += 0.5 * Constants::ANGLE_FORCE_CONST * dtheta * dtheta;
            }
        }
    }

    // Lennard-Jones potential energy
    potential_energy += lj_calculator.compute_lj_energy(particles, box_size);

    // Kinetic energy
    for (const auto& p : particles) {
        if (p.mass <= 0) throw invalid_argument("Particle mass must be positive");
        double m = p.mass * Constants::ATOMIC_MASS_UNIT;
        kinetic_energy += 0.5 * m * p.velocity.magnitude_squared();
    }

    return (potential_energy + kinetic_energy) / Constants::KCAL_PER_MOL_TO_JOULES;
}