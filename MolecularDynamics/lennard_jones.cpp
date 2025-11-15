#include "lennard_jones.h"
#include "constants.h"
#include <cmath>

using namespace std;

LennardJones::LennardJones() {
    lj_epsilon.resize(6);
    lj_sigma.resize(6);

    const double KCAL_TO_J = Constants::KCAL_PER_MOL_TO_JOULES / Constants::AVOGADRO;

    for (int i = 0; i < 6; ++i) {
        lj_epsilon[i] = Constants::LJ_EPSILON[i] * KCAL_TO_J;
        lj_sigma[i] = Constants::LJ_SIGMA[i] * Constants::ANGSTROM_TO_METERS;
    }

    lj_cutoff = Constants::LJ_CUTOFF * Constants::ANGSTROM_TO_METERS;
}

inline int LennardJones::get_atom_type_index(AtomType atom_type) const {
    return static_cast<int>(atom_type);
}

inline Vec3 LennardJones::apply_pbc_vector(const Vec3& v, double box_size) const {
    if (box_size <= 0) return v;
    return Vec3(
        v.x - box_size * round(v.x / box_size),
        v.y - box_size * round(v.y / box_size),
        v.z - box_size * round(v.z / box_size)
    );
}

void LennardJones::calculate_lj_forces(vector<Particle>& particles, double box_size) {
    const size_t n = particles.size();
    if (n == 0) return;

    for (auto& p : particles) {
        p.force = Vec3(0.0, 0.0, 0.0);
    }

    const double cutoff_sq = lj_cutoff * lj_cutoff;

    for (size_t i = 0; i < n - 1; ++i) {
        const int type_i = static_cast<int>(particles[i].type);
        const double eps_i = lj_epsilon[type_i];
        const double sig_i = lj_sigma[type_i];

        for (size_t j = i + 1; j < n; ++j) {
            Vec3 dr = particles[i].position - particles[j].position;

            if (box_size > 0) {
                dr = apply_pbc_vector(dr, box_size);
            }

            const double r_sq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

            if (r_sq > 1e-30 && r_sq < cutoff_sq) {
                const double r = sqrt(r_sq);
                const int type_j = static_cast<int>(particles[j].type);

                const double epsilon = sqrt(eps_i * lj_epsilon[type_j]);
                const double sigma = 0.5 * (sig_i + lj_sigma[type_j]);

                const double sigma_r = sigma / r;
                const double sr6 = pow(sigma_r, 6.0);
                const double sr12 = sr6 * sr6;

                const double force_mag = 24.0 * epsilon / r * (2.0 * sr12 - sr6);
                const Vec3 force_vec = dr * (force_mag / r);

                particles[i].force += force_vec;
                particles[j].force -= force_vec;
            }
        }
    }
}

vector<Vec3> LennardJones::compute_lj_forces(const vector<Particle>& particles, double box_size) const {
    const size_t n = particles.size();
    vector<Vec3> forces(n, Vec3(0.0, 0.0, 0.0));

    const double cutoff_sq = lj_cutoff * lj_cutoff;

    for (size_t i = 0; i < n - 1; ++i) {
        const int type_i = static_cast<int>(particles[i].type);
        const double eps_i = lj_epsilon[type_i];
        const double sig_i = lj_sigma[type_i];

        for (size_t j = i + 1; j < n; ++j) {
            Vec3 dr = particles[i].position - particles[j].position;

            if (box_size > 0) {
                dr = apply_pbc_vector(dr, box_size);
            }

            const double r_sq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

            if (r_sq > 1e-30 && r_sq < cutoff_sq) {
                const double r = sqrt(r_sq);
                const int type_j = static_cast<int>(particles[j].type);

                const double epsilon = sqrt(eps_i * lj_epsilon[type_j]);
                const double sigma = 0.5 * (sig_i + lj_sigma[type_j]);

                const double sigma_r = sigma / r;
                const double sr6 = pow(sigma_r, 6.0);
                const double sr12 = sr6 * sr6;

                const double force_magnitude = 24.0 * epsilon / r * (2.0 * sr12 - sr6);
                const Vec3 force_vector = dr * (force_magnitude / r);

                forces[i] += force_vector;
                forces[j] -= force_vector;
            }
        }
    }

    return forces;
}

double LennardJones::compute_lj_energy(const vector<Particle>& particles, double box_size) const {
    double energy = 0.0;
    const size_t n = particles.size();
    const double cutoff_sq = lj_cutoff * lj_cutoff;

    for (size_t i = 0; i < n - 1; ++i) {
        const int type_i = static_cast<int>(particles[i].type);
        const double eps_i = lj_epsilon[type_i];
        const double sig_i = lj_sigma[type_i];

        for (size_t j = i + 1; j < n; ++j) {
            Vec3 dr = particles[i].position - particles[j].position;

            if (box_size > 0) {
                dr = apply_pbc_vector(dr, box_size);
            }

            const double r_sq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

            if (r_sq > 1e-30 && r_sq < cutoff_sq) {
                const double r = sqrt(r_sq);
                const int type_j = static_cast<int>(particles[j].type);

                const double epsilon = sqrt(eps_i * lj_epsilon[type_j]);
                const double sigma = 0.5 * (sig_i + lj_sigma[type_j]);

                const double sigma_r = sigma / r;
                const double sr6 = pow(sigma_r, 6.0);
                const double sr12 = sr6 * sr6;

                energy += 4.0 * epsilon * (sr12 - sr6);
            }
        }
    }

    return energy;
}