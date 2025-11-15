#ifndef LENNARD_JONES_H
#define LENNARD_JONES_H
#include "particle.h"
#include <vector>

class LennardJones {
private:
    int get_atom_type_index(AtomType atom_type) const;

public:
    std::vector<double> lj_epsilon;
    std::vector<double> lj_sigma;
    double lj_cutoff;

    LennardJones();
    Vec3 apply_pbc_vector(const Vec3& v, double box_size) const;
    void calculate_lj_forces(std::vector<Particle>& particles, double box_size);
    std::vector<Vec3> compute_lj_forces(const std::vector<Particle>& particles, double box_size) const;
    double compute_lj_energy(const std::vector<Particle>& particles, double box_size) const;
};
#endif
