#ifndef LENNARD_JONES_H
#define LENNARD_JONES_H
#include "particle.h"
#include <vector>

class LennardJones {
public:
    std::vector<double> lj_epsilon;
    std::vector<double> lj_sigma;
    double lj_cutoff;

    LennardJones();
    void calculate_lj_forces(std::vector<Particle>& particles, double box_size);
    std::vector<Vec3> compute_lj_forces(const std::vector<Particle>& particles, double box_size) const;
    double compute_lj_energy(const std::vector<Particle>& particles, double box_size) const;
};
#endif
