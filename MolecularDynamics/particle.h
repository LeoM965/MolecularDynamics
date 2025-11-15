#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec3.h"
#include "constants.h"
#include <vector>
#include <string>
#include <array>

// Enum for atom types
enum AtomType { C = 0, H = 1, O = 2, N = 3, Be = 4, Cl = 5 };

// Optimized particle structure
struct Particle {
    Vec3 position, velocity, force;
    double mass;
    AtomType type;

    Particle(double x = 0.0, double y = 0.0, double z = 0.0,
        double mass_ = 1.0, AtomType type_ = C)
        : position(x* Constants::ANGSTROM_TO_METERS,
            y* Constants::ANGSTROM_TO_METERS,
            z* Constants::ANGSTROM_TO_METERS),
        velocity(0.0, 0.0, 0.0), force(0.0, 0.0, 0.0),
        mass(mass_), type(type_) {
    }
};

// Molecule definitions
struct MoleculeDefinition {
    std::string type;
    std::vector<std::pair<AtomType, double>> atoms;
    std::vector<Vec3> positions;
    std::vector<std::pair<int, int>> bonds;
    std::vector<std::array<int, 3>> angles;
    std::vector<double> bond_lengths;
    std::vector<double> angle_degrees;
};

// Molecule in system
struct Molecule {
    std::string type;
    std::vector<int> atom_indices;
    std::vector<std::pair<int, int>> bonds;
    std::vector<std::array<int, 3>> angles;
    std::vector<double> bond_lengths;
    std::vector<double> angle_degrees;
};

#endif // PARTICLE_H