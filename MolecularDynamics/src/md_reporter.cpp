#include "md_reporter.h"
#include <iomanip>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

void MDReporter::ensure_directory(std::string_view path) const {
    try {
        fs::path p{ std::string(path) };
        if (p.has_parent_path()) {
            fs::create_directories(p.parent_path());
        }
    } catch (...) {}
}

void MDReporter::export_topology(const OptimizedMolecularDynamics& engine, const std::string& filename) const {
    ensure_directory(filename);
    std::ofstream file(filename);
    if (!file.is_open()) return;

    const auto& molecules = engine.get_molecules();
    file << "molecule_id,molecule_type,atom_count,bond_count,angle_count\n";
    for (size_t i = 0; i < molecules.size(); ++i) {
        const auto& mol = molecules[i];
        file << i << "," << mol.type << "," << (int)mol.atom_indices.size()
             << "," << (int)mol.bonds.size() << "," << (int)mol.angles.size() << "\n";
    }

    file << "\nBond details (mol_id, atom1, atom2, type)\n";
    for (size_t i = 0; i < molecules.size(); ++i) {
        for (const auto& bond : molecules[i].bonds) {
            file << i << "," << bond.first << "," << bond.second << ",bond\n";
        }
    }

    file << "\nAngle details (mol_id, atom1, atom2, atom3)\n";
    for (size_t i = 0; i < molecules.size(); ++i) {
        for (const auto& angle : molecules[i].angles) {
            file << i << "," << angle[0] << "," << angle[1] << "," << angle[2] << "\n";
        }
    }
}

void MDReporter::export_xyz(const OptimizedMolecularDynamics& engine, const std::string& filename, int step) const {
    ensure_directory(filename);
    std::ofstream file(filename);
    if (!file.is_open()) return;

    const auto& particles = engine.get_particles();
    file << (int)particles.size() << "\nStep " << step << "\n";
    file << std::fixed << std::setprecision(6);

    const double inv_ang = 1.0 / Constants::ANGSTROM_TO_METERS;
    for (const auto& p : particles) {
        file << std::string(Constants::ATOM_SYMBOLS[p.type]) << " " 
             << p.position.x * inv_ang << " " 
             << p.position.y * inv_ang << " " 
             << p.position.z * inv_ang << "\n";
    }
}

void MDReporter::export_energy(const OptimizedMolecularDynamics& engine, const std::string& filename, int step) {
    ensure_directory(filename);
    std::ofstream file(filename, energy_init ? std::ios::app : std::ios::out);
    if (!file.is_open()) return;

    if (!energy_init) {
        file << "step,temp,pot_energy,kin_energy,total_energy\n";
        energy_init = true;
    }

    file << std::fixed << std::setprecision(6);
    double temp = engine.calculate_temperature();
    double kin = engine.calculate_kinetic_energy();
    file << step << "," << temp << ",0.0," << kin << "," << kin << "\n";
}

void MDReporter::export_trajectory(const OptimizedMolecularDynamics& engine, const std::string& filename, int step) {
    ensure_directory(filename);
    std::ofstream file(filename, traj_init ? std::ios::app : std::ios::out);
    if (!file.is_open()) return;

    if (!traj_init) {
        file << "step,atom_id,type,x,y,z\n";
        traj_init = true;
    }

    const auto& particles = engine.get_particles();
    const double inv_ang = 1.0 / Constants::ANGSTROM_TO_METERS;
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        file << step << "," << (int)i << "," << std::string(Constants::ATOM_SYMBOLS[p.type]) << ","
             << p.position.x * inv_ang << "," << p.position.y * inv_ang << "," << p.position.z * inv_ang << "\n";
    }
}

void MDReporter::export_bonds_angles(const OptimizedMolecularDynamics& engine, const std::string& filename, int step) {
    ensure_directory(filename);
    if (!bonds_init) bonds_init = true;
}
