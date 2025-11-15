#include "md_system.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

using namespace std;

void OptimizedMolecularDynamics::export_xyz(const string& filename, int step) const {
    if (particles.empty()) throw runtime_error("No particles to export");

    ofstream file(filename);
    if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
    file << fixed << setprecision(6);

    file << particles.size() << "\n";
    if (step >= 0) {
        file << "Step " << step << " - MD Simulation\n";
    }
    else {
        file << "MD Simulation Export\n";
    }
    static const vector<string> atom_symbols = { "C", "H", "O", "N", "Be", "Cl" };

    const double inv_angstrom = 1.0 / Constants::ANGSTROM_TO_METERS;
    for (const auto& p : particles) {
        if (p.type < 0 || p.type >= static_cast<int>(atom_symbols.size())) {
            throw out_of_range("Invalid atom type: " + to_string(p.type));
        }
        const string& symbol = atom_symbols[p.type];
        const double x = p.position.x * inv_angstrom;
        const double y = p.position.y * inv_angstrom;
        const double z = p.position.z * inv_angstrom;
        file << symbol << " " << x << " " << y << " " << z << "\n";
    }

    file.close();
    cout << "XYZ export saved to: " << filename << "\n";
}

void OptimizedMolecularDynamics::export_csv_trajectory(const string& filename, int step) {
    if (particles.empty()) throw runtime_error("No particles to export");

    if (!trajectory_file_initialized) {
        ofstream file(filename, ios::out);
        if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
        file << "step,atom_id,type,x,y,z,vx,vy,vz,fx,fy,fz,mass\n";
        file.close();
        trajectory_file_initialized = true;
    }

    ofstream file(filename, ios::app);
    if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
    file << fixed << setprecision(6);

    static const vector<string> atom_symbols = { "C", "H", "O", "N", "Be", "Cl" };

    const double inv_angstrom = 1.0 / Constants::ANGSTROM_TO_METERS;
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        if (p.type < 0 || p.type >= static_cast<int>(atom_symbols.size())) {
            throw out_of_range("Invalid atom type: " + to_string(p.type));
        }
        const string& symbol = atom_symbols[p.type];
        const double x = p.position.x * inv_angstrom;
        const double y = p.position.y * inv_angstrom;
        const double z = p.position.z * inv_angstrom;
        file << step << "," << i << "," << symbol << ","
            << x << "," << y << "," << z << ","
            << p.velocity.x << "," << p.velocity.y << "," << p.velocity.z << ","
            << p.force.x << "," << p.force.y << "," << p.force.z << ","
            << p.mass << "\n";
    }

    file.close();
    cout << "Trajectory export saved to: " << filename << "\n";
}

void OptimizedMolecularDynamics::export_energy_data(const string& filename, int step) {
    if (particles.empty() || molecules.empty()) throw runtime_error("No particles or molecules to export energy data");

    if (!energy_file_initialized) {
        ofstream file(filename, ios::out);
        if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
        file << "step,total_energy,temperature,time_fs\n";
        file.close();
        energy_file_initialized = true;
    }

    ofstream file(filename, ios::app);
    if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
    file << fixed << setprecision(6);

    const double energy = integrator.compute_total_energy(particles, box_size, molecules);
    const double temp = calculate_temperature();
    const double time_fs = step * dt * 1e15;

    file << step << "," << energy << "," << temp << "," << time_fs << "\n";

    file.close();
    cout << "Energy data exported to: " << filename << "\n";
}

void OptimizedMolecularDynamics::export_bonds_angles(const string& filename, int step) {
    if (particles.empty() || molecules.empty()) throw runtime_error("No particles or molecules to export bonds and angles");

    if (!bonds_angles_file_initialized) {
        ofstream file(filename, ios::out);
        if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
        file << "step,molecule_id,molecule_type,measurement_type,measurement_id,value,equilibrium,deviation\n";
        file.close();
        bonds_angles_file_initialized = true;
    }

    ofstream file(filename, ios::app);
    if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
    file << fixed << setprecision(6);

    vector<vector<double>> bond_lengths(molecules.size());
    vector<vector<double>> angles(molecules.size());
    for (size_t i = 0; i < molecules.size(); ++i) {
        const Molecule& mol = molecules[i];
        bond_lengths[i].reserve(mol.bonds.size());
        angles[i].reserve(mol.angles.size());

        for (const auto& bond : mol.bonds) {
            const int idx1 = bond.first;
            const int idx2 = bond.second;
            if (idx1 >= static_cast<int>(particles.size()) || idx2 >= static_cast<int>(particles.size())) {
                throw out_of_range("Invalid bond indices: " + to_string(idx1) + ", " + to_string(idx2));
            }
            bond_lengths[i].push_back(calculate_bond_length(particles[idx1], particles[idx2]));
        }

        for (const auto& angle : mol.angles) {
            const int idx1 = angle[0];
            const int idx2 = angle[1];
            const int idx3 = angle[2];
            if (idx1 >= static_cast<int>(particles.size()) || idx2 >= static_cast<int>(particles.size()) ||
                idx3 >= static_cast<int>(particles.size())) {
                throw out_of_range("Invalid angle indices: " + to_string(idx1) + ", " + to_string(idx2) + ", " + to_string(idx3));
            }
            angles[i].push_back(calculate_angle(particles[idx1], particles[idx2], particles[idx3]));
        }
    }

    for (size_t i = 0; i < molecules.size(); ++i) {
        const Molecule& mol = molecules[i];
        for (size_t j = 0; j < mol.bonds.size(); ++j) {
            if (j >= bond_lengths[i].size() || j >= mol.bond_lengths.size()) {
                throw out_of_range("Mismatch in bond data for molecule " + to_string(i));
            }
            const double length = bond_lengths[i][j];
            const double equilibrium = mol.bond_lengths[j];
            const double deviation = abs(length - equilibrium);
            file << step << "," << i << "," << mol.type << ",bond," << j << ","
                << length << "," << equilibrium << "," << deviation << "\n";
        }
        for (size_t j = 0; j < mol.angles.size(); ++j) {
            if (j >= angles[i].size() || j >= mol.angle_degrees.size()) {
                throw out_of_range("Mismatch in angle data for molecule " + to_string(i));
            }
            const double angle = angles[i][j];
            const double equilibrium = mol.angle_degrees[j];
            const double deviation = abs(angle - equilibrium);
            file << step << "," << i << "," << mol.type << ",angle," << j << ","
                << angle << "," << equilibrium << "," << deviation << "\n";
        }
    }

    file.close();
    cout << "Bonds and angles exported to: " << filename << "\n";
}

void OptimizedMolecularDynamics::export_lj_forces(const string& filename, int step) {
    if (particles.empty()) throw runtime_error("No particles to export LJ forces");

    if (!lj_forces_file_initialized) {
        ofstream file(filename, ios::out);
        if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
        file << "step,atom_id,type,fx_lj,fy_lj,fz_lj\n";
        file.close();
        lj_forces_file_initialized = true;
    }

    ofstream file(filename, ios::app);
    if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
    file << fixed << setprecision(6);

    vector<Vec3> lj_forces = integrator.lj_calculator.compute_lj_forces(particles, box_size);

    static const vector<string> atom_symbols = { "C", "H", "O", "N", "Be", "Cl" };

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        if (p.type < 0 || p.type >= static_cast<int>(atom_symbols.size())) {
            throw out_of_range("Invalid atom type: " + to_string(p.type));
        }
        const string& symbol = atom_symbols[p.type];
        file << step << "," << i << "," << symbol << ","
            << lj_forces[i].x << "," << lj_forces[i].y << "," << lj_forces[i].z << "\n";
    }

    file.close();
    cout << "LJ forces exported to: " << filename << "\n";
}

void OptimizedMolecularDynamics::export_molecule_topology(const string& filename) const {
    if (molecules.empty()) throw runtime_error("No molecules to export topology");

    ofstream file(filename);
    if (!file.is_open()) throw runtime_error("Cannot open file for writing: " + filename);
    file << fixed << setprecision(6);

    file << "# Molecule topology - bond and angle information\n";
    file << "molecule_id,molecule_type,atom_count,bond_count,angle_count\n";

    for (size_t i = 0; i < molecules.size(); ++i) {
        const Molecule& mol = molecules[i];
        file << i << "," << mol.type << "," << mol.atom_indices.size()
            << "," << mol.bonds.size() << "," << mol.angles.size() << "\n";
    }

    file << "\n# Bond details\n";
    file << "molecule_id,bond_id,atom1,atom2,equilibrium_length\n";
    for (size_t i = 0; i < molecules.size(); ++i) {
        const Molecule& mol = molecules[i];
        for (size_t j = 0; j < mol.bonds.size(); ++j) {
            if (j >= mol.bond_lengths.size()) {
                throw out_of_range("Mismatch in bond lengths for molecule " + to_string(i));
            }
            file << i << "," << j << "," << mol.bonds[j].first << ","
                << mol.bonds[j].second << "," << mol.bond_lengths[j] << "\n";
        }
    }

    file << "\n# Angle details\n";
    file << "molecule_id,angle_id,atom1,atom2,atom3,equilibrium_angle\n";
    for (size_t i = 0; i < molecules.size(); ++i) {
        const Molecule& mol = molecules[i];
        for (size_t j = 0; j < mol.angles.size(); ++j) {
            if (j >= mol.angle_degrees.size()) {
                throw out_of_range("Mismatch in angle degrees for molecule " + to_string(i));
            }
            file << i << "," << j << "," << mol.angles[j][0] << ","
                << mol.angles[j][1] << "," << mol.angles[j][2] << ","
                << mol.angle_degrees[j] << "\n";
        }
    }

    file.close();
    cout << "Molecule topology exported to: " << filename << "\n";
}