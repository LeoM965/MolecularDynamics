#ifndef MOLECULE_REGISTRY_H
#define MOLECULE_REGISTRY_H

#include "particle.h"
#include <unordered_map>
#include <stdexcept>
#include <vector>
#include <array>
#include <string>
#include <cmath>

using namespace std;

class MoleculeRegistry {
private:
    unordered_map<string, MoleculeDefinition> definitions;

    static const vector<pair<AtomType, double>> h_atoms;
    static const vector<pair<AtomType, double>> h2_atoms;
    static const vector<pair<AtomType, double>> co2_atoms;
    static const vector<pair<AtomType, double>> acetylene_atoms;
    static const vector<pair<AtomType, double>> h3_atoms;
    static const vector<pair<AtomType, double>> hcl_atoms;

    void add_common_atoms(MoleculeDefinition& def, const string& type, AtomType central_type,
        double central_mass, const vector<pair<AtomType, double>>& ligands);

    MoleculeDefinition create_tetrahedral(const string& type, AtomType central_type,
        double central_mass, const vector<pair<AtomType, double>>& ligands);

    MoleculeDefinition create_bent(const string& type, AtomType central_type,
        double central_mass, const vector<pair<AtomType, double>>& ligands);

    MoleculeDefinition create_linear(const string& type, const vector<pair<AtomType, double>>& atoms,
        const vector<double>& bond_lengths);

    MoleculeDefinition create_trigonal_pyramidal(const string& type, AtomType central_type,
        double central_mass, const vector<pair<AtomType, double>>& ligands);


public:
    MoleculeRegistry();
    const MoleculeDefinition& get_definition(const string& type) const;

};

#endif // MOLECULE_REGISTRY_H