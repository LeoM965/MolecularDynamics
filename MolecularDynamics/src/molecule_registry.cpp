#include "molecule_registry.h"

using namespace std;

const vector<pair<AtomType, double>> MoleculeRegistry::h_atoms = { {H, 1.0}, {H, 1.0}, {H, 1.0}, {H, 1.0} };
const vector<pair<AtomType, double>> MoleculeRegistry::h2_atoms = { {H, 1.0}, {H, 1.0} };
const vector<pair<AtomType, double>> MoleculeRegistry::co2_atoms = { {O, 16.0}, {C, 12.0}, {O, 16.0} };
const vector<pair<AtomType, double>> MoleculeRegistry::acetylene_atoms = { {H, 1.0}, {C, 12.0}, {C, 12.0}, {H, 1.0} };
const vector<pair<AtomType, double>> MoleculeRegistry::h3_atoms = { {H, 1.0}, {H, 1.0}, {H, 1.0} };
const vector<pair<AtomType, double>> MoleculeRegistry::hcl_atoms = { {H, 1.0}, {Cl, 35.5} };

void MoleculeRegistry::add_common_atoms(MoleculeDefinition& def, const string& type, AtomType central_type,
    double central_mass, const vector<pair<AtomType, double>>& ligands) {
    if (central_mass <= 0) throw invalid_argument("Central mass must be positive");
    def.type = type;
    def.atoms.reserve(ligands.size() + 1);
    def.atoms.push_back({ central_type, central_mass });
    def.atoms.insert(def.atoms.end(), ligands.begin(), ligands.end());
}

MoleculeDefinition MoleculeRegistry::create_tetrahedral(const string& type, AtomType central_type,
    double central_mass, const vector<pair<AtomType, double>>& ligands) {
    if (ligands.size() != 4) throw invalid_argument("Tetrahedral molecule requires exactly 4 ligands");
    MoleculeDefinition def;
    this->add_common_atoms(def, type, central_type, central_mass, ligands);

    def.positions.reserve(5);
    def.positions.push_back(Vec3(0.0, 0.0, 0.0));
    def.positions.push_back(Vec3(Constants::TETRAHEDRAL_S, Constants::TETRAHEDRAL_S, Constants::TETRAHEDRAL_S));
    def.positions.push_back(Vec3(Constants::TETRAHEDRAL_S, -Constants::TETRAHEDRAL_S, -Constants::TETRAHEDRAL_S));
    def.positions.push_back(Vec3(-Constants::TETRAHEDRAL_S, Constants::TETRAHEDRAL_S, -Constants::TETRAHEDRAL_S));
    def.positions.push_back(Vec3(-Constants::TETRAHEDRAL_S, -Constants::TETRAHEDRAL_S, Constants::TETRAHEDRAL_S));

    def.bonds.reserve(4);
    def.bond_lengths.reserve(4);
    for (int i = 1; i <= 4; ++i) {
        def.bonds.push_back({ 0, i });
        def.bond_lengths.push_back(Constants::TETRAHEDRAL_BOND_LENGTH);
    }

    def.angles = { {1, 0, 2}, {1, 0, 3}, {1, 0, 4}, {2, 0, 3}, {2, 0, 4}, {3, 0, 4} };
    def.angle_degrees = { 109.5, 109.5, 109.5, 109.5, 109.5, 109.5 };

    return def;
}

MoleculeDefinition MoleculeRegistry::create_bent(const string& type, AtomType central_type,
    double central_mass, const vector<pair<AtomType, double>>& ligands) {
    if (ligands.size() != 2) throw invalid_argument("Bent molecule requires exactly 2 ligands");
    MoleculeDefinition def;
    this->add_common_atoms(def, type, central_type, central_mass, ligands);

    def.positions.reserve(3);
    def.positions.push_back(Vec3(0.0, 0.0, 0.0));
    def.positions.push_back(Vec3(Constants::WATER_DX, Constants::WATER_DY, 0.0));
    def.positions.push_back(Vec3(-Constants::WATER_DX, Constants::WATER_DY, 0.0));

    def.bonds = { {0, 1}, {0, 2} };
    def.bond_lengths = { Constants::WATER_BOND_LENGTH, Constants::WATER_BOND_LENGTH };

    def.angles = { {1, 0, 2} };
    def.angle_degrees = { Constants::WATER_ANGLE };

    return def;
}

MoleculeDefinition MoleculeRegistry::create_linear(const string& type, const vector<pair<AtomType, double>>& atoms,
    const vector<double>& bond_lengths) {
    if (atoms.size() < 2) throw invalid_argument("Linear molecule requires at least 2 atoms");
    if (bond_lengths.size() != atoms.size() - 1) throw invalid_argument("Invalid number of bond lengths");
    MoleculeDefinition def;
    def.type = type;
    def.atoms = atoms;

    double total_length = 0.0;
    vector<double> offsets;
    offsets.reserve(atoms.size());
    offsets.push_back(0.0);
    for (double len : bond_lengths) {
        if (len <= 0) throw invalid_argument("Bond length must be positive");
        total_length += len;
        offsets.push_back(total_length);
    }

    double center = total_length * 0.5;
    def.positions.reserve(atoms.size());
    for (size_t i = 0; i < atoms.size(); ++i) {
        def.positions.push_back(Vec3(offsets[i] - center, 0.0, 0.0));
    }

    def.bonds.reserve(atoms.size() - 1);
    def.bond_lengths = bond_lengths;
    for (int i = 0; i < static_cast<int>(atoms.size()) - 1; ++i) { 
        def.bonds.push_back({ i, i + 1 });
    }

    def.angles.reserve(max(0, static_cast<int>(atoms.size()) - 2));
    def.angle_degrees.reserve(max(0, static_cast<int>(atoms.size()) - 2));
    for (int i = 0; i < static_cast<int>(atoms.size()) - 2; ++i) { 
        def.angles.push_back({ i, i + 1, i + 2 });
        def.angle_degrees.push_back(180.0);
    }

    return def;
}

MoleculeDefinition MoleculeRegistry::create_trigonal_pyramidal(const string& type, AtomType central_type,
    double central_mass, const vector<pair<AtomType, double>>& ligands) {
    if (ligands.size() != 3) throw invalid_argument("Trigonal pyramidal molecule requires exactly 3 ligands");
    MoleculeDefinition def;
    this->add_common_atoms(def, type, central_type, central_mass, ligands);

    double h_height = -Constants::NH3_BOND_LENGTH * Constants::NH3_HEIGHT_FACTOR;
    double h_radius = Constants::NH3_BOND_LENGTH * Constants::NH3_RADIUS_FACTOR;

    def.positions.reserve(4);
    def.positions.push_back(Vec3(0.0, 0.0, 0.0));
    for (int i = 0; i < 3; ++i) {
        double theta = i * 2.0 * Constants::PI / 3.0;
        def.positions.push_back(Vec3(h_radius * cos(theta), h_radius * sin(theta), h_height));
    }

    def.bonds = { {0, 1}, {0, 2}, {0, 3} };
    def.bond_lengths = { Constants::NH3_BOND_LENGTH, Constants::NH3_BOND_LENGTH, Constants::NH3_BOND_LENGTH };

    def.angles = { {1, 0, 2}, {1, 0, 3}, {2, 0, 3} };
    def.angle_degrees = { Constants::NH3_ANGLE, Constants::NH3_ANGLE, Constants::NH3_ANGLE };

    return def;
}

MoleculeRegistry::MoleculeRegistry() {
    this->definitions["methane"] = this->create_tetrahedral("methane", C, 12.0, h_atoms);

    this->definitions["water"] = this->create_bent("water", O, 16.0, h2_atoms);

    vector<double> co2_bonds = { Constants::CO2_BOND_LENGTH, Constants::CO2_BOND_LENGTH };
    this->definitions["co2"] = this->create_linear("co2", co2_atoms, co2_bonds);

    vector<double> acetylene_bonds = { Constants::ACETYLENE_BOND_LENGTHS[0], Constants::ACETYLENE_BOND_LENGTHS[1], Constants::ACETYLENE_BOND_LENGTHS[2] };
    this->definitions["acetylene"] = this->create_linear("acetylene", acetylene_atoms, acetylene_bonds);

    this->definitions["ammonia"] = this->create_trigonal_pyramidal("ammonia", N, 14.0, h3_atoms);

    vector<double> hcl_bonds = { Constants::HCL_BOND_LENGTH };
    this->definitions["hcl"] = this->create_linear("hcl", hcl_atoms, hcl_bonds);
}

const MoleculeDefinition& MoleculeRegistry::get_definition(const string& type) const {
    auto it = this->definitions.find(type);
    if (it == this->definitions.end()) {
        throw runtime_error("Unknown molecule type: " + type);
    }
    return it->second;
}