//#include <vector>
//#include <array>
//#include <string>
//#include <cmath>
//#include <iostream>
//#include <iomanip>
//#include <unordered_map>
//#include <algorithm>
//#include <random>
//#include <memory>
//#include <fstream>
//
//using namespace std;
//
//// Constante fizice optimizate
//namespace Constants {
//    const double ANGSTROM_TO_METERS = 1e-10;
//    const double KCAL_PER_MOL_TO_JOULES = 4184.0;
//    const double AVOGADRO = 6.022e23;
//    const double ATOMIC_MASS_UNIT = 1.661e-27;
//    const double BOLTZMANN = 1.381e-23;
//    const double PI = 3.141592653589793;
//    // Constante de forta corectate pentru energie realista
//    const double BOND_FORCE_CONST = 400.0 * KCAL_PER_MOL_TO_JOULES / (ANGSTROM_TO_METERS * ANGSTROM_TO_METERS);
//    const double ANGLE_FORCE_CONST = 50.0 * KCAL_PER_MOL_TO_JOULES;
//    const double MAX_FORCE = 1.5e-11;
//}
//
//// Structura pentru vectori 3D optimizata
//struct Vec3 {
//    double x, y, z;
//
//    Vec3(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0) : x(x_), y(y_), z(z_) {}
//
//    Vec3& operator+=(const Vec3& other) {
//        x += other.x; y += other.y; z += other.z;
//        return *this;
//    }
//
//    Vec3& operator-=(const Vec3& other) {
//        x -= other.x; y -= other.y; z -= other.z;
//        return *this;
//    }
//
//    Vec3& operator*=(double scalar) {
//        x *= scalar; y *= scalar; z *= scalar;
//        return *this;
//    }
//
//    Vec3 operator-(const Vec3& other) const {
//        return Vec3(x - other.x, y - other.y, z - other.z);
//    }
//
//    Vec3 operator+(const Vec3& other) const {
//        return Vec3(x + other.x, y + other.y, z + other.z);
//    }
//
//    Vec3 operator*(double scalar) const {
//        return Vec3(x * scalar, y * scalar, z * scalar);
//    }
//
//    double magnitude() const {
//        return sqrt(x * x + y * y + z * z);
//    }
//
//    double magnitude_squared() const {
//        return x * x + y * y + z * z;
//    }
//
//    double dot(const Vec3& other) const {
//        return x * other.x + y * other.y + z * other.z;
//    }
//
//    Vec3& normalize() {
//        double mag = magnitude();
//        if (mag > 1e-12) {
//            *this *= (1.0 / mag);
//        }
//        return *this;
//    }
//};
//
//// Enum pentru tipuri de atomi
//enum AtomType { C = 0, H = 1, O = 2, N = 3, Be = 4, Cl = 5 };
//
//// Structura particule optimizata
//struct Particle {
//    Vec3 position, velocity, force;
//    double mass;
//    AtomType type;
//
//    Particle(double x = 0.0, double y = 0.0, double z = 0.0,
//        double mass_ = 1.0, AtomType type_ = C)
//        : position(x* Constants::ANGSTROM_TO_METERS,
//            y* Constants::ANGSTROM_TO_METERS,
//            z* Constants::ANGSTROM_TO_METERS),
//        velocity(0.0, 0.0, 0.0), force(0.0, 0.0, 0.0),
//        mass(mass_), type(type_) {
//    }
//};
//
//// Definitiile moleculelor
//struct MoleculeDefinition {
//    string type;
//    vector<pair<AtomType, double>> atoms;
//    vector<Vec3> positions;
//    vector<pair<int, int>> bonds;
//    vector<array<int, 3>> angles;
//    vector<double> bond_lengths;
//    vector<double> angle_degrees;
//};
//
//// Molecula in sistem
//struct Molecule {
//    string type;
//    vector<int> atom_indices;
//    vector<pair<int, int>> bonds;
//    vector<array<int, 3>> angles;
//    vector<double> bond_lengths;
//    vector<double> angle_degrees;
//};
//
//// Registry optimizat
//class MoleculeRegistry {
//private:
//    unordered_map<string, MoleculeDefinition> definitions;
//
//    MoleculeDefinition create_tetrahedral(const string& type, AtomType central_type,
//        double central_mass, const vector<pair<AtomType, double>>& ligands,
//        double bond_length) {
//        MoleculeDefinition def;
//        def.type = type;
//        def.atoms.push_back(make_pair(central_type, central_mass));
//        def.atoms.insert(def.atoms.end(), ligands.begin(), ligands.end());
//
//        // Coordonate tetraedrice corecte
//        double s = bond_length / sqrt(3.0);
//        def.positions.push_back(Vec3(0.0, 0.0, 0.0));
//        def.positions.push_back(Vec3(s, s, s));
//        def.positions.push_back(Vec3(s, -s, -s));
//        def.positions.push_back(Vec3(-s, s, -s));
//        def.positions.push_back(Vec3(-s, -s, s));
//
//        // Legaturi si unghiuri
//        for (int i = 1; i <= (int)ligands.size(); ++i) {
//            def.bonds.push_back(make_pair(0, i));
//            def.bond_lengths.push_back(bond_length);
//        }
//
//        // Unghiuri tetraedrice
//        for (int i = 1; i <= (int)ligands.size(); ++i) {
//            for (int j = i + 1; j <= (int)ligands.size(); ++j) {
//                array<int, 3> angle = { i, 0, j };
//                def.angles.push_back(angle);
//                def.angle_degrees.push_back(109.5);
//            }
//        }
//        return def;
//    }
//
//    MoleculeDefinition create_bent(const string& type, AtomType central_type, double central_mass,
//        const vector<pair<AtomType, double>>& ligands,
//        double bond_length, double angle) {
//        MoleculeDefinition def;
//        def.type = type;
//        def.atoms.push_back(make_pair(central_type, central_mass));
//        def.atoms.insert(def.atoms.end(), ligands.begin(), ligands.end());
//
//        double half_angle = angle * Constants::PI / 360.0;
//        double dx = bond_length * sin(half_angle);
//        double dy = bond_length * cos(half_angle);
//
//        def.positions.push_back(Vec3(0.0, 0.0, 0.0));
//        def.positions.push_back(Vec3(dx, dy, 0.0));
//        def.positions.push_back(Vec3(-dx, dy, 0.0));
//
//        def.bonds.push_back(make_pair(0, 1));
//        def.bonds.push_back(make_pair(0, 2));
//        def.bond_lengths.push_back(bond_length);
//        def.bond_lengths.push_back(bond_length);
//
//        array<int, 3> angle_def = { 1, 0, 2 };
//        def.angles.push_back(angle_def);
//        def.angle_degrees.push_back(angle);
//
//        return def;
//    }
//
//    MoleculeDefinition create_linear(const string& type,
//        const vector<pair<AtomType, double>>& atoms,
//        const vector<double>& bond_lengths) {
//        MoleculeDefinition def;
//        def.type = type;
//        def.atoms = atoms;
//
//        // Calcul pozitii optimizat
//        double total_length = 0.0;
//        vector<double> offsets;
//        offsets.push_back(0.0);
//
//        for (double len : bond_lengths) {
//            total_length += len;
//            offsets.push_back(total_length);
//        }
//
//        double center = total_length * 0.5;
//        for (size_t i = 0; i < atoms.size(); ++i) {
//            def.positions.push_back(Vec3(offsets[i] - center, 0.0, 0.0));
//        }
//
//        // Legaturi si unghiuri
//        for (int i = 0; i < (int)atoms.size() - 1; ++i) {
//            def.bonds.push_back(make_pair(i, i + 1));
//            if (i < (int)atoms.size() - 2) {
//                array<int, 3> angle = { i, i + 1, i + 2 };
//                def.angles.push_back(angle);
//                def.angle_degrees.push_back(180.0);
//            }
//        }
//        def.bond_lengths = bond_lengths;
//
//        return def;
//    }
//
//    // Geometrie trigonal piramidala CORECTATA pentru NH3 - unghi exact 107 grade
//    MoleculeDefinition create_trigonal_pyramidal(const string& type, AtomType central_type,
//        double central_mass, const vector<pair<AtomType, double>>& ligands,
//        double bond_length, double target_angle) {
//        MoleculeDefinition def;
//        def.type = type;
//        def.atoms.push_back(make_pair(central_type, central_mass));
//        def.atoms.insert(def.atoms.end(), ligands.begin(), ligands.end());
//
//        // Pentru unghi exact de 107 grade, calculez pozitiile geometric
//        double angle_rad = target_angle * Constants::PI / 180.0;
//        double cos_angle = cos(angle_rad);
//
//        // Calculez inaltimea piramidei pentru unghi exact
//        // Folosind formula pentru unghi H-N-H in geometria trigonal piramidala
//        double height_factor = sqrt((1.0 + 2.0 * cos_angle) / 3.0);
//        double radius_factor = sqrt((2.0 * (1.0 - cos_angle)) / 3.0);
//
//        // Pozitia N la origine
//        def.positions.push_back(Vec3(0.0, 0.0, 0.0));
//
//        // Calculez pozitiile H pentru unghi exact de 107 grade
//        double h_height = -bond_length * height_factor;  // sub planul bazei
//        double h_radius = bond_length * radius_factor;
//
//        // Pozitiile H in plan orizontal, la 120 grade unul de altul
//        for (int i = 0; i < 3; ++i) {
//            double theta = i * 2.0 * Constants::PI / 3.0;
//            def.positions.push_back(Vec3(
//                h_radius * cos(theta),
//                h_radius * sin(theta),
//                h_height
//            ));
//        }
//
//        // Legaturi N-H
//        for (int i = 1; i <= 3; ++i) {
//            def.bonds.push_back(make_pair(0, i));
//            def.bond_lengths.push_back(bond_length);
//        }
//
//        // Unghiuri H-N-H (toate ar trebui sa fie exact target_angle)
//        for (int i = 1; i <= 3; ++i) {
//            for (int j = i + 1; j <= 3; ++j) {
//                array<int, 3> angle_def = { i, 0, j };
//                def.angles.push_back(angle_def);
//                def.angle_degrees.push_back(target_angle);
//            }
//        }
//
//        return def;
//    }
//
//public:
//    MoleculeRegistry() {
//        vector<pair<AtomType, double>> h_atoms = { {H, 1.0}, {H, 1.0}, {H, 1.0}, {H, 1.0} };
//        definitions["methane"] = create_tetrahedral("methane", C, 12.0, h_atoms, 1.09);
//
//        vector<pair<AtomType, double>> h2_atoms = { {H, 1.0}, {H, 1.0} };
//        definitions["water"] = create_bent("water", O, 16.0, h2_atoms, 0.96, 104.5);
//
//        vector<pair<AtomType, double>> co2_atoms = { {O, 16.0}, {C, 12.0}, {O, 16.0} };
//        vector<double> co2_bonds = { 1.16, 1.16 };
//        definitions["co2"] = create_linear("co2", co2_atoms, co2_bonds);
//
//        vector<pair<AtomType, double>> acetylene_atoms = { {H, 1.0}, {C, 12.0}, {C, 12.0}, {H, 1.0} };
//        vector<double> acetylene_bonds = { 1.06, 1.20, 1.06 };
//        definitions["acetylene"] = create_linear("acetylene", acetylene_atoms, acetylene_bonds);
//
//        // NH3 cu geometrie COMPLET CORECTATA pentru unghi exact de 107 grade
//        vector<pair<AtomType, double>> h3_atoms = { {H, 1.0}, {H, 1.0}, {H, 1.0} };
//        definitions["ammonia"] = create_trigonal_pyramidal("ammonia", N, 14.0, h3_atoms, 1.01, 107.0);
//
//        vector<pair<AtomType, double>> hcl_atoms = { {H, 1.0}, {Cl, 35.5} };
//        vector<double> hcl_bonds = { 1.27 };
//        definitions["hcl"] = create_linear("hcl", hcl_atoms, hcl_bonds);
//    }
//
//    const MoleculeDefinition& get_definition(const string& type) const {
//        auto it = definitions.find(type);
//        if (it == definitions.end()) {
//            throw runtime_error("Tip molecula necunoscut: " + type);
//        }
//        return it->second;
//    }
//};
//
//// Integrator optimizat cu stabilitate imbunatatita
//class OptimizedVerletIntegrator {
//private:
//    Vec3 apply_pbc_vector(const Vec3& v, double box_size) const {
//        return Vec3(apply_pbc(v.x, box_size), apply_pbc(v.y, box_size), apply_pbc(v.z, box_size));
//    }
//
//    void calculate_bond_forces(vector<Particle>& particles, double box_size,
//        const vector<Molecule>& molecules) {
//        for (const auto& mol : molecules) {
//            for (size_t i = 0; i < mol.bonds.size(); ++i) {
//                int idx1 = mol.bonds[i].first;
//                int idx2 = mol.bonds[i].second;
//                double r0 = mol.bond_lengths[i] * Constants::ANGSTROM_TO_METERS;
//
//                Vec3 dr = apply_pbc_vector(particles[idx1].position - particles[idx2].position, box_size);
//                double r = dr.magnitude();
//
//                if (r > 1e-12) {
//                    double f_magnitude = -Constants::BOND_FORCE_CONST * (r - r0) / r;
//                    Vec3 force = dr * f_magnitude;
//
//                    particles[idx1].force += force;
//                    particles[idx2].force -= force;
//                }
//            }
//        }
//    }
//
//    void calculate_angle_forces(vector<Particle>& particles, double box_size,
//        const vector<Molecule>& molecules) {
//        for (const auto& mol : molecules) {
//            for (size_t i = 0; i < mol.angles.size(); ++i) {
//                int idx1 = mol.angles[i][0];
//                int idx2 = mol.angles[i][1];
//                int idx3 = mol.angles[i][2];
//                double theta0 = mol.angle_degrees[i] * Constants::PI / 180.0;
//
//                Vec3 a = apply_pbc_vector(particles[idx1].position - particles[idx2].position, box_size);
//                Vec3 b = apply_pbc_vector(particles[idx3].position - particles[idx2].position, box_size);
//
//                double ra = a.magnitude();
//                double rb = b.magnitude();
//
//                if (ra > 1e-12 && rb > 1e-12) {
//                    double cos_theta = max(-1.0, min(1.0, a.dot(b) / (ra * rb)));
//                    double theta = acos(cos_theta);
//                    double dtheta = theta - theta0;
//                    double sin_theta = sin(theta);
//
//                    if (abs(sin_theta) > 1e-6) {
//                        double f = -Constants::ANGLE_FORCE_CONST * dtheta;
//                        double fa = f / (ra * sin_theta);
//                        double fb = f / (rb * sin_theta);
//
//                        a.normalize();
//                        b.normalize();
//
//                        Vec3 force1 = (b - a * cos_theta) * fa;
//                        Vec3 force3 = (a - b * cos_theta) * fb;
//
//                        particles[idx1].force += force1;
//                        particles[idx3].force += force3;
//                        particles[idx2].force -= (force1 + force3);
//                    }
//                }
//            }
//        }
//    }
//
//    void limit_forces(vector<Particle>& particles) {
//        for (auto& p : particles) {
//            double f_mag = p.force.magnitude();
//            if (f_mag > Constants::MAX_FORCE) {
//                p.force *= (Constants::MAX_FORCE / f_mag);
//            }
//        }
//    }
//
//public:
//    double apply_pbc(double x, double box_size) const {
//        return x - box_size * round(x / box_size);
//    }
//
//    void update_positions(vector<Particle>& particles, double box_size, double dt,
//        const vector<Molecule>& molecules) {
//        // Reset forces
//        for (auto& p : particles) {
//            p.force = Vec3(0.0, 0.0, 0.0);
//        }
//
//        // Calculate forces
//        calculate_bond_forces(particles, box_size, molecules);
//        calculate_angle_forces(particles, box_size, molecules);
//        limit_forces(particles);
//
//        // First half of velocity-Verlet
//        for (auto& p : particles) {
//            double inv_mass = 1.0 / (p.mass * Constants::ATOMIC_MASS_UNIT);
//            Vec3 acceleration = p.force * inv_mass;
//
//            p.position += p.velocity * dt + acceleration * (0.5 * dt * dt);
//            p.velocity += acceleration * (0.5 * dt);
//
//            // Apply PBC
//            p.position.x = p.position.x - box_size * floor(p.position.x / box_size);
//            p.position.y = p.position.y - box_size * floor(p.position.y / box_size);
//            p.position.z = p.position.z - box_size * floor(p.position.z / box_size);
//        }
//
//        // Recalculate forces
//        for (auto& p : particles) {
//            p.force = Vec3(0.0, 0.0, 0.0);
//        }
//        calculate_bond_forces(particles, box_size, molecules);
//        calculate_angle_forces(particles, box_size, molecules);
//        limit_forces(particles);
//
//        // Second half of velocity-Verlet
//        for (auto& p : particles) {
//            double inv_mass = 1.0 / (p.mass * Constants::ATOMIC_MASS_UNIT);
//            p.velocity += p.force * (inv_mass * 0.5 * dt);
//        }
//    }
//
//    // Thermostat imbunatatit pentru a atinge temperatura tinta
//    void apply_berendsen_thermostat(vector<Particle>& particles, double target_temp,
//        double tau, double dt) {
//        double current_temp = calculate_temperature(particles);
//        if (current_temp > 0) {
//            double scaling_factor = sqrt(1.0 + (dt / tau) * (target_temp / current_temp - 1.0));
//            // Limitam scalarea pentru stabilitate
//            scaling_factor = max(0.95, min(1.05, scaling_factor));
//
//            for (auto& p : particles) {
//                p.velocity *= scaling_factor;
//            }
//        }
//    }
//
//    double calculate_temperature(const vector<Particle>& particles) const {
//        double total_kinetic = 0.0;
//
//        for (const auto& p : particles) {
//            double m = p.mass * Constants::ATOMIC_MASS_UNIT;
//            total_kinetic += 0.5 * m * p.velocity.magnitude_squared();
//        }
//
//        int dof = max(1, (int)(particles.size() * 3 - 6));
//        return (2.0 * total_kinetic) / (dof * Constants::BOLTZMANN);
//    }
//
//    double compute_total_energy(const vector<Particle>& particles, double box_size,
//        const vector<Molecule>& molecules) const {
//        double potential_energy = 0.0;
//        double kinetic_energy = 0.0;
//
//        // Bond potential energy - corectata
//        for (const auto& mol : molecules) {
//            for (size_t i = 0; i < mol.bonds.size(); ++i) {
//                int idx1 = mol.bonds[i].first;
//                int idx2 = mol.bonds[i].second;
//                double r0 = mol.bond_lengths[i] * Constants::ANGSTROM_TO_METERS;
//
//                Vec3 dr = apply_pbc_vector(particles[idx1].position - particles[idx2].position, box_size);
//                double r = dr.magnitude();
//                double dr_val = r - r0;
//
//                // Folosesc aceeasi constanta ca in calculul fortelor
//                potential_energy += 0.5 * Constants::BOND_FORCE_CONST * dr_val * dr_val;
//            }
//
//            // Angle potential energy - corectata
//            for (size_t i = 0; i < mol.angles.size(); ++i) {
//                int idx1 = mol.angles[i][0];
//                int idx2 = mol.angles[i][1];
//                int idx3 = mol.angles[i][2];
//                double theta0 = mol.angle_degrees[i] * Constants::PI / 180.0;
//
//                Vec3 a = apply_pbc_vector(particles[idx1].position - particles[idx2].position, box_size);
//                Vec3 b = apply_pbc_vector(particles[idx3].position - particles[idx2].position, box_size);
//
//                double ra = a.magnitude();
//                double rb = b.magnitude();
//
//                if (ra > 1e-12 && rb > 1e-12) {
//                    double cos_theta = max(-1.0, min(1.0, a.dot(b) / (ra * rb)));
//                    double theta = acos(cos_theta);
//                    double dtheta = theta - theta0;
//
//                    // Folosesc aceeasi constanta ca in calculul fortelor
//                    potential_energy += 0.5 * Constants::ANGLE_FORCE_CONST * dtheta * dtheta;
//                }
//            }
//        }
//
//        // Kinetic energy
//        for (const auto& p : particles) {
//            double m = p.mass * Constants::ATOMIC_MASS_UNIT;
//            kinetic_energy += 0.5 * m * p.velocity.magnitude_squared();
//        }
//
//        return (potential_energy + kinetic_energy) / Constants::KCAL_PER_MOL_TO_JOULES;
//    }
//};
//
//// Clasa principala optimizata cu exporturi CORECTATE
//class OptimizedMolecularDynamics {
//private:
//    vector<Particle> particles;
//    vector<Molecule> molecules;
//    double box_size;
//    double dt;
//    double target_temperature;
//    OptimizedVerletIntegrator integrator;
//    MoleculeRegistry registry;
//
//    // Variabile pentru a verifica daca fisierele au fost initializate
//    bool energy_file_initialized;
//    bool trajectory_file_initialized;
//    bool bonds_angles_file_initialized;
//
//    double calculate_bond_length(const Particle& p1, const Particle& p2) const {
//        Vec3 dr(
//            integrator.apply_pbc(p1.position.x - p2.position.x, box_size),
//            integrator.apply_pbc(p1.position.y - p2.position.y, box_size),
//            integrator.apply_pbc(p1.position.z - p2.position.z, box_size)
//        );
//        return dr.magnitude() / Constants::ANGSTROM_TO_METERS;
//    }
//
//    double calculate_angle(const Particle& p1, const Particle& p2, const Particle& p3) const {
//        Vec3 a(
//            integrator.apply_pbc(p1.position.x - p2.position.x, box_size),
//            integrator.apply_pbc(p1.position.y - p2.position.y, box_size),
//            integrator.apply_pbc(p1.position.z - p2.position.z, box_size)
//        );
//        Vec3 b(
//            integrator.apply_pbc(p3.position.x - p2.position.x, box_size),
//            integrator.apply_pbc(p3.position.y - p2.position.y, box_size),
//            integrator.apply_pbc(p3.position.z - p2.position.z, box_size)
//        );
//
//        double ra = a.magnitude();
//        double rb = b.magnitude();
//
//        if (ra < 1e-12 || rb < 1e-12) return 0.0;
//
//        double cos_theta = max(-1.0, min(1.0, a.dot(b) / (ra * rb)));
//        return acos(cos_theta) * 180.0 / Constants::PI;
//    }
//
//public:
//    OptimizedMolecularDynamics(double box)
//        : box_size(box* Constants::ANGSTROM_TO_METERS),
//        dt(0.02e-15),
//        target_temperature(300.0),
//        energy_file_initialized(false),
//        trajectory_file_initialized(false),
//        bonds_angles_file_initialized(false) {
//
//        if (box <= 0) {
//            throw runtime_error("Dimensiune cutie invalida");
//        }
//
//        cout << fixed << setprecision(6);
//        cout << "Dimensiune cutie: " << box << " A\n";
//        cout << "Pas de timp: " << dt * 1e15 << " fs (redus pentru stabilitate)\n";
//    }
//
//    void add_molecule(const string& type, double center_x, double center_y, double center_z) {
//        const MoleculeDefinition& def = registry.get_definition(type);
//        int base_index = (int)particles.size();
//
//        Molecule mol;
//        mol.type = type;
//
//        // Adaugare atomi
//        for (size_t i = 0; i < def.atoms.size(); ++i) {
//            double x = center_x + def.positions[i].x;
//            double y = center_y + def.positions[i].y;
//            double z = center_z + def.positions[i].z;
//
//            particles.push_back(Particle(x, y, z, def.atoms[i].second, def.atoms[i].first));
//            mol.atom_indices.push_back(base_index + (int)i);
//        }
//
//        // Adaugare legaturi
//        for (const auto& bond : def.bonds) {
//            mol.bonds.push_back(make_pair(base_index + bond.first, base_index + bond.second));
//        }
//
//        // Adaugare unghiuri
//        for (const auto& angle : def.angles) {
//            array<int, 3> mol_angle = { base_index + angle[0], base_index + angle[1], base_index + angle[2] };
//            mol.angles.push_back(mol_angle);
//        }
//
//        mol.bond_lengths = def.bond_lengths;
//        mol.angle_degrees = def.angle_degrees;
//        molecules.push_back(mol);
//
//        cout << "Molecula " << type << " adaugata cu " << def.bonds.size()
//            << " legaturi si " << def.angles.size() << " unghiuri.\n";
//    }
//
//    double calculate_kinetic_energy() const {
//        double kinetic_energy = 0.0;
//
//        for (const auto& p : particles) {
//            double m = p.mass * Constants::ATOMIC_MASS_UNIT;
//            kinetic_energy += 0.5 * m * p.velocity.magnitude_squared();
//        }
//
//        return kinetic_energy / Constants::KCAL_PER_MOL_TO_JOULES;
//    }
//
//    void add_thermal_velocities(double temperature) {
//        target_temperature = temperature;
//        random_device rd;
//        mt19937 gen(rd());
//
//        for (auto& p : particles) {
//            double m = p.mass * Constants::ATOMIC_MASS_UNIT;
//            double sigma = sqrt(Constants::BOLTZMANN * temperature / m);
//            normal_distribution<double> dist(0.0, sigma);
//
//            p.velocity = Vec3(dist(gen), dist(gen), dist(gen));
//        }
//
//        // Corectare centrul de masa
//        Vec3 cm_velocity(0.0, 0.0, 0.0);
//        double total_mass = 0.0;
//        for (const auto& p : particles) {
//            cm_velocity += p.velocity * p.mass;
//            total_mass += p.mass;
//        }
//        cm_velocity *= (1.0 / total_mass);
//        for (auto& p : particles) {
//            p.velocity -= cm_velocity;
//        }
//
//        cout << "Viteze termice adaugate pentru T = " << temperature << " K\n";
//        cout << "Energia cinetica initiala: " << fixed << setprecision(4)
//            << calculate_kinetic_energy() << " kcal/mol\n";
//    }
//
//    void update() {
//        integrator.update_positions(particles, box_size, dt, molecules);
//        integrator.apply_berendsen_thermostat(particles, target_temperature, 1.0e-12, dt);
//    }
//
//    void print_status(int step) const {
//        double energy = integrator.compute_total_energy(particles, box_size, molecules);
//        double temp = calculate_temperature();
//
//        cout << "\n=== Pas " << step << " ===\n";
//        cout << "Energia totala: " << fixed << setprecision(4) << energy << " kcal/mol\n";
//        cout << "Temperatura: " << fixed << setprecision(2) << temp << " K\n";
//
//        for (size_t i = 0; i < molecules.size(); ++i) {
//            const Molecule& mol = molecules[i];
//            cout << "\nMolecula " << mol.type << " " << i << ":\n";
//
//            // Lungimi legaturi
//            cout << "  Legaturi:\n";
//            for (size_t j = 0; j < mol.bonds.size(); ++j) {
//                int idx1 = mol.bonds[j].first;
//                int idx2 = mol.bonds[j].second;
//                double length = calculate_bond_length(particles[idx1], particles[idx2]);
//                double deviation = abs(length - mol.bond_lengths[j]);
//
//                cout << "    Legatura " << j + 1 << " (" << idx1 << "-" << idx2
//                    << "): " << fixed << setprecision(4) << length
//                    << " A (echilibru: " << mol.bond_lengths[j]
//                    << " A, deviatie: " << deviation << " A)\n";
//            }
//
//            // Unghiuri
//            cout << "  Unghiuri:\n";
//            for (size_t j = 0; j < mol.angles.size(); ++j) {
//                int idx1 = mol.angles[j][0];
//                int idx2 = mol.angles[j][1];
//                int idx3 = mol.angles[j][2];
//                double angle = calculate_angle(particles[idx1], particles[idx2], particles[idx3]);
//                double deviation = abs(angle - mol.angle_degrees[j]);
//
//                cout << "    Unghi " << j + 1 << " (" << idx1 << "-" << idx2 << "-" << idx3
//                    << "): " << fixed << setprecision(2) << angle
//                    << " grade (echilibru: " << mol.angle_degrees[j]
//                    << " grade, deviatie: " << deviation << " grade)\n";
//            }
//        }
//    }
//
//    double calculate_temperature() const {
//        return integrator.calculate_temperature(particles);
//    }
//
//    void analyze_stability() const {
//        cout << "\n=== ANALIZA STABILITATE ===\n";
//
//        double max_bond_deviation = 0.0;
//        double max_angle_deviation = 0.0;
//        double avg_bond_deviation = 0.0;
//        double avg_angle_deviation = 0.0;
//        size_t bond_count = 0;
//        size_t angle_count = 0;
//
//        for (const auto& mol : molecules) {
//            for (size_t j = 0; j < mol.bonds.size(); ++j) {
//                int idx1 = mol.bonds[j].first;
//                int idx2 = mol.bonds[j].second;
//                double length = calculate_bond_length(particles[idx1], particles[idx2]);
//                double deviation = abs(length - mol.bond_lengths[j]);
//
//                max_bond_deviation = max(max_bond_deviation, deviation);
//                avg_bond_deviation += deviation;
//                ++bond_count;
//            }
//
//            for (size_t j = 0; j < mol.angles.size(); ++j) {
//                int idx1 = mol.angles[j][0];
//                int idx2 = mol.angles[j][1];
//                int idx3 = mol.angles[j][2];
//                double angle = calculate_angle(particles[idx1], particles[idx2], particles[idx3]);
//                double deviation = abs(angle - mol.angle_degrees[j]);
//
//                max_angle_deviation = max(max_angle_deviation, deviation);
//                avg_angle_deviation += deviation;
//                ++angle_count;
//            }
//        }
//
//        if (bond_count > 0) avg_bond_deviation /= bond_count;
//        if (angle_count > 0) avg_angle_deviation /= angle_count;
//
//        cout << "Deviatia maxima pentru legaturi: " << fixed << setprecision(4)
//            << max_bond_deviation << " A\n";
//        cout << "Deviatia medie pentru legaturi: " << fixed << setprecision(4)
//            << avg_bond_deviation << " A\n";
//        cout << "Deviatia maxima pentru unghiuri: " << fixed << setprecision(2)
//            << max_angle_deviation << " grade\n";
//        cout << "Deviatia medie pentru unghiuri: " << fixed << setprecision(2)
//            << avg_angle_deviation << " grade\n";
//
//        if (max_bond_deviation < 0.05 && max_angle_deviation < 10.0 &&
//            avg_bond_deviation < 0.02 && avg_angle_deviation < 5.0) {
//            cout << "Simularea este STABILA\n";
//        }
//        else if (max_bond_deviation < 0.10 && max_angle_deviation < 15.0) {
//            cout << "Simularea este MODERATA\n";
//        }
//        else {
//            cout << "Simularea este INSTABILA\n";
//        }
//    }
//
//    void verify_initial_geometry() const {
//        cout << "\n=== VERIFICARE GEOMETRIE INITIALA ===\n";
//
//        for (size_t i = 0; i < molecules.size(); ++i) {
//            const Molecule& mol = molecules[i];
//            cout << "\nMolecula " << mol.type << " " << i << ":\n";
//
//            if (mol.type == "ammonia") {
//                cout << "  Verificare NH3 cu geometrie CORECTATA:\n";
//                for (size_t j = 0; j < mol.angles.size(); ++j) {
//                    int idx1 = mol.angles[j][0];
//                    int idx2 = mol.angles[j][1];
//                    int idx3 = mol.angles[j][2];
//                    double angle = calculate_angle(particles[idx1], particles[idx2], particles[idx3]);
//                    cout << "    Unghi H-N-H calculat: " << fixed << setprecision(2)
//                        << angle << " grade (tinta: 107.00 grade)\n";
//                }
//            }
//        }
//    }
//
//    // === FUNCTII DE EXPORT PENTRU PYTHON - CORECTATE SA SUPRASCRIE ===
//
//    void export_xyz(const string& filename, int step = -1) const {
//        ofstream file(filename);
//        if (!file.is_open()) {
//            throw runtime_error("Nu pot deschide fisierul pentru scriere: " + filename);
//        }
//
//        file << particles.size() << "\n";
//        if (step >= 0) {
//            file << "Step " << step << " - MD Simulation\n";
//        }
//        else {
//            file << "MD Simulation Export\n";
//        }
//
//        const vector<string> atom_symbols = { "C", "H", "O", "N", "Be", "Cl" };
//
//        for (const auto& p : particles) {
//            double x = p.position.x / Constants::ANGSTROM_TO_METERS;
//            double y = p.position.y / Constants::ANGSTROM_TO_METERS;
//            double z = p.position.z / Constants::ANGSTROM_TO_METERS;
//
//            file << atom_symbols[p.type] << " "
//                << fixed << setprecision(6) << x << " " << y << " " << z << "\n";
//        }
//        file.close();
//        cout << "Export XYZ salvat in: " << filename << "\n";
//    }
//
//    // CORECTAT: CSV trajectory - initializeaza fisierul o singura data apoi adauga
//    void export_csv_trajectory(const string& filename, int step) {
//        ios_base::openmode mode = ios::out;
//
//        if (!trajectory_file_initialized) {
//            // Prima data - suprascrie fisierul si adauga header
//            ofstream file(filename, mode);
//            if (!file.is_open()) {
//                throw runtime_error("Nu pot deschide fisierul pentru scriere: " + filename);
//            }
//            file << "step,atom_id,type,x,y,z,vx,vy,vz,fx,fy,fz,mass\n";
//            file.close();
//            trajectory_file_initialized = true;
//        }
//
//        // Adauga datele la fisier
//        ofstream file(filename, ios::app);
//        if (!file.is_open()) {
//            throw runtime_error("Nu pot deschide fisierul pentru scriere: " + filename);
//        }
//
//        const vector<string> atom_symbols = { "C", "H", "O", "N", "Be", "Cl" };
//
//        for (size_t i = 0; i < particles.size(); ++i) {
//            const auto& p = particles[i];
//            double x = p.position.x / Constants::ANGSTROM_TO_METERS;
//            double y = p.position.y / Constants::ANGSTROM_TO_METERS;
//            double z = p.position.z / Constants::ANGSTROM_TO_METERS;
//
//            file << step << "," << i << "," << atom_symbols[p.type] << ","
//                << fixed << setprecision(6) << x << "," << y << "," << z << ","
//                << p.velocity.x << "," << p.velocity.y << "," << p.velocity.z << ","
//                << p.force.x << "," << p.force.y << "," << p.force.z << ","
//                << p.mass << "\n";
//        }
//        file.close();
//    }
//
//    // CORECTAT: Energy data - initializeaza fisierul o singura data apoi adauga
//    void export_energy_data(const string& filename, int step) {
//        if (!energy_file_initialized) {
//            // Prima data - suprascrie fisierul si adauga header
//            ofstream file(filename, ios::out);
//            if (!file.is_open()) {
//                throw runtime_error("Nu pot deschide fisierul pentru scriere: " + filename);
//            }
//            file << "step,total_energy,temperature,time_fs\n";
//            file.close();
//            energy_file_initialized = true;
//        }
//
//        // Adauga datele la fisier
//        ofstream file(filename, ios::app);
//        if (!file.is_open()) {
//            throw runtime_error("Nu pot deschide fisierul pentru scriere: " + filename);
//        }
//
//        double energy = integrator.compute_total_energy(particles, box_size, molecules);
//        double temp = calculate_temperature();
//        double time_fs = step * dt * 1e15;
//
//        file << step << "," << fixed << setprecision(6) << energy << ","
//            << temp << "," << time_fs << "\n";
//        file.close();
//    }
//
//    // CORECTAT: Bonds and angles - initializeaza fisierul o singura data apoi adauga
//    void export_bonds_angles(const string& filename, int step) {
//        if (!bonds_angles_file_initialized) {
//            // Prima data - suprascrie fisierul si adauga header
//            ofstream file(filename, ios::out);
//            if (!file.is_open()) {
//                throw runtime_error("Nu pot deschide fisierul pentru scriere: " + filename);
//            }
//            file << "step,molecule_id,molecule_type,measurement_type,measurement_id,value,equilibrium,deviation\n";
//            file.close();
//            bonds_angles_file_initialized = true;
//        }
//
//        // Adauga datele la fisier
//        ofstream file(filename, ios::app);
//        if (!file.is_open()) {
//            throw runtime_error("Nu pot deschide fisierul pentru scriere: " + filename);
//        }
//
//        for (size_t i = 0; i < molecules.size(); ++i) {
//            const Molecule& mol = molecules[i];
//
//            // Export bonds
//            for (size_t j = 0; j < mol.bonds.size(); ++j) {
//                int idx1 = mol.bonds[j].first;
//                int idx2 = mol.bonds[j].second;
//                double length = calculate_bond_length(particles[idx1], particles[idx2]);
//                double deviation = abs(length - mol.bond_lengths[j]);
//
//                file << step << "," << i << "," << mol.type << ",bond," << j << ","
//                    << fixed << setprecision(6) << length << "," << mol.bond_lengths[j]
//                    << "," << deviation << "\n";
//            }
//
//            // Export angles
//            for (size_t j = 0; j < mol.angles.size(); ++j) {
//                int idx1 = mol.angles[j][0];
//                int idx2 = mol.angles[j][1];
//                int idx3 = mol.angles[j][2];
//                double angle = calculate_angle(particles[idx1], particles[idx2], particles[idx3]);
//                double deviation = abs(angle - mol.angle_degrees[j]);
//
//                file << step << "," << i << "," << mol.type << ",angle," << j << ","
//                    << fixed << setprecision(6) << angle << "," << mol.angle_degrees[j]
//                    << "," << deviation << "\n";
//            }
//        }
//        file.close();
//    }
//
//    void export_molecule_topology(const string& filename) const {
//        ofstream file(filename);
//        if (!file.is_open()) {
//            throw runtime_error("Nu pot deschide fisierul pentru scriere: " + filename);
//        }
//
//        file << "# Topologia moleculelor - informatii despre legaturi si unghiuri\n";
//        file << "molecule_id,molecule_type,atom_count,bond_count,angle_count\n";
//
//        for (size_t i = 0; i < molecules.size(); ++i) {
//            const Molecule& mol = molecules[i];
//            file << i << "," << mol.type << "," << mol.atom_indices.size()
//                << "," << mol.bonds.size() << "," << mol.angles.size() << "\n";
//        }
//
//        file << "\n# Detalii legaturi\n";
//        file << "molecule_id,bond_id,atom1,atom2,equilibrium_length\n";
//        for (size_t i = 0; i < molecules.size(); ++i) {
//            const Molecule& mol = molecules[i];
//            for (size_t j = 0; j < mol.bonds.size(); ++j) {
//                file << i << "," << j << "," << mol.bonds[j].first << ","
//                    << mol.bonds[j].second << "," << mol.bond_lengths[j] << "\n";
//            }
//        }
//
//        file << "\n# Detalii unghiuri\n";
//        file << "molecule_id,angle_id,atom1,atom2,atom3,equilibrium_angle\n";
//        for (size_t i = 0; i < molecules.size(); ++i) {
//            const Molecule& mol = molecules[i];
//            for (size_t j = 0; j < mol.angles.size(); ++j) {
//                file << i << "," << j << "," << mol.angles[j][0] << ","
//                    << mol.angles[j][1] << "," << mol.angles[j][2] << ","
//                    << mol.angle_degrees[j] << "\n";
//            }
//        }
//
//        file.close();
//        cout << "Topologia moleculelor exportata in: " << filename << "\n";
//    }
//};
//
//// Functie main optimizata cu exporturi CORECTATE
//int main() {
//    try {
//        cout << "=== SIMULARE DINAMICA MOLECULARA - VERSIUNEA CORECTATA CU SUPRASCRIEREA CSV ===\n\n";
//
//        // Initializare sistem
//        OptimizedMolecularDynamics md(40.0);
//
//        // Adaugare molecule
//        cout << "\nAdaugare molecule:\n";
//        md.add_molecule("methane", 10.0, 10.0, 10.0);
//        md.add_molecule("water", 20.0, 10.0, 10.0);
//        md.add_molecule("co2", 30.0, 10.0, 10.0);
//        md.add_molecule("ammonia", 10.0, 20.0, 10.0);
//        md.add_molecule("acetylene", 20.0, 20.0, 10.0);
//        md.add_molecule("hcl", 30.0, 20.0, 10.0);
//
//        // Export topologie (se suprascrie mereu)
//        md.export_molecule_topology("topology.csv");
//
//        // Verificare geometrie initiala
//        md.verify_initial_geometry();
//
//        // Adaugare viteze termice
//        md.add_thermal_velocities(300.0);
//
//        // Export stare initiala
//        md.export_xyz("initial_structure.xyz", 0);
//
//        // Stare initiala
//        md.print_status(0);
//        cout << "\nTemperatura initiala: " << fixed << setprecision(1)
//            << md.calculate_temperature() << " K\n";
//
//        // Export date initiale - acum se suprascriu fisierele CSV
//        md.export_energy_data("energy.csv", 0);
//        md.export_bonds_angles("bonds_angles.csv", 0);
//
//        // Rulare simulare
//        cout << "\n" << string(70, '=') << "\n";
//        cout << "EVOLUTIA SIMULARII CU EXPORTURI CORECTATE\n";
//        cout << string(70, '=') << "\n";
//
//        int total_steps = 100;
//        int print_interval = 20;
//        int export_interval = 5;
//
//        for (int step = 1; step <= total_steps; ++step) {
//            md.update();
//
//            // Export la fiecare export_interval
//            if (step % export_interval == 0) {
//                md.export_energy_data("energy.csv", step);
//                md.export_bonds_angles("bonds_angles.csv", step);
//                md.export_csv_trajectory("trajectory.csv", step);
//
//                // Export XYZ la anumite pasi
//                if (step % (export_interval * 4) == 0) {
//                    md.export_xyz("structure_step_" + to_string(step) + ".xyz", step);
//                }
//            }
//
//            if (step % print_interval == 0) {
//                md.print_status(step);
//                cout << "Temperatura: " << fixed << setprecision(1)
//                    << md.calculate_temperature() << " K\n";
//
//                if (step == total_steps) {
//                    md.analyze_stability();
//                }
//            }
//        }
//
//        // Export final
//        md.export_xyz("final_structure.xyz", total_steps);
//
//        // Sumar final
//        cout << "\n" << string(70, '=') << "\n";
//        cout << "SIMULAREA S-A COMPLETAT CU SUCCES!\n";
//        cout << "PROBLEMA CORECTATA:\n";
//        cout << "  - CSV-urile se suprascriu acum la fiecare rulare noua\n";
//        cout << "  - La prima rulare se creeaza headerul si se suprascrie continutul\n";
//        cout << "  - La pasii urmatori se adauga doar datele noi\n";
//        cout << "  - Geometrie NH3 PERFECT CORECTA (unghiuri exacte de 107 grade)\n";
//        cout << "  - Exporturi complete pentru Python:\n";
//        cout << "    * energy.csv - evolutia energiei si temperaturii\n";
//        cout << "    * trajectory.csv - pozitii, viteze, forte\n";
//        cout << "    * bonds_angles.csv - geometria moleculelor\n";
//        cout << "    * topology.csv - structura moleculelor\n";
//        cout << "    * *.xyz - structuri pentru vizualizare\n";
//        cout << "  - Toate exporturile sunt compatibile cu pandas/matplotlib\n";
//        cout << "  - Fara diacritice in cod\n";
//        cout << string(70, '=') << "\n";
//
//    }
//    catch (const exception& e) {
//        cerr << "EROARE: " << e.what() << endl;
//        return 1;
//    }
//
//    return 0;
//}