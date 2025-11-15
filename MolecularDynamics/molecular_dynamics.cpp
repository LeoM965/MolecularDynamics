#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;

// Constanta lui Boltzmann in J/K
const double k_B = 1.38064852e-23; // Constanta lui Boltzmann in jouli pe kelvin
const double KCAL_PER_MOL_TO_JOULES = 4.184e-21; // Conversie de la kcal/mol la jouli
const double ANGSTROM_TO_METERS = 1e-10; // Conversie de la angstromi la metri
const double M_PI = 3.141592653589793; // Valoare a lui pi

// Structura pentru o particula cu pozitie, viteza, forta, masa si tip atom
struct Particle {
    double x, y, z; // Pozitie in metri
    double vx, vy, vz; // Viteza in m/s
    double fx, fy, fz; // Forta in newtoni
    double mass; // Masa in unitati atomice (u)
    int type; // Tip atom (0 pentru carbon, 1 pentru hidrogen)
};

// Converteste masa din unitati atomice in kilograme
double atomic_mass_unit_to_kg(double mass_u) {
    return mass_u * 1.66053906660e-27; // 1 u = 1.66053906660e-27 kg
}

// Structura pentru o moleculă, care conține atomi legați
struct Molecule {
    vector<int> atom_indices; // Indicii atomilor din vectorul global de particule
    vector<pair<int, int>> bonds; // Perechi de atomi legați (ex. C-H în metan)
};

// Clasa abstractă pentru integratori
class Integrator {
public:
    virtual void update_positions(vector<Particle>& particles, double box_size, double dt, vector<double>& sigma, vector<double>& epsilon, vector<Molecule>& molecules) = 0;
    virtual ~Integrator() {}
};

// Implementarea integratorului Verlet (existent)
class VerletIntegrator : public Integrator {
public:
    void update_positions(vector<Particle>& particles, double box_size, double dt, vector<double>& sigma, vector<double>& epsilon, vector<Molecule>& molecules) override {
        // Calculează forțele
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            it->fx = it->fy = it->fz = 0.0;
        }
        // Forțe intermoleculare (Lennard-Jones)
        for (vector<Particle>::iterator it_i = particles.begin(); it_i != particles.end(); ++it_i) {
            vector<Particle>::iterator it_j = it_i;
            ++it_j;
            for (; it_j != particles.end(); ++it_j) {
                Particle& pi = *it_i;
                Particle& pj = *it_j;
                double dx = pi.x - pj.x;
                double dy = pi.y - pj.y;
                double dz = pi.z - pj.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > 0) {
                    double r = sqrt(r2);
                    double sig = (sigma[pi.type] + sigma[pj.type]) / 2.0;
                    double eps = sqrt(epsilon[pi.type] * epsilon[pj.type]);
                    if (r < 3.0 * sig) {
                        double r6 = pow(sig / r, 6);
                        double r12 = r6 * r6;
                        double force = 24.0 * eps * (2.0 * r12 - r6) / r;
                        double fx = force * dx / r;
                        double fy = force * dy / r;
                        double fz = force * dz / r;
                        pi.fx += fx;
                        pi.fy += fy;
                        pi.fz += fz;
                        pj.fx -= fx;
                        pj.fy -= fy;
                        pj.fz -= fz;
                    }
                }
            }
        }
        // Forțe intramoleculare (armonice)
        const double k = 300.0 * KCAL_PER_MOL_TO_JOULES / (ANGSTROM_TO_METERS * ANGSTROM_TO_METERS);
        const double r0 = 1.09 * ANGSTROM_TO_METERS;
        for (int i = 0; i < molecules.size(); ++i) {
            Molecule& mol = molecules[i];
            for (int j = 0; j < mol.bonds.size(); ++j) {
                int idx1 = mol.bonds[j].first;
                int idx2 = mol.bonds[j].second;
                Particle& p1 = particles[idx1];
                Particle& p2 = particles[idx2];
                double dx = p1.x - p2.x;
                double dy = p1.y - p2.y;
                double dz = p1.z - p2.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                if (r > 0) {
                    double force = -k * (r - r0);
                    double fx = force * dx / r;
                    double fy = force * dy / r;
                    double fz = force * dz / r;
                    p1.fx += fx;
                    p1.fy += fy;
                    p1.fz += fz;
                    p2.fx -= fx;
                    p2.fy -= fy;
                    p2.fz -= fz;
                }
            }
        }
        // Actualizează pozițiile și vitezele (Verlet)
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            double mass_kg = atomic_mass_unit_to_kg(p.mass);
            p.vx += 0.5 * p.fx / mass_kg * dt;
            p.vy += 0.5 * p.fy / mass_kg * dt;
            p.vz += 0.5 * p.fz / mass_kg * dt;
            p.x += p.vx * dt;
            p.y += p.vy * dt;
            p.z += p.vz * dt;
            p.x -= box_size * floor(p.x / box_size);
            p.y -= box_size * floor(p.y / box_size);
            p.z -= box_size * floor(p.z / box_size);
        }
        // Recalculează forțele după actualizarea pozițiilor
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            it->fx = it->fy = it->fz = 0.0;
        }
        // Forțe intermoleculare (Lennard-Jones)
        for (vector<Particle>::iterator it_i = particles.begin(); it_i != particles.end(); ++it_i) {
            vector<Particle>::iterator it_j = it_i;
            ++it_j;
            for (; it_j != particles.end(); ++it_j) {
                Particle& pi = *it_i;
                Particle& pj = *it_j;
                double dx = pi.x - pj.x;
                double dy = pi.y - pj.y;
                double dz = pi.z - pj.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > 0) {
                    double r = sqrt(r2);
                    double sig = (sigma[pi.type] + sigma[pj.type]) / 2.0;
                    double eps = sqrt(epsilon[pi.type] * epsilon[pj.type]);
                    if (r < 3.0 * sig) {
                        double r6 = pow(sig / r, 6);
                        double r12 = r6 * r6;
                        double force = 24.0 * eps * (2.0 * r12 - r6) / r;
                        double fx = force * dx / r;
                        double fy = force * dy / r;
                        double fz = force * dz / r;
                        pi.fx += fx;
                        pi.fy += fy;
                        pi.fz += fz;
                        pj.fx -= fx;
                        pj.fy -= fy;
                        pj.fz -= fz;
                    }
                }
            }
        }
        // Forțe intramoleculare (armonice)
        for (int i = 0; i < molecules.size(); ++i) {
            Molecule& mol = molecules[i];
            for (int j = 0; j < mol.bonds.size(); ++j) {
                int idx1 = mol.bonds[j].first;
                int idx2 = mol.bonds[j].second;
                Particle& p1 = particles[idx1];
                Particle& p2 = particles[idx2];
                double dx = p1.x - p2.x;
                double dy = p1.y - p2.y;
                double dz = p1.z - p2.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                if (r > 0) {
                    double force = -k * (r - r0);
                    double fx = force * dx / r;
                    double fy = force * dy / r;
                    double fz = force * dz / r;
                    p1.fx += fx;
                    p1.fy += fy;
                    p1.fz += fz;
                    p2.fx -= fx;
                    p2.fy -= fy;
                    p2.fz -= fz;
                }
            }
        }
        // Finalizează actualizarea vitezelor
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            double mass_kg = atomic_mass_unit_to_kg(p.mass);
            p.vx += 0.5 * p.fx / mass_kg * dt;
            p.vy += 0.5 * p.fy / mass_kg * dt;
            p.vz += 0.5 * p.fz / mass_kg * dt;
        }
    }
};

// Implementarea integratorului Velocity Verlet
class VelocityVerletIntegrator : public Integrator {
public:
    void update_positions(vector<Particle>& particles, double box_size, double dt, vector<double>& sigma, vector<double>& epsilon, vector<Molecule>& molecules) override {
        // Calculează forțele
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            it->fx = it->fy = it->fz = 0.0;
        }
        // Forțe intermoleculare (Lennard-Jones)
        for (vector<Particle>::iterator it_i = particles.begin(); it_i != particles.end(); ++it_i) {
            vector<Particle>::iterator it_j = it_i;
            ++it_j;
            for (; it_j != particles.end(); ++it_j) {
                Particle& pi = *it_i;
                Particle& pj = *it_j;
                double dx = pi.x - pj.x;
                double dy = pi.y - pj.y;
                double dz = pi.z - pj.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > 0) {
                    double r = sqrt(r2);
                    double sig = (sigma[pi.type] + sigma[pj.type]) / 2.0;
                    double eps = sqrt(epsilon[pi.type] * epsilon[pj.type]);
                    if (r < 3.0 * sig) {
                        double r6 = pow(sig / r, 6);
                        double r12 = r6 * r6;
                        double force = 24.0 * eps * (2.0 * r12 - r6) / r;
                        double fx = force * dx / r;
                        double fy = force * dy / r;
                        double fz = force * dz / r;
                        pi.fx += fx;
                        pi.fy += fy;
                        pi.fz += fz;
                        pj.fx -= fx;
                        pj.fy -= fy;
                        pj.fz -= fz;
                    }
                }
            }
        }
        // Forțe intramoleculare (armonice)
        const double k = 300.0 * KCAL_PER_MOL_TO_JOULES / (ANGSTROM_TO_METERS * ANGSTROM_TO_METERS);
        const double r0 = 1.09 * ANGSTROM_TO_METERS;
        for (int i = 0; i < molecules.size(); ++i) {
            Molecule& mol = molecules[i];
            for (int j = 0; j < mol.bonds.size(); ++j) {
                int idx1 = mol.bonds[j].first;
                int idx2 = mol.bonds[j].second;
                Particle& p1 = particles[idx1];
                Particle& p2 = particles[idx2];
                double dx = p1.x - p2.x;
                double dy = p1.y - p2.y;
                double dz = p1.z - p2.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                if (r > 0) {
                    double force = -k * (r - r0);
                    double fx = force * dx / r;
                    double fy = force * dy / r;
                    double fz = force * dz / r;
                    p1.fx += fx;
                    p1.fy += fy;
                    p1.fz += fz;
                    p2.fx -= fx;
                    p2.fy -= fy;
                    p2.fz -= fz;
                }
            }
        }
        // Velocity Verlet: actualizează pozițiile și vitezele
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            double mass_kg = atomic_mass_unit_to_kg(p.mass);
            double ax = p.fx / mass_kg;
            double ay = p.fy / mass_kg;
            double az = p.fz / mass_kg;
            p.x += p.vx * dt + 0.5 * ax * dt * dt;
            p.y += p.vy * dt + 0.5 * ay * dt * dt;
            p.z += p.vz * dt + 0.5 * az * dt * dt;
            p.x -= box_size * floor(p.x / box_size);
            p.y -= box_size * floor(p.y / box_size);
            p.z -= box_size * floor(p.z / box_size);
            p.vx += 0.5 * ax * dt;
            p.vy += 0.5 * ay * dt;
            p.vz += 0.5 * az * dt;
        }
        // Recalculează forțele
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            it->fx = it->fy = it->fz = 0.0;
        }
        // Forțe intermoleculare (Lennard-Jones)
        for (vector<Particle>::iterator it_i = particles.begin(); it_i != particles.end(); ++it_i) {
            vector<Particle>::iterator it_j = it_i;
            ++it_j;
            for (; it_j != particles.end(); ++it_j) {
                Particle& pi = *it_i;
                Particle& pj = *it_j;
                double dx = pi.x - pj.x;
                double dy = pi.y - pj.y;
                double dz = pi.z - pj.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > 0) {
                    double r = sqrt(r2);
                    double sig = (sigma[pi.type] + sigma[pj.type]) / 2.0;
                    double eps = sqrt(epsilon[pi.type] * epsilon[pj.type]);
                    if (r < 3.0 * sig) {
                        double r6 = pow(sig / r, 6);
                        double r12 = r6 * r6;
                        double force = 24.0 * eps * (2.0 * r12 - r6) / r;
                        double fx = force * dx / r;
                        double fy = force * dy / r;
                        double fz = force * dz / r;
                        pi.fx += fx;
                        pi.fy += fy;
                        pi.fz += fz;
                        pj.fx -= fx;
                        pj.fy -= fy;
                        pj.fz -= fz;
                    }
                }
            }
        }
        // Forțe intramoleculare (armonice Giuliano)
        for (int i = 0; i < molecules.size(); ++i) {
            Molecule& mol = molecules[i];
            for (int j = 0; j < mol.bonds.size(); ++j) {
                int idx1 = mol.bonds[j].first;
                int idx2 = mol.bonds[j].second;
                Particle& p1 = particles[idx1];
                Particle& p2 = particles[idx2];
                double dx = p1.x - p2.x;
                double dy = p1.y - p2.y;
                double dz = p1.z - p2.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                if (r > 0) {
                    double force = -k * (r - r0);
                    double fx = force * dx / r;
                    double fy = force * dy / r;
                    double fz = force * dz / r;
                    p1.fx += fx;
                    p1.fy += fy;
                    p1.fz += fz;
                    p2.fx -= fx;
                    p2.fy -= fy;
                    p2.fz -= fz;
                }
            }
        }
        // Finalizează actualizarea vitezelor
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            double mass_kg = atomic_mass_unit_to_kg(p.mass);
            p.vx += 0.5 * p.fx / mass_kg * dt;
            p.vy += 0.5 * p.fy / mass_kg * dt;
            p.vz += 0.5 * p.fz / mass_kg * dt;
        }
    }
};

// Clasa pentru dinamica moleculară
class MolecularDynamics {
private:
    vector<Particle> particles;
    vector<Molecule> molecules; // Vector de molecule
    vector<Particle> initial_positions; // Pozițiile inițiale pentru calculul MSD
    double box_size; // Dimensiunea cutiei de simulare in metri
    double dt; // Pasul de timp in secunde
    vector<double> sigma; // Parametru Lennard-Jones sigma in metri
    vector<double> epsilon; // Parametru Lennard-Jones epsilon in jouli
    ofstream out_file; // Fisier pentru traiectorii
    ofstream energy_file; // Fisier pentru energie
    double temperature; // Temperatura in kelvini
    double initial_total_energy; // Energie totala initiala pentru verificare
    Integrator* integrator; // Pointer către integrator

public:
    MolecularDynamics(double box, double time_step, double temp, Integrator* integ)
        : box_size(box* ANGSTROM_TO_METERS), dt(time_step * 1e-15), temperature(temp), integrator(integ) {
        out_file.open("trajectory.xyz");
        energy_file.open("energy.csv");
        // Parametrii Lennard-Jones pentru carbon (0) si hidrogen (1)
        sigma = { 3.4 * ANGSTROM_TO_METERS, 2.59 * ANGSTROM_TO_METERS }; // Å la m
        epsilon = { 0.0103 * KCAL_PER_MOL_TO_JOULES, 0.0152 * KCAL_PER_MOL_TO_JOULES }; // kcal/mol la J
        energy_file << "Step,Kinetic_Energy,Potential_Energy,Total_Energy,Energy_Deviation\n";
        srand(time(0));
    }

    ~MolecularDynamics() {
        out_file.close();
        energy_file.close();
        delete integrator;
    }

    // Adauga o particula cu pozitie, masa si tip atom
    void add_particle(double x, double y, double z, double mass, int type) {
        Particle p = { x * ANGSTROM_TO_METERS, y * ANGSTROM_TO_METERS, z * ANGSTROM_TO_METERS,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, mass, type };
        particles.push_back(p);
    }

    // Adaugă o moleculă de metan (CH₄) cu poziția centrului de masă specificată
    void add_methane_molecule(double center_x, double center_y, double center_z) {
        Molecule mol;
        // Adaugă atomul de carbon (centru)
        int carbon_index = particles.size();
        add_particle(center_x, center_y, center_z, 12.0, 0); // Carbon
        mol.atom_indices.push_back(carbon_index);

        // Distanța C-H în metan este ~1.09 Å
        double d = 1.09; // în Å
        // Poziții tetraedrice pentru hidrogeni
        add_particle(center_x + d, center_y + d, center_z + d, 1.0, 1); // H1
        mol.atom_indices.push_back(particles.size() - 1);
        add_particle(center_x - d, center_y - d, center_z + d, 1.0, 1); // H2
        mol.atom_indices.push_back(particles.size() - 1);
        add_particle(center_x + d, center_y - d, center_z - d, 1.0, 1); // H3
        mol.atom_indices.push_back(particles.size() - 1);
        add_particle(center_x - d, center_y + d, center_z - d, 1.0, 1); // H4
        mol.atom_indices.push_back(particles.size() - 1);

        // Definește legăturile C-H
        mol.bonds.push_back(pair<int, int>(carbon_index, carbon_index + 1)); // C-H1
        mol.bonds.push_back(pair<int, int>(carbon_index, carbon_index + 2)); // C-H2
        mol.bonds.push_back(pair<int, int>(carbon_index, carbon_index + 3)); // C-H3
        mol.bonds.push_back(pair<int, int>(carbon_index, carbon_index + 4)); // C-H4

        molecules.push_back(mol);
    }

    // Genereaza o valoare din distributia normala folosind metoda Box-Muller
    double box_muller() {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    }

    // Initializeaza vitezele conform distributiei Maxwell-Boltzmann
    void initialize_velocities() {
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            double mass_kg = atomic_mass_unit_to_kg(p.mass);
            double sigma_v = sqrt(k_B * temperature / mass_kg); // Deviatia standard a vitezei
            p.vx = box_muller() * sigma_v;
            p.vy = box_muller() * sigma_v;
            p.vz = box_muller() * sigma_v;
        }
        // Elimina impulsul net
        double vx_sum = 0.0, vy_sum = 0.0, vz_sum = 0.0;
        double total_mass = 0.0;
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            vx_sum += p.vx * p.mass;
            vy_sum += p.vy * p.mass;
            vz_sum += p.vz * p.mass;
            total_mass += p.mass;
        }
        double inv_mass_total = 1.0 / total_mass;
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            p.vx -= vx_sum * inv_mass_total;
            p.vy -= vy_sum * inv_mass_total;
            p.vz -= vz_sum * inv_mass_total;
        }
        // Afiseaza vitezele initiale pentru depanare
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            cout << "Viteza initiala: (" << p.vx << ", " << p.vy << ", " << p.vz << ") m/s\n";
        }
    }

    // Salvează pozițiile inițiale pentru calculul MSD
    void save_initial_positions() {
        initial_positions = particles;
    }

    // Calculează deplasarea pătratică medie pentru atomii de carbon
    void compute_msd(int step) {
        double msd = 0.0;
        int carbon_count = 0;
        for (int i = 0; i < particles.size(); ++i) {
            if (particles[i].type == 0) { // Doar atomi de carbon
                Particle& p = particles[i];
                Particle& p0 = initial_positions[i];
                double dx = p.x - p0.x;
                double dy = p.y - p0.y;
                double dz = p.z - p0.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                msd += dx * dx + dy * dy + dz * dz;
                carbon_count++;
            }
        }
        msd /= carbon_count; // Media pe atomi de carbon
        ofstream msd_file("msd.csv", ios::app);
        if (step == 0) {
            msd_file << "Step,MSD\n";
        }
        msd_file << step << "," << msd / (ANGSTROM_TO_METERS * ANGSTROM_TO_METERS) << "\n";
        msd_file.close();
    }

    // Calculeaza energia cinetica si potentiala
    void compute_energy(double& kinetic, double& potential) {
        kinetic = 0.0;
        potential = 0.0;
        // Energie cinetică
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            double mass_kg = atomic_mass_unit_to_kg(p.mass);
            kinetic += 0.5 * mass_kg * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
        }
        // Energie potențială intermoleculară (Lennard-Jones)
        for (vector<Particle>::iterator it_i = particles.begin(); it_i != particles.end(); ++it_i) {
            vector<Particle>::iterator it_j = it_i;
            ++it_j;
            for (; it_j != particles.end(); ++it_j) {
                Particle& pi = *it_i;
                Particle& pj = *it_j;
                double dx = pi.x - pj.x;
                double dy = pi.y - pj.y;
                double dz = pi.z - pj.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > 0) {
                    double r = sqrt(r2);
                    double sig = (sigma[pi.type] + sigma[pj.type]) / 2.0;
                    double eps = sqrt(epsilon[pi.type] * epsilon[pj.type]);
                    if (r < 3.0 * sig) {
                        double r6 = pow(sig / r, 6);
                        double r12 = r6 * r6;
                        potential += 4.0 * eps * (r12 - r6);
                    }
                }
            }
        }
        // Energie potențială intramoleculară (armonică)
        const double k = 300.0 * KCAL_PER_MOL_TO_JOULES / (ANGSTROM_TO_METERS * ANGSTROM_TO_METERS); // N/m
        const double r0 = 1.09 * ANGSTROM_TO_METERS; // m
        for (int i = 0; i < molecules.size(); ++i) {
            Molecule& mol = molecules[i];
            for (int j = 0; j < mol.bonds.size(); ++j) {
                int idx1 = mol.bonds[j].first;
                int idx2 = mol.bonds[j].second;
                Particle& p1 = particles[idx1];
                Particle& p2 = particles[idx2];
                double dx = p1.x - p2.x;
                double dy = p1.y - p2.y;
                double dz = p1.z - p2.z;
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                potential += 0.5 * k * (r - r0) * (r - r0);
            }
        }
        // Afiseaza energiile pentru depanare
        cout << "Energie cinetica: " << kinetic << " J, Energie potentiala: " << potential << " J\n";
    }

    // Actualizeaza pozitiile folosind integratorul
    void update_positions() {
        integrator->update_positions(particles, box_size, dt, sigma, epsilon, molecules);
        // Afiseaza fortele pentru depanare
        for (vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
            Particle& p = *it;
            cout << "Forta: (" << p.fx << ", " << p.fy << ", " << p.fz << ") N\n";
        }
    }

    // Salveaza traiectoriile particulelor in fisier XYZ
    void save_trajectory(int step) {
        out_file << particles.size() << "\n";
        out_file << "Pas " << step << "\n";
        for (vector<Particle>::const_iterator it = particles.begin(); it != particles.end(); ++it) {
            const Particle& p = *it;
            if (p.type == 0) out_file << "C ";
            else out_file << "H ";
            out_file << p.x / ANGSTROM_TO_METERS << " " << p.y / ANGSTROM_TO_METERS << " " << p.z / ANGSTROM_TO_METERS << "\n";
        }
    }

    // Salveaza energia in fisier CSV cu deviatia
    void save_energy(int step) {
        double kinetic, potential;
        compute_energy(kinetic, potential);
        double total_energy = kinetic + potential;
        if (step == 0) {
            initial_total_energy = total_energy; // Salveaza energia initiala
        }
        double energy_deviation = total_energy - initial_total_energy;
        energy_file << step << "," << kinetic << "," << potential << "," << total_energy << "," << energy_deviation << "\n";
        // Afiseaza deviatia pentru depanare
        cout << "Deviatia energiei: " << energy_deviation << " J\n";
    }

    // Ruleaza simularea pentru un numar dat de pasi
    void run_simulation(int steps) {
        initialize_velocities();
        save_initial_positions(); // Salvează pozițiile inițiale pentru MSD
        update_positions();
        for (int i = 0; i < steps; ++i) {
            update_positions();
            save_trajectory(i);
            save_energy(i);
            compute_msd(i); // Calculează MSD la fiecare pas
        }
    }
};

int main() {
    // Initializeaza simularea cu o cutie de 20 Å, pas de timp 1 fs, temperatura 300 K
    Integrator* integrator = new VerletIntegrator(); // Poți schimba cu VelocityVerletIntegrator
    MolecularDynamics md(20.0, 1.0, 300.0, integrator);

    // Adaugă trei molecule de metan la poziții diferite
    md.add_methane_molecule(5.0, 5.0, 5.0);   // Moleculă 1
    md.add_methane_molecule(10.0, 10.0, 10.0); // Moleculă 2
    md.add_methane_molecule(15.0, 5.0, 5.0);   // Moleculă 3

    md.run_simulation(100); // Rulează 100 de pași
    return 0;
}