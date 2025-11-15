#include "md_system.h"
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

int main() {
    try {
        cout << "=== MOLECULAR DYNAMICS SIMULATION ===\n\n";

        OptimizedMolecularDynamics md(50.0);

        cout << "\nAdding molecules:\n";
        md.add_molecule("methane", 10.0, 10.0, 10.0);
        md.add_molecule("water", 13.0, 10.0, 10.0);
        md.add_molecule("co2", 16.0, 10.0, 10.0);
        md.add_molecule("ammonia", 10.0, 13.0, 10.0);
        md.add_molecule("acetylene", 13.0, 13.0, 10.0);
        md.add_molecule("hcl", 16.0, 13.0, 10.0);

        md.export_molecule_topology("topology.csv");

        md.add_thermal_velocities(300.0);

        md.export_xyz("initial_structure.xyz", 0);

        md.print_status(0);
        cout << "\nInitial temperature: " << md.calculate_temperature() << " K\n";

        md.export_energy_data("energy.csv", 0);
        md.export_bonds_angles("bonds_angles.csv", 0);
        md.export_csv_trajectory("trajectory.csv", 0);
        md.export_lj_forces("lj_forces.csv", 0);

        cout << "\n" << string(70, '=') << "\n";
        cout << "SIMULATION EVOLUTION\n";
        cout << string(70, '=') << "\n";

        const int total_steps = 250;
        const int print_interval = 5;
        const int export_interval = 5;
        const int xyz_export_interval = export_interval * 2;

        for (int step = 1; step <= total_steps; ++step) {
            md.update();

            if (step % export_interval == 0) {
                md.export_energy_data("energy.csv", step);
                md.export_bonds_angles("bonds_angles.csv", step);
                md.export_csv_trajectory("trajectory.csv", step);
                md.export_lj_forces("lj_forces.csv", step);

                if (step % xyz_export_interval == 0) {
                    const string xyz_filename = "structure_step_" + to_string(step) + ".xyz";
                    md.export_xyz(xyz_filename, step);
                }
            }

            if (step % print_interval == 0) {
                const double temp = md.calculate_temperature();
                md.print_status(step);
                cout << "Temperature: " << temp << " K\n";

                if (step == total_steps) {
                    md.analyze_stability();
                }
            }
        }

        md.export_xyz("final_structure.xyz", total_steps);

    }
    catch (const exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }

    return 0;
}