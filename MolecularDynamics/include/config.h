#ifndef CONFIG_H
#define CONFIG_H

#include <string>

struct SimulationConfig {
    double box_size = 50.0;
    double dt = 0.02e-15;
    double target_temperature = 300.0;
    double langevin_friction = 1.0e10;
    
    int total_steps = 250;
    int export_interval = 5;
    
    std::string output_dir = "output";
    std::string trajectory_file = "output/trajectory.csv";
    std::string energy_file = "output/energy.csv";
    std::string topology_file = "output/topology.csv";
    std::string bonds_angles_file = "output/bonds_angles.csv";
};

#endif
