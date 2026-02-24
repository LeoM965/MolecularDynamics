#include "../include/md_system.h"
#include "../include/config.h"
#include "../include/md_reporter.h"
#include "../include/terminal_ui.h"
#include <vector>
#include <memory>
#include <iostream>

int main() {
    try {
        SimulationConfig config;
        TerminalUI::print_header();

        auto engine = std::make_unique<OptimizedMolecularDynamics>(config.box_size);
        auto reporter = std::make_unique<MDReporter>();

        TerminalUI::log_info("Initializing molecules...");
        engine->add_molecule("methane", 10, 10, 10);
        engine->add_molecule("water", 15, 10, 10);

        engine->add_thermal_velocities(config.target_temperature);
        reporter->export_topology(*engine, config.topology_file);

        TerminalUI::log_info("Running equilibration & production...");
        for (int step = 1; step <= config.total_steps; ++step) {
            engine->update();

            if (step % config.export_interval == 0) {
                reporter->export_energy(*engine, config.energy_file, step);
                reporter->export_trajectory(*engine, config.trajectory_file, step);
                TerminalUI::print_progress(step, config.total_steps, engine->calculate_temperature(), 0.0);
            }
        }
        std::cout << std::endl;

        reporter->export_xyz(*engine, "output/final.xyz", config.total_steps);
        TerminalUI::log_info("Analyzing stability...");
        engine->analyze_stability();
        
        TerminalUI::print_footer();

    } catch (const std::exception& e) {
        TerminalUI::log_error(e.what());
        return 1;
    }
    return 0;
}