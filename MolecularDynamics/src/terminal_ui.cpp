#include "../include/terminal_ui.h"
#include <iostream>
#include <iomanip>

void TerminalUI::print_header() {
    std::cout << "================================================\n"
              << "      MOLECULAR DYNAMICS ENGINE PRO v2.0        \n"
              << "================================================\n";
}

void TerminalUI::print_progress(int step, int total, double temp, double energy) {
    const int bar_width = 30;
    if (total <= 0) return;
    float progress = (float)step / total;
    int pos = (int)(bar_width * progress);

    std::cout << "\r[" << step << "/" << total << "] [";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "#";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << (int)(progress * 100.0) << "% "
              << "Temp: " << std::fixed << std::setprecision(1) << temp << "K "
              << std::flush;
}

void TerminalUI::log_info(const std::string& message) {
    std::cout << "[INFO] " << message << "\n";
}

void TerminalUI::log_warning(const std::string& message) {
    std::cout << "[WARN] " << message << "\n";
}

void TerminalUI::log_error(const std::string& message) {
    std::cerr << "[ERROR] " << message << "\n";
}

void TerminalUI::print_footer() {
    std::cout << "\n================================================\n"
              << "           SIMULATION COMPLETE                 \n"
              << "================================================\n";
}
