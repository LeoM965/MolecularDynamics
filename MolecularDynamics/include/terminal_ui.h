#ifndef TERMINAL_UI_H
#define TERMINAL_UI_H

#include <string>
#include <string_view>

class TerminalUI {
public:
    static void print_header();
    static void print_progress(int step, int total, double temp, double energy);
    static void log_info(const std::string& message);
    static void log_warning(const std::string& message);
    static void log_error(const std::string& message);
    static void print_footer();
};

#endif
