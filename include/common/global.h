#pragma once

#include "common/base.h"
#include <chrono>
#include <thread>
#include <ctime>

namespace Procs{ class Process; }
// [Simulation Settings]
extern int texctrl;
extern bool update_orientation_required, update_shape_required, update_CRSS_required, update_temperature_required;
extern double temp_atmosphere, temperature_ref;

// [Some Material Properties]
extern double rho_material, Cp_material, sigma_e_mat, h_ext, Surface, V_sample, sigma_k;
extern double duty_ratio_J, Amplitude_J, Frequency, ref_current_intensity_0, ref_current_intensity_1, ref_current_intensity_2;
extern double Current_intensity;
extern double deformation_rate;
extern double deformation_max;
extern double bvalue;
extern double shock_int;
extern double shock_fin;
extern double rss_j;
extern int flag_emode;
extern Matrix3d J_tensor;
extern double K_ew;
extern double time_acc;

// [Output fstreams]
extern vector<double> custom_vars;
extern fstream tex_out; //output of the texture
extern fstream density_out; //output of the grain information
extern fstream acc_strain_out;
extern fstream crss_out;
extern fstream ss_out_csv; //output of the macro stress-strain curves
extern fstream ave_ss_out; //output of the average stress-strain curves
extern fstream grain_out; //output of the grain information
extern fstream custom_out; //output of the custom variables

// [Some global objects]
extern Polycs::polycrystal global_polycrys;
extern Procs::Process global_proc;
class Logger;
extern Logger logger;

// Define the Logger class
class Logger {
public:
    Logger() : outfile("EVPSC_log.txt") {}

    // Define the message hierarchy levels inside the Logger class
    enum class MessageLevel {
        DEBUG,
        NOTICE,
        WARN,
        INFO,
        ERROR,
    };

    // Log a message at a specific message level
    void log(MessageLevel level, const std::string& message) {
        // Get the current time
        std::time_t t = std::time(nullptr);
        std::string timestamp = std::asctime(std::localtime(&t));
        timestamp.pop_back();  // Remove the newline character from the end

        // Output the message to the console
        if (console_level_ <= level) {
            std::cout << "[" << timestamp << "] ";
            printPrefix(level);
            std::cout << message << std::endl;
        }

        // Output the message to the log file
        if (file_level_ <= level) {
            outfile << "[" << timestamp << "] ";
            printPrefix(level, outfile);
            outfile << message << std::endl;
        }
    }

    void log(MessageLevel level, const Eigen::MatrixXd& message) {
        // Get the current time
        std::time_t t = std::time(nullptr);
        std::string timestamp = std::asctime(std::localtime(&t));
        timestamp.pop_back();  // Remove the newline character from the end
        Eigen::IOFormat LogFmt(4, 0, ", ", "\n", "["+timestamp+"] "+getPrefix(level)+ "[", "]");
        // Output the message to the console
        if (console_level_ <= level) {
            std::cout << message.format(LogFmt) << std::endl;
        }
        // Output the message to the log file
        if (file_level_ <= level) {
            outfile << message.format(LogFmt) << std::endl;
        }
    }

    void log(MessageLevel level, const std::vector<double>& message) {
        // Get the current time
        std::time_t t = std::time(nullptr);
        std::string timestamp = std::asctime(std::localtime(&t));
        timestamp.pop_back();  // Remove the newline character from the end
        // Output the message to the console
        if (console_level_ <= level) {
            std::cout << "[" << timestamp << "] ";
            printPrefix(level);
            for (int i = 0; i < message.size(); i++) {
                std::cout << message[i] << " ";
            }
            std::cout << std::endl;
        }
        // Output the message to the log file
        if (file_level_ <= level) {
            outfile << "[" << timestamp << "] ";
            printPrefix(level, outfile);
            for (int i = 0; i < message.size(); i++) {
                outfile << message[i] << " ";
            }
            outfile << std::endl;
        }
    }


    // Overloaded functions for different message levels
    void debug(const std::string& message) { log(MessageLevel::DEBUG, message); }
    void info(const std::string& message) { log(MessageLevel::INFO, message); }
    void notice(const std::string& message) { log(MessageLevel::NOTICE, message); }
    void warn(const std::string& message) { log(MessageLevel::WARN, message); }
    void error(const std::string& message) { log(MessageLevel::ERROR, message); }
    void debug(const Eigen::MatrixXd& message) { log(MessageLevel::DEBUG, message); }
    void info(const Eigen::MatrixXd& message) { log(MessageLevel::INFO, message); }
    void notice(const Eigen::MatrixXd& message) { log(MessageLevel::NOTICE, message); }
    void warn(const Eigen::MatrixXd& message) { log(MessageLevel::WARN, message); }
    void error(const Eigen::MatrixXd& message) { log(MessageLevel::ERROR, message); }
    void debug(const std::vector<double>& message) { log(MessageLevel::DEBUG, message); }
    void info(const std::vector<double>& message) { log(MessageLevel::INFO, message); }
    void notice(const std::vector<double>& message) { log(MessageLevel::NOTICE, message); }
    void warn(const std::vector<double>& message) { log(MessageLevel::WARN, message); }
    void error(const std::vector<double>& message) { log(MessageLevel::ERROR, message); }
    // Set the message level filter for console output
    void setConsoleLevel(MessageLevel level) { console_level_ = level; }

    // Set the message level filter for file output
    void setFileLevel(MessageLevel level) { file_level_ = level; }

private:
    std::ofstream outfile;
    MessageLevel console_level_ = MessageLevel::INFO;
    MessageLevel file_level_ = MessageLevel::DEBUG;

    // Helper function to print the message prefix based on the message level
    void printPrefix(MessageLevel level, std::ostream& stream = std::cout) {
        switch (level) {
            case MessageLevel::DEBUG:
                stream << "[ DEBUG] ";
                break;
            case MessageLevel::INFO:
                stream << "[ INFO ] ";
                break;
            case MessageLevel::NOTICE:
                stream << "[NOTICE] ";
                break;
            case MessageLevel::WARN:
                stream << "[ WARN ] ";
                break;
            case MessageLevel::ERROR:
                stream << "[ ERROR] ";
                break;
            default:
                break;
        }
    }
    std::string getPrefix(MessageLevel level) {
        switch (level) {
            case MessageLevel::DEBUG:
                return "[ DEBUG] ";
            case MessageLevel::INFO:
                return "[ INFO ] ";
            case MessageLevel::NOTICE:
                return "[NOTICE] ";
            case MessageLevel::WARN:
                return "[ WARN ] ";
            case MessageLevel::ERROR:
                return "[ ERROR] ";
            default:
                return "";
        }
    }
};

