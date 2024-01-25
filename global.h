#ifndef GLOBALS_H
#define GLOBALS_H
#include "Polycrystals.h"
#include "Processes.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>
#include <ctime>
#include <Eigen/Dense>

extern double temp_atmosphere, temperature_ref;
extern vector<double> custom_vars;
extern fstream tex_out; //output of the texture
extern fstream density_out; //output of the grain information
extern fstream acc_strain_out;
extern fstream crss_out;
extern fstream ss_out_csv; //output of the macro stress-strain curves
extern fstream ave_ss_out; //output of the average stress-strain curves
extern fstream grain_out; //output of the grain information
extern fstream custom_out; //output of the custom variables

void initial_output_files();
void output_info();
void output_grain_info(int i);

extern Polycs::polycrystal global_polycrys;
extern Procs::Process global_proc;
class Logger;
extern Logger logger;

void update_progress(double progress_f);

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


#endif

