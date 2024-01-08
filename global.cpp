#include "global.h"
#include "Polycrystals.h"
#include "Processes.h"
using namespace std;

double temp_atmosphere = 0.0; // temperature of the atmosphere
double temperature_ref = 0.0; // reference temperature: free of thermal stress
double rho_material; // the density of the mat
double Cp_material; // the specific heat of the mat
double sigma_e_mat; //the electricity conductivity
double h_ext; // the convection constant between the mat and the atmos
double S; //surface of the sample
double sigma_k = 5.67e-8; // the stefan-boltzmann const W/(K^4*m^2)
Logger logger;
Polycs::polycrystal global_polycrys;
Procs::Process global_proc;

void update_progress(double progress_f)
{
    const int bar_width = 70;
    int bar_position = (int)(bar_width * progress_f);

    std::cout << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < bar_position) {
            std::cout << "=";
        } else if (i == bar_position) {
            std::cout << ">";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << (int)(progress_f * 100) << "%\r";
    std::cout.flush();
}
