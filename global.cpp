#include "global.h"
#include "Polycrystals.h"
#include "Processes.h"
using namespace std;

double temp_atmosphere = 0.0; // temperature of the atmosphere
double temperature_ref = 0.0; // reference temperature: free of thermal stress used to characterize the thermal expansion strain quantity
double rho_material; // the density of the mat //在程序任意地方调用，然后在
double Cp_material; // the specific heat of the mat
double sigma_e_mat; //the electricity conductivity
double h_ext; // the convection constant between the mat and the atmos
double Surface; //surface of the sample
double sigma_k = 5.67e-8; // the stefan-boltzmann const W/(K^4*m^2)
double V_sample; // the volume of the sample
double duty_ratio_J = 0.0;//
double Amplitude_J = 0.0;//
double Frequency = 0.0; 
Logger logger;
Polycs::polycrystal global_polycrys;
Procs::Process global_proc;

vector<double> custom_vars(10,0.0); //custom variables
fstream tex_out("Tex.out", ios::out); //output of the texture
fstream density_out("Density.csv",ios::out); 
fstream acc_strain_out("Acc_Strain.csv",ios::out);
fstream crss_out("CRSS.csv",ios::out);
fstream ss_out_csv("str_str.csv",ios::out); //output of the macro stress-strain curves
fstream ave_ss_out("ave_str_str.csv",ios::out); //output of the average stress-strain curves
fstream grain_out("grain_info.csv",ios::out); //output of the grain information"
fstream custom_out("custom.csv",ios::out); //output of the custom variables

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

