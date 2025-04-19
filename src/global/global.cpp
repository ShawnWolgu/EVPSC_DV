#include "common/common.h"

// [Simulation Settings]
int texctrl = 0; //print the texture every n steps(0 means only print at the end)
bool update_orientation_required = true; //update the orientation 1:yes 0:no
bool update_shape_required = false; //update the ellipsoid shape 1:yes 0:nobool
bool update_CRSS_required = true; //update the CRSS 1:yes 0:no
bool update_temperature_required = false; //update the temperature 1:yes 0:no

// [Temperature]
double temp_atmosphere = 0.0; // temperature of the atmosphere
double temperature_ref = 0.0; // reference temperature: free of thermal stress used to characterize the thermal expansion strain quantity

// [Some Material Properties]
double rho_material; // the density of the mat
double Cp_material; // the specific heat of the mat
double sigma_e_mat; //the electricity conductivity
double h_ext; // the convection constant between the mat and the atmos
double Surface; //surface of the sample
double sigma_k = 5.67e-8; // the stefan-boltzmann const W/(K^4*m^2)
double V_sample; // the volume of the sample
// [Electricity current control]
double duty_ratio_J = 0.0;//
double Amplitude_J = 0.0;//
double Frequency = 1.0; 
double ref_current_intensity_0 = 0.0;
double ref_current_intensity_1 = 0.0; 
double ref_current_intensity_2 = 0.0; //reference current intensity
double rss_j = 0.0;
double Current_intensity = 0.0;
double bvalue = 0.0;
double shock_int = 0.0;
double shock_fin = 0.0;
double time_acc = 0.0;
double K_ew = 0.0;
int flag_emode = 0;
Matrix3d J_tensor = Matrix3d::Zero();

// [Some global objects]
Logger logger;
Polycs::polycrystal global_polycrys;
Procs::Process global_proc;

// [Output fstreams]
vector<double> custom_vars(12,0.0); //custom variables
fstream tex_out("Tex.out", ios::out); //output of the texture
fstream density_out("Density.csv",ios::out); 
fstream acc_strain_out("Acc_Strain.csv",ios::out);
fstream crss_out("CRSS.csv",ios::out);
fstream ss_out_csv("str_str.csv",ios::out); //output of the macro stress-strain curves
fstream ave_ss_out("ave_str_str.csv",ios::out); //output of the average stress-strain curves
fstream grain_out("grain_info.csv",ios::out); //output of the grain information"
fstream custom_out("custom.csv",ios::out); //output of the custom variables

