#pragma once

#include "common/base.h"
#include <iomanip>
#include <csignal>

namespace Procs{

class Process{
private:
    Matrix3d UDWdot_input; //the BC gradient in process file
    Matrix3d Ddot_input; //the strain rate calculated by the velocity gradient in process file
    Matrix3d Sdot_input; //the stress rate tensor in process file
    Vector3d Efdot_input; //the electric field vector (dot) in process files
    Vector3d Eddot_input; //the electric displacement vector (dot) in process files

    Matrix3i IUDWdot; //the flag of known (control by Udot_input) and unknown (calculated by EVPSC) velocity components
    Vector6i IDdot; //the flag of known and unknown strain rate components
    Vector6i ISdot; //the flag of known and unknown stress rate componets
    Vector3i IEfdot; //the flag of known and unknown electric field vector components
    Vector3i IEddot; //the flag of known and unknown electric displacement vector components

    double Eincr; //the increment of strain in every istep
    int Ictrl;
    int Nsteps; // total steps
    int istep; // current step
    double temperature_input; // /(K) temperature 
    double tempK_rate = 0.0; // /(K/s) temperature rate
    double tempK_end = 0.0; // /(K) end temperature
    double max_timestep;

    int texctrl; //print the texture every n steps(0 means only print at the end)

public:
    Process();
    ~Process();

    //get the total steps and increment of a process from files
    //Vector4d 0: Nsteps; 1: Ictrl; 2: Eincr; 3: Temperature;
    void load_ctrl(Vector4d);

    //get the velocity gradient in a process file
    void get_Udot(Matrix3d);
    //get the Cauchy stress tensor in a process file
    void get_Sdot(Matrix3d);
    //get the electric field vector in a process file
    void get_Efdot(Vector3d Min);
    //get the electric displacement vector in a process file
    void get_Eddot(Vector3d Min);

    //get the known and unknown flag tensor in a process file
    void get_IUdot(Matrix3i);
    void get_ISdot(Vector6i);
    void get_IEfdot(Vector3i Vin);
    void get_IEddot(Vector3i Vin);

    void set_tempK_control(double rate, double end_temp);
    double calculate_current_intensity(double time) const;
    void timestep_control();
    void loading(Polycs::polycrystal &);

    /////
    //Output functions:
    //output of stress&strain curves  
    void Out_sscurves(Polycs::polycrystal &);

    //output of texture
    void Out_texture(Polycs::polycrystal &, int);
    void Out_texset(int);

    //output of grain information
    void init_grain_info(Polycs::polycrystal &, int);
    void Out_grain_info(Polycs::polycrystal &, int);

};

}
