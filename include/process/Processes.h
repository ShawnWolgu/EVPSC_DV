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

    Matrix3i IUDWdot; //the flag of known (control by Udot_input) and unknown (calculated by EVPSC) velocity components
    Vector6i IDdot; //the flag of known and unknown strain rate components
    Vector6i ISdot; //the flag of known and unknown stress rate componets

    double Eincr; //the increment of strain in every istep
    int Ictrl;
    int Nsteps; // total steps
    int istep; // current step
    double temperature_input; // /(K) temperature 
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

    //get the known and unknown flag tensor in a process file
    void get_IUdot(Matrix3i);
    void get_ISdot(Vector6i);

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
