#pragma once

#include "common/base.h"
#include "solver/Grains.h"
#include <omp.h>
#include <iomanip>
#include <memory>

class Gausspoint;
struct materialPhase;

namespace Polycs{
class polycrystal{
    private:
        double aveGrainSize = 30; //grain size (um)
        Gausspoint* Gpsets  = NULL; 
    
        //the shape of ellipsoid
        Vector3d ell_axis_o; //original ellipsoid axis
        Vector3d ell_axis; //axial length of ellipsoid
        Matrix3d ell_axisb;
        double ell_crit_shape = 25;
        Vector3d ellip_ang; //the rotate angle of the ellipsoid
        bool Iflat = 0; //flag of stretch the ellipsoid (0: yes; 1: no)
        bool Ishape = 0; //flag of individual ellipsoid for grain (1: yes; 0: no)
    
        //Elastic consistent
        Matrix6d elasticModu_SC; //The elastic consistent stiffness
        Matrix6d elasticModu_old;//The initial guessed elastic consistent stiffness
        Matrix6d elasticCompli; //CSC^-1
        Matrix3d therm_expansion_ave; //The average thermal expansion tensor
        Matrix3d therm_expans_ave_old; //The average thermal expansion tensor in last increment
        //visco-plastic consistent
        Matrix5d vPlasticModu_SC; //The visco-plastic stiffness 
        Matrix5d vPlasticModu_old; 
        /* Matrix5d vPlasticCompli_SC; //to Public: The visco-plastic compliance C_VP_SC^-1 */
        Vector6d D_vp_0; //the macro back-extrapolated term (follow Equ[5-41b])
    
        Matrix3d Fij_m; //the macro deformation tensor in grain
        Matrix3d UDWdot_m; //the BC velocity gradient;
        Matrix3d Udot_m; // L matrix
        Matrix3i IUdot; //the flag of known (control by Udot_m) and unknown (calculated by EVPSC) velocity components
        Vector6i IDdot; //the flag of known and unknown strain rate components
        Vector6i ISdot; //the flag of known and unknown stress componets
        /* Matrix3d Dij_m; //to Public: the macro strain rate tensor */
        Matrix3d Wij_m; //the macro rotation rate tensor
        // Udot_m = Dij_m + Wij_m
        Matrix3d Udot_AV;
        Matrix3d Dij_AV; //the average strain rate tensor of all garins
        Matrix3d Dije_AV; //the elastic part
        Matrix3d Dijp_AV; //the vp part
    
        Matrix3d Eps_m = Matrix3d::Zero(); //macro strain
        Matrix3d Sig_m = Matrix3d::Zero(); //macro stress 
        Matrix3d Sig_rate; // stress increment
        Matrix3d Sig_AV; //the average stress tensor of all garins
        Matrix3d Sig_m_old; //macro stress in last increment
    
        //some parameters for error control and iteration control
        double SC_err_m = 0.01; //the error limit of Self-consistent compliance or stiffness
        int SC_iter_m = 30; //the max iteration number of SC
        double errD_m = 0.01; //error limit between the input macro strain rate and output at each iteration 
        double errS_m = 0.01; //error limit of the macro stress
        double err_g_AV = 0.01; //error limit of the average grain stress
    
        Matrix3d Sig_in;
        Matrix3d Dij_in;
        Matrix3d sig_in_AV;
        // Matrices for restoration
        Matrix3d Wij_m_old, Dij_m_old, Dije_AV_old, Dijp_AV_old;
    
    public:
        polycrystal();
    
        vector<grain> g;
        int grains_num = 0;
        int family_count = 0;
        vector<int> modes_count_by_family;
        vector<double> density_by_family, acc_strain_by_family, crss_by_family;
        double error_SC = 0., twin_threshold = 1., temperature_poly = 0., temp_poly_old = 0.;

        Matrix3d Dij_m; //the macro strain rate tensor
        Matrix3d thermal_strain_m;//the macro thermal strain tensor
        Matrix3d ther_strain_m_old;//the macro thermal strain tensor in last increment
        Matrix3d elastic_strain_m = Matrix3d::Zero();//the macro elastic strain tensor
        Matrix3d plastic_strain_m = Matrix3d::Zero();//the macro plastic strain tensor
        Vector5d DVP_AV;
        Matrix5d vPlasticCompli_SC; //The visco-plastic compliance C_VP_SC^-1
        
        vector<Matrix3d> strain_phases;
        vector<Matrix3d> stress_phases;
    
        void add_grains(int phase_id, int grain_count, const vector<Vector4d> &eulerData, materialPhase *mat);
        void ini_Udot_m(Matrix3d);
        void ini_Sig_m(Matrix3d);
        void set_IUdot(Matrix3i);
        void set_ISdot(Vector6i);
        void set_boundary_conditions(Matrix3d, Matrix3d, Matrix3i, Vector6i);
    
        //input the number of grains
        int check_grains_n();
        int add_grain(Vector4d, int tp_id, int mode_id);
        void weightNormalization();
    
        int crysCharacterInitialization(const materialPhase &mat);
        int check_gmode();
    
        //input the normal and Burgers vector of deformation system in a mode
        //input parameters:
        //(MatrixXd) sn, (int) flag of twin or slip, (int) number of systems,(int) mode label
        //needs loop over grains
        int ini_sn(MatrixXd, int, int, int);
        int check_sn();
    
        int cal_aveGrainSize();
        //input the grain size;
    
        int ini_hardening(double, VectorXd, VectorXd, int, int);
        //input the hardening parameters
        //input parameters:
        //double nrsx_in; VectorXd CRSS_p_in; VectorXd hst_in; int modei
        int check_hardening();
    
        int Update_Fij(double);
        int Update_shape();
        Vector3d get_ell_axis(); 
        Vector3d get_ellip_ang(); 
    
        //the singular step according to a certain process
        int EVPSC(int, double);
    
        //calculate the macro&grain elastic compliance
        int Selfconsistent_E(int, double, int); 
    
        //calculate the macro&grain VP compliance
        int Selfconsistent_P(int, double, int); 
    
        // status management
        void save_status();
        void restore_status(bool reset);
        void update_status(double time_incre);
        void update_twin_control();
    
        // calculate the Sig_m and Dij_m
        void Cal_Sig_m(double);  
        // calculate the Sig_g and Dij_g including the elastic and vp part
        double Cal_Sig_g(double);
        void Update_AV(); //update the volume average value
    
        //output
        Vector6d get_Sig_m();
        Vector6d get_Sig_ave();
        Vector6d get_Eps_m();
        void printEuler(fstream &, int phase_id);
        void update_info_by_family();
        void update_phase_info();
        // Set the temperature in polyX and all grains to the given value
        void set_temperature(double);
        // Update the temperature in polyX and all grains by the given increment
        void update_temperature(double time_incre);
    };
}
