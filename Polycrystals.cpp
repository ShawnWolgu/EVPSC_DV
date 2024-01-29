#include "Polycrystals.h"
#include "Toolbox.h"
#include "global.h"
#include <memory>

using namespace Polycs;
using namespace std;

polycrystal::polycrystal()
{

    //initial the macro stress&strain
    Eps_m = Matrix3d::Zero();
    Sig_m = Matrix3d::Zero();

    //initial the shape of ellipsoid
    ell_axis = Vector3d::Ones();
    ellip_ang << 90,90,90;
    ell_axisb = Euler_trans(ellip_ang);
    Fij_m = Matrix3d::Identity();

    //initial the VP consistent
    M_VP_SC = 1e-40 * Matrix5d::Identity();
    C_VP_SC = M_VP_SC.inverse();
    D0 = Vector6d::Zero();

    Msup<<1,0,0,0,0,0,
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,2,0,0,
    0,0,0,0,2,0,
    0,0,0,0,0,2;

    Vector10d xth,xph,wth,wph;
    //integral points and weights
    Integralpoint3 alpha, aww;
    Integralpoint6 aa6, aaww6; //coordinate and weigts in Fourier space 
    Integralpoint1 ww;

    Gpsets = new Gausspoint[11];

    for(int Gpcase = 0; Gpcase < 11; Gpcase++ ){
        switch(Gpcase){
            case 0:
            xth << 4.71236594e-02,0.241774723e0,0.565131843e0,0.968887568e0,1.37937832e0,
            1.76221442e0,2.17270517e0,2.57646084e0,2.89981818e0,3.09446883e0;
            wth << 0.120191820e0,0.264987558e0,0.373805553e0,0.420841277e0,0.390970200e0,
            0.390970260e0,0.420841366e0,0.373805553e0,0.264987499e0,0.120192111e0;
            break;
            case 1:
            xth << 1.57080423e-02,0.144995824e0,0.425559640e0,0.829968274e0,1.31460333e0,
            1.82698941e0,2.31162453e0,2.71603298e0,2.99659705e0,3.12588477e0;
            wth << 5.41692823e-02,0.207461149e0,0.348739326e0,0.452716887e0,0.507709801e0,
            0.507709682e0,0.452716798e0,0.348738998e0,0.207461327e0,5.41692935e-02;
            break; 
            case 2:
            xth << 3.76990959e-02,0.198626831e0,0.483041346e0,0.871647120e0,1.32964790e0,
            1.81194484e0,2.26994562e0,2.65855122e0,2.94296598e0,3.10389376e0;
            wth << 9.68142375e-02,0.224478707e0,0.341134071e0,0.430180043e0,0.478189558e0,
            0.478189170e0, 0.430180043e0, 0.341134191e0, 0.224478647e0, 9.68143344e-02;
            break;
            case 3:
            xth << 3.45576368e-02,0.187556863e0,0.468425453e0,0.859980166e0,1.32527423e0,
            1.81631863e0,2.28161263e0,2.67316723e0,2.95403576e0,3.10703516e0;
            wth << 8.95763785e-02,0.217725381e0,0.341026783e0,0.435772508e0,0.486694932e0,
            0.486695170e0,0.435772508e0,0.341026902e0,0.217725128e0,8.95764604e-02;
            break;
            case 4:
            xth << 3.14158052e-02,0.177928671e0,0.457155794e0,0.851592362e0,1.32222414e0,
            1.81936860e0,2.29000044e0,2.68443704e0,2.96366405e0,3.1101768e0;
            wth << 8.26927349e-02,0.213228315e0,0.342008322e0,0.440196186e0,0.492670894e0,
            0.492670983e0,0.440195888e0,0.342008322e0, 0.213227972e0, 8.26930404e-02;
            break;
            case 5:
            xth << 2.98452154e-02,0.173592165e0,0.452448040e0,0.848216832e0,1.32101476e0,
            1.82057810e0,2.29337597e0,2.68914461e0,2.96800065e0,3.11174774e0;
            wth << 7.93928578e-02,0.211627841e0,0.342669785e0,0.442057431e0,0.495048553e0,
            0.495048642e0,0.442057490e0,0.342670023e0,0.211627468e0,7.93929026e-02;
            break;
            case 6:
            xth << 2.67036632e-02,0.165752888e0,0.444431901e0,0.842614472e0,1.31902647e0,
            1.82256627e0,2.29897833e0,2.69716072e0,2.97583985e0,3.11488938e0;
            wth << 7.30879456e-02,0.209402516e0,0.344104946e0,0.445234656e0,0.498966068e0,
            0.498966306e0,0.445234746e0, 0.344104946e0,0.209402665e0,7.30878562e-02;
            break;
            case 7:
            xth << 2.67036632e-02,0.165752888e0,0.444431901e0,0.842614472e0,1.31902647e0,
            1.82256627e0,2.29897833e0,2.69716072e0,2.97583985e0,3.11488938e0;
            wth << 7.30879456e-02,0.209402516e0,0.344104946e0,0.445234656e0,0.498966068e0,
            0.498966306e0,0.445234746e0,0.344104946e0,0.209402665e0,7.30878562e-02;
            break;
            case 8:
            xth <<2.43473575e-02,0.160516247e0,0.439386278e0,0.839168847e0,1.31781363e0,
            1.82377899e0,2.30242372e0,2.70220637e0,2.98107672e0,3.11724544e0;
            wth << 6.86219111e-02,0.208388865e0,0.345189095e0,0.447236270e0,0.501360059e0,
            0.501359940e0,0.447236151e0,0.345189214e0,0.208388969e0,6.86219335e-02;
            break;
            case 9:
            xth << 2.19910536e-02,0.155757755e0,0.434985727e0,0.836206555e0,1.31677616e0,
            1.82481658e0,2.30538607e0,2.70660710e0,2.98583508e0,3.11960149e0;
            wth << 6.43825606e-02,0.207786217e0,0.346235514e0,0.448981822e0,0.503410578e0,
            0.503410578e0,0.448981792e0,0.346235693e0,0.207785636e0,6.43827692e-02;
            break;
            case 11:
            xth << 2.04204638e-02,0.152822554e0,0.432348520e0,0.834448099e0,1.31616223e0,
            1.82543063e0,2.30714464e0,2.70924401e0,2.98877001e0,3.12117243e0;
            wth << 6.16818815e-02,0.207559645e0,0.346902698e0,0.450027168e0,0.504624724e0,
            0.504624426e0,0.450027317e0,0.346902847e0,0.207559645e0,6.16819337e-02;
            break;
        }
        wph = wth;
        xph = xth;	
        double sinth, costh, simbtet;
        Matrix3d aa, aaww;
        int ny;
        for (int ith = 0; ith < Intn; ith++){
            sinth = sin(xth(ith));
            costh = cos(xth(ith));
            simbtet = wth(ith) * sinth / (2.0*M_PI);

            for (int iph = 0; iph < Intn; iph++){
                ny = iph + ith * Intn;
                ww(ny) = simbtet*wph(iph);
                alpha(0,ny) = sinth*cos(xph(iph));
                alpha(1,ny) = sinth*sin(xph(iph));
                alpha(2,ny) = costh;
                //
                for(int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++){
                        aa(i,j) = alpha(i,ny) * alpha(j,ny);
                        aaww(i,j) = aa(i,j) * ww(ny);
                    }
                aa6.col(ny) = voigt(aa);
                aaww6.col(ny) = voigt(aaww);
                for(int i = 0; i < 3; i++)
                    aww(i,ny) = alpha(i,ny) * ww(ny);
            }
        }
        Gpsets[Gpcase].Gpaa6 = aa6;
        Gpsets[Gpcase].Gpaaww6 = aaww6;
        Gpsets[Gpcase].Gpalpha = alpha;
        Gpsets[Gpcase].Gpaww = aww;
        Gpsets[Gpcase].Gpww = ww;
    }
}

void polycrystal::set_BC_const(Matrix3d udot_input, Matrix3d sig_input, Matrix3i iudot_input, Vector6i isig_input){
    Sig_m = sig_input;
    ISdot = isig_input;
    set_IUdot(iudot_input);
    ini_Udot_m(udot_input);
}

void polycrystal::ini_Udot_m(Matrix3d Udot_input)
{   
    UDWdot_m = Udot_input;
    Dij_m = 0.5*(UDWdot_m + UDWdot_m.transpose());
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
	if (IUdot(i,j) == 2) Dij_m(i,j) = UDWdot_m(i,j);
	if (IUdot(i,j) == 3) Dij_m(i,j) = 0.0;
    }
    Wij_m = UDWdot_m - Dij_m;
    
    Dij_AV = Dij_m;
    Dije_AV = Dij_m;
    Dijp_AV = Matrix3d::Zero();
}

void polycrystal::ini_Sig_m(Matrix3d Min){Sig_m = Min;}

void polycrystal::set_IUdot(Matrix3i Min)
{   
    IUdot = Min;
    for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++)
	if(IUdot(i,j) > 1 && IUdot(j,i) != 0){
	    logger.error("IUdot configuration not support!");
	    exit(1);
	}
    IDdot(0)=IUdot(0,0);
    IDdot(1)=IUdot(1,1);
    IDdot(2)=IUdot(2,2);
    IDdot(3)=int((IUdot(1,2)+IUdot(2,1)) == 2);
    IDdot(4)=int((IUdot(0,2)+IUdot(2,0)) == 2);
    IDdot(5)=int((IUdot(0,1)+IUdot(1,0)) == 2);
}
void polycrystal::set_ISdot(Vector6i Min){ISdot = Min;}

int polycrystal::grains_n(int n)
{
    grains_num = n;
    g = new grain[n*20];
    //change the number of grains
    for(int i = 0; i < n; i++)  g[i].grain_i = i;
    return 0;
}

// New grain born based on the parent grain. Under tested
int polycrystal::add_grain(Vector4d vin, int tp_id, int mode_id){
    logger.debug("Add a new grain because of the twin sys " + to_string(mode_id) + " in grain " + to_string(tp_id));
    /* grain* g_tmp = new grain(g[tp_id]); */
    grains_num++;
    int id = grains_num-1;
    g[id] = *new grain(g[tp_id]);
    g[id].grain_i = id;
    g[id].ini_euler_g(vin);
    g[id].ini_gmode_g(g[tp_id]);
    /* g[id].lat_hard_mat = g[tp_id].lat_hard_mat; */
    g[id].set_stress_g(g[tp_id].get_stress_g());
    /* g[id].print_latent_matrix(); */
    if (mode_id != -1){
        if (Twin* twinPtr = dynamic_cast<Twin*>(g[id].gmode[mode_id])){
            twinPtr->set_parent(tp_id);
        }
    }
    /* for(int j = 0; j < g[id].modes_num; ++j){ */
    /*     g[id].gmode[j]->update_status(g[id],0.0); */
        /* g[id].gmode[j]->print(); */
    /* } */
    return id;
}

int polycrystal::check_grains_n()
{
    logger.debug("the number of grains: "+std::to_string(grains_num));
    return grains_num; //print the number of grains
}

void polycrystal::ini_euler(Vector4d vin, int i){g[i].ini_euler_g(vin);} 

void polycrystal::Norm_weight()
{
    double total_w = 0;
    for(int i = 0; i < grains_num; i++)
        total_w += g[i].get_weight_g();
    
    for(int i = 0; i < grains_num; i++){
        g[i].set_weight_g(g[i].get_weight_g()/total_w);
        g[i].weight_ref = g[i].get_weight_g();
    }

}

int polycrystal::ini_cry(json &j){
    crysym = j["crysym"];   Miller_n = j["Miller_n"]; 
    Cdim = to_vector(j, "Cdim", 3);  Cang = to_vector(j, "Cang", 3);
    Trans_Miller = to_matrix(j, "Trans_Miller", 3, Miller_n);
    Mabc = to_matrix(j, "Mabc", 3, 3);
    return 0;
}

int polycrystal::ini_cry(string strin, VectorXd vin)
{
    //transform str to lower case
    for (int i = 0; i < strin.size(); i++)
		strin[i] = tolower(strin[i]);
    crysym = strin;    
    Cdim = vin(seq(0,2));
    Cang = vin(seq(3,5)) / 180 * M_PI; //Converting degrees to radians
    
    //calculate conversion matrix of Miller indices according to the crysym
    if(!crysym.compare("hexag"))
    {
        Miller_n = 4;
        MatrixXd Mtemp(3,4);
        Trans_Miller = Mtemp;
        Trans_Miller <<
        1, 0, -1, 0,
        0, 1, -1, 0,
        0, 0,  0, 1;
    }
    else if(!crysym.compare("cubic")) 
    {
        Miller_n = 3;
        MatrixXd Mtemp(3,3);
        Trans_Miller = Mtemp;
        Trans_Miller <<
        1, 0,  0,
        0, 1,  0,
        0, 0,  1;
    }

    Mabc(0,0)=sin(Cang(1));
    Mabc(1,0)=0.;
    Mabc(2,0)=cos(Cang(1));
    Mabc(0,1)=(cos(Cang(2))-cos(Cang(0))*cos(Cang(1)))/sin(Cang(1));
    Mabc(2,1)=cos(Cang(0));
    Mabc(1,1)=sqrt(1.0-pow(Mabc(0,1),2)-pow(Mabc(2,1),2));
    Mabc(0,2)=0.;
    Mabc(1,2)=0.;
    Mabc(2,2)=1.;

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            Mabc(i,j) = Cdim(j) * Mabc(i,j);

    return 0;
}

int polycrystal::get_Millern(){return Miller_n;}

void polycrystal::check_cry()
{
    logger.debug("the crysym: "+crysym);
    logger.debug("the cdim and cang:");
    logger.debug(Cdim.transpose());
    logger.debug(Cang.transpose());
}

void polycrystal::ini_Cij6(MatrixXd Min)
{   
    Cij6 = Min;
    voigt(Cij6,Cijkl);
}

int polycrystal::check_Cij6()
{
    logger.debug("elastic constant:");
    logger.debug(Cij6);
    return 0;
}

int polycrystal::ini_therm(Vector6d vin)
{
    therm = vin;
    for(int i = 0; i < grains_num; i++)
    {
        g[i].therm_expansion_config(therm);
    }
    return 0;
}

int polycrystal::check_therm()
{
    logger.debug("Thermal coefficient:");
    logger.debug(therm.transpose());
    return 0;
}

int polycrystal::ini_gmode(int n)
{
    for(int i = 0; i < grains_num; i++)
    {
        g[i].ini_gmode_g(n);
    }
    return 0;
}

int polycrystal::ini_gmode(json &j)
{
    for(int i = 0; i < grains_num; i++)
    {
        g[i].ini_gmode_g(j);
    }
    return 0;
}

int polycrystal::check_gmode()
{
    for(int i = 0; i < grains_num; i++)
    {
	logger.debug("the number of modes in Grain "+std::to_string(i)+":");
	logger.debug(std::to_string(g[i].check_gmode_g()));
    }
    return 0;
}

void polycrystal::ini_from_json(json &sx_json){
    ini_cry(sx_json);
    ini_Cij6(to_matrix(sx_json, "Cij6", 6, 6));
    ini_therm((Vector6d)(to_vector(sx_json, "therm", 6)));
    ini_GZ(sx_json["GZ"]);
    //ini_heateffect//自己写
    ini_gmode(sx_json);
    family_count = sx_json["family_num"];
    VectorXd temp_vec = to_vector(sx_json, "modes_count_by_family", family_count);
    for (int i = 0; i < family_count; i++) {
        modes_count_by_family.push_back(temp_vec(i));
        density_by_family.push_back(0.0);
        acc_strain_by_family.push_back(0.0);
        crss_by_family.push_back(0.0);
    }
    for(int i = 0; i < grains_num; ++i) g[i].set_lat_hard_mat();
    g[0].print_latent_matrix();
	for(int i = 0; i < grains_num; ++i){
	   for(int j = 0; j < g[i].modes_num; ++j){
	       g[i].gmode[j]->update_status(g[i],0.0);
	   }
	}
}


int polycrystal::ini_GZ(double x)
{
    GZ = x;
    return 0;
}

int polycrystal::check_hardening()
{
    for(int i = 0; i < grains_num; i++)
    {
	logger.debug("the hardening parameters in Grain "+std::to_string(i)+":");
        g[i].check_hardening_g();
    }    
    return 0;
}

int polycrystal::Selfconsistent_E(int Istep, double ERRM, int ITMAX)
{
    // Calculate the Elastic stiffness in Jaumann rate in all grains
    //(Rotate from crystal to Sample axes according to the euler angle)
    // and sum with the weight
    // the result is CUB
    Matrix6d CUB; // CUB is the volume average Elastic stiffness of all grains
    CUB = Matrix6d::Zero();
    double C4SA[3][3][3][3];  // Elastic stiffness Rotate from crystal to Sample axes
    double C4SAS[3][3][3][3]; // ...in Jaumann rate
    Matrix3d sig_g;

    // Calculate CUB
    for(int G_n = 0; G_n < grains_num; G_n++)
    {
        // Rotate the tensor Cijkl in grain to Sample Axes
        Matrix3d euler_g = g[G_n].get_Euler_M_g();
        voigt(rotate_C66(Cij6, euler_g.transpose()), C4SA);
        g[G_n].Update_Cij6_SA_g(voigt(C4SA));
        sig_g = g[G_n].get_stress_g();
        // the elastic stiffness invovling Jaumann rate
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < 3; k++)
                    for(int l = 0; l < 3; l++)
                        C4SAS[i][j][k][l] = C4SA[i][j][k][l] - sig_g(i,j)*Iij(k,l);
        g[G_n].Update_Mij6_J_g(Chg_basis6(C4SAS).inverse());
        //store the Jaumann rate elastic stiffness in grains
        CUB += Chg_basis6(C4SAS) * g[G_n].get_weight_g();
    }

    if(Istep == 0)  CSC = CUB;  //first step, use the volume average
    SSC = CSC.inverse();

    //loop to make the guessed elastic stiffness CSC to the Eshelby calculated CNEW 
    Matrix6d SSC_new;

    int IT = 0; //loop flag
    double RER = 2*ERRM;
    while((RER >= ERRM) && (IT < ITMAX))
    {
        IT++;
        SSC_new = Matrix6d::Zero();
        Chg_basis(CSC, C4SA);

        Matrix6d Ctilde;
        Matrix6d S66; //the Sijkl Equ[5-33] in sample axes in Manual 7d
        double R4_SA[3][3][3][3]; //PIijkl
        double RSinv_SA[3][3][3][3];
        Matrix6d S66inv; 
        Matrix3d axisb_t;
        Vector3d axis_t;
        Matrix6d Me_g; 

        //solve the eshelby tensor in the common ellipsoid
        if(Ishape == 0)
        {
            /* Rotate the macro Elastic stiffness from Sample axes to the ellipsoid axes; */
            /* if Ishape = 0, which means the grains share the same ellipsoid axes  */
            axisb_t = ell_axisb; 
            axis_t = ell_axis; 
            /* Rotate the self consistent Elastic stiffness (CSC) */
            /* from sample axes to ellipsoid axes */
            Matrix6d C66 = rotate_C66(voigt(C4SA), axisb_t.transpose());
            //calculate the Eshelby tensor
            double S4_EA[3][3][3][3] = {0}; 
            double R4_EA[3][3][3][3] = {0};
            int case_c = Eshelby_case(axis_t);
            g[0].Eshelby_E(S4_EA,R4_EA,axis_t,C66,Gpsets[case_c].Gpaa6,Gpsets[case_c].Gpaaww6,Gpsets[case_c].Gpalpha);
            //rotate the eshelby tensor back to the sample axes
            S66 = voigttoB6(rotate_C66(voigt(S4_EA), axisb_t));
            rot_4th(R4_EA, axisb_t, R4_SA);
            // Calculate Ctilde refer to Equ[5-33] in Manual 7d
            // M~=(I-S)^-1 * S * M
            // C~ = (M~)^-1 = C * (S^-1 - I)
            S66inv = S66.inverse(); // 6x6 of S^-1
            Matrix6d S66inv_I = S66inv - Matrix6d::Identity(); // 6x6 of (S^-1 - I)
            Ctilde = CSC * S66inv_I;
        }

        Matrix6d B_g_ave = Matrix6d::Zero(); //<B_g>

        for(int G_n = 0; G_n < grains_num; G_n++)
        {
            //solve the eshelby tensor in all the ellipsoid of grain 
            if(Ishape == 1)
            {
                //rotate the macro Elastic stiffness from sample axes to the ellipsoid axes;
                // if Ishape = 1, which means the ellipsoid axes varies in grains,
                axis_t = g[G_n].get_ell_axis_g();
                axisb_t = g[G_n].get_ell_axisb_g();
                Matrix6d C66 = rotate_C66(CSC, axisb_t.transpose());
                //calculate the Eshelby tensor
                double S4_EA[3][3][3][3] = {0}; 
                double R4_EA[3][3][3][3] = {0};
                int case_c = Eshelby_case(axis_t);
                g[G_n].Eshelby_E(S4_EA,R4_EA,axis_t,C66,Gpsets[case_c].Gpaa6,Gpsets[case_c].Gpaaww6,Gpsets[case_c].Gpalpha);
                //rotate the eshelby tensor back to the sample axes
                S66 = rotate_C66(voigt(S4_EA), axisb_t);
                rot_4th(R4_EA, axisb_t, R4_SA);
                // Calculate Ctilde refer to Equ[5-33] in Manual 7d
                // M~=(I-S)^-1 * S * M
                // C~ = (M~)^-1 = C * (S^-1 - I)
                S66inv = S66.inverse(); // 6x6 of S^-1
                Matrix6d S66inv_I = S66inv - Matrix6d::Identity();
                Ctilde = CSC * S66inv_I;
            }

            //store some matrix into grain
            g[G_n].Update_Metilde_g(Ctilde.inverse()); //store the C~
            double S66inv4th[3][3][3][3] = {0};
            Chg_basis(S66inv,S66inv4th);
            mult_4th(R4_SA,S66inv4th,RSinv_SA);
            g[G_n].Update_RSinv_C_g(RSinv_SA); //store the PI*(S^-1)

            //Calculate the localization tensor B_g, _g means the value depends on the grain
            //refer to Equ[5-35] in Manual 7d
            // B_g = (M_g + M~)^-1 * (M_ + M~)
            Matrix6d Metilde = Ctilde.inverse();
            Me_g = g[G_n].get_Mij6_J_g();
            Matrix6d Part1 = Me_g + Metilde;
            Matrix6d Part1_inv = Part1.inverse();
            Matrix6d Part2 = SSC + Metilde; 
            Matrix6d B_g = Part1_inv * Part2;
            B_g_ave += B_g * g[G_n].get_weight_g();
            //Calculate the New elastic consistent stiffness CNEW
            //refer to Equ[5-40a] in Manual 7d
            SSC_new += Me_g * B_g * g[G_n].get_weight_g();
        } //loop over grains
        SSC_new = SSC_new * B_g_ave.inverse();
        //
        //error between CSC and CNEW
        RER=Errorcal(SSC,SSC_new);
        SSC = 0.5*(SSC_new+SSC_new.transpose());
        CSC = SSC.inverse();
        logger.notice("**Error in ESC iteration "+to_string(IT)+": \t"+to_string(RER));
        if(isnan(RER)) return 1;
    } //while loop
    // Thermal Expansion Part
    Matrix6d B_g_ret_ave = Matrix6d::Zero(); //<B_g_retro>
    Vector6d B_g_ret_thrmexp = Vector6d::Zero(); //<B_g_retro*alpha_g>
    for(int G_n = 0; G_n < grains_num; G_n++){
        Matrix6d Me_tilde_g = g[G_n].get_Metilde_g();
        Matrix6d Me_g = g[G_n].get_Mij6_J_g();
        Matrix6d B_g_retro = (Me_tilde_g + SSC) * (Me_tilde_g + Me_g).inverse();
        B_g_ret_ave += B_g_retro * g[G_n].get_weight_g();
        Vector6d therm_exp_g_basis = Chg_basis6(g[G_n].get_therm_expansion());
        B_g_ret_thrmexp += B_g_retro * therm_exp_g_basis * g[G_n].get_weight_g();
    }
    therm_expansion_ave = Chg_basis(B_g_ret_ave.inverse() * B_g_ret_thrmexp);
    //SSC = Msup * Btovoigt(CSC.inverse());
    return 0;
}

int polycrystal::Selfconsistent_P(int Istep, double ERRM, int ITMAX)
{
    Matrix5d MNEW; //VP compliance updated in every do-while
    Vector5d D0_new; //the back-extrapolated term updated in every do-while
    Matrix5d B_g_ave; // <B_g> in Equ[5-41a] average of B_g
    Matrix5d B_g_switch_ave;
    Vector5d B_switch_d0_ave;
    int IT = 0; //loop flag
    double RER = 2*ERRM;
    while((RER >= ERRM) && (IT < ITMAX))
    {   
        IT++;
        MNEW = Matrix5d::Zero();
        B_g_ave = Matrix5d::Zero();
        B_g_switch_ave = Matrix5d::Zero();
        B_switch_d0_ave = Vector5d::Zero();
        D0_new = Vector5d::Zero();
        Matrix5d Mtilde;
        Matrix5d S55;//the Sijkl Equ[5-33] in sample axes in Manual 7d
        double R4_SA[3][3][3][3]; //PIijkl
        double RSinv_SA[3][3][3][3];
        Matrix5d S55_inv;
        Matrix5d R55;
        Matrix3d axisb_t;
        Vector3d axis_t;
        //solve the eshelby tensor in the common ellipsoid
        if(Ishape == 0)
        {
            /* Rotate the macro VP stiffness (C_VP_SC) from sample axes to the ellipsoid axes; */
            /* if Ishape = 0, which means the ellipsoid axes keep unchanged in grains, */
            /* the transform matrix should be taken from the polycrystal  */
            axisb_t = ell_axisb; 
            axis_t = ell_axis;
            Matrix6d C66 = rotate_C66(Btovoigt(C_VP_SC), axisb_t.transpose());
            //Calculate Eshelby tensor
            double S4_EA[3][3][3][3] = {0}; 
            double R4_EA[3][3][3][3] = {0};
            int case_c = Eshelby_case(axis_t);
            g[0].Eshelby_P(S4_EA,R4_EA,axis_t,C66,Gpsets[case_c].Gpaa6,Gpsets[case_c].Gpaaww6,\
                           Gpsets[case_c].Gpalpha,Gpsets[case_c].Gpaww,Gpsets[case_c].Gpww);
            //rotate the eshelby tensor back to the sample axes
            S55 = voigttoB5(rotate_C66(voigt(S4_EA), axisb_t));
            rot_4th(R4_EA, axisb_t, R4_SA);
            //Calculate Mtilde refer to Equ[5-33] in Manual 7d
            //M~=(I-S)^-1 * S * M
            Matrix5d I_S55, I_S55_inv; // (I-S) and (I-S)^-1
            I_S55 = Matrix5d::Identity() - S55;
            I_S55_inv = I_S55.inverse();
            Mtilde = I_S55_inv * S55 * M_VP_SC;
            S55_inv = S55.inverse();
        }

        DVP_AV = Vector5d::Zero();
        for(int G_n = 0; G_n < grains_num; G_n++)
        {
            /* Rotate the macro VP stiffness (C_VP_SC) from sample axes to the ellipsoid axes; */
            /* if Ishape = 1, which means the ellipsoid axes varies in grains, */
            /* the transform matrix should be taken out of each grain   */
            if(Ishape == 1)
            {
                axis_t = g[G_n].get_ell_axis_g();
                axisb_t = g[G_n].get_ell_axisb_g();
                Matrix6d C66 = rotate_C66(Btovoigt(C_VP_SC), axisb_t.transpose());
                //Calculate Eshelby tensor
                double S4_EA[3][3][3][3] = {0}; 
                double R4_EA[3][3][3][3] = {0};
                int case_c = Eshelby_case(axis_t);
                g[G_n].Eshelby_P(S4_EA,R4_EA,axis_t,C66,Gpsets[case_c].Gpaa6,Gpsets[case_c].Gpaaww6,\
                                 Gpsets[case_c].Gpalpha,Gpsets[case_c].Gpaww,Gpsets[case_c].Gpww);
                //rotate the eshelby tensor back to the sample axes
                S55 = voigttoB5(rotate_C66(voigt(S4_EA), axisb_t));
                rot_4th(R4_EA, axisb_t, R4_SA);
                //Calculate Mtilde refer to Equ[5-33] in Manual 7d
                //M~=(I-S)^-1 * S * M
                Matrix5d I_S55, I_S55_inv; // (I-S) and (I-S)^-1
                I_S55 = Matrix5d::Identity() - S55;
                I_S55_inv = I_S55.inverse();
                Mtilde = I_S55_inv * S55 * M_VP_SC;
                S55_inv = S55.inverse();
            }
            g[G_n].Update_Mptilde_g(Mtilde);  //store the M~
            //only affine case (interaction = 1) 
            double S66inv4th[3][3][3][3] = {0};
            Chg_basis(S55.inverse(),S66inv4th);
            mult_4th(R4_SA,S66inv4th,RSinv_SA);
            g[G_n].Update_RSinv_VP_g(RSinv_SA); //store the PI*(S^-1)

            Matrix5d M_g = g[G_n].get_Mpij6_g();
            Vector5d d0_g = g[G_n].get_d0_g(); 
            double weight_g = g[G_n].get_weight_g();

            //Calculate the localization tensor B_g, _g means the value depends on the grain
            // B_g = (M_g + M~)^-1 * (M_ + M~) 
            Matrix5d term_M = (M_g + Mtilde).inverse();
            Matrix5d B_g = term_M * (M_VP_SC + Mtilde);
            Matrix5d B_g_switch = (M_VP_SC + Mtilde) * term_M;
            // b_g = (M_g + M~)^-1 * (d- - d-_g)
            Vector5d b_gv =  term_M * (Chg_basis5(voigt(D0))-d0_g);
            // MNEW = <M_g * B_g> Equ[5-40a]
            MNEW += M_g * B_g * weight_g;
            //<B_g> Equ[5-40b]
            B_g_ave += B_g * weight_g;
            B_g_switch_ave += B_g_switch * weight_g;
            B_switch_d0_ave += B_g_switch * d0_g * weight_g;
            DVP_AV +=  Chg_basis5(g[G_n].get_Dijp_g()) * weight_g;
        } //loop over grains

        // <M_g * B_g> * <B_g>^-1
        MNEW = MNEW * B_g_ave.inverse();
        Matrix5d MNEW2 = 0.5*(MNEW+MNEW.transpose());
        MNEW = MNEW2;
        D0_new = B_g_switch_ave.inverse() * B_switch_d0_ave;
        /* calculate the error between the input M_VP_SC of do-while loop and the output MNEW of the loop  */
        RER = Errorcal(M_VP_SC, MNEW);
        M_VP_SC = MNEW;
        D0 = voigt(Chg_basis(D0_new));
        C_VP_SC = M_VP_SC.inverse();

        logger.notice("**Error in VPSC iteration :" + to_string(IT) + "\t" + to_string(RER));
        if(isnan(RER)) return 1;
    }//while loop
    return 0;
}


int polycrystal::Update_Fij(double Tincr)
{
    Matrix3d Fnew;
    Fnew = Matrix3d::Zero();
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
    {
        Fnew(i,j) += (Tincr*Udot_AV(i,k)+Iij(i,k))*Fij_m(k,j);
    }
    Fij_m = Fnew;
    return 0;
}

int polycrystal::Update_shape()
{
    Matrix3d FT = Fij_m.transpose();
    Matrix3d BX = Fij_m * FT;

    //Solve the eigen vector of BX and sort the value from largest to smallest
    Matrix3d BX_vectors;
    Vector3d BX_value;
    Jacobi(BX,BX_value,BX_vectors);
    Eigsrt(BX_vectors, BX_value);
    
    Matrix3d B = BX_vectors;
    Vector3d W = BX_value;

    /* Redefine Axis(1) (the second) to be the largest to improve the accuracy in calculation of Eshelby tensor */
    /* IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT HANDED BY EXCHANGING 1 AND 2. */
    double Sign = -1, temp = 0.0;
    if(B.determinant() <= 0) Sign = 1;
    for(int i = 0; i < 3; i++)
    {
        temp = B(i,0);
        B(i,0) = B(i,1);
        B(i,1) = temp * Sign;
    }
    temp = W(0); W(0) = W(1); W(1) = temp;

    // Update the stretching of ellipsoid
    double Ratmax=sqrt(W(1)/W(2));
    double Ratmin=sqrt(W(0)/W(2));
    ell_axisb = B;
    if(!Iflat)
        for(int i = 0; i < 3; i++)
            ell_axis(i) = sqrt(W(i));
    // if Iflat = 1, the axis of ellipsoid keeps unchange
    
    // Update the ellipsoid orientation
    Matrix3d BT = B.transpose();
    ellip_ang = Euler_trans(BT);
    
    // Update the Iflat_g according to the Max axes ratio of ellipsoid
    if((!Iflat)&&(Ratmax >= ell_crit_shape)) Iflat = 1;

    // Iflat_g = 1; recaculates the Fij_g in grain
    else if((Iflat)&&(Ratmax >= ell_crit_shape))
    {
        W(1) = W(1)/4;
        if(Ratmin >= 0.5 * ell_crit_shape) W(0) = W(0)/4;
        for(int i = 0; i < 3; i++) ell_axis(i) = sqrt(W(i));
        Matrix3d Fijx(ell_axis.asDiagonal());
        Fij_m = B * Fijx * BT;
    }
    return 0;
}

int polycrystal::EVPSC(int istep, double Tincr){
    double errd, errs, err_g, err_av;
    int max_iter = 30, break_th = 5, break_counter = 0;
    save_status();

    for(int i = 0; i <= max_iter; ++i)
    {
        //save the input for error calculation
        Sig_in = Sig_m;
        Dij_in = Dij_m;
        sig_in_AV = Matrix3d::Zero();

        for(int G_n = 0; G_n < grains_num; ++G_n) 
            sig_in_AV += g[G_n].get_stress_g() * g[G_n].get_weight_g();

        ///////////
        logger.notice("        \tIteration " + std::to_string(i+1) + "\t        ");
        int return_scE = Selfconsistent_E(istep, SC_err_m, SC_iter_m);
        if (return_scE == 1){ error_SC = 2; return 1;}
        int return_scP = Selfconsistent_P(istep, SC_err_m, SC_iter_m);
        if (return_scP == 1){ error_SC = 2; return 1;}
        Cal_Sig_m(Tincr); 
        double sig_frac = Cal_Sig_g(Tincr);
        if (sig_frac < 0.5){ error_SC = 2; return 1;}
        Update_AV();
        ///////////

        errs = Errorcal(Sig_m, Sig_in);
        errd = Errorcal(Dij_m, Dij_in);
        err_g = Errorcal(Sig_AV, sig_in_AV);
        double err_stress = Errorcal(Sig_AV, Sig_m);
        /* logger.debug("Sig_AV: "); logger.debug(Sig_AV); */
        /* logger.debug("Sig_m: "); logger.debug(Sig_m); */
        err_av = std::max({errs, errd, err_g})/(errS_m + errD_m + err_g_AV) * 4;
        if ((err_av / error_SC) > 1.01) break_counter += 1;
        error_SC = err_av;

        logger.notice("Error in macro stress tensor:\t" + std::to_string(errs));
        logger.notice("Error in strain rate tensor:\t" + std::to_string(errd));
        logger.notice("Error in average grain stress:\t" + std::to_string(err_g));
        logger.notice("Error between stress tensors:\t" + std::to_string(err_stress));
        logger.notice("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=");
        bool err_stress_flag = (err_stress <= err_g_AV*1.)||(err_stress*calc_equivalent_value(Sig_m)<0.5);
        if((errs<errS_m)&&(errd<errD_m)&&(err_g<err_g_AV)&&err_stress_flag) break;
        if((errs<0.01* errS_m)&&(errd<0.01*errD_m)&&(err_g<0.01* err_g_AV)) break;
        if ((i == max_iter) || (break_counter == break_th) )return 1;
    }

    Eps_m += Dij_m * Tincr; //update the macro strain tensor

<<<<<<< HEAD
    //update the shape of ellipsoid
    if(Ishape == 0)
    {
        Update_Fij(Tincr);
        if(Iupdate_shp) Update_shape();
    }
    // Update temperature in grains here
    temperature_poly = 0.0;
    for(int G_n = 0; G_n < grains_num; ++G_n){
        g[G_n].update_temperature(Tincr);
        temperature_poly += g[G_n].temperature * g[G_n].get_weight_g();
    }
    /* logger.debug("Temperature in atmosphere: " + std::to_string(temp_atmosphere) + " K"); */
    /* logger.debug("Temperature in polycrystal: " + std::to_string(temperature_poly) + " K"); */
    // double Delta_T = (Tincr*Surface/(V_sample*rho_material*Cp_material))*(h_ext*(temp_atmosphere-temperature_poly)+sigma_k*(pow(temp_atmosphere,4)-pow(temperature_poly)));
    // temperature_poly += Delta_T; 
    // // we have to redistribute the temperature?
    // for(int G_n = 0; G_n < grains_num; ++G_n){
    //     g[G_n].temperature = temperature_poly;
    //    }
//test mode 1.14

    // this is to calculate the temperature after the thermal conduction.

    //update the state in deformation systems and 
    // crystalline orientation 
     for(int G_n = 0; G_n < grains_num; ++G_n)
    {
        /* logger.debug("Stress in grain " + std::to_string(G_n) + ":"); */
        /* logger.debug(g[G_n].get_stress_g()); */
        g[G_n].Update_shear_strain(Tincr);//seems this function is no more needed
        g[G_n].update_strain(Tincr);
        if(Iupdate_ori) g[G_n].update_orientation(Tincr, Wij_m, Dije_AV, Dijp_AV);
        if(Iupdate_CRSS) g[G_n].update_modes(Tincr);
        if(Ishape == 1)
        {
            g[G_n].Update_Fij_g(Tincr);
            if(Iupdate_shp) g[G_n].Update_shape_g();
        }
        if(Iupdate_CRSS) update_twin_control();
    }
=======
    if(update_temperature_required) update_temperature(Tincr);
    update_status(Tincr);
    update_info_by_family();
>>>>>>> 6f8c8fa27d07fdc544418580ee2acaef7ff1449d
    return 0;
}

void polycrystal::save_status(){
    // Save the current status after a successful iteration
    temp_poly_old = temperature_poly;
    Wij_m_old = Wij_m;
    Dij_m_old = Dij_m;
    Dije_AV_old = Dije_AV;
    Dijp_AV_old = Dijp_AV;
    Sig_m_old = Sig_m;
    COLD = CSC;
    therm_expans_ave_old = therm_expansion_ave;
    ther_strain_m_old = thermal_strain_m;
    C_VP_SC_old = C_VP_SC;
    for(int G_n = 0; G_n < grains_num; ++G_n)
        g[G_n].save_status_g();
}

void polycrystal::restore_status(bool flag){
    // Restore the status before the iteration
    temperature_poly = temp_poly_old;
    Wij_m = Wij_m_old;
    Dij_m = Dij_m_old;
    Dije_AV = Dije_AV_old;
    Dijp_AV = Dijp_AV_old;
    Sig_m = Sig_m_old;
    CSC = COLD;  C_VP_SC = C_VP_SC_old;
    therm_expansion_ave = therm_expans_ave_old;
    thermal_strain_m = ther_strain_m_old;
    SSC = CSC.inverse();  M_VP_SC = C_VP_SC.inverse();
    for(int G_n = 0; G_n < grains_num; ++G_n)
        g[G_n].restore_status_g();
}

// update the state in deformation systems and crystalline orientation and shapes
void polycrystal::update_status(double time_incre){
    if(Ishape == 0)
    {
        Update_Fij(time_incre);
        if(update_shape_required) Update_shape();
    }
    #pragma omp parallel for num_threads(Mtr)
    for(int G_n = 0; G_n < grains_num; ++G_n){
        g[G_n].update_strain(time_incre);
        if(update_orientation_required) g[G_n].update_orientation(time_incre, Wij_m, Dije_AV, Dijp_AV);
        if(update_CRSS_required) g[G_n].update_modes(time_incre);
        if(update_CRSS_required) update_twin_control();
        if(Ishape == 0) continue;
        g[G_n].Update_Fij_g(time_incre);
        if(update_shape_required) g[G_n].Update_shape_g();
    }
    #pragma omp barrier
    update_info_by_family();
}

void polycrystal::Cal_Sig_m(double Tincr)
{
    double temperature_diff = temperature_poly - temperature_ref;
    Wij_m = Matrix3d::Zero();
    //why not Wij = Udot - Dij ?
    //because the Udot need be update
    //some components are not imposed
    //IUdot: 0: not imposed; 1: imposed; 2: control Dij; 3: control Wij
    for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++){
        if(IUdot(i,j) > 1 && IUdot(j,i) != 0){
            logger.error("IUdot configuration not support!");
            exit(1);
        }
        if(IUdot(i,j) == 1 && IUdot(j,i) == 1){
            Wij_m(i,j) = 0.5*(UDWdot_m(i,j) - UDWdot_m(j,i));
            Wij_m(j,i) = -Wij_m(i,j);
            }
        else if(IUdot(i,j) == 1){
            Wij_m(i,j) = UDWdot_m(i,j) - Dij_m(i,j);
            Wij_m(j,i) = -Wij_m(i,j);
            }
        else if(IUdot(i,j) == 3){
            Wij_m(i,j) = UDWdot_m(i,j);
            Wij_m(j,i) = -Wij_m(i,j);
        }
        Udot_m(i,j) = Dij_m(i,j) + Wij_m(i,j);
    }

    Matrix3d Sig_J = Wij_m * Sig_m - Sig_m * Wij_m; //Jaumann rate part
    Matrix3d therm_strain_rate = (therm_expansion_ave * temperature_diff - ther_strain_m_old)/Tincr;
    Vector6d De = Bbasisadd(Chg_basis6(Dij_m), -DVP_AV) - Chg_basis6(therm_strain_rate); //De = D - Dp - D0 - Dt
    //Calculate AX = B, X is the macro stress
    Matrix6d A = SSC/Tincr;
    Matrix3d Mtemp = Sig_m_old / Tincr;
    Vector6d B = De + SSC * ( Chg_basis6(Mtemp) + Chg_basis6(Sig_J));
    /* According to the IUdot and ISDOT to solve AX = B, transform to solve At Xt = Bt */
    /* Xt in the unkown set of Udot and Sig */
    Vector6i Is = ISdot;
    Vector6i Id = IDdot;
    Vector6d profac;
    profac << 1,1,1,2,2,2;

    Vector6d BC_D = voigt(Chg_basis(B));
    Vector6d BC_S = voigt(Sig_m);
    Matrix6d AUX2 = Btovoigt(A);

    Vector6d AUX11;
    Matrix6d AUX21;

    for(int i = 0; i != 6; i++){
        AUX11(i) = -Id(i)*BC_D(i);
        for(int j = 0; j != 6; j++){
            AUX11(i) = AUX11(i) + AUX2(i,j)*Is(j)*BC_S(j)*profac(j);
            AUX21(i,j) = Is(j)*int(i==j) - Id(j)*AUX2(i,j)*profac(j);
        }
    }
    Vector6d AUX6 = AUX21.inverse()*AUX11;

    Vector6d Sig_x, B_x;
    B_x = mult_dot(BC_D,Id) + mult_dot(AUX6,Is);
    Sig_x = mult_dot(BC_S,Is) + mult_dot(AUX6,Id);
    De = Chg_basis6(voigt(B_x)) - SSC * ( Chg_basis6(Mtemp) + Chg_basis6(Sig_J));
    Vector6d Dij_m_v = Bbasisadd(De, DVP_AV) + Chg_basis6(therm_strain_rate); //D = De + Dp + D0 + Dt
    thermal_strain_m = therm_expansion_ave * temperature_diff;
    Dij_m = Chg_basis(Dij_m_v);
    Sig_m = voigt(Sig_x);
}

double polycrystal::Cal_Sig_g(double Tincr)
{
    #pragma omp parallel for num_threads(Mtr)
    for(int G_n = 0; G_n < grains_num; G_n++){
        g[G_n].grain_stress(Tincr, Wij_m, Dij_m, Dije_AV, Dijp_AV, Sig_m, Sig_m_old, therm_expansion_ave);
    }
    #pragma omp barrier
    int stress_count = 0;
    for (int G_n = 0; G_n < grains_num; G_n++){
        stress_count += g[G_n].if_stress;
    }
    return double(stress_count)/grains_num;
}


void polycrystal::Update_AV()
{
    Dije_AV = Matrix3d::Zero();
    Dijp_AV = Matrix3d::Zero();
    Sig_AV = Matrix3d::Zero();
    Udot_AV = Matrix3d::Zero();
    double wei;
    for(int G_n = 0; G_n < grains_num; G_n++)
    {
        wei = g[G_n].get_weight_g();
        Dije_AV += g[G_n].get_Dije_g() * wei;
        Dijp_AV += g[G_n].get_Dijp_g() * wei;
        Sig_AV += g[G_n].get_stress_g() * wei;
        Udot_AV += g[G_n].get_Udot_g() * wei;
        Dij_AV += g[G_n].get_Dij_g() * wei;
    }
}

Vector6d polycrystal::get_Sig_m(){return voigt(Sig_m);}
Vector6d polycrystal::get_Sig_ave(){return voigt(Sig_AV);}
Vector6d polycrystal::get_Eps_m(){return voigt(Eps_m);}

void polycrystal::get_euler(fstream &texfile)
{
    IOFormat Outformat(StreamPrecision);
    for(int i = 0; i < grains_num; i++)
    {
        texfile << setprecision(4) << scientific << g[i].get_euler_g().transpose().format(Outformat) << "  ";
        texfile << setprecision(4) << scientific << g[i].get_weight_g_eff() << endl;
        for (int j = 0; j < g[i].modes_num; j++)
        {
            if (g[i].gmode[j]->type != mode_type::twin) continue;
            if (g[i].gmode[j]->disloc_density < 0.001) continue;//disloc_density is child_frac for twinG mode;
            texfile << setprecision(4) << scientific << g[i].get_euler_g(j).transpose().format(Outformat) << "  ";
            texfile << setprecision(4) << scientific << g[i].get_weight_g(j) << endl;
        }
    }
}

Vector3d polycrystal::get_ell_axis(){return ell_axis;}
Vector3d polycrystal::get_ellip_ang(){return ellip_ang;}

void polycrystal::update_twin_control(){
    double V_eff = 0., V_acc = 0., A_1 = -1., A_2 = -1.;
    for (int G_n = 0; G_n < grains_num; G_n++){
        V_eff += g[G_n].child_frac * g[G_n].get_weight_g();
        V_acc += (g[G_n].twin_term_flag == true) ? g[G_n].get_weight_g() : 0;
        if (A_1 > -0.01 && A_2 > -0.01) continue;
        for (int imode = 0; imode < g[G_n].modes_num; imode++){
            if (g[G_n].gmode[imode]->type != mode_type::twin) continue;
            if (A_1 < 0.) A_1 = g[G_n].gmode[imode]->harden_params[5];
            if (A_2 < 0.) A_2 = g[G_n].gmode[imode]->harden_params[6];
            break;
        }
    }
    twin_threshold = min(1.0, A_1 + A_2 * (V_eff / V_acc));
}

void polycrystal::update_info_by_family(){
    PMode* modePtr = nullptr;
    int modes_from = 0, modes_to = 0;
    for (int family_id = 0; family_id < family_count; family_id++){
        if (family_id != 0) modes_from += modes_count_by_family[family_id - 1];
        modes_to += modes_count_by_family[family_id];
        double density_in_family = 0., acc_strain_in_family = 0., crss_in_family = 0.;
        for (int grain_i = 0; grain_i < grains_num; grain_i++){
            double weight_g = g[grain_i].get_weight_g();
            for (int mode_i = modes_from; mode_i < modes_to; mode_i++){
                modePtr = g[grain_i].gmode[mode_i];
                density_in_family += modePtr->disloc_density * weight_g;
                acc_strain_in_family += modePtr->acc_strain * weight_g;
                crss_in_family += (modePtr->crss/((double)(modes_to-modes_from))) * weight_g;
            }
        }
        density_by_family[family_id] = density_in_family;
        acc_strain_by_family[family_id] = acc_strain_in_family;
        crss_by_family[family_id] = crss_in_family;
    }
}

void polycrystal::set_temperature(double temperature_in){
    temperature_poly = temperature_in;
    for (int i = 0; i < grains_num; ++i) {
        g[i].temperature = temperature_in;
    }
}

// Calculate the temperature after the thermal conduction.
void polycrystal::update_temperature(double time_incre){
    temperature_poly = 0.0;
    for(int G_n = 0; G_n < grains_num; ++G_n){
        g[G_n].update_temperature(time_incre);
        temperature_poly += g[G_n].temperature * g[G_n].get_weight_g();
    }
    double scale = time_incre * Surface / (V_sample * rho_material * Cp_material);
    double convection_term = h_ext * (temp_atmosphere - temperature_poly);
    double radiation_term = sigma_k * (pow(temp_atmosphere, 4) - pow(temperature_poly, 4));
    double Delta_T = (convection_term + radiation_term) * scale;
    temperature_poly += Delta_T; 
    // we have to redistribute the temperature?
    set_temperature(temperature_poly);
}
