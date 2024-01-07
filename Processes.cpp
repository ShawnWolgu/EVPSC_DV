#include "Processes.h"
#include "global.h"
using namespace Procs;

Process::Process()
{
    ss_out.open("str_str.out",ios::out); 
    ss_out << setw(10) << "E11" << setw(11)<< "E22" << setw(11)<< "E33";
    ss_out << setw(11) << "E32" << setw(11)<< "E13" << setw(11)<< "E12";
    ss_out << setw(11) << "S11" << setw(11)<< "S22" << setw(11)<< "S33";
    ss_out << setw(11) << "S32" << setw(11)<< "S13" << setw(11)<< "S12" << endl;
    tex_out.open("Tex.out",ios::out);
    disloc_out.open("Disloc.csv",ios::out);
    /* time_out.open("RunTimeFrac.csv",ios::out); */
    rate_out.open("Rate.csv",ios::out);
    crss_out.open("CRSS.csv",ios::out);
    rss_out.open("RSS.csv",ios::out);
    euler_out.open("Euler.csv",ios::out);
    ss_out_csv.open("str_str.csv",ios::out);
    ss_grain_out.open("str_str_grain.csv",ios::out);
    ave_ss_out.open("ave_str_str.csv",ios::out);
    de_out.open("de.csv",ios::out);
    Mvp_out.open("Mvp.csv",ios::out);
}

Process::~Process()
{
    ss_out.close();
    tex_out.close();
    disloc_out.close();
    /* time_out.close(); */
    rate_out.close();
    crss_out.close();
    euler_out.close();
    ss_out_csv.close();
    ss_grain_out.close();
    ave_ss_out.close();
    de_out.close();
    Mvp_out.close();
}

void Process::load_ctrl(Vector4d Vin)
{
    Nsteps = int(Vin(0));
    Ictrl = int(Vin(1)) - 1;
    Eincr = Vin(2);
    Temp = Vin(3);

    if(!texctrl) texctrl = Nsteps;
}

void Process::get_Udot(Matrix3d Min)
{   
    UDWdot_input = Min;
    Ddot_input = 0.5*(Min + Min.transpose());
    
    //calculate Time increment Tincr
    Vector6d Vtemp = voigt(Ddot_input);
    if(Vtemp(Ictrl) == 0 || Eincr == 0){
        Tincr = 1;
    }
    else{
        Tincr = Eincr / Vtemp(Ictrl);
    }
}
void Process::get_Sdot(Matrix3d Min){Sdot_input = Min;}
void Process::get_IUdot(Matrix3i Min){IUDWdot = Min;}
void Process::get_ISdot(Vector6i Vin){ISdot = Vin;}

void Process::loading(Polycs::polycrystal &pcrys){
    temp_atmosphere = Temp;
    if (pcrys.temperature_poly < 1e-3){
        // Temperature is not set, this is the first loading step, free of thermal stress
        pcrys.temperature_poly = Temp;
        temperature_ref = Temp;
        for (int i = 0; i < pcrys.grains_num; ++i) {
            pcrys.g[i].temperature = Temp;
        }
    }
    pcrys.set_BC_const(UDWdot_input, Sdot_input, IUDWdot, ISdot);
    Out_sscurves(pcrys);
    double coeff_step = 1, current_step = 1.;
    for(int istep = 0; istep < Nsteps; ++istep)
    {
	logger.info("");
	logger.info("**********\tSTEP\t"+ to_string(istep) + "\t**********");
	logger.info("");
	double pct_step = 0; 
	update_progress(pct_step);
	do{
        current_step = min(1.0 - pct_step, coeff_step);
	    logger.notice("Step " + to_string(istep) + ":\t" + to_string(pct_step) + " to " + to_string(pct_step + current_step));
        int return_SC = pcrys.EVPSC(istep, current_step * Tincr, Iupdate_ori, Iupdate_shp, Iupdate_CRSS);
	    if (return_SC != 0) {
            pcrys.restore_status(false);
            logger.warn("Not convergent... Retry with a smaller increment.");
            logger.notice("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=");
            if(isnan(pcrys.error_SC)) coeff_step *= 0.1;
            else coeff_step *= 0.5;
	    	if (coeff_step < 1e-10) { logger.error("Not convergent... Abort."); exit(1);}
		    continue;
	    }
	    pct_step += current_step;
	    update_progress(pct_step);
        coeff_step *= 2;
        coeff_step = min(coeff_step, 1.0);
	} while (pct_step < 1-1e-10);
	cout.flush();
        Out_sscurves(pcrys);
        if(!((istep+1)%texctrl))
            Out_texture(pcrys,istep);
	if(istep == 0) init_grain_info(pcrys, 0);
	Out_grain_info(pcrys, 0);
    }
    Out_texture(pcrys,Nsteps);
}

void Process::Out_sscurves(Polycs::polycrystal &pcrys)
{
    IOFormat Outformat(StreamPrecision);
    ss_out << setprecision(4) << scientific << pcrys.get_Eps_m().transpose().format(Outformat)<< " ";
    ss_out << setprecision(4) << scientific << pcrys.get_Sig_m().transpose().format(Outformat) << endl;
}

void Process::Out_texture(Polycs::polycrystal &pcrys, int istep)
{
    IOFormat Outformat(StreamPrecision);
    logger.notice("Output texture at step " + to_string(istep+1));
    tex_out << "TEXTURE AT STEP = " << istep+1 << endl;
    tex_out << setprecision(4) << pcrys.get_ell_axis().transpose().format(Outformat)<< endl; 
    tex_out << setprecision(4) << pcrys.get_ellip_ang().transpose().format(Outformat) << endl << endl;
    pcrys.get_euler(tex_out);
    tex_out << endl;
}

void Process::init_grain_info(Polycs::polycrystal &pcrys, int num){
    IOFormat Outformat(StreamPrecision);
    logger.notice("Output grain info at step 0");
    //output dislocation density
    grain *g_this = &pcrys.g[num];
    int mode_num = g_this->modes_num;
    
    disloc_out << "EVM";
    for(int i = 0; i < mode_num; ++i) disloc_out << "," << "Mode" << i+1;
    disloc_out << endl;
    disloc_out << calc_equivalent_value(g_this->get_strain_g()) << ", ";
    /* disloc_out << g_this->gmode[0].rho_sat; */
    for(int i = 0; i < mode_num; ++i) disloc_out << "," << g_this->gmode[i]->disloc_density;
    disloc_out << endl;
    
    /* time_out << "EVM"; */
    /* for(int i = 0; i < mode_num; ++i) time_out << "," << "Mode" << i+1; */
    /* time_out << endl; */
    /* time_out << calc_equivalent_value(g_this->get_strain_g()); */
    /* for(int i = 0; i < mode_num; ++i) time_out << "," << (g_this->gmode[i].t_run)/(g_this->gmode[i].t_wait + g_this->gmode[i].t_run); */
    /* time_out << endl; */
    
    rate_out << "EVM";
    for(int i = 0; i < mode_num; ++i) rate_out << "," << "Mode" << i+1;
    rate_out << endl;
    rate_out << calc_equivalent_value(g_this->get_strain_g());
    for(int i = 0; i < mode_num; ++i) rate_out << "," << abs(g_this->gmode[i]->shear_rate);
    rate_out << endl;

    crss_out << "EVM";
    for(int i = 0; i < mode_num; ++i) crss_out << "," << "Mode" << i+1;
    crss_out << endl;
    crss_out << calc_equivalent_value(g_this->get_strain_g());
    for(int i = 0; i < mode_num; ++i) crss_out << "," << g_this->gmode[i]->crss;
    crss_out << endl;

    rss_out << "EVM";
    for(int i = 0; i < mode_num; ++i) rss_out << "," << "Mode" << i+1;
    rss_out << endl;
    rss_out << calc_equivalent_value(g_this->get_strain_g());
    for(int i = 0; i < mode_num; ++i) rss_out << "," << g_this->gmode[i]->rss;
    rss_out << endl;
    
    euler_out << "EVM,phi1,PHI,phi2\n";
    euler_out << calc_equivalent_value(g_this->get_strain_g());
    Vector3d v_out = g_this->get_euler_g();
    for(int i = 0; i < 3; ++i) euler_out << "," << v_out(i);
    euler_out << endl;

    ss_out_csv << "E11,E22,E33,E23,E13,E12,S11,S22,S33,S23,S13,S12,TempK\n";
    for(int i = 0; i < 6; ++i) ss_out_csv << pcrys.get_Eps_m()(i) << ",";
    for(int i = 0; i < 5; ++i) ss_out_csv << pcrys.get_Sig_m()(i) << ",";
    ss_out_csv << pcrys.get_Sig_m()(5) << "," << pcrys.temperature_poly << endl;

    ave_ss_out << "E11,E22,E33,E23,E13,E12,S11,S22,S33,S23,S13,S12\n";
    for(int i = 0; i < 6; ++i) ave_ss_out << pcrys.get_Eps_m()(i) << ",";
    for(int i = 0; i < 5; ++i) ave_ss_out << pcrys.get_Sig_ave()(i) << ",";
    ave_ss_out << pcrys.get_Sig_ave()(5) << endl;

    ss_grain_out << "EVM,SVM,E11,E22,E33,E23,E13,E12,S11,S22,S33,S23,S13,S12\n";
    ss_grain_out << calc_equivalent_value(g_this->get_strain_g());
    ss_grain_out << "," << calc_equivalent_value(g_this->get_stress_g());
    for(int i = 0; i < 6; ++i) ss_grain_out << "," << voigt(g_this->get_strain_g())(i);
    for(int i = 0; i < 6; ++i) ss_grain_out << "," << voigt(g_this->get_stress_g())(i);
    ss_grain_out << endl;

    de_out << "EVM,de11,de22,de33,de23,de13,de12\n";
    de_out << calc_equivalent_value(voigt(pcrys.get_Eps_m()));
    Vector6d De = Bbasisadd(Chg_basis6(pcrys.Dij_m), -pcrys.DVP_AV);
    for(int i = 0; i < 6; ++i) de_out << "," << De(i);
    de_out << endl;

    Mvp_out << "EVM,Mvp_norm\n";
    Mvp_out << calc_equivalent_value(voigt(pcrys.get_Eps_m()));
    Mvp_out << "," << (pcrys.M_VP_SC).norm() << endl;
}

void Process::Out_grain_info(Polycs::polycrystal &pcrys, int num){
    IOFormat Outformat(StreamPrecision);
    //output dislocation density
    grain *g_this = &pcrys.g[num];
    int mode_num = g_this->modes_num;

    disloc_out << calc_equivalent_value(g_this->get_strain_g()) << ", ";
    /* disloc_out << g_this->gmode[0].rho_sat; */
    for(int i = 0; i < mode_num; ++i) disloc_out << "," << g_this->gmode[i]->disloc_density;
    disloc_out << endl;
    
    /* time_out << calc_equivalent_value(g_this->get_strain_g()); */
    /* for(int i = 0; i < mode_num; ++i) time_out << "," << (g_this->gmode[i].t_run)/(g_this->gmode[i].t_wait + g_this->gmode[i].t_run); */
    /* time_out << endl; */
    
    rate_out << calc_equivalent_value(g_this->get_strain_g());
    for(int i = 0; i < mode_num; ++i) rate_out << "," << abs(g_this->gmode[i]->shear_rate);
    rate_out << endl;
    
    crss_out << calc_equivalent_value(g_this->get_strain_g());
    for(int i = 0; i < mode_num; ++i) crss_out << "," << g_this->gmode[i]->crss;
    crss_out << endl;
    
    rss_out << calc_equivalent_value(g_this->get_strain_g());
    for(int i = 0; i < mode_num; ++i) rss_out << "," << g_this->gmode[i]->rss;
    rss_out << endl;
    
    euler_out << calc_equivalent_value(g_this->get_strain_g());
    Vector3d v_out = g_this->get_euler_g();
    for(int i = 0; i < 3; ++i) euler_out << "," << v_out(i);
    euler_out << endl;

    for(int i = 0; i < 6; ++i) ss_out_csv << pcrys.get_Eps_m()(i) << ",";
    for(int i = 0; i < 5; ++i) ss_out_csv << pcrys.get_Sig_m()(i) << ",";
    ss_out_csv << pcrys.get_Sig_m()(5) << "," << pcrys.temperature_poly << endl;
    
    for(int i = 0; i < 6; ++i) ave_ss_out << pcrys.get_Eps_m()(i) << ",";
    for(int i = 0; i < 5; ++i) ave_ss_out << pcrys.get_Sig_ave()(i) << ",";
    ave_ss_out << pcrys.get_Sig_ave()(5) << endl;

    ss_grain_out << calc_equivalent_value(g_this->get_strain_g());
    ss_grain_out << "," << calc_equivalent_value(g_this->get_stress_g());
    for(int i = 0; i < 6; ++i) ss_grain_out << "," << voigt(g_this->get_strain_g())(i);
    for(int i = 0; i < 6; ++i) ss_grain_out << "," << voigt(g_this->get_stress_g())(i);
    ss_grain_out << endl;

    de_out << calc_equivalent_value(voigt(pcrys.get_Eps_m()));
    Vector6d De = Bbasisadd(Chg_basis6(pcrys.Dij_m), -pcrys.DVP_AV);
    for(int i = 0; i < 6; ++i) de_out << "," << De(i);
    de_out << endl;

    Mvp_out << calc_equivalent_value(voigt(pcrys.get_Eps_m()));
    Mvp_out << "," << (pcrys.M_VP_SC).norm() << endl;
}

void Process::Out_texset(int input){texctrl = input;}

void Process::Update_ctrl(Vector3i input){
    Iupdate_ori = input(0);
    Iupdate_shp = input(1);
    Iupdate_CRSS = input(2);
    }
