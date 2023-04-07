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
    grain_out.open("Grain.csv",ios::out);
}

Process::~Process()
{
    ss_out.close();
    tex_out.close();
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
    Udot_input = Min;
    Ddot_input = 0.5*(Min + Min.transpose());
    
    //calculate Time increment Tincr
    Vector6d Vtemp = voigt(Ddot_input);
    Tincr = Eincr / Vtemp(Ictrl);
}
void Process::get_Sdot(Matrix3d Min){Sdot_input = Min;}
void Process::get_IUdot(Matrix3i Min){IUdot = Min;}
void Process::get_ISdot(Vector6i Vin){ISdot = Vin;}

void Process::loading(Polycs::polycrystal &pcrys)
{
    pcrys.ini_Udot_m(Udot_input);
    pcrys.ini_Sig_m(Sdot_input);
    pcrys.set_IUdot(IUdot);
    pcrys.set_ISdot(ISdot);
    double coeff_step = 1;
    for(int istep = 0; istep < Nsteps; ++istep)
    {
	logger.info("");
	logger.info("**********\tSTEP\t"+ to_string(istep) + "\t**********");
	logger.info("");
	double pct_step = 0; 
	update_progress(pct_step);
	do{
	    coeff_step = min(1.0 - pct_step, coeff_step);
	    logger.notice("Step " + to_string(istep) + ":\t" + to_string(pct_step) + " to " + to_string(pct_step + coeff_step));
            int return_SC = pcrys.EVPSC(istep, coeff_step * Tincr, Iupdate_ori, Iupdate_shp, Iupdate_CRSS);
	    if (return_SC == 1) {
		pcrys.restore_status();
		logger.warn("Not convergent... Retry with a smaller increment.");
		logger.notice("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=");
		if(isnan(pcrys.error_SC)) coeff_step *= 0.1;
		else coeff_step *= sqrt(sqrt(1/pcrys.error_SC));
	    	if (coeff_step < 1e-6) { logger.error("Not convergent... Abort."); exit(1);}
		continue;
	    }
	    pct_step += coeff_step;
	    update_progress(pct_step);
	    coeff_step *= sqrt(sqrt(1/pcrys.error_SC));
	} while (pct_step < 1-1e-6);
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
    grain_out << "EVM, ";
    for(int i = 0; i < mode_num; ++i) grain_out << "Mode" << i+1 << ", ";
    grain_out << endl;
}

void Process::Out_grain_info(Polycs::polycrystal &pcrys, int num){
    IOFormat Outformat(StreamPrecision);
    //output dislocation density
    grain *g_this = &pcrys.g[num];
    int mode_num = g_this->modes_num;
    grain_out << calc_equivalent_value(g_this->get_strain_g()) << ", ";
    for(int i = 0; i < mode_num; ++i) grain_out << g_this->gmode[i].disloc_density << ", ";
    grain_out << endl;
}

void Process::Out_texset(int input){texctrl = input;}

void Process::Update_ctrl(Vector3i input){
    Iupdate_ori = input(0);
    Iupdate_shp = input(1);
    Iupdate_CRSS = input(2);
    }
