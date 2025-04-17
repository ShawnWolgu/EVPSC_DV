#include "common/common.h"
#include "solver/Grains.h"
#include "mechanism/PMode.h"
#include "io/Output.h"

using namespace Procs;

Process::Process(){};

Process::~Process()
{
    /* Out_texture(global_polycrys, istep); */
    tex_out.close();
    density_out.close();
    ss_out_csv.close();
    ave_ss_out.close();
    custom_out.close();
    crss_out.close();
}

void Process::load_ctrl(Vector4d Vin)
{
    Nsteps = int(Vin(0));
    Ictrl = int(Vin(1));
    Eincr = Vin(2);
    temperature_input = Vin(3);

    if(!texctrl) texctrl = Nsteps;
}

void Process::get_Udot(Matrix3d Min)
{   
    UDWdot_input = Min;
    Ddot_input = 0.5*(Min + Min.transpose());
    
}
void Process::get_Sdot(Matrix3d Min){Sdot_input = Min;}
void Process::get_IUdot(Matrix3i Min){IUDWdot = Min;}
void Process::get_ISdot(Vector6i Vin){ISdot = Vin;}
void Process::timestep_control(){
    //calculate Time increment Tincr
    //Ictrl <= 6: strain rate control
    const int IJV[6][2]= {0,0,1,1,2,2,1,2,0,2,0,1};
    if(Ictrl <= 6){
        int I1 = IJV[Ictrl-1][0];
        int I2 = IJV[Ictrl-1][1];
        if(IUDWdot(I1,I2) == 0 || IUDWdot(I2,I1) == 0){
            logger.error("Error. The velocity gradient is not set.");
            logger.error("Please check the process file.");
            throw std::runtime_error("Error. The velocity gradient is not set.");
            return;
        }
        Vector6d Vtemp = voigt(Ddot_input);
        double d_control = Vtemp(Ictrl-1);
        if (d_control == 0){
            logger.error("Error. The control strain rate is zero.");
            logger.error("Please check the process file.");
            throw std::runtime_error("Error. The control strain rate is zero.");
            return;
        }
        max_timestep = abs(Eincr) / abs(d_control);
        logger.info("max_timestep = " + to_string(max_timestep));
        return;
    }
    //Ictrl == 7: stress rate control
    else if(Ictrl == 7){
        max_timestep = abs(Eincr);
        return;
    }
    //Else: throw error
    else{
        logger.error("Error. The control type is not set.");
        logger.error("Please check the process file.");
        throw std::runtime_error("Error. The control type is not set.");
        return;
    }
}

void Process::loading(Polycs::polycrystal &pcrys){
    timestep_control();
    time_acc = 0;
    temp_atmosphere = temperature_input; //loading temp_init from the txt.in 
    // Temperature is not set, this is the first loading step, free of thermal stress
    if (temperature_ref < 1e-3)  temperature_ref = temperature_input;
    if (pcrys.temperature_poly < 1e-3)  pcrys.set_temperature(temperature_input);
    pcrys.set_boundary_conditions(UDWdot_input, Sdot_input, IUDWdot, ISdot);
    initial_output_files();
    double coeff_step = 1, current_step = 1.;
    int success_count = 0;
    for(istep = 0; istep < Nsteps; ++istep)
    {
        logger.info("");
        logger.info("**********\tSTEP\t"+ to_string(istep) + "\t**********");
        logger.info("");
        bool is_convergent = true;
        double pct_step = 0; 
        update_progress(pct_step);
        do{
            current_step = min(1.0 - pct_step, coeff_step);
            logger.notice("Step " + to_string(istep) + ":\t" + to_string(pct_step) + " to " + to_string(pct_step + current_step));
            int return_SC = pcrys.EVPSC(istep, current_step * max_timestep);
            if (return_SC != 0) {
                success_count = 0;
                pcrys.restore_status(false);
                logger.warn("Not convergent... Retry with a smaller increment.");
                logger.notice("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=");
                if(isnan(pcrys.error_SC)) coeff_step *= 0.1;
                else coeff_step *= 0.66;
                if (coeff_step < 1e-10) {
                    logger.error("Not convergent... Abort.");
                    is_convergent = false;
                    break;
                }
                continue;
            }
            time_acc += current_step * max_timestep;
            if (flag_emode == 0){
                Current_intensity = 0.0;
            }else if (flag_emode == 1){
                Current_intensity = J_intensity_pulse(time_acc, duty_ratio_J, Amplitude_J, Frequency);
            }else if(flag_emode == 2){
                Current_intensity = J_shock_sim(time_acc, deformation_max, deformation_rate, Amplitude_J, shock_int, shock_fin); 
            }else{
                logger.warn("Error. Please check the emode.");
            } // 辨别电流模式
            // Current_intensity = J_shock_sim(time_acc, deformation_max, deformation_rate, Amplitude_J, shock_int, shock_fin);
            custom_vars[5] = Current_intensity; //输出电流to csv
            pct_step += current_step;
            update_progress(pct_step);
            success_count++;
            if (success_count > 3) coeff_step *= 1.5;
            else if (success_count > 6) coeff_step *= 2;
            coeff_step = min(coeff_step, 1.0);
        } while (pct_step < 1-1e-10);
        cout.flush();
        if(!((istep+1)%texctrl)) Out_texture(pcrys,istep);
        output_info();
        output_grain_info(0);
        if(!is_convergent) {
            Out_texture(pcrys, istep);
            break;
        }
    }
    /* Out_texture(pcrys,Nsteps); */
}

void Process::Out_texture(Polycs::polycrystal &pcrys, int istep)
{
    IOFormat Outformat(StreamPrecision);
    logger.notice("Output texture at step " + to_string(istep));
    tex_out << "TEXTURE AT STEP = " << istep+1 << endl;
    tex_out << setprecision(4) << pcrys.get_ell_axis().transpose().format(Outformat)<< endl; 
    tex_out << setprecision(4) << pcrys.get_ellip_ang().transpose().format(Outformat) << endl << endl;
    pcrys.get_euler(tex_out);
    tex_out << endl;
}

void Process::Out_texset(int input){texctrl = input;}

