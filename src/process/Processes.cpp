#include "common/common.h"
#include "solver/Grains.h"
#include "mechanism/PMode.h"
#include "io/Output.h"

using namespace Procs;

Process::Process(){};

Process::~Process(){};

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
    logger.info("Input velocity gradient tensor: ");
    logger.info(UDWdot_input);
    logger.info("Input strain rate tensor: ");
    logger.info(Ddot_input);
}
void Process::get_Sdot(Matrix3d Min){
    Sdot_input = Min;
    logger.info("Input stress rate tensor: ");
    logger.info(Sdot_input);
}
void Process::get_Efdot(Vector3d Min){
    Efdot_input = Min;
    logger.info("Input electric field vector: ");
    logger.info(Efdot_input.transpose());
}
void Process::get_Eddot(Vector3d Min){
    Eddot_input = Min;
    logger.info("Input electric displacement vector: ");
    logger.info(Eddot_input.transpose());
}

void Process::get_IUdot(Matrix3i Min){
    IUDWdot = Min;
    logger.info("Input velocity gradient flag: ");
    logger.info(IUDWdot);
}
void Process::get_ISdot(Vector6i Vin){
    ISdot = Vin;
    logger.info("Input stress rate flag: ");
    logger.info(ISdot.transpose());
}
void Process::get_IEfdot(Vector3i Vin){
    IEfdot = Vin;
    logger.info("Input electric field flag: ");
    logger.info(IEfdot.transpose());
}
void Process::get_IEddot(Vector3i Vin){
    IEddot = Vin;
    logger.info("Input electric displacement flag: ");
    logger.info(IEddot.transpose());
}

void Process::timestep_control(){
    //calculate Time increment Tincr
    //Ictrl <= 6: strain rate control
    const int IJV[6][2]= {0,0,1,1,2,2,1,2,0,2,0,1};
    if(Ictrl <= 6 && Ictrl != 0){
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
    //Ictrl == 8: electric field & displacement control
    else if(Ictrl == 7 || Ictrl == 8){
        max_timestep = abs(Eincr);
        return;
    }
    //Ictrl == 0: temperature control
    else if(Ictrl == 0){
        max_timestep = abs(Eincr);
        logger.info("Temperature loading control");
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
    prepare_loading(pcrys);
    for(istep = 0; istep < Nsteps; ++istep)
    {
        if(!process_step(pcrys, istep)) {
            Out_texture(pcrys, this_step);
            break;
        }
        output_step_result(pcrys);
        this_step++;
    }
    /* Out_texture(pcrys,Nsteps); */
}

void Process::prepare_loading(Polycs::polycrystal &pcrys){
    timestep_control();
    time_acc = 0;
    temp_atmosphere = temperature_input; //loading temp_init from the txt.in 
    // Temperature is not set, this is the first loading step, free of thermal stress
    if (temperature_ref < 1e-3)  temperature_ref = temperature_input;
    if (pcrys.temperature_poly < 1e-3)  pcrys.set_temperature(temperature_input);
    pcrys.set_boundary_conditions(UDWdot_input, Sdot_input, IUDWdot, ISdot);
}

bool Process::process_step(Polycs::polycrystal &pcrys, int istep) {
    logger.info("");
    logger.info("**********\tSTEP\t" + to_string(istep) + "\t**********");
    logger.info("");
    bool is_convergent = true;
    double pct_step = 0;
    double coeff_step = 1, current_step = 1.;
    int success_count = 0;
    update_progress(pct_step);
    do {
        current_step = min(1.0 - pct_step, coeff_step);
        logger.notice("Step " + to_string(istep) + ":\t" + to_string(pct_step) + " to " + to_string(pct_step + current_step));
        int return_SC = pcrys.EVPSC(istep, current_step * max_timestep);
        if (return_SC != 0) {
            handle_nonconvergence(pcrys, coeff_step, is_convergent, success_count);
            if (!is_convergent) break;
            continue;
        }
        update_after_success(pcrys, current_step, success_count);
        pct_step += current_step;
        if (success_count > 3) coeff_step *= 1.5;
        else if (success_count > 6) coeff_step *= 2;
        coeff_step = min(coeff_step, 1.0);
        update_progress(pct_step);
    } while (pct_step < 1 - 1e-10);
    cout.flush();
    return is_convergent;
}

void Process::handle_nonconvergence(Polycs::polycrystal &pcrys, double &coeff_step, bool &is_convergent, int &success_count) {
    success_count = 0;
    pcrys.restore_status(false);
    logger.warn("Not convergent... Retry with a smaller increment.");
    logger.notice("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=");
    if (isnan(pcrys.error_SC)) coeff_step *= 0.1;
    else coeff_step *= 0.66;
    if (coeff_step < 1e-10) {
        logger.error("Not convergent... Abort.");
        is_convergent = false;
    }
}

void Process::update_after_success(Polycs::polycrystal &pcrys, double current_step, int &success_count) {
    time_acc += current_step * max_timestep;
    update_current_intensity(time_acc);    
    update_temperature(pcrys, current_step);
    custom_vars[5] = Current_intensity;     
    success_count++;
}

void Process::update_current_intensity(double time_acc) {
    if (flag_emode == 0) {
        Current_intensity = 0.0;
    } else if (flag_emode == 1) {
        Current_intensity = J_intensity_pulse(time_acc, duty_ratio_J, Amplitude_J, Frequency);
    } else if (flag_emode == 2) {
        Current_intensity = J_shock_sim(time_acc, max_timestep * Nsteps, Amplitude_J, shock_int, shock_fin);
    } else {
        logger.warn("Error. Please check the emode.");
    }
}

void Process::update_temperature(Polycs::polycrystal &pcrys, double current_step) {
    if (Ictrl == 0) {
        double dtempK = tempK_rate * current_step * max_timestep;
        if (temp_atmosphere < tempK_end && tempK_rate > 0.0) {
            temp_atmosphere += dtempK;
            if (temp_atmosphere > tempK_end) temp_atmosphere = tempK_end;
        } else if (temp_atmosphere > tempK_end && tempK_rate < 0.0) {
            temp_atmosphere += dtempK;
            if (temp_atmosphere < tempK_end) temp_atmosphere = tempK_end;
        }
        pcrys.set_temperature(temp_atmosphere);
    }
}

void Process::output_step_result(Polycs::polycrystal &pcrys) {
    if (!((this_step + 1) % texctrl))
        Out_texture(pcrys, this_step);
    output_info();
    output_phase_info();
    output_grain_info(0);
}

void Process::Out_texture(Polycs::polycrystal &pcrys, int this_step)
{
    IOFormat Outformat(StreamPrecision);
    logger.notice("Output texture at step " + to_string(this_step));
    for(int phase_id = 0; phase_id < phaseCount; ++phase_id){
        fstream &tex_out_i = tex_out[phase_id];
        if (!tex_out_i.is_open()) {
            logger.error("Failed to open output file for phase " + to_string(phase_id));
            continue;  // 跳过该相
        }
        tex_out_i << "TEXTURE AT STEP = " << this_step+1 << endl;
        tex_out_i << setprecision(4) << pcrys.get_ell_axis().transpose().format(Outformat)<< endl; 
        tex_out_i << setprecision(4) << pcrys.get_ellip_ang().transpose().format(Outformat) << endl << endl;
        pcrys.printEuler(tex_out_i, phase_id);
        tex_out_i << endl;
    }
}

void Process::Out_texset(int input){texctrl = input;}

void Process::set_tempK_control(double rate, double end_temp)
{
    tempK_rate = rate;
    tempK_end = end_temp;
    if (tempK_rate > 0) {
        logger.notice("Temperature is increasing at a rate of " + to_string(tempK_rate) + " K/s");
    } else if (tempK_rate < 0) {
        logger.notice("Temperature is decreasing at a rate of " + to_string(-tempK_rate) + " K/s");
    } else {
        logger.notice("Temperature is constant.");
    }
}

double Process::calculate_current_intensity(double time) const {
    switch (flag_emode) {
        case 0:  // 无电流模式
            return 0.0;
        case 1:  // 脉冲电流模式
            return J_intensity_pulse(time, duty_ratio_J, Amplitude_J, Frequency);
        case 2: {  // 冲击电流模式
            double total_time = max_timestep * static_cast<double>(Nsteps);
            return J_shock_sim(time, total_time, Amplitude_J, shock_int, shock_fin);
        }
        default:
            logger.warn("Invalid electric current mode. Defaulting to no current.");
            return 0.0;
    }
}

