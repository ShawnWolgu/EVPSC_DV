#include "Processes.h"
#include "global.h"

void initial_output_files(){
    IOFormat Outformat(StreamPrecision);
    Polycs::polycrystal* pcrys = &global_polycrys;

    ss_out_csv << "EVM,SVM,E11,E22,E33,E23,E13,E12,S11,S22,S33,S23,S13,S12,TempK\n";
    ss_out_csv << 0.0 << "," << 0.0 << ",";
    for(int i = 0; i < 6; ++i) ss_out_csv << 0.0 << ",";
    for(int i = 0; i < 6; ++i) ss_out_csv << 0.0 << ",";
    ss_out_csv << pcrys->temperature_poly << endl;

    ave_ss_out << "EVM,SVM,E11,E22,E33,E23,E13,E12,S11,S22,S33,S23,S13,S12\n";
    ave_ss_out << 0.0 << "," << 0.0;
    for(int i = 0; i < 6; ++i) ave_ss_out << "," << 0.0;
    for(int i = 0; i < 6; ++i) ave_ss_out << "," << 0.0;
    ave_ss_out << endl;

    //output dislocation density
    int family_num = pcrys->family_count;
    density_out << "EVM";
    for(int i = 0; i < family_num; ++i) density_out << "," << "Mode " << i+1;
    density_out << endl;
    density_out << 0.0;
    for(int i = 0; i < family_num; ++i) density_out << "," << 0.0;
    density_out << endl;

    acc_strain_out << "EVM";
    for(int i = 0; i < family_num; ++i) acc_strain_out << "," << "Mode " << i+1;
    acc_strain_out << endl;
    acc_strain_out << 0.0;
    for(int i = 0; i < family_num; ++i) acc_strain_out << "," << 0.0;
    acc_strain_out << endl;

    crss_out << "EVM";
    for(int i = 0; i < family_num; ++i) crss_out << "," << "Mode " << i+1;
    crss_out << endl;
    crss_out << 0.0;
    for(int i = 0; i < family_num; ++i) crss_out << "," << 0.0;
    crss_out << endl;

    int custom_length = custom_vars.size();
    custom_out << "EVM";
    for(int i = 0; i < custom_length; ++i) custom_out << "," << "Custom Var " << i+1;
    custom_out << endl;
    custom_out << 0.0;
    for(int i = 0; i < custom_length; ++i) custom_out << "," << 0.0;
    custom_out << endl;

    grain_out << "Grain ID,EVM_g,SVM_g,Temp";
    int modes_count = pcrys->g[0].modes_num;
    for(int i = 0; i < modes_count; ++i) grain_out << ",Mode_density " << i+1;
    grain_out << endl;
    grain_out << 0 << "," << 0.0 << "," << 0.0 << "," << pcrys->temperature_poly;
    for(int i = 0; i < modes_count; ++i) grain_out << "," << 0.0;
    grain_out << endl;

    logger.notice("Output files are initialized.");
}

void output_info(){
    IOFormat Outformat(StreamPrecision);
    Polycs::polycrystal* pcrys = &global_polycrys;
    double equi_strain = calc_equivalent_value(pcrys->get_Eps_m());

    ss_out_csv << equi_strain << ",";
    ss_out_csv << calc_equivalent_value(pcrys->get_Sig_m()) << ",";
    for(int i = 0; i < 6; ++i) ss_out_csv << pcrys->get_Eps_m()(i) << ",";
    for(int i = 0; i < 6; ++i) ss_out_csv << pcrys->get_Sig_m()(i) << ",";
    ss_out_csv << pcrys->temperature_poly << endl;

    ave_ss_out << equi_strain << "," ;
    ave_ss_out << calc_equivalent_value(pcrys->get_Sig_ave());
    for(int i = 0; i < 6; ++i) ave_ss_out << "," << pcrys->get_Eps_m()(i);
    for(int i = 0; i < 6; ++i) ave_ss_out << "," << pcrys->get_Sig_ave()(i);
    ave_ss_out << endl;

    //output dislocation density
    int family_num = pcrys->family_count;
    density_out << equi_strain;
    for(int i = 0; i < family_num; ++i) density_out << "," << pcrys->density_by_family[i];
    density_out << endl;

    acc_strain_out << equi_strain;
    for(int i = 0; i < family_num; ++i) acc_strain_out << "," << pcrys->acc_strain_by_family[i];
    acc_strain_out << endl;

    crss_out << equi_strain;
    for(int i = 0; i < family_num; ++i) crss_out << "," << pcrys->crss_by_family[i];
    crss_out << endl;
    
    int custom_length = custom_vars.size();
    custom_out << equi_strain;
    for(int i = 0; i < custom_length; ++i) custom_out << "," << custom_vars[i];
    custom_out << endl;
}

void output_grain_info(int i){
    IOFormat Outformat(StreamPrecision);
    Polycs::polycrystal* pcrys = &global_polycrys;
    if (i >= pcrys->grains_num) return;
    grain* grain = &pcrys->g[i];
    int modes_count = grain->modes_num;
    grain_out << grain->grain_i << "," << calc_equivalent_value(grain->get_strain_g()) << ",";
    grain_out << calc_equivalent_value(grain->get_stress_g()) << "," << pcrys->temperature_poly;
    for(int j = 0; j < modes_count; ++j) grain_out << "," << grain->gmode[j]->disloc_density;
    grain_out << endl;
}
