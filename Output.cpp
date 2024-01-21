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

    int custom_length = custom_vars.size();
    custom_out << "EVM";
    for(int i = 0; i < custom_length; ++i) custom_out << "," << "Custom Var " << i+1;
    custom_out << endl;
    custom_out << 0.0;
    for(int i = 0; i < custom_length; ++i) custom_out << "," << 0.0;
    custom_out << endl;

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

    int custom_length = custom_vars.size();
    custom_out << equi_strain;
    for(int i = 0; i < custom_length; ++i) custom_out << "," << custom_vars[i];
    custom_out << endl;
}
