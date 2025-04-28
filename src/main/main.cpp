#include "common/base.h"
#include "common/common.h"
#include <time.h>
#include <chrono>
#include <regex>
#include "io/Input.h"
#include "io/Output.h"

int main()
{
    auto start_time = std::chrono::high_resolution_clock::now();
    logger.info("EVPSC_CPP Start!");

    string tex_path, sx_path, load_path;
    if(EVPSCinput(tex_path, sx_path, load_path, global_proc)) exit(0);

    global_materials.resize(phaseCount);
    for (int phase_id = 0; phase_id < phaseCount; phase_id++){
        phase_out.push_back(fstream("phase_out_" + std::to_string(phase_id) + ".csv", ios::out));
        tex_out.push_back(fstream("tex_out_" + std::to_string(phase_id) + ".txt", ios::out));
        int grain_count = 0; 
        string tex_path_i = extractPath(tex_path, phase_id);
        vector<Vector4d> eulerData;
        if(texinput(tex_path_i, grain_count, eulerData)) exit(0);
        logger.info("The number of grains in phase " + std::to_string(phase_id) + " is: " + std::to_string(grain_count));

        materialPhase mat;
        string sx_path_i = extractPath(sx_path, phase_id);
        if(sxinput(sx_path_i, mat)) exit(0);
        global_materials[phase_id] = mat;

        global_polycrys.add_grains(phase_id, grain_count, eulerData, &global_materials[phase_id]);
    }
    global_polycrys.weightNormalization();
    global_polycrys.cal_aveGrainSize();

    initial_output_files();
    for (int i = 0; i < processCount; i++){
        string load_path_i = extractPath(load_path, i);
        if(loadinput(load_path_i, global_proc)) exit(0);
        global_proc.loading(global_polycrys);
    }

    density_out.close();
    acc_strain_out.close();
    crss_out.close();
    ss_out_csv.close();
    ave_ss_out.close();
    custom_out.close();
    grain_out.close();
    for(fstream &tex_out_i : tex_out) tex_out_i.close();
    for(fstream &phase_out_i : phase_out) phase_out_i.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    std::string end_message = "EVPSC_CPP End! The run time is: " + 
                             std::to_string(elapsed_seconds.count()) + " sec";
    logger.info(end_message);
    return 0;
}
