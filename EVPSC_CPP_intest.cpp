#include <time.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <Eigen/Dense>

#include "Input.h"
#include "Toolbox.h"
#include "global.h"
#include "Processes.h"
#include "Polycrystals.h"

int main()
{
    logger.info("EVPSC_CPP Start!");

    string tex_path, sx_path, load_path;
    if(EVPSCinput(tex_path, sx_path, load_path, global_proc)) exit(0);
    if(texinput(tex_path, global_polycrys)) exit(0);
    if(sxinput(sx_path, global_polycrys)) exit(0);
    if(loadinput(load_path, global_proc)) exit(0);
    
    double Start = clock();
    global_proc.loading(global_polycrys);
    double End = clock();
    string end_message = "EVPSC_CPP End! The run time is: " + std::to_string((double)(End - Start) / (Mtr * CLOCKS_PER_SEC)) + " sec";
    logger.info(end_message);
    return 0;
}
