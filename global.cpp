#include "global.h"
#include "Polycrystals.h"
#include "Processes.h"
using namespace std;

double temp_atmosphere = 0.0;
Logger logger;
Polycs::polycrystal global_polycrys;
Procs::Process global_proc;

void update_progress(double progress_f)
{
    const int bar_width = 70;
    int bar_position = (int)(bar_width * progress_f);

    std::cout << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < bar_position) {
            std::cout << "=";
        } else if (i == bar_position) {
            std::cout << ">";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << (int)(progress_f * 100) << "%\r";
    std::cout.flush();
}
