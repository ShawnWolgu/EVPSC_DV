#pragma once

#include "common/base.h"
#include <regex>
#include "Eigen/src/Core/Matrix.h"

struct materialPhase;
namespace Procs{ class Process; }

void skipLine(std::ifstream& file, int n = 1);
int EVPSCinput(string &,string &,string &, Procs::Process &); //read .in file
int texinput(string, Polycs::polycrystal &); //read .tex file
int texinput(const string& fname, int &grain_count, vector<Vector4d> &euler_info);
int sxinput(const string& fname, materialPhase &mat);
int loadinput(string, Procs::Process &Proc);
VectorXd getnum(string, int);
// add more input functions here
vector<double> getnum_vec(string, int);
vector<double> get_vector(MatrixXd &matrix);
vector<double> get_vector(VectorXd &matrix);
string extractPath(const std::string& ftex, int a);
