#pragma once

#include "common/base.h"
#include <regex>
#include "Eigen/src/Core/Matrix.h"

namespace Procs{ class Process; }

void skipLine(std::fstream& file, int n = 1);
int EVPSCinput(string &,string &,string &, Procs::Process &); //read .in file
int sxinput(string, Polycs::polycrystal &); //read .sx file
int texinput(string, Polycs::polycrystal &); //read .tex file
int loadinput(string, Procs::Process &Proc);
VectorXd getnum(string, int);
// add more input functions here
vector<double> getnum_vec(string, int);
vector<double> get_vector(MatrixXd &matrix);
vector<double> get_vector(VectorXd &matrix);
MatrixXd cal_sn_info(MatrixXd &Min, vector<double> Mabc, vector<double> Trans_Miller, int Miller_n, int system_n);
