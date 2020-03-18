/*
 **********************************************************
 
				 Power Grid Simulator
		 (Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: mna.h,v $
 *    Authors: Duo Li and Han Zhou
 *    Functions: header file for MNA direct solver
 *    $Date: 2020/03/17 $
 *
 */


#ifndef MNA_H
#define MNA_H

#include <itpp/base/vec.h>
#include <itpp/base/smat.h>
#include <itpp/base/mat.h>
#include "cs.h"

using namespace itpp;
using namespace std;

typedef struct{
    vec time;
    vec value;
} Source;


void form_vec(vec &v, double start, double step, double end);

void mna_solve(cs_dl *G, cs_dl *C, cs_dl *B, Source *VS, int nVS, Source *IS, int nIS, double tstep, double tstop, const ivec &port, mat &sim_port_value, vector<int> &tc_node, vector<string> &tc_name, int num, int ir_info, char *ir_name);

#endif
