/*
 **********************************************************
 
 				Power Grid Simulator
 		(Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: mna.cpp,v $
 *    Authors: Duo Li and Han Zhou
 *    Functions: MNA direct solver using CXSparse
 *    $Date: 2020/03/17 $
 *
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <itpp/base/timing.h>
#include <itpp/base/smat.h>
#include <itpp/base/mat.h>
#include <itpp/base/vec.h>
#include <itpp/base/specmat.h>
#include <itpp/base/algebra/lapack.h>
#include <itpp/base/algebra/lu.h>
#include <itpp/base/algebra/ls_solve.h>
#include <itpp/base/algebra/svd.h>
#include <itpp/signal/transforms.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/math/min_max.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/sort.h>
#include "umfpack.h"
#include "cs.h"
#include "mna.h"

using namespace itpp;
using namespace std;

void form_vec(vec &v, double start, double step, double end) {
	if (step != 0) {
		v.set_size(floor_i((end-start)/step) + 1 + 1);
		for (int i = 0; i < v.size(); i++) {
			v(i) = start + i*step;
        }
        if (v.get(v.length()-1) > end)
			v.set(v.length()-1, end);
	} else {
        v.set_size(1);
        v(0) = start;
    }
}

void mna_solve(cs_dl *G, cs_dl *C, cs_dl *B, 
			   Source *VS, int nVS, Source *IS, int nIS, 
			   double tstep, double tstop, const ivec &port, mat &sim_port_value, 
			   vector<int> &tc_node, vector<string> &tc_name, int num, int ir_info,
			   char *ir_name)
{
  Real_Timer interp2_run_time;
  Real_Timer ir_run_time;

  vec max_value, min_value, avg_value, ir_value;
  double max_ir, min_ir, avg_ir;
  int max_ir_idx, min_ir_idx;
  ivec sorted_max_value_idx, sorted_min_value_idx, 
	sorted_avg_value_idx, sorted_ir_value_idx;
  int nNodes = tc_node.size();
  int display_num = num<tc_node.size()?num:tc_node.size();
  max_value.set_size(nNodes);
  min_value.set_size(nNodes);
  avg_value.set_size(nNodes);
  sorted_max_value_idx.set_size(nNodes);
  sorted_min_value_idx.set_size(nNodes);
  sorted_avg_value_idx.set_size(nNodes);
  sorted_ir_value_idx.set_size(nNodes);
  UF_long n = G->n;
  vec u_col(nVS+nIS);
  u_col.zeros();
  vec w(n);
  w.zeros();
  vec ts;
//   cout<<"reach 1"<<endl;
  form_vec(ts, 0, tstep, tstop);
  sim_port_value.set_size(port.size(), ts.size());
  double temp;
  int* cur = new int[nVS+nIS];
  for(int i = 0; i < nVS+nIS; i++){
	cur[i] = 0;
  }
//  cout << "reach 2"<<endl;
  vector<int> const_v, const_i, var_v, var_i;
  for(int j = 0; j < nVS; j++){
	if (VS[j].time.size() == 1)
	  const_v.push_back(j);
	else
	  var_v.push_back(j);
  }
  for(int j = 0; j < nIS; j++){
	if (IS[j].time.size() == 1)
	  const_i.push_back(j);
	else
	  var_i.push_back(j);
  }
//   cout<<"reach 3"<<endl;
  /* DC simulation */
  for(vector<int>::iterator it = const_v.begin(); it != const_v.end(); ++it){
	u_col(*it) = VS[*it].value(0);
  }
  for(vector<int>::iterator it = const_i.begin(); it != const_i.end(); ++it){
	u_col(nVS+(*it)) = IS[*it].value(0);
  }
  // cout << "reach here"<<endl;
  for (int i = 0; i < 1; i++){
  	for(vector<int>::iterator it = var_v.begin(); it != var_v.end(); ++it){
	  u_col(*it) = VS[*it].value(0); // XXLiu
  	}
	//cout << "reach here"<<endl;
  	for(vector<int>::iterator it = var_i.begin(); it != var_i.end(); ++it){
  	  u_col(nVS+(*it)) = IS[*it].value(0); // XXLiu
  	}
  	cs_dl_gaxpy(B, u_col._data(), w._data());
  }
//  cout << "reach 4"<<endl;
  vec xres(n);
  xres.zeros();
  vec x(n);
  x.zeros();
  cs_dls *Symbolic;
  cs_dln *Numeric;
  int order = 2;
  double tol = 1e-10;
  Symbolic = cs_dl_sqr(order, G, 0);
  Numeric = cs_dl_lu(G, Symbolic, tol);
  cs_dl_ipvec(Numeric->pinv, w._data(), x._data(), n);
  cs_dl_lsolve(Numeric->L, x._data());
  cs_dl_usolve(Numeric->U, x._data());
  cs_dl_ipvec(Symbolic->q, x._data(), xres._data(), n);  
  cs_dl_sfree(Symbolic);
  cs_dl_nfree(Numeric);
//  cout << "reach 5"<<endl;
  for (int j = 0; j < port.size(); j++){
	sim_port_value.set(j, 0, xres(port(j)));
  }
  if (ir_info){
	ir_run_time.start();
	for (int j = 0; j < nNodes; j++){
	  max_value(j) = xres(tc_node[j]);
	  min_value(j) = xres(tc_node[j]);
	  avg_value(j) = xres(tc_node[j]);
	}
	ir_run_time.stop();
  }
  
  std::cout.setf(std::ios::fixed,std::ios::floatfield); 
  std::cout.precision(2);
  std::cout << "interpolation2  \t: " << interp2_run_time.get_time() << std::endl;
  std::cout << "IR analysis     \t: " << ir_run_time.get_time() << std::endl;
}
