/*
 **********************************************************
 
 				Power Grid Simulator
 		(Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: em_cmd_ibm.cpp,v $
 *    Author: Xin Huang and Han Zhou
 *    Functions: main function (IBM Benchmarks)
 *    $Date: 2020/03/17$
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <itpp/base/timing.h>
#include <itpp/base/smat.h>
#include "cs.h"
#include <ctime>
#include "element.h"
#include "mna.h"
#include "parser.h"
#include "function_ibm.h"
using namespace itpp;
using namespace std;

/*//shyu
void EM_check (res_map &rmap, map<string, double> &vmap, vector<res_geo_vec> &vtree) {

  cout << "Check the tree EM constraints.\n" << endl;
  double A;  //branch surface area;
  double Vcirt = 0.07;  //EM critical voltage
  int treeIdx = 0;
  int treeNum = vtree.size();  //total number of trees

  cout << treeNum << endl;

}*///shyu end

int main(int argc, char* argv[]){
  clock_t start;
 // double duration;
  start = clock();

    if (argc == 1 || argc > 3){
        cout << "usage: em_cmd circuit_name\n";
		exit(-1);
	}
  
	Real_Timer parser_run_time;
	CPU_Timer parser_cpu_time;
	parser_run_time.start();
	parser_cpu_time.start();
	
	char cktname[100];
	strcpy(cktname, argv[1]);

	int ir_info = 0;
	for (int i = 2; i < argc;){
	  if (strcmp(argv[i],"-ir") == 0){
		ir_info = 1;
		i++;
	  }else{
        cout << "usage: etbr_cmd circuit_name [-ir]\n";
		exit(-1);
	  }
	}

	int display_ir_num = 20;

	Source* VS, * IS;
	NodeList *nodePool;
	int nVS,nIS,nport,nNodes,nL,nsubnode;
	double tstep = 0, tstop = 0;
	int dc_sign = 0;
	ivec port;
	vector<string> port_name;
	vector<int> tc_node;
	vector<string> tc_name;
	matrix* G, *C, *B; 
	int size_G, size_C, row_B, col_B;
	NodeList::Node* node_tmp;

	if((nodePool = new NodeList)==NULL){
	  printf("Out of memory!\n"); exit(1);
	}

	double tcur = 0;
	double max_stress;
	double max_IRdrop;
    string min_tnuc_rname;
	bool isEMimmune = false;
	double max_IR_drop;
	string first_failure;
	string minvol_nname;
	res_map resmap; // excludes vias and _X_n** external resistors
	str_vec res_grow;
	vector<res_geo_vec> trees; // vdd trees
	map<string, str_vec> node_res_map;
	map<string, double> node_vol;
	map<string, double_vec> vhistory;
	map<string, double> node_stress;
	map<int, treeInfo> tree_info; // vdd trees, tree number, <tree volume, tree length>
	vector<map<string, node_vec> > node_node_maps;

	cout << "/****************************************/" << endl;
	cout << "/*                                      */" << endl;
	cout << "/*          STEP1: INITIALIZATION       */" << endl;
	cout << "/*                                      */" << endl;
	cout << "/****************************************/" << endl;

	// 1. parse circut, get resmap & stamps
	printf("Parser ...\n");
	nsubnode = 0;
	parser_sub(cktname, tstep, tstop, nsubnode, nVS, nIS, nL, nodePool);
	printf("get port information.\n");
	
	/* get port information */
	nport = nodePool->numPort();
	port.set_size(nport);
	for (int i = 0; i < nport; i++){
	  node_tmp = nodePool->getPort(i);
	  string pn = nodePool->getPortName(i);
	  port(i) = node_tmp->row_no;
	  port_name.push_back(pn);
	}
	for (int i = 0; i < nodePool->getTCNum(); i++){
	  node_tmp = nodePool->getTCNode(i);
	  string tcn = nodePool->getTCName(i);
	  tc_node.push_back(node_tmp->row_no);
	  tc_name.push_back(tcn);
	}
	if (tstep == 0 && tstop == 0){
	  dc_sign = 1;
	}

	/* initialize G,C,B,VS,IS */
	nNodes = nodePool->numNode();
	printf("node number: %d\n", nNodes);
	printf("voltage source number: %d\n", nVS);
	printf("current source number: %d\n", nIS);
	size_G = nNodes + nL + nVS + nsubnode;
	size_C = size_G;
	row_B = size_G;
	col_B = nVS + nIS;
	
    if((G = new matrix(size_G, size_G)) == NULL){
	  printf("Out of memory!\n"); exit(1);
	}
	
	if((C = new matrix(size_C, size_C)) == NULL){
	  printf("Out of memory!\n"); exit(1);
	}

	if((B = new matrix(row_B, col_B)) == NULL){
	  printf("Out of memory!\n"); exit(1);
	}
	if((VS = new Source[nVS])==NULL){ 
	  printf("Out of memory.\n"); exit(1); 
	}
	if((IS = new Source[nIS])==NULL){ 
	  printf("Out of memory.\n"); exit(1); 
	}
  
	printf("stamp circuit...\n");
	cs_dl* Gs, *Cs, *Bs;
	if(nsubnode != 0){
	  stamp_sub(cktname, nL, nIS, nVS, nNodes, tstep, tstop, VS, IS, G,  C,  B, nodePool); // XXLiu
	  Bs = B->mat2csdl();
	  delete B;
	  printf("B cs done.\n");
	  Cs = C->mat2csdl();
	  printf("C cs done.\n");
	  delete C;
	  Gs = G->mat2csdl();
	  delete G;
	  printf("G cs done.\n");
	}
	else{
	  stampB(cktname, nL, nIS, nVS, nNodes, tstop, VS, IS, B, nodePool);
	  Bs = B->mat2csdl();
	  delete B;
	  printf("B cs done.\n");
	  stampC(cktname, nL, nVS, nNodes, C, nodePool);
	  Cs = C->mat2csdl();
	  delete C;
	  printf("C cs done.\n");
	  stampG(cktname, nL, nVS, nNodes, G, nodePool);
	  Gs = G->mat2csdl();
	  printf("G cs done.\n");
	  printf("finish stamp.\n");
	}

	/* get resmap */
	printf("parse netlist, get resmap & eff_resmap......\n");
	ParseNetlist(resmap, cktname);

	/* get trees */
	printf("get trees...\n");
	GetNodeResMap(resmap, node_res_map);
    map<string, str_vec>::iterator it;
    for(it=node_res_map.begin(); it!=node_res_map.end(); ++it){
        if(it->second.size()==0){
            cout << "WARNING: NO BRANCH CONNECTED TO NODE:: " << it->first << endl; 
        }
    }

	GetTree(resmap, node_res_map, trees);

	parser_run_time.stop();
	parser_cpu_time.stop();
	
	// 2. some info needed
	printf("get tree info...\n");
	GetTreeInfo(trees, tree_info);

	//printf("get node_node_map......\n");
	//GetNodeNodeMap(trees, node_res_map, node_node_maps);

	// 3. get initial node v and branch j
	cout << "**** MNA solver starts ****" << endl;
	mat sim_port_value;
	char ir_name[100];
	strcpy(ir_name, cktname);
	strcat(ir_name, ".ir.mna");
	mna_solve(Gs, Cs, Bs, VS, nVS, IS, nIS, tstep, tstop, port, sim_port_value, tc_node, tc_name, display_ir_num, ir_info, ir_name);
	cout << "**** MNA solver ends ****" << endl;

	/* initial node voltage */
	vec pv = sim_port_value.get_col(0);
	SaveNodeVol(nport, pv, port_name, node_vol);
	SaveNodeVolHistory(node_vol, vhistory); // save node voltage history :: v(t)

	/* initial current density */
	GetCurDen(resmap, node_vol); // node_vol needs to be:: .empty(u)

//SHyu
	EM_check(resmap, node_vol, trees, tree_info);
//SHyu end
}
