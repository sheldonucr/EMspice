/*
 **********************************************************
 
 				Power Grid Simulator
 		(Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: function_ibm.h,v $
 *    Authors: Xin Huang, Han Zhou, and Shuyuan Yu
 *    Functions: header file
 *    $Date: 2020/03/17 $
 *
 */


#ifndef FUNCTION_IBM_XH_H_
#define FUNCTION_IBM_XH_H_

#define MAX 100
#define MIN 0.1
#include "parser.h"
#include "parameter.h"
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <stack>
#include <sstream>
#include <iostream>
#include <fstream>
#include <climits>
#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <itpp/base/smat.h>
#include <itpp/base/mat.h>

using namespace itpp;
using namespace std;

struct node {
    string name;
    int net;
    int x;
    int y;
};

struct resInfo {
    node left_node;
    node right_node;
    double value;
    double h;
    double w;
    double l;
    int rc;
    int rcnum;
    double curden;
    int ntree;
    double gradius;
    double tnuc;
    string voidnode;
    double T;
    resInfo()
    {
        l=100e-6;
        w=6e-7;
        h=6.5e-7;
        T=373;
        tnuc = INT_MAX;
    }
};

struct geoInfo {
    node left_node;
    node right_node;
    double h;
    double w;
    int rc;
    int rcnum;
    double gradius;
    double l;
    int ntree;
    geoInfo()
    {
        l=100e-6;
        w=7.5e-7;
        h=6.5e-7;
    }
};

struct treeInfo {
    double volume;
    double len;
    double jlSquare;
    double stressL;
};

typedef map<string, resInfo> res_map;
typedef map<string, string> str_map;
typedef map<string, double> u_map;
typedef set<string> str_set;
typedef vector<double> double_vec;
typedef vector<string> str_vec;
typedef stack<string> str_stack;
typedef vector<node> node_vec;
typedef pair<string, geoInfo> res_geo;
//typedef pair<string, resInfo> res_geo;
typedef vector<res_geo> res_geo_vec;
typedef map<string,double_vec> vHistory;

struct find_node{
    string _name;
    find_node(const string &name):_name(name){}
    bool operator()(const node &n) const{
        return n.name == _name;
    }
};

struct find_res{
    string _rname;
    find_res(const string &rname):_rname(rname){}
    bool operator()(const res_geo &r) const{
        return r.first == _rname;
    }
};

////////////////////////////////////////////////

void GetNodeLoc(string s, int &net, int &x, int &y);

void ParseNetlist(res_map &resmap, char* cktname);

void GetNodeResMap(res_map &rmap, map<string, str_vec> &nmap);

void GetTree(res_map &rmap, map<string, str_vec> &nmap,
		vector<res_geo_vec> &vtree);

void GetTreeInfo(vector<res_geo_vec> &vtree, map<int,treeInfo> &tinfo);

void GetNodeNodeMap(vector<res_geo_vec> &vtree, map<string, str_vec> &nmap,
		vector<map<string, node_vec> > &nnmaps);


////////////////////////////////////////////////

void SaveNodeVolHistory(map<string, double> &vmap, map<string, double_vec> &vhistory);

void SaveNodeVol(int nport, vec pv, str_vec port_name, map<string, double> &vmap);

void GetCurDen(res_map &resmap, map<string, double> &node_vol);

void EM_check(res_map &rmap, map<string,double> &vmap, vector<res_geo_vec> &vtree, map<int, treeInfo> &tImap);

void UpdateRes(str_vec &rgrow, res_map &rmap, vector<res_geo_vec> &vtree, map<int,treeInfo> &tinfo, map<string, double> &smap, double tcur);

void UpdateBackG(matrix* G, NodeList* nodePool, res_map rmap);

void UpdateG(matrix* G, NodeList* nodePool, res_map rmap);

////////////////////////////////////////////////

bool CheckIRDrop(map<string, double> &node_vol, double &max_IR_drop,
		string &minvol_nname, string &first_failure);

///////////////////////////////////////////////

void StressSolve(vector<res_geo_vec> &vtree, map<string, double> &smap,
		res_map &rmap); 

///////////////////////////////////////////////
void CheckNucleation(res_map &rmap, vector<res_geo_vec> &vtree, str_vec &rgrow,
		map<string, double> &smap, map<string, str_vec> &nmap, vector<map<
				string, node_vec> > &nnmaps, double tcur, str_vec &void_nodes);

double MinNucleationTime(res_map &rmap, vector<res_geo_vec> &vtree,
		str_vec &rgrow, map<string, double> &smap, map<string, str_vec> &nmap,
		vector<map<string, node_vec> > &nnmaps, string &min_tnuc_rname,
		bool &isEMimmune, str_vec &void_nodes);


///////////////////////////////
void StresstoFile(map<string, double> &smap, const char* filename);

void RgrowtoFile(str_vec &rgrow, res_map &rmap, str_vec &void_nodes, const char* filename);

void IRdroptoFile(map<string, double> &node_vol, const char* filename);

#endif /* FUNCTION_IBM_XH_H_ */
