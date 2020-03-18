/*
 **********************************************************
 
				 Power Grid Simulator
		 (Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: function_ibm.cpp,v $
 *    Authors: Xin Huang and Han Zhou
 *    Functions: First step of EM analysis for IBM-format
 *    power grid benchmarks
 *    $Date: 2020/03/17 $
 *
 */

#include "function_ibm.h"
#include <string.h>
#include <stdio.h>
#include <typeinfo>
#include <iostream>
using namespace std;

/*****************************************************
 * Parser
 ****************************************************/

void GetNodeLoc(string s, int &net, int &x, int &y) {
	int xx = 0;
	int yy = 0;
	int nn = 0;
	int xbegin;
	int ybegin;

	for (int i = 1; i < s.length(); i++) {
		if (s[i] == '_') {
			xbegin = i + 1;
			break;
		}
		nn = nn * 10 + (s[i] - '0');
	}
	for (int i = xbegin; i < s.length(); i++) {
		if (s[i] == '_') {
			ybegin = i + 1;
			break;
		}
		xx = xx * 10 + (s[i] - '0');
	}
	for (int i = ybegin; i < s.length(); i++) {
		yy = yy * 10 + (s[i] - '0');
	}
	x = xx;
	y = yy;
	net = nn;
}

// Function #1
void ParseNetlist(res_map &rmap, char* cktname) {
	string str1;
	string dump;
	string sline;
	string name;
	string lnode, rnode;
	double value;
	int wflag = 0;

	fstream file1;
	char ch;
	file1.open("tree_width.txt", ios::in);
	ch = file1.get();
	if (file1.eof()) {
		wflag = 1;
	}
	file1.close();

	ifstream file(cktname);

	cout << "Read the netlist..." << endl;
	if (!file.is_open()) {
		cout << "Could not open the file" << endl;
		exit(1);
	}
	while (!file.eof()) {
		getline(file, sline);
		if (sline.empty())
			continue;

		if (sline[0] == 'R' | sline[0] == 'r') {
//debug
//		cout << sline << endl;
//
			istringstream str(sline);
			resInfo info;
			str >> name >> lnode >> rnode >> value;

			bool xflag = (lnode[1] == 'X') || (rnode[1] == 'X'); // external resistor
			if (!xflag) {
				int x1, y1;
				int x2, y2;
				int net1, net2;

				GetNodeLoc(lnode, net1, x1, y1);
				GetNodeLoc(rnode, net2, x2, y2);
				bool vflag = (net1 != net2); // via

				if (!vflag) {
					info.value = value;
					info.left_node.name = lnode;
					info.right_node.name = rnode;
					info.left_node.x = x1;
					info.right_node.x = x2;
					info.left_node.y = y1;
					info.right_node.y = y2;
					info.left_node.net = net1;
					info.right_node.net = net2;

					if ((x1 != x2) && (y1 == y2)) {
						info.rc = 0;
						info.rcnum = y1;
						info.l = abs(x1 - x2) * (L_UNIT/1000);
					} else if ((x1 == x2) && (y1 != y2)) {
						info.rc = 1;
						info.rcnum = x1;
						info.l = abs(y1 - y2) * (L_UNIT/1000);
					} else {
						cout << "error: found nodes::x1!=x2 && y1!=y2!" << endl;
						exit(1);
					}
					switch (net1) {
					case 0:
					case 1:
					case 2:
					case 3:
						info.h = H_M34;
                    //    info.w = W_M34;
						break;
					case 4:
					case 5:
					case 6:
					case 7:
						info.h = H_M56;
                    //    info.w = W_M56;
						break;
					default:
						break;
					}
					if (wflag == 1) {
						info.w = (info.l * G_RESIS_CU) / info.value / info.h;
						//cout << info.w << endl;
					}
					pair<string, resInfo> res_pair(name, info);
					rmap.insert(res_pair);
// debug
//		cout << name << info.left_node.name << info.right_node.name << endl;
//
				}
			}
		}
	}
	file.close();
//SHyu
/*	res_map::iterator resit;
	for(resit = rmap.begin(); resit != rmap.end(); ++resit){
		cout << resit->first << " " << resit->second.left_node.name << " resistance: "<< resit->second.value << endl;
	}
*///SHyu end
}

void PrintTree(vector<res_geo_vec> &vtree) {
	vector<res_geo_vec>::iterator vtit;
	int num = 0;
	for (vtit = vtree.begin(); vtit != vtree.end(); ++vtit) {
		cout << endl << "TREE  " << num << "   :" << endl;
		for (res_geo_vec::iterator it = (*vtit).begin(); it != (*vtit).end(); ++it) {
			cout << it->first << "    ";
		}
		cout << endl;
		num++;
	}
}

void PrintNodeMap(map<string, str_vec> &nmap) {
	map<string, str_vec>::iterator it;
	vector<string>::iterator vit;
	cout << endl;
	cout << "print the node map... " << endl;
	for (it = nmap.begin(); it != nmap.end(); ++it) {
		cout << it->first << "  ";
		for (vit = (it->second).begin(); vit != (it->second).end(); ++vit) {
			cout << *vit << "  ";
		}
		cout << endl;
	}
	cout << endl;
}

bool CheckNumberofNodes(str_set &ckres, res_map &rmap) {
	str_set::iterator it;
	int nnodes = 0;
	str_set cktnodes;
	res_map::iterator rit;
	for (it = ckres.begin(); it != ckres.end(); ++it) {
		rit = rmap.find(*it);
		if (rit == rmap.end()) {
			cout << "find res failed!" << endl;
			exit(1);
		}
		string lnode = rit->second.left_node.name;
		string rnode = rit->second.right_node.name;
		if (cktnodes.find(lnode) == cktnodes.end()) {
			cktnodes.insert(lnode);
			nnodes++;
		}
		if (cktnodes.find(rnode) == cktnodes.end()) {
			cktnodes.insert(rnode);
			nnodes++;
		}
	}
	if (nnodes == (ckres.size() + 1)) {
		//	cout << "correct tree..." << endl;
		return true;
	} else if (nnodes < (ckres.size() + 1)) {
		//	cout << "exists circles:: #branch - #nodes = " << ckres.size() - nnodes << endl;
		return false;
	} else {
		//	cout << " ERROR: #nodes > #branche + 1... " << endl;
		return false;
	}
}

void SaveTree(str_set &ckres, res_map &rmap, vector<res_geo_vec> &vtree,
		int ntree) {
   // cout << endl << "%%%% tree::" << ntree << endl << endl;
	res_geo_vec tree;
	for (str_set::iterator it = ckres.begin(); it != ckres.end(); ++it) {
		resInfo info = rmap.find(*it)->second;
		geoInfo ginfo;
		info.ntree = ntree;
		rmap[*it].ntree = ntree;
		// GET GEO INFO FROM RESMAP
		ginfo.left_node = info.left_node;
		ginfo.right_node = info.right_node;
		ginfo.ntree = info.ntree;
		ginfo.rc = info.rc;
		ginfo.rcnum = info.rcnum;
		ginfo.h = info.h;
		ginfo.w = info.w;
		ginfo.l = info.l;
		//ginfo.R_b = info.R_b;
		tree.push_back(make_pair(*it, ginfo));
	}
	vtree.push_back(tree);
}


//Function #2
void GetNodeResMap(res_map &rmap, map<string, str_vec> &nmap) {

	res_map::iterator rit;
	map<string, str_vec>::iterator nit;

	for (rit = rmap.begin(); rit != rmap.end(); ++rit) {
		nit = nmap.find(rit->second.left_node.name);
		if (nit == nmap.end()) {
			vector<string> resVec;
			resVec.push_back(rit->first);
			nmap.insert(pair<string, vector<string> > (
					rit->second.left_node.name, resVec));
		} else {
			nit->second.push_back(rit->first);
		}

		nit = nmap.find(rit->second.right_node.name);
		if (nit == nmap.end()) {
			vector<string> resVec;
			resVec.push_back(rit->first);
			nmap.insert(pair<string, vector<string> > (
					rit->second.right_node.name, resVec));
		} else {
			nit->second.push_back(rit->first);
		}
	}
	//PrintNodeMap(nmap);
     map<string, str_vec>::iterator it;
    for(it=nmap.begin(); it!=nmap.end(); ++it){
        if(it->second.size()==0){
            cout << "WARNING: NO BRANCH CONNECTED TO NODE:: " << it->first << endl; 
        }
//SHyu
/*
	for(it=nmap.begin(); it!=nmap.end(); ++it) {
		str_vec::iterator itt;
		cout << it->first;
		for (int i=0; i<it->second.size(); i++) {
			cout << " " << it->second[i] ;
		}
		cout << endl;
	}
*/
//SHyu end
    }
}

//Function #3
void GetTree(res_map &rmap, map<string, str_vec> &nmap,
		vector<res_geo_vec> &vtree) {
//SHyu
	res_map::iterator rit;
	int ntree;  // tree number
	int lntree;  // tree number of last branch
	int nNet;  // net number
	int ntbegin;  //tree number begin
	int tree_name_l;  // tree name length
	string sname;  // branch name
	map<string, int> tree_map;  //relation between branch name and its tree number
	vector<string> tree_name_vec;  // tree name vector
	res_geo_vec tree;
	geoInfo ginfo;

	vector<string>::iterator tnvit;  // tree name vector iterator

	FILE *fid;  // open file for writing relationship between tree name and tree number
	fid = fopen("tree_number_r.txt", "a+");

	for (rit=rmap.begin(); rit!=rmap.end(); ++rit) {
		sname = rit->first;
		nNet = 0;
		ntree = 0;
		for (int i = 1; i < sname.length(); i++) {
			if (sname[i] == '-') {
				ntbegin = i+1;
				break;
			}
			nNet = nNet*10 + (sname[i]-'0');
		}
		for (int i=ntbegin; i<sname.length(); i++) {
			if (sname[i] == '-') {
				tree_name_l = (i-1);  // get the length of tree name
				break;
			}
			ntree = ntree*10 + (sname[i]-'0');
		}
		string treename(sname.substr(1,tree_name_l));
		tree_name_vec.push_back(treename);
		pair<string, int> brch_tree_pair(sname, ntree);  // get the relation of branch and its pair
		tree_map.insert(brch_tree_pair);
//debug
//		cout << sname << "  " << ntree << endl;
//
	}
	tree_name_vec.erase(unique(tree_name_vec.begin(), tree_name_vec.end()), tree_name_vec.end());
	for (int i=0; i<tree_name_vec.size(); i++) {
//		cout << "Tree name: " << tree_name_vec[i] << " Tree Index: " << i << endl;
		fprintf(fid, "Tree name: %s Tree Index: %d\r\n", tree_name_vec[i].c_str(), i);
	}
//	cout << "tree name vector size: " << tree_name_vec.size() << endl;

	map<string, int>::iterator trit;
	for (trit=tree_map.begin(); trit!=tree_map.end(); ++trit) {
		ntree = trit->second;  //get tree number
		resInfo info = rmap.find(trit->first)->second;  // find the branch information
		if (trit == tree_map.begin()) {
			lntree = trit->second;
		}
		rmap[trit->first].ntree = ntree;
		info.ntree = ntree;
		ginfo.left_node = info.left_node;
		ginfo.right_node = info.right_node;
		ginfo.ntree = info.ntree;
		ginfo.rc = info.rc;
		ginfo.rcnum = info.rcnum;
		ginfo.h = info.h;
		ginfo.w = info.w;
		ginfo.l = info.l;
		if (ntree == lntree) {
			tree.push_back(make_pair(trit->first, ginfo));
		} else {
			vtree.push_back(tree);
			tree.clear();
			tree.push_back(make_pair(trit->first, ginfo));
		}
		lntree = ntree;

//		cout << "Branch name: " << trit->first << " Tree Index: " << trit->second << endl;
	}
	vtree.push_back(tree);
//	cout << "Tree size: " << vtree.size() << endl;

/*	for (vector<res_geo_vec>::iterator vecit = vtree.begin(); vecit != vtree.end(); ++vecit) {
		tree = *vecit;
		res_geo_vec::iterator brchit;
		for (brchit=tree.begin(); brchit!=tree.end(); ++brchit) {
			cout << "Tree Index: " << brchit->second.ntree << " ";
			cout << "Branch name: " << brchit->first << " ";
			cout << "length: " << brchit->second.l*1e6 << " ";
			cout << "width: " << brchit->second.w*1e6 << " ";
			cout << "thickness: " << brchit->second.h*1e6 << " \n";  
		}
		cout << endl;
	}  */
	fclose(fid);

//rewrite
/*	str_stack nstack; // when a stack is empty(except the initial condition)-> a tree is found
	str_set lres;
	str_set cknode; // when a tree is found, clean it
	str_set ckres; // when a tree is found, put it in treeGroup and clean it
	str_vec::iterator vecit;
	str_set::iterator sit;
	res_map::iterator rit;
	map<string, str_vec>::iterator nit;
	//debug
	int ntree = 0;

	for (rit = rmap.begin(); rit != rmap.end(); ++rit) {
		lres.insert(rit->first);
	}

	while (!lres.empty()) {
		string rname = *(lres.begin()); // randomly get one segment, then do DFS
		ckres.insert(rname);
		lres.erase(lres.begin());
		map<string, resInfo>::iterator it;
		it = rmap.find(rname);
		// if node is unchecked, put in stack & checkedNode
		if (cknode.find(it->second.left_node.name) == cknode.end()) {
			nstack.push(it->second.left_node.name);
			cknode.insert(it->second.left_node.name);
		}
		if (cknode.find(it->second.right_node.name) == cknode.end()) {
			nstack.push(it->second.right_node.name);
			cknode.insert(it->second.right_node.name);
		}
		while (!nstack.empty()) {
			string nname = nstack.top();
			nstack.pop();
			nit = nmap.find(nname);
			// put all nodes (unchecked) of all connected R (unchecked) in stack
			for (vecit = nit->second.begin(); vecit != nit->second.end(); ++vecit) {
				if (ckres.find(*vecit) == ckres.end()) {
					// if node is unchecked, put in stack & cknode
					string lnode;
					string rnode;
					lnode = rmap.find(*vecit)->second.left_node.name;
					rnode = rmap.find(*vecit)->second.right_node.name;

					if (cknode.find(lnode) == cknode.end()) {
						nstack.push(lnode);
						cknode.insert(lnode);
					}
					if (cknode.find(rnode) == cknode.end()) {
						nstack.push(rnode);
						cknode.insert(rnode);
					}
					ckres.insert(*vecit);
					lres.erase(lres.find(*vecit));
				}
			}
		}
		if (nstack.empty()) {
			bool flag = CheckNumberofNodes(ckres, rmap); // FOR DEBUGING; ONLY SAVE CORRECT TREES!!!
			if (flag) {
				SaveTree(ckres, rmap, vtree, ntree);
				ntree++; // number of tree
			}
			ckres.clear();
			cknode.clear();
		}
	}
*/
//	cout << "The # of tree: " << vtree.size() << endl;
//	PrintTree(vtree);
}


//Function #4
void GetTreeInfo(vector<res_geo_vec> &vtree, map<int, treeInfo> &tinfo) {
	if (vtree.empty()) {
		cout << "tree vector is empty!" << endl;
		exit(1);
	}
	treeInfo info;
	vector<res_geo_vec>::iterator vtit;
	res_geo_vec::iterator it;
	res_geo_vec tree;
	double Vtree = 0;
	double ltree = 0;
	int ntree = 0;
	for (vtit = vtree.begin(); vtit != vtree.end(); ++vtit) {
		tree = *vtit;
		for (it = tree.begin(); it != tree.end(); ++it) {
			Vtree += it->second.h * it->second.w * it->second.l;
			ltree += it->second.l;
		}
		info.volume = Vtree;
		info.len = ltree;
		tinfo.insert(make_pair(ntree, info));
		Vtree = 0;
		ltree = 0;
		ntree++;
	}
}


//Function #5
void GetNodeNodeMap(vector<res_geo_vec> &vtree, map<string, str_vec> &nmap,
		vector<map<string, node_vec> > &nnmaps) {
	node Lnode;
	node Rnode;
	node lnode;
	node rnode;
	str_vec adj_res;
	node_vec adj_node;
	res_geo_vec tree;
	res_geo_vec::iterator it;
	str_vec::iterator rit; // RES
	vector<res_geo_vec>::iterator tit;

	for (tit = vtree.begin(); tit != vtree.end(); ++tit) {
		tree = *tit;
		map<string, node_vec> tnmap;
		for (it = tree.begin(); it != tree.end(); ++it) {
			Lnode = it->second.left_node;
			Rnode = it->second.right_node;
			// If Lnode is not in map (KEY).
			if (tnmap.find(Lnode.name) == tnmap.end()) {
				adj_res = nmap[Lnode.name];
				for (rit = adj_res.begin(); rit != adj_res.end(); ++rit) {
					res_geo_vec::iterator tmp_it;
					tmp_it = find_if(tree.begin(), tree.end(), find_res(*rit));
					if (tmp_it == tree.end()) {
						cout << "find adj res in tree failed.." << endl;
						exit(1);
					}
					lnode = tmp_it->second.left_node;
					rnode = tmp_it->second.right_node;
					if ((lnode.name != Lnode.name) && (find_if(
							adj_node.begin(), adj_node.end(), find_node(
									lnode.name)) == adj_node.end())) {
						adj_node.push_back(lnode);
					}
					if ((rnode.name != Lnode.name) && (find_if(
							adj_node.begin(), adj_node.end(), find_node(
									rnode.name)) == adj_node.end())) {
						adj_node.push_back(rnode);
					}
				}
				tnmap[Lnode.name] = adj_node;
				adj_node.clear();
			}
			// If Rnode is not in map (KEY).
			if (tnmap.find(Rnode.name) == tnmap.end()) {
				adj_res = nmap[Rnode.name];
				for (rit = adj_res.begin(); rit != adj_res.end(); ++rit) {
					res_geo_vec::iterator tmp_it;
					tmp_it = find_if(tree.begin(), tree.end(), find_res(*rit));
					if (tmp_it == tree.end()) {
						cout << "find adj res in tree failed.." << endl;
						exit(1);
					}
					lnode = tmp_it->second.left_node;
					rnode = tmp_it->second.right_node;
					if ((lnode.name != Rnode.name) && (find_if(
							adj_node.begin(), adj_node.end(), find_node(
									lnode.name)) == adj_node.end())) {
						adj_node.push_back(lnode);
					}
					if ((rnode.name != Rnode.name) && (find_if(
							adj_node.begin(), adj_node.end(), find_node(
									rnode.name)) == adj_node.end())) {
						adj_node.push_back(rnode);
					}
				}
				tnmap[Rnode.name] = adj_node;
				adj_node.clear();
			}
		}
		nnmaps.push_back(tnmap);
	}
}


/*****************************************************
 *
 * Electrical part:: save nodal voltage history:: v(t)
 *
 ****************************************************/
//Function #7
void SaveNodeVolHistory(map<string, double> &vmap,
		map<string, double_vec> &vhistory) {
	map<string, double>::iterator it;
	for (it = vmap.begin(); it != vmap.end(); it++) {
		if (vhistory.find(it->first) == vhistory.end()) {
			double_vec tmp;
			tmp.push_back(it->second);
			vhistory[it->first] = tmp;
		} else {
            // for debug
            int len = vhistory[it->first].size();
            double vold = vhistory[it->first][len-1];
            if((vold - it->second)>1e-2 || (vold - it->second)<-1e-2)
                cout << "**Vol Change Node:" << it->first << endl;

			vhistory[it->first].push_back(it->second);
		}
	}
}

//Function #8
void SaveNodeVol(int nport, vec pv, str_vec port_name, map<string,double> &vmap) {
	for (int i = 0; i < nport; i++) {
		vmap.insert(pair<string, double> (string(port_name[i]), pv(i)));
	}
/*	map<string, double>::iterator vit;
	for(vit=vmap.begin(); vit!=vmap.end(); ++vit) {
		cout << vit->first << "  " << vit->second << endl;
	}  //output voltage map
*/
}

//Function #9
void GetCurDen(res_map &rmap, map<string, double> &node_vol) {
	res_map::iterator it;
//SHyu output curden
//
	string branchname;
	int lnet, rnet, lx, rx, ly, ry;  //define node infomation which can be obtainedfrom map 
//
	double brchleng, brchwidth;  //branch length, width
	int noderecx = 0;
	int noderecy = 0;
//SHyu end
	if (rmap.empty()) {
		cout << "the rmap is empty" << endl;
		exit(1);
	}
	for (it = rmap.begin(); it != rmap.end(); it++) {
		resInfo info = it->second;
//SHyu output node information
//
		lnet = info.left_node.net;
		lx = info.left_node.x;
		ly = info.left_node.y;
		rnet = info.right_node.net;
		rx = info.right_node.x;
		ry = info.right_node.y;
		brchleng = info.l;
		brchwidth = info.w;
//
//SHyu end

		map<string, double>::iterator it1 = node_vol.find(info.left_node.name);
		map<string, double>::iterator it2 = node_vol.find(info.right_node.name);
		if (it1 == node_vol.end()) {
			cout << "Node:: " << info.left_node.name << " is not found!"
					<< endl;
			exit(1);
		}
		if (it2 == node_vol.end()) {
			cout << "Node:: " << info.right_node.name << " is not found!" << endl;
			exit(1);
		}
		double vl = it1->second;
		double vr = it2->second;
		double curden = (vl - vr) / info.value / info.w / info.h;// decrease the current source
		(it->second).curden = curden;
// deal with the x and y location of left node is smaller than the location of right node
		if ((lx-rx+ly-ry)>0) {
			(it->second).curden *= -1;
		}
//
		//!Test_hzhou
//		cout << it->first << " current density = " << curden << ", width = " << info.w*1e9 << endl;
//SHyu output curden
/*
		FILE *tree00;	
		tree00 = fopen("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_curden.txt","a+");
		branchname = it->first;
		fprintf(tree00,"\n%s, curden = %f",branchname.c_str(),curden);
		fclose(tree00);
		FILE *tree000;
		tree000 = fopen("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_leng.txt","a+");
		fprintf(tree000, "A%f\r\n", brchleng);
		fclose(tree000);
		FILE *tree0;
		tree0 = fopen("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_info.txt","a+");
		if((noderecx - lx)!=0 && (noderecy - ly)!=0 && noderecx != 0 && noderecy != 0)
		{
		  fprintf(tree0, "\nxx");
		}
		fprintf(tree0, "\n%s n%d_%d_%d n%d_%d_%d %f", branchname.c_str(), lnet, lx, ly, rnet, rx, ry, info.value);
		fclose(tree0);

		noderecx = rx;
		noderecy = ry;
*/
//SHyu end
		//Shyu output width
		
		//FILE *treew;
		//treew = fopen("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_width.txt","a");
		//std::ofstream file("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_width.txt");
		//file << it->first;
		//fprintf(treew,"%S tree_width = %f\n",it->first,info.w*1e9);

		//cout << typeid(it->first).name() << endl;

		//fclose(treew);
		
		//Shyu output width end

	}
/*	FILE *tree00;
	tree00 = fopen("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_info.txt","a+");
	fprintf(tree00, "\nxx");
	fclose(tree00);
*/
}


//fucntion 10: EM check

void EM_check (res_map &rmap, map<string, double> &vmap, vector<res_geo_vec> &vtree, map<int, treeInfo> &tImap) {

  cout << "Check the tree EM constraints." << endl;	
  double A = 0;  //branch surface area
  double A_t = 0;  //total area
  double VA = 0;  //the product of voltage and one of the node adjacent area
  double w = 0;  //the width of the branch
  double l = 0;  //the length of the branch
  double l_cath;
  double l_cath_l;
  double l_cath_r;
  double h = 0;  //the thickness of the branch
  double Vcrit = 3.69*1e-3;  //EM critical voltage
  double Ve = 0;  //EM voltage
  int treeIdx = 0; 
  int treeNum = vtree.size();  //total number of trees
  int numFT = 0;  //# of failed tree
  res_geo_vec brch;  //All the branch information in one tree
  string cathode;
  double lnode_vol;  //find left node voltage based on vmap
  double rnode_vol;  //find right node voltage based on vmap
  double cath_vol;  //cathode voltage
  double curden;  //branch current density
  double curden_cath;
  double curden_cath_l;
  double curden_cath_r;
  double Vcross;  //cross section area volume, also means critical volume
  double Vsat = 0;  //saturation volume
  int cath_flag = 0;  //find the branch that before or after the cathode 
  int cath_left_flag = 0;  // if the cathode is on the left side of the branch
  int cath_right_flag = 0;  // if the cathode is on the right
  int half_fail = 0;
  char ch1;  // used to judge whether the tree_width.txt file is empty
  char ch2;

  vector<res_geo_vec>::iterator vtit;  //tree vector iterator
  res_geo_vec::iterator brchit;  //branch information iterator
  map<string, double>::iterator vit;  //node voltage iterator
  res_map::iterator rit;  //branch resistance and curden information

  cout << treeNum << endl;
  
/*  for(vit=vmap.begin(); vit!=vmap.end(); ++vit) {
	  cout << vit->first << endl;
	  cout << vit->second << endl;
  }
*/ //hzhou Jan
  //FILE *nname; //open file for writing node information
  //nname = fopen("pg_nodes.txt", "a+");

  FILE *tinfo;  //open file for writing tree information
//  tinfo = fopen("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_info.txt","a+");
  tinfo = fopen("tree_info.txt", "a+");
  FILE *curd;  //open file for writing branch currrent density
  curd = fopen("tree_curden.txt", "a+");
//  curd = fopen("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_curden.txt","a+");
  fstream file1;
  fstream file2;
  int wflag = 0;
  int ftflag = 0;
  file1.open("tree_width.txt", ios::in);
  file2.open("tree_ftid.txt", ios::in);
  ch1 = file1.get();
  ch2 = file2.get();
  if(file1.eof()){
    wflag = 1;
  }
  file1.close();
  if(file2.eof()) {
    ftflag = 1;
  }
  file2.close();
  FILE *width;  //open file for writing branch width
  width = fopen("tree_width.txt", "a+");
  FILE *leng;  //open file for wirting branch length
  leng = fopen("tree_leng.txt", "a+");
  FILE *FTID;  //failed tree Index
  FTID = fopen("tree_fid.txt", "a+");
  FILE *ftid0;  //failed tree Index for record in C++
  ftid0 = fopen("tree_ftid.txt", "a+");
  FILE *voltage;  // open file for writing node voltage
  voltage = fopen("tree_node_voltage.txt", "a+");
  FILE *node; // open file for writing node axis
  node = fopen("tree_axis.txt", "a+");
  FILE *cath;
  cath = fopen("tree_cathode.txt", "a+");
  
//  FILE *filter_flag;
//  filter_flag = fopen("tree_fail_flag.txt", "r+");

  for(vtit=vtree.begin(); vtit!=vtree.end(); ++vtit) {
    int brchIdx = 0;  //branch ID, used to judge whether it is the first node of the tree
    int brchIdx_max;
//    cout << "Tree " << treeIdx << " ";
    brch = *vtit;

//find cathode
    for(brchit=brch.begin(), brchIdx=0; brchit!=brch.end(); ++brchit) {
      int lx = brchit->second.left_node.x;
      int ly = brchit->second.left_node.y;
      int rx = brchit->second.right_node.x;
      int ry = brchit->second.right_node.y;
      if(brchIdx == 0) {
        cathode = brchit->second.left_node.name;
        vit = vmap.find(cathode);
	cath_vol = vit->second;
      }
      vit = vmap.find(brchit->second.left_node.name);  //get left node voltage
      lnode_vol = vit->second;
      vit = vmap.find(brchit->second.right_node.name);  //get right node voltage
      rnode_vol = vit->second;
      if(lnode_vol <= cath_vol) {
	cathode = brchit->second.left_node.name;
	cath_vol = lnode_vol;
	cath_left_flag = 1;
	cath_right_flag = 0;
	cath_flag = brchIdx;
      }
      if(rnode_vol <= cath_vol) {
	cathode = brchit->second.right_node.name;
	cath_vol = rnode_vol;
	cath_left_flag = 0;
	cath_right_flag = 1;
	cath_flag = brchIdx+1;
	if (lx-rx+ly-ry>0) {
	  cath_flag = brchIdx;
	}
      }
      brchIdx++;
// Vsat minus debug
//	if (treeIdx==66 || treeIdx==30) {
//		cout << brchit->first << " " << brchit->second.left_node.name << " " << lnode_vol << " " << brchit->second.right_node.name << " " << rnode_vol << " Cathode voltage: " << cath_vol << endl;
//	}
//
    }

//    if (treeIdx == 25 || treeIdx == 55 || treeIdx == 56 || treeIdx == 64 || treeIdx == 49 || treeIdx == 50) {
//	fprintf(cath, "%d\t%s\t%f\r\n", treeIdx, cathode.c_str(), cath_vol);
//	fprintf(cath, "%f\t", cath_vol);  
//    }  

// Vsat minus debug
//	if (treeIdx==66 || treeIdx==30) {
//		cout << "cathode: " << cathode << "  Branch Index: " << cath_flag << " Cathode voltage: " << cath_vol << endl;
//	}
//

//    for (vit=vmap.begin(); vit!=vmap.end(); ++vit) {
//      fprintf(voltage, "%s\t%.6f\r\n", vit->first.c_str(), vit->second);
//    }
//    output voltage drop
///*  
    brchIdx_max = brchIdx;  
    brchIdx = 0;
    for(brchit=brch.begin(); brchit!=brch.end(); ++brchit) {
      double lnode_vol;
      double rnode_vol;
      int lx = brchit->second.left_node.x;
      int ly = brchit->second.left_node.y;
      int rx = brchit->second.right_node.x;
      int ry = brchit->second.right_node.y;
      double delta;
      double volt;
      int delta_x;
      int delta_y;
      int temp_x;
      int temp_y;
// debug

      if (lx-rx+ly-ry>0) {
	if (brchIdx == 0) {
	  vit = vmap.find(brchit->second.right_node.name);
	  fprintf(voltage, "%s\t%.6f\r\n", vit->first.c_str(), vit->second);
	  fprintf(node, "%d\t%d\r\n", rx, ry);
	  vit = vmap.find(brchit->second.left_node.name);
	  fprintf(voltage, "%s\t%.6f\r\n", vit->first.c_str(), vit->second);
	  fprintf(node, "%d\t%d\r\n", lx, ly);
	} else {
	  vit = vmap.find(brchit->second.left_node.name);
	  fprintf(voltage, "%s\t%.6f\r\n", vit->first.c_str(), vit->second);
	  fprintf(node, "%d\t%d\r\n", lx, ly);
	}
/*	fprintf(voltage, "%s\t%.6f\r\n", brchit->first.c_str(), rnode_vol);
	fprintf(node, "%d\t%d\r\n", rx, ry);
	delta = (lnode_vol-rnode_vol)/20;
	delta_x = (lx-rx)/20;
	delta_y = (ly-ry)/20;
	for (int k=1; k<20; k++) {
	  volt = rnode_vol + delta*k;
	  temp_x = rx + delta_x*k;
	  temp_y = ry + delta_y*k;
	  fprintf(node, "%d\t%d\r\n", temp_x, temp_y);
	  fprintf(voltage, "\t\t%.6f\r\n", volt);
	}  */
      } else {
	if (brchIdx == 0) {
	  vit = vmap.find(brchit->second.left_node.name);
	  fprintf(voltage, "%s\t%.6f\r\n", vit->first.c_str(), vit->second);
	  fprintf(node, "%d\t%d\r\n", lx, ly);
	  vit = vmap.find(brchit->second.right_node.name);
	  fprintf(voltage, "%s\t%.6f\r\n", vit->first.c_str(), vit->second);
	  fprintf(node, "%d\t%d\r\n", rx, ry);
	} else {
	  vit = vmap.find(brchit->second.right_node.name);
	  fprintf(voltage, "%s\t%.6f\r\n", vit->first.c_str(), vit->second);
	  fprintf(node, "%d\t%d\r\n", rx, ry);
	}
/*	fprintf(voltage, "%s\t%.6f\r\n", brchit->first.c_str(), lnode_vol);
	fprintf(node, "%d\t%d\r\n", lx, ly);
	delta = (rnode_vol-lnode_vol)/20;
	delta_x = (rx-lx)/20;
	delta_y = (ry-ly)/20;
	for (int k=1; k<20; k++) {
	  volt = lnode_vol + delta*k;
	  temp_x = rx + delta_x*k;
	  temp_y = ry + delta_y*k;
	  fprintf(node, "%d\t%d\r\n", temp_x, temp_y);
	  fprintf(voltage, "\t\t%.6f\r\n", volt);
	}   */
//	if (treeIdx == 32) {
//	  cout << brchit->second.left_node.name << " " << lnode_vol << " " << brchit->second.right_node.name << " " << rnode_vol << endl;
//	}
      }
/*	if (treeIdx==87) {
		if (brchIdx == cath_flag) {
			fprintf(vol_6, "%.6f\r\n", vit->second);
		}
	}
	if (treeIdx==92) {
		if (brchIdx == cath_flag) {
			fprintf(vol_7, "%.6f\r\n", vit->second);
		}
	}
	if (treeIdx==883) {
		if (brchIdx == cath_flag) {
			fprintf(vol_38, "%.6f\r\n", vit->second);
		}
	}
	if (treeIdx==884) {
		if (brchIdx == cath_flag) {
			fprintf(vol_39, "%.6f\r\n", vit->second);
		}
	}
	if (treeIdx==902) {
		if (brchIdx == cath_flag) {
			fprintf(vol_44, "%.6f\r\n", vit->second);
		}
	}
	if (treeIdx==903) {
		if (brchIdx == cath_flag) {
			fprintf(vol_45, "%.6f\r\n", vit->second);
		}
	}*/
      brchIdx++;
    }
    brchIdx = 0;

//    cout << cathode << " " << cath_vol << endl;

    for(brchit=brch.begin(); brchit!=brch.end(); ++brchit) {
//      cout << brchit->first << " ";  //output branch name
      vit = vmap.find(brchit->second.left_node.name);  //get relative left node voltage
      lnode_vol = vit->second - cath_vol;
//      cout << vit->second << " " << lnode_vol << " ";  // output left node voltage
      vit = vmap.find(brchit->second.right_node.name);  //get relative right node voltage
      rnode_vol = vit->second - cath_vol;
//      cout << vit->second << " " << rnode_vol << "\n";  //output right node voltage
      w = brchit->second.w;
      l = brchit->second.l;
      h = brchit->second.h;
      A = l*w;
      A_t += A;
      VA += A*(lnode_vol + rnode_vol);
//zeyu curden output
//	rit = rmap.find(brchit->first);
//	curden = rit->second.curden;
//	fprintf(curd, "%s, curden = %f\r\n", rit->first.c_str(), curden);
//
    }
    Ve = VA/A_t;
    Vcross = w * w * h;  //calculate cross section area volume

	if (treeIdx > 854 && treeIdx < 909) {
		Vcross = Vcross/25/2.5;
	}

/*    if (treeIdx == 49 || treeIdx == 50 || treeIdx == 66 || treeIdx == 67) {
	Vcross = Vcross/5;  
    }  */

//zeyu test
//    cout << "Ve: " << Ve << endl;
//
   // int fail_flag;
   // fscanf(filter_flag, "%d", &fail_flag);
   //

    int id_fit=0;

    if(ftflag == 0) {
	ifstream inf("tree_ftid.txt");
	int ftree_id[1000];
	int ii=0;
	while(inf>>ftree_id[ii]) {
	  //cout << "ii =" << ii << "tree_id =" << ftree_id[ii] << endl;
	  ++ii;
	}
	inf.close();
	for(int jj=0; jj<ii; jj++) {
	 // cout << "jj =" << jj << " " << "tree_id =" << ftree_id[jj] << " " << "TreeIdx: " << treeIdx << endl;
	  if(ftree_id[jj] == treeIdx) {
	    id_fit = 1;
	    //cout << treeIdx << ftree_id[jj] << id_fit << endl;
	    break;
	  }
	}
    }

    if(wflag == 0 && id_fit == 1) {
      for(brchit=brch.begin(); brchit!=brch.end(); ++brchit) {
	rit = rmap.find(brchit->first);
        curden = rit->second.curden;
	w = brchit->second.w;
        l = brchit->second.l;
	h = brchit->second.h;
	resInfo info = rit->second;
	fprintf(tinfo, "%s n%d_%d_%d n%d_%d_%d %f\r\n", rit->first.c_str(), info.left_node.net, info.left_node.x, info.left_node.y, info.right_node.net, info.right_node.x, info.right_node.y, info.value);
	fprintf(curd, "%s, curden = %f\r\n", rit->first.c_str(), curden);
	fprintf(leng, "L%.10f\r\n", l);
      }
      fprintf(tinfo, "xx\r\n");
      treeIdx++;
      VA = 0;
      Vsat = 0;
      A_t = 0;
      cath_flag = 0;
      continue;
    }
    else if(wflag == 0 && id_fit == 0) {
	treeIdx++;
	VA = 0;
	Vsat = 0;
	A_t = 0;
	cath_flag = 0;
	continue;
    }

    if(Ve > Vcrit) {
//    if(treeIdx == 42 || treeIdx == 52 || treeIdx == 62 || treeIdx == 72 || treeIdx == 77 || treeIdx == 87 || (treeIdx>90 && treeIdx<145 && (treeIdx%10==2 || treeIdx%10==7)) || (treeIdx>715 && treeIdx<810 && (treeIdx%10==2 || treeIdx%10==7)) || (treeIdx>881 && treeIdx<889) || (treeIdx>901 && treeIdx<909)) {
//    if(treeIdx == 14 || treeIdx == 30 || treeIdx == 32 || treeIdx == 52 || (treeIdx>53 && treeIdx<66)) {
// calculate the saturation volume, judge whether the void will larger than the cross section
      half_fail++;
      brchIdx = 0;
      for(brchit=brch.begin(); brchit!=brch.end(); ++brchit) {
	rit = rmap.find(brchit->first);
        curden = rit->second.curden;  // get branch current density
	
	if (cath_flag == 0) {
	  if (brchIdx == 0) {
	    curden_cath = curden;
	  }
	} else if (cath_flag == brchIdx_max) {
	  if (brchIdx == brchIdx_max) {
	    curden_cath = curden;
	  }
	} else {
	  if (brchIdx == cath_flag-1) {  // since the cathode is always on the left
	    curden_cath_l = curden;
	  } else if (brchIdx == cath_flag) {
	    curden_cath_r = curden;
	  }
	}
	
	w = brchit->second.w;
        l = brchit->second.l;
	h = brchit->second.h;

	if (cath_flag == 0) {
	  if (brchIdx == 0) {
	    l_cath = l;
	  }
	} else if (cath_flag == brchIdx_max) {
	  if (brchIdx == brchIdx_max) {
	    l_cath = l;
	  }
	} else {
	  if (brchIdx == cath_flag-1) {  // since the cathode is always on the left
	    l_cath_l = l;
	  } else if (brchIdx == cath_flag) {
	    l_cath_r = l;
	  }
	}
//Vsat minus debug
//	if (treeIdx==66 || treeIdx==30) {
//		cout << brchit->first << "   curden: "  << curden << endl;
//	}
// The direction of current of the branch after the cathode should add "-"
	if (cath_flag == 0) {
	  if (brchIdx == 0) {
	    Vsat -= (curden_cath*l_cath*G_RESIS_CU*E*G_Z/G_OMEGA)*(l_cath*w*h/2/G_B);
	  } else {
	    Vsat -= (curden_cath*l_cath + curden*l*2)*G_RESIS_CU*E*G_Z/G_OMEGA*(l*w*h/2/G_B);
	  }
	} else if (cath_flag == brchIdx_max) {
	  if (brchIdx == brchIdx_max) {
	    Vsat += (curden_cath*l_cath*G_RESIS_CU*E*G_Z/G_OMEGA)*(l_cath*w*h/2/G_B);
	  } else {
	    Vsat += (curden_cath*l_cath + curden*l*2)*G_RESIS_CU*E*G_Z/G_OMEGA*(l*w*h/2/G_B);
	  }
	} else {
	  if (brchIdx < cath_flag-1) {
	    Vsat += (curden_cath_l*l_cath_l + curden*l*2)*G_RESIS_CU*E*G_Z/G_OMEGA*(l*w*h/2/G_B);
	  } else if (brchIdx == cath_flag-1) {
	    Vsat += (curden_cath_l*l_cath_l*G_RESIS_CU*E*G_Z/G_OMEGA)*(l_cath_l*w*h/2/G_B);
	  } else if (brchIdx == cath_flag) {
	    Vsat -= (curden_cath_r*l_cath_r*G_RESIS_CU*E*G_Z/G_OMEGA)*(l_cath_r*w*h/2/G_B);
	  } else {
	    Vsat -= (curden_cath_r*l_cath_r + curden*l*2)*G_RESIS_CU*E*G_Z/G_OMEGA*(l*w*h/2/G_B);
	  }
	}

// debug
/*	if(treeIdx == 211) {
	  cout << brchIdx << ": " << brchit->first << " l: " << l*1e6 << " w: " << w*1e6 << " h: " << h*1e6 << " current density: " << curden << endl;
	  if(brchIdx == cath_flag) {
	    cout << brchit->second.left_node.name.c_str() << " " << cathode << endl;
	  }
	}  */
	brchIdx++;
      }
//	cout << l*1e6 << " " << w*1e6 << " " << h*1e6 << " " << G_RESIS_CU*1e6 << " " << E*1e6 << " " << G_Z << " " << G_OMEGA << " " << G_B << endl;
      brchIdx = 0;
      if(Vsat>Vcross) {
//      if(1) {
//      if(treeIdx == 42 || treeIdx == 52 || treeIdx == 62 || treeIdx == 72 || treeIdx == 77 || treeIdx == 87 || (treeIdx>90 && treeIdx<145 && (treeIdx%10==2 || treeIdx%10==7)) || (treeIdx>715 && treeIdx<810 && (treeIdx%10==2 || treeIdx%10==7)) || (treeIdx>881 && treeIdx<889) || (treeIdx>901 && treeIdx<909)) {
//      if(treeIdx == 14 || treeIdx == 30 || treeIdx == 32 || treeIdx == 52 || (treeIdx>53 && treeIdx<66)) {
	half_fail--;
	numFT++;
	fprintf(FTID, "F%d\r\n", treeIdx);
	fprintf(ftid0, "%d\r\n", treeIdx);
	for(brchit=brch.begin(); brchit!=brch.end(); ++brchit) {
	  rit = rmap.find(brchit->first);
// debug
//	vit = vmap.find(brchit->second.left_node.name);
//	cout << brchit->first << "  " << vit->first << " " << vit->second << " ";
//	vit = vmap.find(brchit->second.right_node.name);
//	cout << vit->first << " " << vit->second << endl;
//
          curden = rit->second.curden;
	  resInfo info = rit->second;
	  w = brchit->second.w;
	  l = brchit->second.l;
	  fprintf(tinfo, "%s n%d_%d_%d n%d_%d_%d %f\r\n", rit->first.c_str(), info.left_node.net, info.left_node.x, info.left_node.y, info.right_node.net, info.right_node.x, info.right_node.y, info.value);
	  fprintf(curd, "%s, curden = %f\r\n", rit->first.c_str(), curden);
	  if (wflag == 1) {
	    fprintf(width, "W%.10f\r\n", w);
	  }
	  fprintf(leng, "L%.10f\r\n", l);
	}
	fprintf(tinfo, "xx\r\n");
//	fprintf(curd, "xx\r\n");
//	fprintf(width, "xx\r\n");
      }
    }
//    cout << endl;
//    cout << "Tree " << treeIdx << " Ve: " << Ve << " Failed tree #: " << numFT << " Vsat: " << Vsat*1e18 << "um^3 " << "Vcross: " << Vcross*1e18 << "um^3" << " thickness: " << h*1e6 << "um" << endl;
    cout.precision(5);
    cout << "Tree " << treeIdx << " Ve: " << Ve << " Failed Tree #: " << numFT << " Vsat: " << Vsat*1e18 << " Vcross: " << Vcross*1e18 << endl;
//    cout << cath_flag << endl;
    treeIdx++;
    VA = 0;
    Vsat = 0;
    A_t = 0;
    cath_flag = 0;

  }  // tree iteration end

  fprintf(cath, "\r\n");

  fclose(tinfo);
  fclose(curd); 
  fclose(width);
  fclose(leng);
  fclose(FTID);
  fclose(ftid0);
  fclose(voltage);
  fclose(node);
  fclose(cath);

	cout << "The # of trees that exceed the critical voltage but immortal: " << half_fail << endl;  

}
