/*
 **********************************************************
 	
				Power Grid Simulator
 		(Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: element.cpp,v $
 *    Authors: Ning Mi, Xin Huang, and Han Zhou
 *    Functions: node class and matrix class
 *    $Date: 2020/03/17 $
 *
 */


#include "element.h"
#include "UFconfig.h"
#include <iostream>
#include <stdlib.h>
#include <algorithm>

NodeList::NodeList(){
	nodeSize = 0;
}

NodeList::~NodeList(){
  nodeData.clear();
  portData.clear();
  nodeTable.clear();
}

int NodeList::findNode2(const char* name){
  int addr;
  string sname = name;
  map<string, int>::iterator it;

  it = nodeTable.find(sname);
  
  if(it != nodeTable.end()){ // name in the nodelist
    return it->second;
  }
  else
    return -1;
}

int NodeList::findorPushNode(const char* name){
	int addr;
	string sname = name;
	map<string, int>::iterator it;
	
	it = nodeTable.find(sname);

	if (it != nodeTable.end()){
		addr = it->second;
		return addr;
	}
	else{
		addr = nodeData.size();
		Node tempnode;
		if (!strcmp(name,"0") || !strcmp(name,"gnd"))
			tempnode.row_no = GNDNODE;
		else 
			tempnode.row_no = nodeSize++;

		nodeData.push_back(tempnode);
		nodeTable[sname] = addr;
	}
	return addr;

}		 

void NodeList::setRowNum(int i, int index){
  nodeData[i].row_no = index;
}

void NodeList::incRowNum(int i, int inc_index){
  nodeData[i].row_no = nodeData[i].row_no + inc_index;
}

void NodeList::pushPort(const char* name)
{
  int i = findNode2(name);
  if(nodeData[i].row_no >= 0){
    portData.push_back(i);
    string sname(name);
    portName.push_back(sname);
  }else if(nodeData[i].row_no == GNDNODE){
	printf("Print GND node. \n");
  }
  else
    printf("Print port node %s does not exist. \n", name);
}

void NodeList::pushTCNode(const char* name)
{
  int i = findNode2(name);
  if(nodeData[i].row_no >= 0){
    tcData.push_back(i);
    string sname(name);
    tcName.push_back(sname);
  }else if(nodeData[i].row_no == GNDNODE){
  }
  else
    printf("tap current node %s does not exist. \n", name);
}

bool entry_comp(Entry a, Entry b) {
    return a.i < b.i;
}

matrix::matrix(int m, int n) {
    int i;
    
    rowsize = m;
    colsize = n;
    nnz = 0;
    colIndex = (Entry**) malloc(n * sizeof(Entry*));
    num_per_col = (int*) malloc(n * sizeof(int));
    for (i = 0; i < colsize; i++) {
        num_per_col[i] = 0;
        colIndex[i] = NULL;
    }
    uncompact_data = NULL;
    compact_data = NULL;
}

matrix::~matrix() {
    for (int i = 0; i < colsize; i++) {
        if (colIndex[i] != NULL)
            free(colIndex[i]);
    }
    
    free(colIndex);
    free(num_per_col);
    if (compact_data != NULL)
        free(compact_data);
    if (uncompact_data != NULL)
        delete uncompact_data;
}

void matrix::pushEntry(int i, int j, double value) {
    int k, num;
    int mark = 0;
    int nblock;
    if (num_per_col[j] == 0) {
        colIndex[j] = (Entry*) malloc(ROW_NUM * sizeof(Entry));
        colIndex[j][0].i = i;
        colIndex[j][0].value = value;
        num_per_col[j]++;
        nnz++;
    } else {
        for (k = 0; k < num_per_col[j]; k++) {
            if (colIndex[j][k].i == i) {
                colIndex[j][k].value += value;
                mark = 1;
                break;
            }
        }
        if (mark != 1) {
            num = num_per_col[j];
            if (num % ROW_NUM == 0) {
                nblock = num / ROW_NUM + 1;
                colIndex[j] = (Entry*) realloc(colIndex[j], nblock * ROW_NUM
                                               * sizeof(Entry));
            }
            colIndex[j][num].i = i;
            colIndex[j][num].value = value;
            num_per_col[j]++;
            nnz++;
        }
    }
}

void matrix::sort() {
    int j;
    for (j = 0; j < colsize; j++) {
        colIndex[j] = (Entry*) realloc(colIndex[j], num_per_col[j]
                                       * sizeof(Entry));
        std::sort(colIndex[j], colIndex[j] + num_per_col[j], entry_comp);
    }
}

cs_dl* matrix::mat2csdl() {
    if (nnz == 0)
        return NULL;
    cs_dl *T = cs_dl_spalloc((UF_long) rowsize, (UF_long) colsize,
                             (UF_long) nnz, 1, 0);
    //Entry* temp_entry;
    UF_long i, j, nnz_tmp;
    
    nnz_tmp = 0;
    for (i = 0; i < colsize; i++) {
        T->p[i] = nnz_tmp;
        for (j = 0; j < num_per_col[i]; j++) {
            T->i[nnz_tmp] = colIndex[i][j].i;
            T->x[nnz_tmp] = colIndex[i][j].value;
            nnz_tmp++;
        }
    }
    T->p[colsize] = nnz;
    return T;
}
