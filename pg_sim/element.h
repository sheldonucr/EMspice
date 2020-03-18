/*
 **********************************************************
 
				 Power Grid Simulator
		 (Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: element.h,v $
 *    Authors: Ning Mi and Han Zhou
 *    Functions: header file for node class and matrix class
 *    $Date: 2020/03/17 $
 *
 */


#ifndef __ELEMENT_H
#define __ELEMENT_H

#include <stdio.h>
#include <vector>
#include <string>
#include <map>
#include <string.h>

using namespace std;

const int GNDNODE = -1;

/**********************************
* class NodeList
*********************************/
class NodeList
{
public:
	struct Node
	{
		int row_no;
	};

private:
	vector<Node> nodeData;
	vector<int>  portData;
	vector<string> portName;
	vector<int> tcData;
	vector<string> tcName;
	map<string, int> nodeTable;
	int nodeSize;

public:
	NodeList(void);
	~NodeList(void);
	
	int findorPushNode_map(const char* name);
	int findorPushNode(const char* name);
	int findNode2(const char* name);
	void pushPort(const char* name);
	void pushTCNode(const char* name);
	void setRowNum(int i, int index); //set the row_no of ith 
	void incRowNum(int i, int inc_index); // increase by inc_index

	inline Node* getNode(int i)        //access the ith entry
	{  return &(nodeData[i]);  } 
	inline Node* getPort(int i)        //access the ith entry
	{  return &(nodeData[portData[i]]);  } 
	inline string getPortName(int i)        //access the ith entry
	{  return portName[i];  } 
	inline Node* getTCNode(int i)        //access the ith entry
	{  return &(nodeData[tcData[i]]);  } 
	inline string getTCName(int i)        //access the ith entry
	{  return tcName[i];  } 
	inline int getTCNum() {return tcData.size();}
	int numNode(void){  return nodeSize; }
	int numPort(void){  return portData.size(); }
};


/**************************************
 * class matrix
 **************************************/
#include <vector>
#include "cs.h"

#define ROW_NUM 5

using namespace std;

struct Entry
{
    int i;
    double value;
};

class matrix
{
public:
    
private:
    Entry* compact_data;
    Entry* uncompact_data;
    
public:
    int colsize,rowsize;
    int nnz;
    
    Entry** colIndex;
    int* num_per_col;
    
    matrix(int m, int n);
    ~matrix();
    
    void pushEntry(int i, int j, double value);
    
    void sort();
    
    cs_dl* mat2csdl();
};

bool entry_comp (Entry a, Entry b);

#endif
