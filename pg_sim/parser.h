/*
 **********************************************************
 
 				Power Grid Simulator
		 (Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: parser.cpp,v $
 *    Authors: Ning Mi and Han Zhou
 *    Functions: header file for parser
 *    $Date: 2020/03/17 $
 *
 */


#ifndef __PARSER_H
#define __PARSER_H

#include "element.h"
#include "mna.h"
#include <itpp/base/vec.h>
#include <itpp/base/smat.h>
#include <itpp/base/mat.h>
#include <itpp/base/math/elem_math.h>

// XXLiu: change READ_BLOCK_SIZE to 1000. Before, it was 400, set by Duo.
#define READ_BLOCK_SIZE 300000000 
#define NAME_BLOCK_SIZE 2000
#define VALUE_BLOCK_SIZE 400

using namespace itpp;

typedef struct{
    char type;
    char* node1;
    char* node2;
    double value;
} subelement;

void parser_sub(const char* filename, double& tstep, double& tstop, int& num_subnode, int& nVS, int& nIS, int& nL, NodeList* nodePool);

void stamp_sub(const char* filename, int nL, int nIS, int nVS, int& nNodes,
	       double& tstep, double& tstop, Source *VS, Source *IS,
	       matrix* G, matrix* C, matrix* B, NodeList* nodePool);

void stampG(const char* filename, int nL, int nVS, int nNodes, matrix* G, NodeList* nodePool);

void stampC(const char* filename, int nL, int nVS, int nNodes, matrix* C, NodeList* nodePool);

void stampB(const char* filename, int nL, int nIS, int nVS, int nNodes, double tstop, Source *VS, Source *IS, matrix* B, NodeList* nodePool);

#endif
