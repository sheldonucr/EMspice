/*
 **********************************************************
 
				 Power Grid Simulator
		 (Netlist Parser and Nodal Voltage Solver)
 
 **********************************************************
 */

/*
 *    $RCSfile: parser.cpp,v $
 *    Authors: Ning Mi and Han Zhou
 *    Functions: parser--read circuit information from spice file
 *    $Date: 2020/03/17 $
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "parser.h"
#include "parameter.h"

double StrToNum(char* strnum)
{
  double value;
  char* endptr;

  value = strtod(strnum, &endptr);

  switch(endptr[0])
    {
    case 'T': case 't':
      return value*pow(10.0,12);
    case 'G': case 'g':
      return value*pow(10.0, 9);
    case 'K': case 'k':
      return value*pow(10.0, 3);
    case 'M': case 'm':
      switch(endptr[1])
	{
	case 'E': case 'e':
	  return (value*pow(10.0, 6));
	default:
	  return (value*pow(10.0, -3));
	}
    case 'U': case 'u':
      return value*pow(10.0,(int)-6);
    case 'n': case 'N':
      return value*pow(10.0,(int)-9);

    case 'p': case 'P':
      return value*pow(10.0,(int)-12);
    case 'f': case 'F':
      return value*pow(10.0,(int)-15);
    case 'e': case 'E':
      return value;
    default:
      return value;
    }
}

void parser_sub(const char* filename, double& tstep, double& tstop, int& num_subnode, int& nVS, int& nIS, int& nL, NodeList* nodePool)
{
  char* strline;
  char* node1, *node2;
  char *tstepstr, *tstopstr;
  char *portstr;
  char *incfilename, *tmp_incfilename;
  char *subcktportstr, *subcktportname;

  int num_R,num_C,num_L,num_V,num_I,nport;
  int n1, n2;
  int pos;
  int spacemark;

  int subportnum =0;
  int subcktportnum = 0;	
  int i,j,p;

  bool subckt=false, mainckt=false, subctn=false;
  NodeList *subnode;

  FILE* fid_tmp = NULL;

  FILE* fid = fopen(filename, "r");

  if(fid == NULL){
    printf("Open file Error!\n");
    exit(-1);
  }

  strline = (char*)malloc(READ_BLOCK_SIZE*sizeof(char));
  incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  tmp_incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node1 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node2 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  portstr = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  subcktportstr = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  subcktportname = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));

  tstepstr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  tstopstr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));

  // get the director of file
  i = 0;
  pos = 0;
  while(filename[i]!='\0'){
    if(filename[i] == '/')
      pos = i;
    i++;
  }

  if(pos != 0){
    pos++;
    strncpy(incfilename,filename,pos);
    incfilename[pos] = '\0';
  }

  //get number of R,C,L,I,V, construct nodelist
  num_R = 0;
  num_C = 0;
  num_L = 0;
  num_V = 0;
  num_I = 0;
 

  while(!feof(fid) || (fid_tmp!=NULL && !feof(fid_tmp))){
    if(feof(fid)==0){
      fgets(strline, READ_BLOCK_SIZE,fid);
    }
    else{
      fgets(strline, READ_BLOCK_SIZE, fid_tmp);
    }

    switch(strline[0])
      {
      case'R': case'r':
	num_R++;
	if(subckt == false){
	  if(sscanf(strline, "%*s %s %s %*s", node1, node2) == 2){
	    n1 = nodePool->findorPushNode(node1);
	    n2 = nodePool->findorPushNode(node2);
	  }
	}
	else{
	  if(subctn == true) subctn = false;
	  if(sscanf(strline, "%*s %s %s", node1, node2) == 2){
	    n1 = subnode->findorPushNode(node1);
	    n2 = subnode->findorPushNode(node2);
	  }
	}
	//printf("In R branch with %s and %s, total node number is %d.\n", node1, node2, nodePool->numNode());
	break;
	
      case'C': case'c':
	num_C++;
	if(subckt == false){
	  if(sscanf(strline, "%*s %s %s %*s", node1, node2) == 2){
	    n1 = nodePool->findorPushNode(node1);
	    n2 = nodePool->findorPushNode(node2);
	  }
	}
	else{
	  if(subctn == true) subctn = false;
	  if(sscanf(strline, "%*s %s %s", node1, node2) == 2){
	    n1 = subnode->findorPushNode(node1);
	    n2 = subnode->findorPushNode(node2);
	  }
	}
	//printf("In C branch, total node number is %d.\n", nodePool->numNode());
	break;
      
	  case'V': case'v':
	num_V++;
	if(subckt == false){
	  if(sscanf(strline, "%*s %s %s %*s", node1, node2) == 2){
	    n1 = nodePool->findorPushNode(node1);
	    n2 = nodePool->findorPushNode(node2);
	  }
	}
	else{
	  if(subctn == true) subctn = false;
	  if(sscanf(strline, "%*s %s %s", node1, node2) == 2){
	    n1 = subnode->findorPushNode(node1);
	    n2 = subnode->findorPushNode(node2);
	  }
	}
	//printf("In V branch, total node number is %d.\n", nodePool->numNode());
	break;

      case'I': case'i':
	num_I++;
	if(subckt == false){
	  if(sscanf(strline, "%*s %s %s %*s", node1, node2) == 2){
	    n1 = nodePool->findorPushNode(node1);
		nodePool->pushTCNode(node1);
	    n2 = nodePool->findorPushNode(node2);
		nodePool->pushTCNode(node2);
	  }
	}
	else{
	  if(subctn == true) subctn = false;
	  if(sscanf(strline, "%*s %s %s", node1, node2) == 2){
	    n1 = subnode->findorPushNode(node1);
		nodePool->pushTCNode(node1);
	    n2 = subnode->findorPushNode(node2);
		nodePool->pushTCNode(node2);
	  }
	}
	//printf("In I branch, total node number is %d.\n", nodePool->numNode());
	break;
	
      case'X': case'x':
	mainckt = true;
	subcktportnum = 0;
	if(sscanf(strline, "%*s %s", subcktportstr) == 1){
	  i = 0;
	  while(strline[i]!= ' ')
	    i++;
	  while(strline[i]!= '\0'){
	    j = 0;
	    while(strline[i] != ' ' && strline[i]!='\n') //copy a port name
	      subcktportname[j++] = strline[i++]; 
	    subcktportname[j] = '\0';
	    if(j != 0){
	      n1 = nodePool->findorPushNode(subcktportname);
              subcktportnum++;
	    }
	    i++;
	  }
	  subcktportname[0] = '\0';
	}
	break;
      
	  case'.':
	if(strncmp(strline, ".tran", 5)==0 || strncmp(strline, ".TRAN", 5)==0){//get simulation start and stop time
	  if(sscanf(strline, "%*s %s %s", tstepstr, tstopstr)==2){
	    tstep = StrToNum(tstepstr);
	    tstop = StrToNum(tstopstr);
	  }
	}
	else if (strncmp(strline, ".print", 6)==0 ){ //get port name
	  i = 0;
	  if (strline[12] != 'v' && strline[12] != 'V' && strline[12] != 'i' && strline[12] != 'I')
	    printf("Invalid command: %s\n",strline);
	  else{
	    while(strline[i]!='\0'){
	      if(strline[i]=='('){
		i++;
		j = 0;
		while(strline[i]!=')'){
		  portstr[j++] = strline[i++];
		}
		portstr[j] = '\0';
		nodePool->pushPort(portstr);
	      }
	      i++;
	    }
	  }
	}
	else if (strncmp(strline, ".include", 8)==0 || strncmp(strline, ".INCLUDE", 8)==0){ //.include
	  if(sscanf(strline, "%*s %s", tmp_incfilename) == 1){
	    i = 0;
	    j = pos;
	    while(tmp_incfilename[i] != '\0'){
	      if(tmp_incfilename[i]!='\"'){
		incfilename[j] = tmp_incfilename[i];
		j++;
	      }
	      i++;
	    }
	    incfilename[j] = '\0';

	    fid_tmp = fid;

	    fid = fopen(incfilename,"r");
	    if(fid == NULL){
	      printf("Open file Error!\n");
	      exit(-1);
	    }	    
	  }
	}
	else if(strncmp(strline, ".SUBCKT", 7)==0 || strncmp(strline, ".subckt", 7)==0)
	  {
	    subckt = true;
	    subctn = true;
	    mainckt = false;
	    subportnum = 0;
	    if((subnode = new NodeList) == NULL){
	      printf("Out of memory!\n");exit(1);
	    }
	    i = 0;
	    spacemark = 0;
	    while(spacemark < 2 && strline[i]!='\n'){
	      if(strline[i] == ' '){
		spacemark++;
	      }
	      i++;
	    }
	    p = 0;
	    while(strline[i]!='\0'){
	      while(strline[i]!=' ' && strline[i]!='\n'){
		subcktportname[p++] = strline[i++];
	      }
	      if(p != 0){
		subcktportname[p] = '\0';
		subnode->findorPushNode(subcktportname);
		subportnum++;
		p = 0;
	      }
	      i++;
	    }
	  }
	else if((subckt == true) && (strncmp(strline, ".ends", 5) == 0 || strncmp(strline, ".ENDS", 5) == 0))
	  {
	    subckt = false;
	    num_subnode = num_subnode + subnode->numNode()-subportnum;
	    delete subnode;
	    subportnum = 0;
	  }
	break;

      default:
	break;
      }
  }
  fclose(fid);

  if(fid_tmp != NULL){
    fclose(fid_tmp);
  }

  nIS = num_I;
  nVS = num_V;
  nL = num_L;


  free(strline);   
  free(tmp_incfilename);
  free(incfilename);
  free(node1);
  free(node2);
  free(portstr);
  free(tstepstr);
  free(tstopstr);
  free(subcktportstr);
  free(subcktportname);

}


void stamp_sub(const char* filename, int nL, int nIS, int nVS, int& nNodes,
	       double& tstep, double& tstop, Source *VS, Source *IS,
	       matrix* G, matrix* C, matrix* B, NodeList* nodePool)
{
  char* strline;
  char* node1, *node2;
  char *timestr, *stimestr, *valuestr;
  char *incfilename, *tmp_incfilename;
  char *istr, *istr_tmp;
  char *subcktportname, *subcktportstr;

  double ptime, value, Rvalue;

  int num_point;
  int n1, n2, n1_tmp, n2_tmp;
  int num_L_tmp, num_I_tmp, num_V_tmp;
  int index_i, index_j;
  int i,j,k,p;
  int tmp_portnum,  subportnum, subnodenum;
  int mark_value;
  int pos;
  int num_mainckt, num_subckt;
  int inc_index;
  int spacemark;
  int substart, subend;
  int last_npoint_I, last_npoint_V;

  bool isI, isV;
  bool subckt, mainckt, subctn=false;
  bool subgnd, subportgnd;

  NodeList* subNode;
  vector<subelement> subData;

  FILE* fid_tmp=NULL;

  FILE* fid = fopen(filename, "r");
  if(fid == NULL){
    printf("Open file Error!\n");
    exit(-1);
  }

  strline = (char*)malloc(READ_BLOCK_SIZE*sizeof(char));
  incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  tmp_incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node1 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node2 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  istr = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  istr_tmp = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  subcktportname = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  subcktportstr = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  timestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  stimestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  valuestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));

  // get the director of file
  i = 0;
  pos = 0;
  while(filename[i]!='\0'){
    if(filename[i] == '/')
      pos = i;
    i++;	  
  }

  if(pos != 0){
    pos++;
    strncpy(incfilename,filename,pos);
    incfilename[pos] = '\0';
  }


  /* stamp circuit to G,C,B and waveform */   
  num_mainckt = nodePool->numNode(); // node number in main ckt
  num_subckt = 0;    //node number in subckt
  num_L_tmp = 0;
  num_I_tmp = -1;
  num_V_tmp = -1;
  while(!feof(fid) || (fid_tmp!=NULL && !feof(fid_tmp))){
    if(feof(fid)==0){
      fgets(strline, READ_BLOCK_SIZE,fid);
    }
    else{
      if (fgets(strline, READ_BLOCK_SIZE, fid_tmp) == NULL) break;
    }

    switch(strline[0])
      {
	//  ************ Resistor ****************
      case'R': case'r':
	if(subckt == false){
	  if(sscanf(strline, "%*s %s %s %s", node1, node2, valuestr) == 3){
	    Rvalue = StrToNum(valuestr);
	    value = 1.0/Rvalue;
	    n1_tmp = nodePool->findorPushNode(node1);
	    n2_tmp = nodePool->findorPushNode(node2);
	    n1 = nodePool->getNode(n1_tmp)->row_no;
	    n2 = nodePool->getNode(n2_tmp)->row_no;
	    if (n1 != GNDNODE) G->pushEntry(n1, n1, value);
	    if (n2 != GNDNODE) G->pushEntry(n2, n2, value);
	    if (n1 != GNDNODE && n2 != GNDNODE)
	      {
		G->pushEntry(n1,n2,-value);
		G->pushEntry(n2,n1,-value);
	      }
	  }
	  else
	    printf("Fail in obtaining resistence value.\n");
	}
	else{
	  if(subctn == true)
	    subctn = false;
	  if(sscanf(strline, "%*s %s %s %s", node1, node2, valuestr) == 3){
	    Rvalue = StrToNum(valuestr);
	    value = 1/Rvalue;
	    subNode->findorPushNode(node1);
	    subNode->findorPushNode(node2);
	    if((subportgnd == false) && (!strcmp(node1, "0") || !strcmp(node1, "gnd") || !strcmp(node2, "0") || !strcmp(node2, "gnd")) )
	      subgnd = true;
	    subelement subele;
	    subele.type = 'R';
	    subele.node1 = (char*)malloc(strlen(node1)*sizeof(char));
	    subele.node2 = (char*)malloc(strlen(node2)*sizeof(char));
	    strcpy(subele.node1, node1);
	    strcpy(subele.node2, node2);
	    subele.value = value;
	    subData.push_back(subele);
	  }
	}
	break;

	//  ************ Capacitor ***************
      case'C': case'c':
	if(subckt == false){
	  if(sscanf(strline, "%*s %s %s %s", node1, node2, valuestr) == 3){
	    value = StrToNum(valuestr);
	    n1_tmp = nodePool->findorPushNode(node1);
	    n2_tmp = nodePool->findorPushNode(node2);
	    n1 = nodePool->getNode(n1_tmp)->row_no;
	    n2 = nodePool->getNode(n2_tmp)->row_no;
	    if (n1 != GNDNODE) C->pushEntry(n1, n1, value);
	    if (n2 != GNDNODE) C->pushEntry(n2, n2, value);
	    if (n1 != GNDNODE && n2 != GNDNODE)
	      {
		C->pushEntry(n1,n2,-value);
		C->pushEntry(n2,n1,-value);
	      }
	  }
	  else
	    printf("Fail in obtaining capacitor value.\n");
	}
	else{
	  if(subctn == true)
	    subctn = false;
	  if(sscanf(strline, "%*s %s %s %s", node1, node2, valuestr) == 3){
	    value = StrToNum(valuestr);
	    subNode->findorPushNode(node1);
	    subNode->findorPushNode(node2);
	    if(subportgnd == false && (!strcmp(node1, "0") || !strcmp(node1, "gnd") || !strcmp(node2, "0") || !strcmp(node2, "gnd")) )
	      subgnd = true;
	    subelement subele;
	    subele.type = 'C';
	    subele.node1 = (char*)malloc(strlen(node1)*sizeof(char));
	    subele.node2 = (char*)malloc(strlen(node2)*sizeof(char));
	    strcpy(subele.node1, node1);
	    strcpy(subele.node2, node2);
	    subele.value = value;
	    subData.push_back(subele);
	  }
	}
	break;

	// ************ Voltage *******************
      case'V': case'v':
	if(subckt == false){
	  num_V_tmp++;
	  isV = true; isI = false;
	  if(sscanf(strline, "%*s %s %s %s",  node1, node2, istr) == 3){
	    n1_tmp = nodePool->findorPushNode(node1);
	    n2_tmp = nodePool->findorPushNode(node2);
	    n1 = nodePool->getNode(n1_tmp)->row_no;
	    n2 = nodePool->getNode(n2_tmp)->row_no;
	    
	    index_i = nNodes + nL + num_V_tmp;
	    index_j = num_V_tmp;
	    
	    if(n1 != GNDNODE) {
	      G->pushEntry(n1, index_i, 1);
	      G->pushEntry(index_i, n1, -1);
	    }
	    if(n2 != GNDNODE) {
	      G->pushEntry(n2, index_i, -1);
	      G->pushEntry(index_i, n2, 1);
	    }
	    B->pushEntry(index_i, index_j, -1);
	    
	    //DC
	    value = StrToNum(istr);
	      
	    VS[num_V_tmp].time.set_size(2, false);
	    VS[num_V_tmp].value.set_size(2, false);
	      
	    VS[num_V_tmp].time(0) = 0;
	    VS[num_V_tmp].time(1) = tstop;
	    VS[num_V_tmp].value(0) = value;
	    VS[num_V_tmp].value(1) = value;
	  }
	}
	break;

	// ******************** Current source ******************
      case'I': case'i':
	if(subckt == false){
	  num_I_tmp++;
	  isV = false; isI = true;
	  if(sscanf(strline, "%*s %s %s %s", node1, node2, istr) == 3){
	    n1_tmp = nodePool->findorPushNode(node1);
	    n2_tmp = nodePool->findorPushNode(node2);
	    n1 = nodePool->getNode(n1_tmp)->row_no;
	    n2 = nodePool->getNode(n2_tmp)->row_no;
	    
	    index_j = nVS + num_I_tmp; 
	    
	    if(n1 != GNDNODE) {
	      B->pushEntry(n1,index_j,-1);
	    }
	    if(n2 != GNDNODE) {
	      B->pushEntry(n2,index_j,1);
	    }
	    
	    //DC
	      value = StrToNum(istr)*I_SCALE;//special
	      IS[num_I_tmp].time.set_size(2, false);
	      IS[num_I_tmp].value.set_size(2, false);
	      
	      IS[num_I_tmp].time(0) = 0;
	      IS[num_I_tmp].time(1) = tstop;
	      IS[num_I_tmp].value(0) = value;
	      IS[num_I_tmp].value(1) = value;
	  }
	}
	break;

	// ************ subckt used in mainckt
      case'X':
	mainckt = true;
	tmp_portnum = 0;
	if(sscanf(strline, "%*s %s", subcktportstr) == 1){
	  i = 0;
	  while(strline[i]!= ' ')
	    i++;
	  while(strline[i]!= '\0'){
	    p = 0;
	    while(strline[i] != ' ' && strline[i]!= '\n') //copy a port name
	      subcktportname[p++] = strline[i++]; 
	    subcktportname[p] = '\0';
	    if(p != 0){
	      if(tmp_portnum < subportnum){
		if(subcktportname[0] == 'P'){ // end of subckt instance
		  mainckt = false;
		  tmp_portnum++;
		  break;
		}
		n1_tmp = nodePool->findorPushNode(subcktportname);
		n1 = nodePool->getNode(n1_tmp)->row_no;
		subNode->setRowNum(tmp_portnum, n1);
		tmp_portnum++;
	      }
	      else if(tmp_portnum == subportnum){
		inc_index = num_mainckt + num_subckt;
		if(subportgnd == true){//gnd in the subckt port
		  substart = tmp_portnum+1;
		  subend = subNode->numNode()+1;
		}
		else if(subgnd == true){//gnd inside subckt
		  substart = tmp_portnum;
		  subend = subNode->numNode()+1;
		}
		else{//no gnd in subckt or in subckt port
		  substart = tmp_portnum;
		  subend = subNode->numNode();
		}
		for(j = substart; j<subend; j++){
		  if(subNode->getNode(j)->row_no != GNDNODE)
		    subNode->incRowNum(j,inc_index);
		}
		for(j = 0; j<subData.size(); j++){
		  n1_tmp = subNode->findorPushNode(subData[j].node1);
		  n2_tmp = subNode->findorPushNode(subData[j].node2);
		  n1 = subNode->getNode(n1_tmp)->row_no;
		  n2 = subNode->getNode(n2_tmp)->row_no;
		  if(subData[j].type == 'R'){
		    if(n1 != GNDNODE) G->pushEntry(n1, n1, value);
		    if (n2 != GNDNODE) G->pushEntry(n2, n2, value);
		    if (n1 != GNDNODE && n2 != GNDNODE)
		      {
			G->pushEntry(n1,n2,-value);
			G->pushEntry(n2,n1,-value); 
		      }
		  }
		  else if(subData[j].type == 'C'){
		    if (n1 != GNDNODE) C->pushEntry(n1, n1, value);
		    if (n2 != GNDNODE) C->pushEntry(n2, n2, value);
		    if (n1 != GNDNODE && n2 != GNDNODE)
		      {
			C->pushEntry(n1,n2,value);
			C->pushEntry(n2,n1,value);
		      }
		  }
		}
		num_subckt = num_subckt + subNode->numNode() - subportnum;
		delete subNode;
		for (j = 0; j<subData.size(); j++){
		  free(subData[j].node1);
		  free(subData[j].node2);
		}
		subData.clear();
		if(subcktportname[0]=='P')
		  mainckt = false;
		if(subgnd == true)
		  subgnd = false;
		if(subportgnd == true)
		  subportgnd = false;
		tmp_portnum++;
		//break;
	      }
	      else{
		if(subcktportname[0]=='P')
		  mainckt = false;

	      }
	    }
	    i++;
	  }
	}
	break;

        //  ********** .include ********************
      case'.':
	if (strline[1]=='i' && strline[2]=='n'){ //.include
	  if(sscanf(strline, "%*s %s", tmp_incfilename) == 1){
	    
	    i = 0;
	    j = pos;
	    while(tmp_incfilename[i] != '\0'){
	      if(tmp_incfilename[i]!='\"'){
		incfilename[j] = tmp_incfilename[i];
		j++;
	      }
	      i++;
	    }
	    incfilename[j] = '\0';


	    fid_tmp = fid;

	    fid = fopen(incfilename,"r");
	    if(fid == NULL){
	      printf("Open file Error!\n");
	      exit(-1);
	    }	    
	  }
	}
	else if(strncmp(strline, ".SUBCKT", 7)==0 || strncmp(strline, ".subckt", 7)==0)
	  {
	    subckt = true;
	    subportnum = 0;
	    //mainckt = false;
	    subctn = true; //mark +
	    if((subNode = new NodeList) == NULL){
	      printf("Out of memory!\n");exit(1);
	    }
	    i = 0;
	    spacemark = 0;
	    while(spacemark < 2 && strline[i]!='\0'){
	      if(strline[i] == ' '){
		spacemark++;
	      }
	      i++;
	    }
	    p = 0;
	    while(strline[i]!='\0'){
	      while(strline[i]!=' ' && strline[i]!='\n'){
		subcktportname[p++] = strline[i++];
	      }
	      subcktportname[p] = '\0';
	      if(p != 0){
		subNode->findorPushNode(subcktportname);
		subportnum++;
		p = 0;
		if(!strcmp(subcktportname,"0")||!strcmp(subcktportname,"gnd"))
		  subportgnd = true;
	      }
	      i++;
	    }
	  }

	else if((subckt == true) && (strncmp(strline, ".ends", 5) == 0 || strncmp(strline, ".ENDS", 5) == 0))
	  {
	    subckt = false;
	  }
	else if(strncmp(strline, ".end", 4) == 0 || strncmp(strline, ".END", 4) == 0)
	  {
	    if(last_npoint_I != 0){
	      IS[num_I_tmp].time.set_size(last_npoint_I,true);
	      IS[num_I_tmp].value.set_size(last_npoint_I,true);
	    }
	    if(last_npoint_V != 0){
	      VS[num_V_tmp].time.set_size(last_npoint_V,true);
	      VS[num_V_tmp].value.set_size(last_npoint_V,true);
	    }
	  }
	break;

      default:
	break;
      }
  }
  fclose(fid);

  if(fid_tmp != NULL){
    fclose(fid_tmp);
  }

  G->sort();
  C->sort();
  B->sort();

  free(strline);   
  free(tmp_incfilename);
  free(incfilename);
  free(node1);
  free(node2);
  free(istr);
  free(istr_tmp);
  free(timestr);
  free(stimestr);
  free(valuestr);
  free(subcktportstr);
  free(subcktportname);

}


void stampG(const char* filename, int nL, int nVS, int nNodes, matrix* G, NodeList* nodePool)
{
  char* strline;
  char* node1, *node2;
  char *timestr, *stimestr, *valuestr;
  char *incfilename, *tmp_incfilename;

  double ptime, value, Rvalue;

  int num_point;
  int n1, n2, n1_tmp, n2_tmp;
  int num_L_tmp, num_I_tmp, num_V_tmp;
  int index_i, index_j;
  int i,j;
  int pos;

  bool isI, isV;
//SHyu
//  char noderec[] = "123";
//SHyu end

  //NodeList::Node* node_tmp;

  FILE* fid_tmp=NULL;

  FILE* fid = fopen(filename, "r");
  if(fid == NULL){
    printf("Open file Error!\n");
    exit(-1);
  }

  strline = (char*)malloc(READ_BLOCK_SIZE*sizeof(char));
  incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  tmp_incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node1 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node2 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  timestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  stimestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  valuestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));

  /* get the director of file */
  i = 0;
  pos = 0;
  while(filename[i]!='\0'){
    if(filename[i] == '/')
      pos = i;
    i++;
  }

  if(pos != 0){
    pos++;
    strncpy(incfilename,filename,pos);
    incfilename[pos] = '\0';
  }

  /* stamp circuit to G,C,B and waveform */   
  num_L_tmp = 0;
  num_I_tmp = -1;
  num_V_tmp = -1;

  while(!feof(fid)||(fid_tmp != NULL && !feof(fid_tmp))){

    if(feof(fid)==0){
      fgets(strline, READ_BLOCK_SIZE,fid);
    }
    else{
      if (fgets(strline, READ_BLOCK_SIZE, fid_tmp) == NULL) break;
    }

    switch(strline[0])
      {
	//  ************ Resistor ****************
      case'R': case'r':
	if(sscanf(strline, "%*s %s %s %s", node1, node2, valuestr) == 3){
	  Rvalue = StrToNum(valuestr);
	  value = 1.0/Rvalue;
	  n1_tmp = nodePool->findorPushNode(node1);
	  n2_tmp = nodePool->findorPushNode(node2);
	  n1 = nodePool->getNode(n1_tmp)->row_no;
	  n2 = nodePool->getNode(n2_tmp)->row_no;
//SHyu output brch info
//
//	FILE * tree0;
//	tree0 = fopen("/home/eegrad/shyu/dac19_EMspice_workrepo/tree_info.txt","a+");

//	if(strcmp(noderec, node1)!= 0)
//	{
//	  fprintf(tree0, "\nxx");
//	}

//	fprintf(tree0, "\n%s %s %s", node1, node2, valuestr);
//need to output branchname
//	fclose(tree0);
//
//	strcpy(noderec, node2);
//
//SHYU end
	  if (n1 != GNDNODE) G->pushEntry(n1, n1, value);
	  if (n2 != GNDNODE) G->pushEntry(n2, n2, value);
	  if (n1 != GNDNODE && n2 != GNDNODE)
	    {
	      G->pushEntry(n1,n2,-value);
	      G->pushEntry(n2,n1,-value);
	    }
	}
	else{
	  printf("%s \n", strline);
	  printf("Fail in obtaining resistence value.\n");
	}
	break;

	// ************ Voltage *******************
      case'V': case'v':
	num_V_tmp++;
	isV = true; isI = false;
	if(sscanf(strline, "%*s %s %s",  node1, node2) == 2){
	  n1_tmp = nodePool->findorPushNode(node1);
	  n2_tmp = nodePool->findorPushNode(node2);
	  n1 = nodePool->getNode(n1_tmp)->row_no;
	  n2 = nodePool->getNode(n2_tmp)->row_no;
	  
	  index_i = nNodes + nL + num_V_tmp;
	  index_j = num_V_tmp;

	  if(n1 != GNDNODE) {
	    G->pushEntry(n1, index_i, 1);
	    G->pushEntry(index_i, n1, -1);
	  }
	  if(n2 != GNDNODE) {
	    G->pushEntry(n2, index_i, -1);
	    G->pushEntry(index_i, n2, 1);
	  }
	}
	break;

        //  ********** .include ********************
      case'.':
	if (strline[1]=='i' && strline[2]=='n'){ //.include
	  if(sscanf(strline, "%*s %s", tmp_incfilename) == 1){
	    i = 0;
	    j = pos;
	    while(tmp_incfilename[i] != '\0'){
	      if(tmp_incfilename[i]!='\"'){
		incfilename[j] = tmp_incfilename[i];
		j++;
	      }
	      i++;
	    }
	    incfilename[j] = '\0';


	    fid_tmp = fid;

	    fid = fopen(incfilename,"r");
	    if(fid == NULL){
	      printf("Open file Error!\n");
	      exit(-1);
	    }	    
	  }
	}
	break;

      default:
	break;
      }
  }
  fclose(fid);

  if(fid_tmp != NULL){
    fclose(fid_tmp);
  }

  G->sort();

  free(strline);   
  free(tmp_incfilename);
  free(incfilename);
  free(node1);
  free(node2);
  free(timestr);
  free(stimestr);
  free(valuestr);

}

void stampC(const char* filename, int nL, int nVS, int nNodes, matrix* C, NodeList* nodePool)
{
  char* strline;
  char* node1, *node2;
  char *timestr, *stimestr, *valuestr;
  char *incfilename, *tmp_incfilename;

  double ptime, value, Rvalue;

  int num_point;
  int n1, n2, n1_tmp, n2_tmp;
  int num_L_tmp, num_I_tmp, num_V_tmp;
  int index_i, index_j;
  int i,j;
  int pos;

  bool isI, isV;

  FILE* fid_tmp=NULL;

  FILE* fid = fopen(filename, "r");
  if(fid == NULL){
    printf("Open file Error!\n");
    exit(-1);
  }

  strline = (char*)malloc(READ_BLOCK_SIZE*sizeof(char));
  incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  tmp_incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node1 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node2 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  timestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  stimestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  valuestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));

  /* get the director of file */
  i = 0;
  pos = 0;
  while(filename[i]!='\0'){
    if(filename[i] == '/')
      pos = i;
    i++;
  }

  if(pos != 0){
    pos++;
    strncpy(incfilename,filename,pos);
    incfilename[pos] = '\0';
  }

  /* stamp circuit to G,C,B and waveform */   
  num_L_tmp = 0;
  num_I_tmp = -1;
  num_V_tmp = -1;

  while(!feof(fid)||(fid_tmp!=NULL && !feof(fid_tmp))){
    if(feof(fid)==0){
      fgets(strline, READ_BLOCK_SIZE,fid);
    }
    else{
      if (fgets(strline, READ_BLOCK_SIZE, fid_tmp) == NULL) break;
    }


    switch(strline[0])
      {

	//  ************ Capacitor ***************
      case'C': case'c':
	if(sscanf(strline, "%*s %s %s %s", node1, node2, valuestr) == 3){
	  value = StrToNum(valuestr);
	  n1_tmp = nodePool->findorPushNode(node1);
	  n2_tmp = nodePool->findorPushNode(node2);
	  n1 = nodePool->getNode(n1_tmp)->row_no;
	  n2 = nodePool->getNode(n2_tmp)->row_no;
	  if (n1 != GNDNODE) C->pushEntry(n1, n1, value);
	  if (n2 != GNDNODE) C->pushEntry(n2, n2, value);
	  if (n1 != GNDNODE && n2 != GNDNODE)
	    {
	      C->pushEntry(n1,n2,-value);
	      C->pushEntry(n2,n1,-value);
	    }
	}
	else
	  printf("Fail in obtaining capacitor value.\n");

	break;

        //  ********** .include ********************
      case'.':
	if (strline[1]=='i' && strline[2]=='n'){ //.include
	  if(sscanf(strline, "%*s %s", tmp_incfilename) == 1){
	    i = 0;
	    j = pos;
	    while(tmp_incfilename[i] != '\0'){
	      if(tmp_incfilename[i]!='\"'){
		incfilename[j] = tmp_incfilename[i];
		j++;
	      }
	      i++;
	    }
	    incfilename[j] = '\0';


	    fid_tmp = fid;

	    fid = fopen(incfilename,"r");
	    if(fid == NULL){
	      printf("Open file Error!\n");
	      exit(-1);
	    }	    
	  }
	}
	break;


      default:
	break;
      }
  }
  fclose(fid);

  if(fid_tmp != NULL){
    fclose(fid_tmp);
  }

  C->sort();

  free(strline);   
  free(tmp_incfilename);
  free(incfilename);
  free(node1);
  free(node2);
  free(timestr);
  free(stimestr);
  free(valuestr);

}

void stampB(const char* filename, int nL, int nIS, int nVS, int nNodes, double tstop, Source *VS, Source *IS, matrix* B, NodeList* nodePool)
{
  char* strline;
  char* node1, *node2;
  char *timestr, *stimestr, *valuestr;
  char *incfilename, *tmp_incfilename;
  char *istr;

  char isrcName[128];
  
  double ptime, value, Rvalue;

  int num_point;
  int n1, n2, n1_tmp, n2_tmp;
  int num_L_tmp, num_I_tmp, num_V_tmp;
  int index_i, index_j;
  int i,j,k,mark_value, pos;
  int nblock;
  int last_npoint_I, last_npoint_V;

  bool isI, isV;

  FILE* fid_tmp=NULL;

  FILE* fid = fopen(filename, "r");
  if(fid == NULL){
    printf("Open file Error!\n");
    exit(-1);
  }

  strline = (char*)malloc(READ_BLOCK_SIZE*sizeof(char));
  incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  tmp_incfilename = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node1 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  node2 = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  istr = (char*)malloc(NAME_BLOCK_SIZE*sizeof(char));
  timestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  stimestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));
  valuestr = (char*)malloc(VALUE_BLOCK_SIZE*sizeof(char));

  /* get the director of file */
  i = 0;
  pos = 0;
  while(filename[i]!='\0'){
    if(filename[i] == '/')
      pos = i;
    i++;
  }

  if(pos != 0){
    pos++;
    strncpy(incfilename,filename,pos);
    incfilename[pos] = '\0';
  }

  /* stamp circuit to G,C,B and waveform */   
  num_L_tmp = 0;
  num_I_tmp = -1;
  num_V_tmp = -1;
  last_npoint_I = 0;
  last_npoint_V = 0;

  while(!feof(fid) || (fid_tmp != NULL && !feof(fid_tmp))){

    if(feof(fid)==0){
      fgets(strline, READ_BLOCK_SIZE,fid);
    }
    else{
      if (fgets(strline, READ_BLOCK_SIZE, fid_tmp) == NULL) break;
    }

    switch(strline[0])
      {


	// ************ Voltage *******************
      case'V': case'v':
	num_V_tmp++;
	isV = true; isI = false;
	if(sscanf(strline, "%*s %s %s %s",  node1, node2, istr) == 3){
	  n1_tmp = nodePool->findorPushNode(node1);
	  n2_tmp = nodePool->findorPushNode(node2);
	  n1 = nodePool->getNode(n1_tmp)->row_no;
	  n2 = nodePool->getNode(n2_tmp)->row_no;
	  
	  index_i = nNodes + nL + num_V_tmp;
	  index_j = num_V_tmp;

	  B->pushEntry(index_i, index_j, -1);

	  //DC
	  value = StrToNum(istr);
	  VS[num_V_tmp].time.set_size(2, false);
	  VS[num_V_tmp].value.set_size(2, false);
	  VS[num_V_tmp].time(0) = 0;
	  VS[num_V_tmp].time(1) = tstop;
	  VS[num_V_tmp].value(0) = value;
	  VS[num_V_tmp].value(1) = value;
	}
	break;

	// ******************** Current source ******************
      case'I': case'i':
	num_I_tmp++;
	isV = false; isI = true;
	if(sscanf(strline, "%s %s %s %s", isrcName, node1, node2, istr) == 4){ // XXLiu
	  n1_tmp = nodePool->findorPushNode(node1);
	  n2_tmp = nodePool->findorPushNode(node2);
	  n1 = nodePool->getNode(n1_tmp)->row_no;
	  n2 = nodePool->getNode(n2_tmp)->row_no;
	  
	  index_j = nVS + num_I_tmp; 

	  if(n1 != GNDNODE) {
	    B->pushEntry(n1,index_j,-1);
	  }
	  if(n2 != GNDNODE) {
	    B->pushEntry(n2,index_j,1);
	  }

	  //DC
	  value = StrToNum(istr)*I_SCALE;//special
	  IS[num_I_tmp].time.set_size(2, false);
	  IS[num_I_tmp].value.set_size(2, false);  
	  IS[num_I_tmp].time(0) = 0;
	  IS[num_I_tmp].time(1) = tstop;
	  IS[num_I_tmp].value(0) = value;
	  IS[num_I_tmp].value(1) = value;
	}

	break;

      default:
	break;
      }
  }
  fclose(fid);

  if(fid_tmp != NULL){
    fclose(fid_tmp);
  }

  B->sort();

  free(strline);   
  free(tmp_incfilename);
  free(incfilename);
  free(node1);
  free(node2);
  free(istr);
  free(timestr);
  free(stimestr);
  free(valuestr);
}
