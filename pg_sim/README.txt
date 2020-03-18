

The power grid simulator (PGSim) is a netlist parser, IR drop simulator (nodal 
voltage solver) and EM filter.


*******************************************************************************
*
*						Introduction
*
*******************************************************************************

Input: an IBM-format spice file (.sp) of a power grid.
Output: 
	1. The EM condition of the power grid, including total interconnect tree
	number, EM voltage and saturation volume (if exists) of each tree, failed
	tree number and index. (On screen)
	2. Nodal voltage of the whole power grid. (File: tree_node_voltage.txt)
	3. Current density of each branch. (File: tree_curden.txt) 
	4. Width of each branch. (File: tree_width.txt)

*******************************************************************************
*
*						Useage
*
*******************************************************************************

**Source code:
	element.cpp		element.h		mna.cpp				mna.h		
	parser.cpp		parser.h		function_ibm.cpp	function_ibm.h
	em_cmd_ibm.cpp	parameter.h

**Library depended: 
	1. IT++: IT++ is a C++ library of mathematical, signal processing and
	communication classes and functions. (Version: itpp-4.3.1)
	2. UMFPACK: multifrontal LU factorization.
	3. CSparse and CXSparse: a concise sparse Cholesky factorization package.


**Commands to compile the source code: 
	make clean;make


**Executable file:
	em_cmd


**Commands to run with one spice file (.sp):
	rm tree_*
	touch tree_ftid.txt
	touch tree_width.txt
	./em_cmd {filename}


*******************************************************************************
*
*						PGSim in EMSpice
*
*******************************************************************************
PGSim is the first step of EMSpice, it helps parse the power grid and determine 
its electrical condition, and filter out all mortal wires.

If mortal wire exists, it will cause either early failure (leads to open circuit)
or late failure (leads to resistance increase). Therefore, iterative coupled 
simulation will be performed. 

In each iteration, the EMSpice (EM_spice.py) will first call PGSim, with the 
input of either original or modified spice file where the modification was
performed through EMSpice during the last iteration, then output new voltage, 
current and width information for stress simulation.

