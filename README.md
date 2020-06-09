# EMspice v1.0 

Last update: June 2020

EMspice is a Coupled EM-IR Analysis Tool for Full-Chip Power Grid EM and IR Check and Sign-off. EMspice was developed by Zeyu Sun, Han Zhou, Yibo Liu and Sheldon Tan at UC Riverside in 2020. EMspice performs the mult-physics electrical and stress analysis of multi-segment interconnect wires. The repot consists of phython version and matlab version impementation of EMspice. 

EMspice v1.0 is released under the terms of BSD 3-Clause License (see LICENSE for details). 

## References 

C. Cook, Z. Sun, E. Demircan, M. Shroff,  S. X-.D. Tan, “Fast electromigration stress evolution analysis for interconnect trees using Krylov subspace method”, IEEE Transactions on Very Large Scale Integrated Systems  (TVLSI),  Vol. 26, no. 5, pp. 969-980, May 2018.

Z. Sun, S. Sadiqbatcha, H. Zhao, and S. X.-D. Tan, “Saturation volume estimation for  multi-segment copper interconnect wire”, IEEE Transactions on Very Large Scale Integrated Systems  (TVLSI), vol. 27, no. 7, pp.1063-8210, July, 2019, DOI: https://doi.org/10.1109/TVLSI.2019.2901824

Z. Sun, S. Yu, H. Zhou, Y. Liu and S. X.-D. Tan, “EMSpice: physics-based electromigration check using coupled electronic and stress simulation”, IEEE Transaction on Device and Materials  Reliability (T-DMR), Mach 2020. 10.1109/TDMR.2020.2981628

Sheldon X.-D. Tan, Mehdi Tahoori, Taeyoung Kim, Shengcheng Wang, Zeyu Sun and Saman Kiamehr, “VLSI Systems Long-Term Reliability -- Modeling, Simulation and Optimization”,  Springer Publisher, 2019. DOI: 10.1007/978-3-030-26172-6, ISBN: 978-3-030-26171-9 (https://www.springer.com/gp/book/9783030261719). 

This project was funded by National Science Fundation (CCF-1816361 and CCF-1527324). See details at https://intra.ece.ucr.edu/~stan/project/vsclab_wiki_new/index.php/VLSI_reliability,_resilience,_fault-tolerant_computing_and_dynamic_reliability_management

## The major components of the program

**ICC to spice code in "EMspice_python/icc2spice/"**

**PG_solver code in EMspice_python/pg_sim/"**

**Coupled simulation code and framework in "EMspice_python/emspice_python/"**

**GUI code in "EMspice_python/em_gui/"**

Here are details on how to run the code

## Step 1. Get the spice file from ICC result

This step we convert the output from ICC to spice file which can be read by PG_solver

Find the {pnafilename}.pna file from icc simulation. Run the syn2spice.py to get the spice file. 

Input is {pnafilename}.pna file

Output is {spicefilename}.sp file

Command is python syn2spice.py

## Step 2. Coupled simulation using EM solver and PG solver for EM and IR drop checks
 
This step we do the coupled simulation with PG_solver and EM_solver

Use the spice file ({spicefilename}.sp) from the step 1, and the constant for simulation in the file EM_spice_constant.txt 
as the input for this step. EM_spice.py  is the EM solver and em_cmd is the binary code of PG solver. 
Do the filter after we get the current density from PG solver and then do coupled simulation with EM solver 
and PG solver. The stress, void and resistance are updated after each simulation step. 

Input is Spice({spicefilename}.sp) file(from step 1) EM_spice_constant.txt (constant for simulation).

Output is Void.txt: position of void. Lvoid.txt: Length of void. Length of void at corresponding position. 
u_stress.txt: stress on each segment tree_node_voltage.txt nodal voltage on all points.

Command is python EM_spice.py {spicefilename}.sp

## Step 3 Display the information in GUI

This step display the simulation result from the GUI. 

The GUI part take structure information current density stress and void location and size information 
and display that to the user graphally. 

Input is {pnafilename}.pna file and {spicefilename}.sp file for geomotery information, u_stress for stress information, u_curden.txt for current density information, Void and LVoid for void information and tree_node_voltage.txt for voltage information. 

Output is displayed to user with MatplotLib. 

More information please refer to detailed guidlines in "EMspice_python/docs/"
