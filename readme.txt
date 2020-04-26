This directory include a few modules:
emspice_python: the emspice codes in python
emspice_matlab: the emspice codes in the matlab
pg_sim: the linear DC solver for power grid networks
icc2spce: convert the Synopysis icc p/g networks format into spice format (from *.pna to *.sp)
docs: readme files and documents for different modules

Basically flows work like this. 
(1) First you need to convert the *.pna files from ICC into the spice *.sp file. 
(2) Run the EM spice file with output .sp file
Please see the readme files in each moduel
