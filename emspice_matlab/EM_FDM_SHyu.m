clc;
clear;
clear all;
warning off Matlab:nearlySingularMatrix;

tic;

% EM FDM solver
%% ---------- physical parameters
T  = 373; % Temperature in Kelvin
D0 = 7.5e-8;% -8  %maximum diffusion coefficient (at infinite temperature; m2 / s)
E=1.609e-19;
Ea = 0.84*E;%activation energy for diffusion in dimensions of (J atom�?)
kB = 1.3806e-23;%Boltzmann constant.
Da = D0*exp(-Ea/(kB*T));%diffusion coefficient (m2 / s)
B0  =  1e11;%effective bulk elasticity modulus
Omega=1.182e-29;%atomic lattice volume

cu_res = 3e-8;
hCu = 9e-8;
Ta_res = 1.76e-07;
hTa = 2e-8;
Z=10;

kappa = (Da*B0*Omega)/(kB*T);
% j1=2e9;
% g1 = (E*Z*cu_res*j1)/(Omega);

nx=20; % number of nodes
wire_length=1; %normalized length

T0 = zeros(1,nx);

t_total = 3600*24*30;  % one month time
t_num = 10;
tstep = t_total/t_num;
% t_num=int32(t_total/tstep);
time = ones(t_num);
time = time*tstep;
tstamp=[0:tstep:t_total-tstep];

tstamp_flag = 1;
istep0 = 1;
stop_flag = 0;

% --- Creates peicewise linear current density
cdnum = 10 % number of different current densities
cdList = ones(cdnum) % creates the PWL waveform ( just constant now )
pwl = pwl_create(t_num,cdList); % cretes current density trace ( current density at each time step )

init_flag = 1;  % initial status flag, if = 1, means the process is in the initial status; else, means not
iterIdx = 0;  % The EM violation FDTD method iteration #, every time step is the demorlization of t_total

% fid_resis_6 = fopen('u_resis_6.txt', 'a+');
% fid_resis_7 = fopen('u_resis_7.txt', 'a+');
% fid_resis_38 = fopen('u_resis_38.txt', 'a+');
% fid_resis_39 = fopen('u_resis_39.txt', 'a+');
% fid_resis_44 = fopen('u_resis_44.txt', 'a+');
% fid_resis_45 = fopen('u_resis_45.txt', 'a+');

%% iteration part
for istep = 1:1201
    
    iterIdx = iterIdx + 1;  % to name the output stress, current density, and Lvoid file
    if stop_flag == 1
        break;
    end
    
    %% run c++ part to get branch, current information
    ExeFileName = 'em_cmd';
    ExeFilePath = fullfile('./', ExeFileName);
    %% Param = [' ', 'armcore.sp'];
    Param = [' ', 'armcore.sp'];
    Cmd = [ExeFilePath, Param];
    system(Cmd);
    
    %% input branch name, node information, resistance
    wires0 = '';    %store branch information of the same tree once temporarily
    wires = {};     %tree information
    fid = fopen('tree_info.txt');
    treenum = 1;    %tree #
    while ~feof(fid)
        while 1
            str = fgetl(fid);    %only read one line, in order to tell the difference from branch information & flag "xx"
            if strcmp(str,'xx')
                break;              %get out of while circulation when meet flag "xx"
            end
            str = strrep(str, '_', ' ');      %get branch information temporarily
            wires0 = [wires0,str];        %continuously get branch information: branch 1branch 2branch 3...
        end
        wires0 = textscan(wires0,'R%s n%s %d %d n%s %d %d %f','CommentStyle','*');      %list the branch information branch by branch based on a certain pattern, change the wires0 from string variable into cell variable
        wires{treenum} = wires0;        %store the tree information into the certain location of cell variable "wires"
        wires0 = '';           %reset wires0
        treenum = treenum + 1;
    end
    fclose(fid);   
    treenum = treenum - 1;  % get the correct total tree number
  
    %% get branch # of every tree
    brchnum = ones(1,treenum);       %record branch number of every tree
    for i = 1:treenum
        brchnum(i) = numel(wires{i}{1});
    end
    
    fname = ['tree_info_', num2str(iterIdx), '.txt'];
    fid_info = fopen(fname, 'a+');
    for i = 1:treenum
        for j = 1:brchnum(i)
            fprintf(fid_info, 'R%s n%s %d %d n%s %d %d %f\r\n', char(wires{i}{1}(j)), char(wires{i}{2}(j)), wires{i}{3}(j), wires{i}{4}(j), char(wires{i}{5}(j)), wires{i}{6}(j), wires{i}{7}(j), wires{i}{8}(j));
        end
    end
    fclose(fid_info);
    
    %% input branch current density
    J_list = cell(1,treenum);        %current density information, but we only need the last element
    J_list_norm = cell(1,treenum);
    J_list0 = '';
    fid = fopen('tree_curden.txt');
    for i = 1:treenum
        for j = 1:brchnum(i)
            str = fgetl(fid);
            J_list0 = [J_list0, str];        
        end
        J_list{i} = textscan(J_list0, '%s %s %s %f');
        J_list{i} = J_list{i}{end};
        J_list_norm{i} = J_list{i}/norm(J_list{i},inf);
        J_list0 = '';
    end
    fclose(fid);
    
    fname = ['tree_curden_', num2str(iterIdx), '.txt'];
    fid_curd = fopen(fname, 'a+');
    for i = 1:treenum
        for j = 1:brchnum(i)
            fprintf(fid_curd, '%f\r\n', J_list{i}(j));
        end
    end
    fclose(fid_curd);
    
    %% input branch width and get tree width
    if init_flag == 1
        fid = fopen('tree_width.txt');
        brchwidth = cell(1,treenum);
        brchwidth0 = '';
        for i = 1:treenum
            for j = 1:brchnum(i)
                str = fgetl(fid);
                brchwidth0 = [brchwidth0, str];
            end
            brchwidth{i} = textscan(brchwidth0, 'W%f');     % get the length of branch
            brchwidth0 = '';
        end
        fclose(fid);
        twidth = zeros(1, treenum);    % the width of the cross section
        twidth1 = zeros(1, treenum);    % the width of the tree
        for i=1:treenum
            twidth(i) = brchwidth{i}{1}(1);
        end
        for i=1:treenum
            twidth1(i) = twidth(i);
        end
%         for i=37:50
%             twidth(i) = twidth(i)/25*4*1.5;
%         end
%         for i=3:4
%             twidth(i) = twidth(i)/25*4;
%         end
    end
    
    %% input branch length and get tree length
    if init_flag == 1
        fid = fopen('tree_leng.txt');
    else
        fid = fopen('tree_leng1.txt');
    end
    brchleng = cell(1,treenum);
    brchleng0 = '';
    for i = 1:treenum
        for j = 1:brchnum(i)
            str = fgetl(fid);
            brchleng0 = [brchleng0, str];
        end
        brchleng{i} = textscan(brchleng0, 'L%f');     % get the length of branch
        brchleng0 = '';
    end
    fclose(fid);

    treeleng = zeros(1,treenum);                      % get the length of tree
    for i = 1:treenum
        for j = 1:brchnum(i)
            treeleng(i) = treeleng(i)+brchleng{i}{1}(j);
        end
    end
    
    fname = ['tree_leng_', num2str(iterIdx), '.txt'];
    fid_leng = fopen(fname, 'a+');
    for i = 1:treenum
        for j = 1:brchnum(i)
            fprintf(fid_leng, '%f\r\n', brchleng{i}{1}(j));
        end
    end
    fclose(fid_leng);
    
    %% tree node voltage
    fid = fopen('tree_node_voltage.txt');
    volt0 = '';
    node_num = 0;
    while ~feof(fid)
        str = fgetl(fid);   
        volt0 = [volt0,str];
        node_num = node_num+1;
    end
    volt = textscan(volt0,'n%s %f','CommentStyle','*');
    fclose(fid);
    
    fname = ['tree_node_voltage_', num2str(iterIdx), '.txt'];
    fid_nodevolt = fopen(fname, 'a+');
    for i = 1:node_num
            fprintf(fid_nodevolt, '%f\r\n', volt{2}(i));
    end
    fclose(fid_nodevolt);
    
    %% input spice information
    wire_spice = {};
    %% fid = fopen('armcore.sp');
    fid = fopen('armcore.sp');
    line = 0;                                    %   scan each line in spice file
    while ~feof(fid)
        line = line + 1;
        str = fgetl(fid);
        wire_spice{line} = str;
    end
    fclose(fid);
    line = max(size(wire_spice));                   %   total # of lines of the spice file
    
    fname = ['armcore_', num2str(iterIdx), '.sp'];
    fid_sp = fopen(fname, 'a+');
    for i = 1:line
        fprintf(fid_sp, '%s\r\n', char(wire_spice{i}));
    end
    fclose(fid_sp);
    
    %% input failed tree index and compared it with the formal ones
    if init_flag == 1
        fail_Idx = zeros(1,treenum);
        fcell0 = cell(1,1);
        fstr0 = '';
        fid = fopen('tree_fid.txt');
        for i = 1:treenum
            str = fgetl(fid);
            fstr0 = [fstr0, str];
        end
        fcell0 = textscan(fstr0, 'F%d');
        for i = 1:treenum
            fail_Idx(i) = fcell0{1}(i);
        end
        fail_Idx_last = fail_Idx;
        fclose(fid);
    else
        fail_Idx = zeros(1,treenum);
        fcell0 = cell(1,1);
        fstr0 = '';
        fid = fopen('tree_fid.txt');
        for i = 1:treenum
            str = fgetl(fid);
            fstr0 = [fstr0, str];
        end
        fcell0 = textscan(fstr0, 'F%d');
        for i = 1:treenum
            fail_Idx(i) = fcell0{1}(i);
        end
        fclose(fid);
        if ~isequal(fail_Idx, fail_Idx_last)
            fid = fopen('tree_fid.txt', 'wt');
            fprintf(fid, '');
            fclose(fid);
%             break;  % means failed tree index changed
        end
    end
    
    %%  tree_info, tree_curden, tree_width, tree_leng, tree_fid reset to empty file

    fid  = fopen('tree_info.txt', 'wt');
    fprintf(fid, '');
    fclose(fid);
    fid  = fopen('tree_curden.txt', 'wt');
    fprintf(fid, '');
    fclose(fid);
%     fid  = fopen('tree_width.txt', 'wt');
%     fprintf(fid, '');
%     fclose(fid);
    fid = fopen('tree_leng.txt', 'wt');
    fprintf(fid, '');
    fclose(fid);
%     fid = fopen('tree_fid.txt', 'wt');
%     fprintf(fid, '');
%     fclose(fid);
    fid  = fopen('tree_number_r.txt', 'wt');
    fprintf(fid, '');
    fclose(fid);
    fid  = fopen('tree_node_voltage.txt', 'wt');
    fprintf(fid, '');
    fclose(fid);
    fid  = fopen('tree_axis.txt', 'wt');
    fprintf(fid, '');
    fclose(fid);
    
    out_spice_flag = 0;
    g1 = (E*Z*cu_res*norm(J_list{i},inf))/(Omega);
    
    
    %% Process
    
    if init_flag == 1  % set the initial status (stress, void flag, void length) of the tree
        initC = cell(1,treenum);
        initC_last = cell(1,treenum);
        void = zeros(1,treenum);                     %  flag "void: if the void starts growing(0/1)"
        Lvoid = zeros(1,treenum);                    %  record void length, and judge whether it is over cross section width
        Lvoid_last = Lvoid;
        Lvoid1 = zeros(1,treenum);
        Lvoid2 = zeros(1,treenum);
        Lvoid_last1 = Lvoid1;
        Lvoid_last2 = Lvoid2;
        max_stress = zeros(1,treenum);               %  "max_stress: the biggest stress on the wire (cathode stress)"
        max_stress_location = zeros(1,treenum);      %  "max_stress_location: cathode"
        nx_total = ones(1,treenum);
        for i = 1:treenum
            nx_total(i) = (nx*brchnum(i))-(brchnum(i)-1);
            initC{i} = zeros(nx_total(i),1);
        end
        init_flag = 0;
    end
    
    % open file for stess, Lvoid, current density output everytime
    fname1 = ['u_stress_', num2str(iterIdx), '.txt'];
    fid_stress = fopen(fname1, 'a+');
    fname2 = ['u_curden_', num2str(iterIdx), '.txt'];
    fid_curden = fopen(fname2, 'a+');
    fname3 = ['u_Lvoid_', num2str(iterIdx), '.txt'];
    fid_Lvoid = fopen(fname3, 'a+');
%     fname4 = ['u_resis_', num2str(iterIdx), '.txt'];
    fid_resis = fopen('u_resis.txt', 'a+');
    
    for i = 1:treenum
        Lvoid0 = 0;
        Lvoid01 = 0;
        Lvoid02 = 0;
        Lvoid_temp = 0;
        % output curden
        for k = 1:brchnum(i)
            fprintf(fid_curden, 'R%s %f\r\n', char(wires{i}{1}(k)), J_list{i}(k));
        end
        %
        
        for j = 1:1   % how many step the matlab FEM part runs
            
            if void(i) == 0
                [G,C,B] = matrix_formation(nx,wires{i},J_list_norm{i});
                initC{i} = initC{i}/g1/treeleng(i);    %  normalized stress
                a1 = length(initC{i});                   %   total point to be analyzed
                a2 = length(tstamp);              %   total time steps
                max_stress_point = zeros(a1,1);      %   record all the nodes whose stress is over 500MPa
                [sol] = backEuler_fdm(G,diag(C,0),B,tstamp/treeleng(i)^2*kappa,initC{i},pwl);
%                 tstamp = tstamp*treeleng(i)^2/kappa;        %   denormalized time scale
                initC{i} = sol(a2,:)'*g1*treeleng(i);   %  denormalized stress
                
                dx = treeleng(i)/(a1-1);
                [v,p] = max(initC{i});              %   v = "value", p = "position"
                max_stress(i) = v;
                max_stress_location(i) = p;
%                 plot(sol(a2,:));
                if max_stress(i) > 5e8
                    void(i) = 1;
                end
            else
                initC{i} = initC{i}/g1/treeleng(i);    %  normalized stress
                wires_temp = cell(1,length(wires{i}));
                if max_stress_location(i) == 1
                    
                    [G,C,B] = growth_matrix_formation(nx,wires{i},J_list_norm{i});
                    [sol] = backEuler_fdm(G,diag(C,0),B,tstamp/treeleng(i)^2*kappa,initC{i},pwl);
                    initC{i} = sol(a2,:)'*g1*treeleng(i);   %  denormalized stress
                    
                elseif max_stress_location(i) == nx_total(i)
                    
                    initC_temp = initC{i}(end:-1:1);  %  inverse
                    for k = 1:length(wires{i})
                        wires_temp{k} = wires{i}{k}(end:-1:1);     % truncate and inverse the left part
                    end
                    J_list_norm_temp = -J_list_norm{i}(end:-1:1);
                    [G,C,B] = growth_matrix_formation(nx,wires_temp,J_list_norm_temp);            % growth phase, can only deal with the case that void grows on the left side of the tree
                    [sol] = backEuler_fdm(G,diag(C,0),B,tstamp/treeleng(i)^2*kappa,initC_temp,pwl);
                    initC_temp = sol(a2,:)'*g1*treeleng(i);   %  denormalized stress
                    initC{i} = initC_temp(end:-1:1);  %  inverse again 
                    
                else
                    
                    for k = 1:length(wires{i})
                        wires_temp{k} = wires{i}{k}((max_stress_location(i)-1)/(nx-1):-1:1);     % truncate and inverse the left part
                    end
                    J_list_norm_temp = -J_list_norm{i}((max_stress_location(i)-1)/(nx-1):-1:1);  % notice the direction of current density should also inverse with the inversion of the left part of the tree
                    initC_temp = initC{i}(max_stress_location(i):-1:1);
                    [G,C,B] = growth_matrix_formation(nx,wires_temp,J_list_norm_temp);            % growth phase, can only deal with the case that void grows on the left side of the tree
                    [sol] = backEuler_fdm(G,diag(C,0),B,tstamp/treeleng(i)^2*kappa,initC_temp,pwl);
                    initC_temp = sol(a2,:)'*g1*treeleng(i)*((max_stress_location(i)-1)/(nx-1)/brchnum(i));   %  denormalized stress
                    initC_temp1 = initC_temp(end:-1:1);
                    % deal with the right part of the tree
                    for k = 1:length(wires{i})
                        wires_temp{k} = wires{i}{k}((max_stress_location(i)-1)/(nx-1)+1:1:end);     % truncate and inverse the left part
                    end
                    J_list_norm_temp = J_list_norm{i}((max_stress_location(i)-1)/(nx-1)+1:1:end);
                    initC_temp = initC{i}(max_stress_location(i):1:end);
                    [G,C,B] = growth_matrix_formation(nx,wires_temp,J_list_norm_temp);            % growth phase, can only deal with the case that void grows on the left side of the tree
                    [sol] = backEuler_fdm(G,diag(C,0),B,tstamp/treeleng(i)^2*kappa,initC_temp,pwl);
                    initC_temp = sol(a2,:)'*g1*treeleng(i)*((brchnum(i)-(max_stress_location(i)-1)/(nx-1))/brchnum(i));   %  denormalized stress
                    
                    initC_temp(1) = (initC_temp(1) + initC_temp1(end))/2;
                    initC_temp1 = initC_temp1(1:1:end-1);
                    initC{i} = [initC_temp1', initC_temp']';
                    
                end
                
%                 [G,C,B] = growth_matrix_formation(nx,wires{i},J_list_norm{i},max_stress_location(i));           %growth phase
%                 initC{i} = initC{i}/g1/treeleng(i);    %  normalized stress
%                 [sol] = backEuler_fdm(G,diag(C,0),B,tstamp,initC{i},pwl);
%                 initC{i} = sol(a2,:)'*g1*treeleng(i);   %  denormalized stress
                a1 = length(initC{i});                   %   total point to be analyzed
                a2 = length(tstamp);              %   total time steps
                dx = treeleng(i)/(a1-1);
                
                if max(abs(initC{i}-initC_last{i})) < 1e6
                    stop_flag = 1;
                end
                
                Lvoid_temp = 0;
                Lvoid_temp1 = 0;
                Lvoid_temp2 = 0;
                    
                %% calculate Lvoid, and adjust the stress of the next location
                if max_stress_location(i) == a1                     
                    for k=1:a1
                        Lvoid_temp = Lvoid_temp-initC{i}(k)/B0*dx;
                    end
                    Lvoid(i) = Lvoid_temp;
                    
%                     treeleng(i) = treeleng(i)-Lvoid(i);
                    dx = treeleng(i)/(a1-1);
%                     for k=1:(a1-2)
%                         Lvoid0 = Lvoid0-initC{i}(k)/B0*dx;
%                     end
%                     Lvoid0 = Lvoid(i)-Lvoid0;
%                     initC{i}(a1-1) = Lvoid0*B0/dx;
                    brchleng{i}{1}(brchnum(i)) = brchleng{i}{1}(brchnum(i))-(Lvoid(i)-Lvoid_last(i));        % update the branch length;
                    
                    if Lvoid(i)>twidth(i)
                        
                        out_spice_flag = 1;
                        
                        wires{i}{8}(brchnum(i)) = Ta_res*(Lvoid(i)-Lvoid_last(i))/(hTa*twidth1(i)) + wires{i}{8}(brchnum(i));
                        fprintf(fid_resis, 'R%s\t%.6f\r\n', char(wires{i}{1}(brchnum(i))), wires{i}{8}(brchnum(i)));
                        
%                         switch i
%                             case 6 
%                                 fprintf(fid_resis_6, '%.6f\r\n', sum(wires{i}{8}));
%                             case 7
%                                 fprintf(fid_resis_7, '%.6f\r\n', sum(wires{i}{8}));
%                             case 38
%                                 fprintf(fid_resis_38, '%.6f\r\n', sum(wires{i}{8}));
%                             case 39
%                                 fprintf(fid_resis_39, '%.6f\r\n', sum(wires{i}{8}));
%                             case 44
%                                 fprintf(fid_resis_44, '%.6f\r\n', sum(wires{i}{8}));
%                             case 45
%                                 fprintf(fid_resis_45, '%.6f\r\n', sum(wires{i}{8}));
%                         end
                        
%                         if tstamp_flag == 1
%                             tstamp = tstamp/10;
%                             tstamp_flag = 0;
%                         end
                        %% update wire_spice
                        brchname = {};
                        for k = 1:line
                            if strcmp(wire_spice{k}, '')
                                continue;                      % avoid textscan function meet error
                            end
                            brchname = textscan(wire_spice{k}, 'R%s %s %s %s', 'CommentStyle','*');
                            if strcmp(brchname{1}, wires{i}{1}(brchnum(i)))       % compare the branch name whose resistance requires update
                                brchname0 = strrep(wires{i}{1}(brchnum(i)), '-', ' ');
                                netID = textscan(char(brchname0), '%d %d %d');
                                netID = netID{1};
                                if (netID == 25 && J_list_norm{i}(brchnum(i))>0)
                                    if wires{i}{3}(brchnum(i)) == wires{i}{6}(brchnum(i))
                                        new_node_name = ['n', '1', '_', num2str(wires{i}{6}(brchnum(i))), '_', num2str(wires{i}{7}(brchnum(i))+1)];
                                    else
                                        new_node_name = ['n', '1', '_', num2str(wires{i}{6}(brchnum(i))+1), '_', num2str(wires{i}{7}(brchnum(i)))];
                                    end
                                    wire_spice{k} = strrep(wire_spice{k}, brchname{3}, new_node_name);
                                    wire_spice{k} = char(wire_spice{k});
%                                     break;    %  early failure
                                end
                                
                                str = num2str(wires{i}{8}(brchnum(i)));
                                wire_spice{k} = strrep(wire_spice{k}, brchname{4}, str);
                                wire_spice{k} = char(wire_spice{k});
                                if (netID == 25 && J_list_norm{i}(brchnum(i))>0)
                                    wire_spice{end-1} = [wire_spice{end-1}, ' v(', new_node_name, ')'];
                                end
                                % wire_spice{end-1} = [wire_spice{end-1}, ' v(', new_node_name, ')'];
                                break;                         % since only need to update one of the wire_spice cells
                            end
                        end
                    end
                    
                elseif max_stress_location(i) == 1
                    for k=1:a1
                        Lvoid_temp = Lvoid_temp-initC{i}(k)/B0*dx;
                    end
                    Lvoid(i) = Lvoid_temp;
%                     treeleng(i) = treeleng(i)-Lvoid(i);
                    dx = treeleng(i)/(a1-1);
%                     Lvoid_last(i) = Lvoid(i);
%                     for k=3:a1
%                         Lvoid0 = Lvoid0-initC{i}(k)/B0*dx;
%                     end
%                     Lvoid0 = Lvoid(i)-Lvoid0;
%                     initC{i}(2) = Lvoid0*B0/dx;
                    brchleng{i}{1}(1) = brchleng{i}{1}(1)-(Lvoid(i)-Lvoid_last(i));       % update the branch length;
                    
                    if Lvoid(i)>twidth(i)   % update the resistance
%                         if tstamp_flag == 1
%                             tstamp = tstamp/10;
%                             tstamp_flag = 0;
%                         end
                        
                        out_spice_flag = 1;
                        
                        wires{i}{8}(1) = Ta_res*(Lvoid(i)-Lvoid_last(i))/(hTa*twidth1(i)) + wires{i}{8}(1);
                        fprintf(fid_resis, 'R%s\t%.6f\r\n', char(wires{i}{1}(1)), wires{i}{8}(1));
                        
%                         switch i
%                             case 6 
%                                 fprintf(fid_resis_6, '%.6f\r\n', sum(wires{i}{8}));
%                             case 7
%                                 fprintf(fid_resis_7, '%.6f\r\n', sum(wires{i}{8}));
%                             case 38
%                                 fprintf(fid_resis_38, '%.6f\r\n', sum(wires{i}{8}));
%                             case 39
%                                 fprintf(fid_resis_39, '%.6f\r\n', sum(wires{i}{8}));
%                             case 44
%                                 fprintf(fid_resis_44, '%.6f\r\n', sum(wires{i}{8}));
%                             case 45
%                                 fprintf(fid_resis_45, '%.6f\r\n', sum(wires{i}{8}));
%                         end
                        
                        %% update wire_spice
                        brchname = {};
                        for k = 1:line
                            if strcmp(wire_spice{k}, '')
                                continue;                      % avoid textscan function meet error
                            end
                            brchname = textscan(wire_spice{k}, 'R%s %s %s %s', 'CommentStyle','*');
                            if strcmp(brchname{1}, wires{i}{1}(1))       % compare the branch name whose resistance requires update
                                brchname0 = strrep(wires{i}{1}(1), '-', ' ');
                                netID = textscan(char(brchname0), '%d %d %d');
                                netID = netID{1};
                                if (netID == 25 && J_list_norm{i}(1)<0)
                                    if wires{i}{3}(1) == wires{i}{6}(1)
                                        new_node_name = ['n', '1', '_', num2str(wires{i}{3}(1)), '_', num2str(wires{i}{4}(1)+1)];
                                    else
                                        new_node_name = ['n', '1', '_', num2str(wires{i}{3}(1)+1), '_', num2str(wires{i}{4}(1))];
                                    end
                                    wire_spice{k} = strrep(wire_spice{k}, brchname{2}, new_node_name);
                                    wire_spice{k} = char(wire_spice{k});
%                                     break;    %  early failure
                                end
                                
                                str = num2str(wires{i}{8}(1));
                                wire_spice{k} = strrep(wire_spice{k}, brchname{4}, str);
                                wire_spice{k} = char(wire_spice{k});
                                if (netID == 25 && J_list_norm{i}(1)<0)
                                    wire_spice{end-1} = [wire_spice{end-1}, ' v(', new_node_name, ')'];
                                end
                                % wire_spice{end-1} = [wire_spice{end-1}, ' v(', new_node_name, ')'];
                                break;                         % since only need to update one of the wire_spice cells
                            end
                        end
                    end
                    
                else
                    for k=1:max_stress_location(i)-1
                        Lvoid_temp1 = Lvoid_temp1-initC{i}(k)/B0*dx;
                    end
                    for k=max_stress_location(i)+1:a1
                        Lvoid_temp2 = Lvoid_temp2-initC{i}(k)/B0*dx;
                    end
                    Lvoid1(i) = Lvoid_temp1;
                    Lvoid2(i) = Lvoid_temp2;
%                     treeleng(i) = treeleng(i)-Lvoid1(i)-Lvoid2(i);
                    dx = treeleng(i)/(a1-1);
                    Lvoid(i) = Lvoid1(i)+Lvoid2(i);
                    if Lvoid1(i) < 0
                        Lvoid2(i) = Lvoid(i);
                        Lvoid1(i) = 0;
                    end
                    if Lvoid2(i) < 0
                        Lvoid1(i) = Lvoid(i);
                        Lvoid2(i) = 0;
                    end
                    
%                     for k=1:max_stress_location(i)-2
%                         Lvoid01 = Lvoid01-initC{i}(k)/B0*dx;
%                     end
%                     for k=max_stress_location(i)+2:a1
%                         Lvoid02 = Lvoid02-initC{i}(k)/B0*dx;
%                     end
%                     Lvoid01 = Lvoid1(i)-Lvoid01;
%                     Lvoid02 = Lvoid2(i)-Lvoid02;
%                     initC{i}(max_stress_location(i)-1) = Lvoid01*B0/dx;
%                     initC{i}(max_stress_location(i)+1) = Lvoid02*B0/dx;
%                     Lvoid01(i) = 0;      %  otherwise, the Lvoid01 will accumulate
%                     Lvoid02(i) = 0;
                    brchleng{i}{1}((max_stress_location(i) - 1)/(nx-1)) = brchleng{i}{1}((max_stress_location(i) - 1)/(nx-1))-(Lvoid1(i)-Lvoid_last1(i));
                    brchleng{i}{1}((max_stress_location(i) - 1)/(nx-1)+1) = brchleng{i}{1}((max_stress_location(i) - 1)/(nx-1)+1)-(Lvoid2(i)-Lvoid_last2(i));
                    
                    if Lvoid(i)>twidth(i)   % update the resistance
                        
                        out_spice_flag = 1;
                        
                        wires{i}{8}((max_stress_location(i) - 1)/(nx-1)) = Ta_res*(Lvoid1(i)-Lvoid_last1(i))/2/(hTa*twidth1(i)) + wires{i}{8}((max_stress_location(i) - 1)/(nx-1));
%                         fprintf(fid_resis, 'R%s\t%.6f\r\n', char(wires{i}{1}((max_stress_location(i) - 1)/(nx-1))), wires{i}{8}((max_stress_location(i) - 1)/(nx-1)));
                        wires{i}{8}((max_stress_location(i) - 1)/(nx-1)+1) = Ta_res*(Lvoid2(i)-Lvoid_last2(i))/2/(hTa*twidth1(i)) + wires{i}{8}((max_stress_location(i) - 1)/(nx-1)+1);
%                         fprintf(fid_resis, 'R%s\t%.6f\r\n', char(wires{i}{1}((max_stress_location(i) - 1)/(nx-1)+1)), wires{i}{8}((max_stress_location(i) - 1)/(nx-1)+1));
                        
                        fprintf(fid_resis, '%d\t%.6f\t%.6f\r\n', istep, wires{i}{8}((max_stress_location(i) - 1)/(nx-1)), wires{i}{8}((max_stress_location(i) - 1)/(nx-1)+1));
                        
%                         switch i
%                             case 6 
%                                 fprintf(fid_resis_6, '%.6f\r\n', sum(wires{i}{8}));
%                             case 7
%                                 fprintf(fid_resis_7, '%.6f\r\n', sum(wires{i}{8}));
%                             case 38
%                                 fprintf(fid_resis_38, '%.6f\r\n', sum(wires{i}{8}));
%                             case 39
%                                 fprintf(fid_resis_39, '%.6f\r\n', sum(wires{i}{8}));
%                             case 44
%                                 fprintf(fid_resis_44, '%.6f\r\n', sum(wires{i}{8}));
%                             case 45
%                                 fprintf(fid_resis_45, '%.6f\r\n', sum(wires{i}{8}));
%                         end

%                         if tstamp_flag == 1
%                             tstamp = tstamp/10;
%                             tstamp_flag = 0;
%                         end
                        %% update wire_spice
                        brchname = {};
                        process = 0;
                        for k = 1:line
                            if strcmp(wire_spice{k}, '')
                                continue;                      % avoid textscan function meet error
                            end
                            brchname = textscan(wire_spice{k}, 'R%s %s %s %s', 'CommentStyle','*');
                            if strcmp(brchname{1}, wires{i}{1}((max_stress_location(i) - 1)/(nx-1)))       % compare the branch name whose resistance requires update
                                brchname0 = strrep(wires{i}{1}((max_stress_location(i) - 1)/(nx-1)), '-', ' ');
                                netID = textscan(char(brchname0), '%d %d %d');
                                netID = netID{1};
                                cath_brch_flag = (max_stress_location(i) - 1)/(nx-1);
                                if (netID == 25 && J_list_norm{i}((max_stress_location(i) - 1)/(nx-1))>0)
                                    
                                    if wires{i}{3}(cath_brch_flag)+wires{i}{4}(cath_brch_flag)-wires{i}{6}(cath_brch_flag)-wires{i}{7}(cath_brch_flag)>0
                                        if wires{i}{3}(cath_brch_flag) == wires{i}{6}(cath_brch_flag)   %  x axis equal
                                            new_node_name = ['n', '1', '_', num2str(wires{i}{3}((max_stress_location(i) - 1)/(nx-1))), '_', num2str(wires{i}{4}((max_stress_location(i) - 1)/(nx-1))+1)];
                                        else
                                            new_node_name = ['n', '1', '_', num2str(wires{i}{3}((max_stress_location(i) - 1)/(nx-1))+1), '_', num2str(wires{i}{4}((max_stress_location(i) - 1)/(nx-1)))];
                                        end
                                        wire_spice{k} = strrep(wire_spice{k}, brchname{2}, new_node_name);
                                    else
                                        if wires{i}{3}(cath_brch_flag) == wires{i}{6}(cath_brch_flag)
                                            new_node_name = ['n', '1', '_', num2str(wires{i}{6}((max_stress_location(i) - 1)/(nx-1))), '_', num2str(wires{i}{7}((max_stress_location(i) - 1)/(nx-1))+1)];
                                        else
                                            new_node_name = ['n', '1', '_', num2str(wires{i}{6}((max_stress_location(i) - 1)/(nx-1))+1), '_', num2str(wires{i}{7}((max_stress_location(i) - 1)/(nx-1)))];
                                        end
                                        wire_spice{k} = strrep(wire_spice{k}, brchname{3}, new_node_name);
                                    end
                                    wire_spice{k} = char(wire_spice{k});
%                                     break;    %  early failure
                                end
                                
                                str = num2str(wires{i}{8}((max_stress_location(i) - 1)/(nx-1)));
                                wire_spice{k} = strrep(wire_spice{k}, brchname{4}, str);
                                wire_spice{k} = char(wire_spice{k});
                                process = process + 1;
                            end
                            if strcmp(brchname{1}, wires{i}{1}((max_stress_location(i) - 1)/(nx-1) + 1))       % compare the branch name whose resistance requires update
                                brchname0 = strrep(wires{i}{1}(1), '-', ' ');
                                netID = textscan(char(brchname0), '%d %d %d');
                                netID = netID{1};
                                cath_brch_flag = (max_stress_location(i) - 1)/(nx-1)+1;
                                if (netID == 25 && J_list_norm{i}((max_stress_location(i) - 1)/(nx-1))>0)
                                    if wires{i}{3}(cath_brch_flag)+wires{i}{4}(cath_brch_flag)-wires{i}{6}(cath_brch_flag)-wires{i}{7}(cath_brch_flag)>0
                                        if wires{i}{3}(cath_brch_flag) == wires{i}{6}(cath_brch_flag)
                                            new_node_name = ['n', '1', '_', num2str(wires{i}{6}((max_stress_location(i) - 1)/(nx-1) + 1)), '_', num2str(wires{i}{7}((max_stress_location(i) - 1)/(nx-1) + 1)+1)];
                                        else
                                            new_node_name = ['n', '1', '_', num2str(wires{i}{6}((max_stress_location(i) - 1)/(nx-1) + 1)+1), '_', num2str(wires{i}{7}((max_stress_location(i) - 1)/(nx-1) + 1))];
                                        end
                                        wire_spice{k} = strrep(wire_spice{k}, brchname{3}, new_node_name);
                                    else
                                        if wires{i}{3}(cath_brch_flag) == wires{i}{6}(cath_brch_flag)
                                            new_node_name = ['n', '1', '_', num2str(wires{i}{3}((max_stress_location(i) - 1)/(nx-1) + 1)), '_', num2str(wires{i}{4}((max_stress_location(i) - 1)/(nx-1) + 1)+1)];
                                        else
                                            new_node_name = ['n', '1', '_', num2str(wires{i}{3}((max_stress_location(i) - 1)/(nx-1) + 1)+1), '_', num2str(wires{i}{4}((max_stress_location(i) - 1)/(nx-1) + 1))];
                                        end
                                        wire_spice{k} = strrep(wire_spice{k}, brchname{2}, new_node_name);
                                    end
                                    wire_spice{k} = char(wire_spice{k});
%                                     break;    %  early failure
                                end
                                
                                str = num2str(wires{i}{8}((max_stress_location(i) - 1)/(nx-1) + 1));
                                wire_spice{k} = strrep(wire_spice{k}, brchname{4}, str);
                                wire_spice{k} = char(wire_spice{k});
                                process = process + 1; 
                                % add node
                                if (netID == 25 && J_list_norm{i}((max_stress_location(i) - 1)/(nx-1))>0)
                                    wire_spice{end-1} = [wire_spice{end-1}, ' v(', new_node_name, ')'];
                                end
                                
                            end
                            if process == 2
                                break;            
                            end
                        end
                    end
                    
                end
                Lvoid_last1(i) = Lvoid1(i);
                Lvoid_last2(i) = Lvoid2(i);
                Lvoid_last(i) = Lvoid(i);
                
            end
             initC_last{i} = initC{i};
        end
        
%         switch i
% %             case 1
% %                 fprintf(fid_resis_6, '%.6f\r\n', sum(wires{i}{8}));
% %             case 2
% %                 fprintf(fid_resis_7, '%.6f\r\n', sum(wires{i}{8}));
% %             case 3
% %                 fprintf(fid_resis_44, '%.6f\r\n', sum(wires{i}{8}));
% %             case 4
% %                 fprintf(fid_resis_45, '%.6f\r\n', sum(wires{i}{8}));
%             case 6 
%                 fprintf(fid_resis_6, '%.6f\r\n', sum(wires{i}{8}));
%             case 7
%                 fprintf(fid_resis_7, '%.6f\r\n', sum(wires{i}{8}));
%             case 38
%                 fprintf(fid_resis_38, '%.6f\r\n', sum(wires{i}{8}));
%             case 39
%                 fprintf(fid_resis_39, '%.6f\r\n', sum(wires{i}{8}));
%             case 44
%                 fprintf(fid_resis_44, '%.6f\r\n', sum(wires{i}{8}));
%             case 45
%                 fprintf(fid_resis_45, '%.6f\r\n', sum(wires{i}{8}));
%         end
        
        %% output stress
        for i0 = 1:brchnum(i)
            fprintf(fid_stress, 'R%s %f\r\n', char(wires{i}{1}(i0)), initC{i}((i0-1)*(nx-1)+1));
            for j0 = 2:(nx-1)
                fprintf(fid_stress, '       %f\r\n', initC{i}((i0-1)*(nx-1)+j0));
            end
        end
        fprintf(fid_stress, '       %f\r\n', initC{i}(a1));
        %
        
        %% output Lvoid
        if Lvoid1(i)+Lvoid2(i)~=0
            fprintf(fid_Lvoid, 'R%s n%s_%d_%d Lvoid: %fum\r\n', char(wires{i}{1}((max_stress_location(i) - 1)/(nx-1))), char(wires{i}{5}((max_stress_location(i) - 1)/(nx-1))), wires{i}{6}((max_stress_location(i) - 1)/(nx-1)), wires{i}{7}((max_stress_location(i) - 1)/(nx-1)), Lvoid1(i)*1e6);
            fprintf(fid_Lvoid, 'R%s n%s_%d_%d Lvoid: %fum\r\n', char(wires{i}{1}((max_stress_location(i) - 1)/(nx-1)+1)), char(wires{i}{2}((max_stress_location(i) - 1)/(nx-1)+1)), wires{i}{3}((max_stress_location(i) - 1)/(nx-1)+1), wires{i}{4}((max_stress_location(i) - 1)/(nx-1)+1), Lvoid2(i)*1e6);
        else
            if max_stress_location(i) == 1
                fprintf(fid_Lvoid, 'R%s n%s_%d_%d Lvoid: %fum\r\n', char(wires{i}{1}(1)), char(wires{i}{2}(1)), wires{i}{3}(1), wires{i}{4}(1), Lvoid(i)*1e6);
            else
                fprintf(fid_Lvoid, 'R%s n%s_%d_%d Lvoid: %fum\r\n', char(wires{i}{1}(end)), char(wires{i}{5}(end)), wires{i}{6}(end), wires{i}{7}(end), Lvoid(i)*1e6);
            end
            
        end
        
    end
    
    %%  output spice

    %% fid = fopen('armcore.sp', 'wt');
    
    if out_spice_flag == 1
        fid2 = fopen('armcore.sp', 'wt');
        for i = 1:line
            fprintf(fid2, '%s\n', wire_spice{i});
        end
        fclose(fid2);
        out_spice_flag = 0;
    end
    
    %% output branch length
    
    fid1 = fopen('tree_leng1.txt', 'wt');
    for i = 1:treenum
        for j = 1:brchnum(i)
            fprintf(fid1, 'L%.10f\r\n', brchleng{i}{1}(j));
        end
    end
    fclose(fid1);
    
    %% close stress, curden, Lvoid file
    fclose(fid_stress);
    fclose(fid_curden);
    fclose(fid_Lvoid);
    fclose(fid_resis);

end

% fclose(fid_resis_6);
% fclose(fid_resis_7);
% fclose(fid_resis_38);
% fclose(fid_resis_39);
% fclose(fid_resis_44);
% fclose(fid_resis_45);

toc