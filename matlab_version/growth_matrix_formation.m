function [G,C,B,nx_total,w_lengs,J_list] = growth_matrix_formation(nx0,wires0,J_list_norm0)


% % Function for EM Discrete Matrix formation 
% % parameters: 
% % L : length of wire
% % nx : Number of sample points along L
% % ini_condition: Vector of initional conditions
% % t : time to simulate
% % nt : size of time steps
% % T  : Temperature
% % j  : current density @ boundary 

% %Physical parameters
T = 373;% Temperature in Kelvin
D0=7.8e-5;%maximum diffusion coefficient (at infinite temperature; m2 / s),
Ea=1.609e-19;%activation energy for diffusion in dimensions of (J atom�?)
kB=1.3806e-23;%Boltzmann constant.
Da=D0*exp(-Ea/(kB*T));%diffusion coefficient (m2 / s)
B=1e11;%effective bulk elasticity modulus
% %Omega = .95;
Omega=1.182e-29;%atomic lattice volume
E=1.609e-19;
cu_res=1.7e-8;
Z=10;
% %L = 80e-6;
% 
% %kappa=(Da*B*Omega)/(kB*T);
% %G=(E*Z*cu_res*j)/(Omega);

wires=wires0;
J_list_norm=J_list_norm0;

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % wires{1} holds the name of a wire
% % % wires{2} holds the name of the net the wire is associated with
% % % wires{3} holds the x value of starting nodes
% % % wires{4} holds the y value of starting nodes
% % % wires{5} holds the names of the net the wires are assocaited with
% % % wires{6} holds the x valye of destination nodes
% % % wires{7} holds the y value of destination nodes
% % % wires{8} Resistance of wires
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_wires = numel(wires{1}); %number of wire segments in tree
w_lengs = zeros(n_wires,2); %will hold lengths of each segment

for n = 1:n_wires
    xL = abs(wires{6}(n) - wires{3}(n));
    if (xL > 0)
        w_lengs(n,1) = xL;
        w_lengs(n,2) = 1;
    else
        yL = abs(wires{7}(n) - wires{4}(n));
        if(yL > 0)
            w_lengs(n,1) = yL;
            w_lengs(n,2) = 2;
        end
    end
end

ys = wires{7};
min(ys);
yn=ys-min(ys);
w_lengs;
num_Junctions = n_wires-1;

nx = nx0;
nx_total = (nx*n_wires)-(num_Junctions);
del = 3e-7;

%%***********************
dx=1/(nx_total-1);

G = zeros(nx_total);
B = zeros(nx_total,n_wires);

for i = 2:nx_total-1
    G(i,i-1)=G(i,i-1)+1;
    G(i,i) = G(i,i) -2;
    G(i,i+1) = G(i,i+1)+1;
end

%%*************************************
% G(1,1) = G(1,1)*(-1-(dx/del));

G(1,1) = G(1,1) -1;
G(1,2) = G(1,2) +1;

G(end,end) = G(end,end) -1;
G(end,end-1) = G(end,end-1) +1;

G(1,1) = G(1,1)*(-1-(dx/del));

% msl=max_stress_location0;
% G(msl,msl)=G(msl,msl)*(-1-(dx/del));

index = 1;
t=zeros(1,n_wires);
% %tic
for i = 1:n_wires
    B(index,i)=1;
    B(index+nx-1,i)=-1;
    index = index + nx -1;
%     %t(i)=toc;
end


C=ones([nx_total,1]);
B=B*J_list_norm/dx;
B(1) = 0;
% B(msl) = 0;
G=-G;
G=G/(dx^2);

% %plot(t)
% %t(end)
% %%% OLD IMPLEMENTATION%%%%%
% % %%%% file IO and parameter extraction %%%%
% % fileID = fopen(filename,'r');
% % nets = textscan(fileID, '%s %f %f %f %f %f')
% % %matrix_out = null;
% % n_nets = numel(nets{1});
% % 
% 
% % %%  Grid generation %%
% % %min_x1 = min(nets{2});
% % %min_x2 = min(nets{4});
% % x_vals = [nets{2};nets{4}];
% % %min_x = min(min_x1,min_x2);
% % xL = range(x_vals)+1 % row size
% % %min_y1 = min(nets{3});
% % %min_y2 = min(nets{5});
% % y_vals = [nets{3};nets{5}]; 
% % %min_y = min(min_y1,min_y2);
% % yL = range(y_vals)+1;% Col size, or number of rows
% % area = xL*yL; % Total grid size
% % 
% % 
% % G = zeros(area,area)
% % % for net = 1:n_nets
% % %     start_node = [nets{2}(net),nets{4}(net)];
% % %     desti_node = [nets{2}(net),nets{4}(net)];
% % %     G(
% % %     

     
end
