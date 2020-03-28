function[results] = backEuler_fdm(G,C,B,tstamps,ini,u)
% Solver for the FDM matrix
%dim=size(Gh,1);
x = ini;

dt=tstamps(3)-tstamps(2);
tmax = length(tstamps);


for t=2:tmax
   x(:,t) = (G + C/(dt)) \ (B*u(t) + (C/(dt))*x(:,t-1)); 
end

results = x';

end

