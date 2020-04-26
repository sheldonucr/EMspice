function [ pwl ] = pwl_create( nt, cdList )
% Creates the PWL function for input current density
% takes input as number of time steps (nt) and the
% current density profile (cdList)
%cdList
tmp = size(cdList);
nj = tmp(2); % number of different current densities

jdur = nt/nj;% duration of each current density
pwl = zeros(1,nt); % Peicewise Linear function of current densities

for t = 0:nj-1
    for d = 1:jdur
        pwl(1,t*jdur+d) = cdList(t+1);
    end
end
%figure(5)
%plot(pwl)


end

