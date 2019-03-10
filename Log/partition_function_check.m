k = 0.002510;

dz = (0:0.001:2000);
[foo N_dz] = size(dz);
E = zeros(N_dz,1);


for(i=1:1:N_dz)
  E(i,1) = 0.5*k*dz(i)*dz(i);
end
% plot(dz, E)
p_E = exp(-E/4.1);
%plot(dz, p_E)
trapz(dz, p_E)

%ans =

%   46.7004
 
