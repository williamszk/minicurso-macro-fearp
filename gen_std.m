function [F] = gen_std(x,sh)
% auxiliary code to Estimation_Simulation_NKPC by Fernando Barros

beta = x(1);
eta = x(2);
global targets pi_0 mc  W

[T,S] = size(sh);

inflation = zeros(T,S);

inflation(1,:) = pi_0*ones(1,S);

for t = 2:T
    for s = 1:S
        inflation(t,s) = 1/(beta)*inflation(t-1,s) - (1-eta)*(1-eta*beta)/(eta*beta)*mc(t-1) + sh(t,s);
    end
end

m1 = 0;
for s = 1:S
    
    m1 = m1 + inflation(2:T,s)'*inflation(1:T-1,s);
end

m1 = m1/(S*T);

m2 = 0;

for s = 1:S
    
    m2 = m2 + inflation(3:T,s)'*inflation(1:T-2,s);
end

m2 = m2/(S*T);

m3 = 0;

for s = 1:S
    
    m3 = m3 + inflation(4:T,s)'*inflation(1:T-3,s);
end

m3 = m3/(S*T);

moments = [m1 m2 m3];

diff = targets - moments;

F = diff*W*diff';


end

