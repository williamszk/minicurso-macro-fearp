function [vf] = turbo(V,pol,Grid)

global beta delta alpha theta

V_0 = V;

C = Grid.^alpha + ((1-delta).*Grid) - pol;

for i=1:10
    C = smooth(C)';
end

for i=1:10
    vf = C.^(1-theta)./(1-theta) + beta*V_0;
    V_0 = vf;
end

end

