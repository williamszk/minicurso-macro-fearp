function [V1,Pol,C] = VF_iteration(Grid,V)

global beta delta alpha theta

N = length(Grid);

tol = 0.0001;
itmax = 10000;
it = 0;
dist = 10;
TV = V;
Pol = V;

while (dist>tol && it<itmax)
    it = it + 1;
    for i=1:N;
        C = Grid(i)^alpha+((1-delta)*Grid(i))-Grid;
        ind = find(C<0);
        C(ind) = eps; 
        aux = C.^(1-theta)./(1-theta) + beta*V;
        [TV(i),pos] = max(aux);
        Pol(i) = Grid(pos);
    end
    dist = abs(max(TV-V));
    V = TV;
    
%     if it>15 && N > 100
%         V = turbo(V,Pol,Grid);
%     end
end

C = Grid.^alpha + ((1-delta).*Grid) - Pol;

it

V1 = V;


end

