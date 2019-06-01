function [ F ] = simulation_moments(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global targets shocks z_0 K_0 W

alpha = x(1);
beta  = x(2);
rho   = x(3);

[T,S] = size(shocks);

Z = zeros(T,S); K = Z; C = Z;

Z(1,:) = ones(1,S)*z_0;
K(1,:) = ones(1,S)*K_0;
C(1,:) = (1-alpha*beta).*(K(1,:).^alpha).*exp(Z(1,:));

for t=2:T
    Z(t,:) = rho.*Z(t-1,:) + shocks(t,:);
    K(t,:) = beta.*alpha.*exp(Z(t-1,:)).*(K(t-1,:).^alpha);
    C(t,:) = (1-alpha*beta).*(K(t,:).^alpha).*exp(Z(t,:));
end

Y = K.^alpha.*exp(Z);

DC = C(2:T,:)./C(1:T-1,:) - 1;

m1 = sum(sum(DC))/((T-1)*S);

DY = C(2:T,:)./C(1:T-1,:) - 1;

m2 = sum(sum(DY))/((T-1)*S);

Ratio = C./Y;

m3 = sum(sum(Ratio))/(T*S);

moments = [m1;m2;m3];

g = targets-moments;

F = g'*W*g;

end

