% This code run a SMM, OLS and Indirect Inference estimation for a New
% Keynesian Phillips Curve based on simulated data
% Author: Fernando Barros
% June 07, 2017 

%% Initializing

% cleaning the workspace

clc
clear
close all

% some objects that are common in all m-files

global shocks targets pi_0 mc T S W

% load the data

load('dados_simulation_NKPC.mat')

T = length(mc); % get the length of the data

time = 1:T;

% see the data on scream

figure(1)
plot(time,mc); hold on; grid on;
plot(time,inflation_data);
xlabel('time')
title('inflation and mc data')
legend('mc','inflation')


%% SMM

S = 100; % the number of simulations

% computing moments from data

m1 = inflation_data(2:T)*inflation_data(1:T-1)';
m1 = m1/T;

m2 = inflation_data(3:T)*inflation_data(1:T-2)';
m2 = m2/T;

m3 = inflation_data(4:T)*inflation_data(1:T-3)';
m3 = m3/T;
    
targets = [m1 m2 m3];

% select inital point for simulated data

pi_0 = inflation_data(1);

% generate S sequences of shocks

shocks = normrnd(0,1,[T,S]);

figure(2)
plot(time,shocks(:,1)); hold on; grid on;
%plot(time,shocks(:,floor(S/4)));
plot(time,shocks(:,floor(S/2)));
%plot(time,shocks(:,floor(3*S/4)));
plot(time,shocks(:,S));
xlabel('time')
title('Shocks')
% legend('shock_1','shock_2','shock_3','shock_4','shock_5')
legend('shock_1','shock_2','shock_3')

% weighting matrix

W = [1 0 0;
     0 1 0;
     0 0 1];

% initial guess for optimization 
 
x0 = [3,0.4];

% options for optimzation

options = optimset('Display','iter','TolFun',1e-8,'TolX',1e-8,'MaxIter',1000000,'MaxFunEvals',1000000);

% call the solver

[x, fval] = fminsearch(@obj,x0,options);

beta_SMM = x(1);

eta_SMM = x(2);

min_SMM = fval;

% let's taka a look at the obj function

b_grid = x(1)-0.01:0.001:x(1)+0.01;

e_grid = x(2)-0.01:0.001:x(2)+0.01;

% get standard deviations

N = 100;

[BB,EE] = get_standard_deviation(N,2);

for i=1:length(BB)
    if BB(i) > 10;
        BB(i) = 10;
    end
end

std_beta = std(BB);
std_eta  = std(EE);


figure(9)
histogram(BB,N/10);

figure(10)
histogram(EE,N/10);


% ploting objective function

for b=1:length(b_grid)
        val_obj_SMM_b(b) = obj([b_grid(b),x(2)]);
end

figure(3)
plot(b_grid,val_obj_SMM_b); grid on
xlabel('\beta')
title('SMM - Obj function - varing \beta')

for e=1:length(e_grid)
        val_obj_SMM_e(e) = obj([x(1),e_grid(e)]);
end

figure(4)
plot(e_grid,val_obj_SMM_e); grid on
xlabel('\eta')
title('SMM - Obj function - varing \eta')

for b=1:length(b_grid)
    for e=1:length(e_grid)
        val_obj_SMM(b,e) = obj([b_grid(b),e_grid(e)]);
    end
end

figure(5)
surf(b_grid,e_grid,val_obj_SMM)
xlabel('\beta')
ylabel('\eta')
title('SMM - Surface of Obj function - varing \beta and \eta')

% alternative solver

% [x, fval] = ga(@obj,2,[],[],[],[],[],[],[],options)

%% OLS

Y = inflation_data(2:T)'; X = [inflation_data(1:T-1)' mc(1:T-1)'];

OLS(Y,X,0)

x = OLS(Y,X,0);

estimators_ols = x;

beta_ols = 1/x(1);

aux_f = @(e)(e^(-1)*x(1)*(1-e)*(1-e/x(1)) - x(2));

eta_ols = fsolve(aux_f,0.2);


%% Indirect Inference

% weighting matrix

W_if = eye(2);

obj_if = @(x)([estimators_ols(1)-1/x(1) estimators_ols(2)-((1-x(2))*(1-x(1)*x(2))/(x(1)*x(2)))]*W_if*[estimators_ols(1)-1/x(1); estimators_ols(2)-((1-x(2))*(1-x(1)*x(2))/(x(1)*x(2)))]);

x0 = [1,1];

[x, fval] = fminsearch(obj_if,x0,options);

beta_if = x(1);

eta_if = x(2);

min_if = fval;

% ploting objective function

b_grid = x(1)-0.01:0.001:x(1)+0.01;

e_grid = x(2)-0.01:0.001:x(2)+0.01;

for b=1:length(b_grid)
        val_obj_if_b(b) = obj_if([b_grid(b),x(2)]);
end

figure(6)
plot(b_grid,val_obj_if_b); grid on
xlabel('\beta')
title('I Inf - Obj function - varing \beta')

for e=1:length(e_grid)
        val_obj_if_e(e) = obj_if([x(1),e_grid(e)]);
end

figure(7)
plot(e_grid,val_obj_if_e); grid on
xlabel('\eta')
title('I Inf - Obj function - varing \eta')

for b=1:length(b_grid)
    for e=1:length(e_grid)
        val_obj_if(b,e) = obj_if([b_grid(b),e_grid(e)]);
    end
end

figure(8)
surf(b_grid,e_grid,val_obj_if)
xlabel('\beta')
ylabel('\eta')
title('I Inf - Surface of Obj function - varing \beta and \eta')

%% Display results

   disp(' ')
   disp(' ')
   disp(' ')
   disp(' ')
   
   disp('           ***************************************')
   disp('                           Results        ') 
   disp(' ')
   disp('                    OLS      SMM         IF')
   disp('           ***************************************')
fprintf('           beta  %7.4f   %7.4f    %7.4f ',beta_ols, beta_SMM, beta_if)
disp(' ')
fprintf('            eta  %7.4f   %7.4f    %7.4f ',eta_ols, eta_SMM, eta_if)
disp(' ')
disp('           ***************************************')
 disp(' ')
   disp(' ')
   disp(' ')
   disp(' ')
    disp(' ')
   disp(' ')
   disp(' ')
   disp(' ')
