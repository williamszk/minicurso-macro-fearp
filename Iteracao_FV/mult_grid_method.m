clear
close all
clc

%% This code runs a mult grid method to solve a non-stochasct growth model

%% Parameters

global beta delta alpha theta

beta  = 0.98;
delta = 0.10;
alpha = 0.30;
theta = 2.00;

k_max = 1000;
k_min = 1;

tic

%% grid 1

N_1 = 100;

Grid_1 = linspace(k_min,k_max,N_1);

V_0 = zeros(1,N_1);

[V_1, G_1, C_1] = VF_iteration(Grid_1,V_0);

%% grid 2

N_2 = 1000;

Grid_2 = linspace(k_min,k_max,N_2);

V_0 = spline(Grid_1,V_1,Grid_2);

[V_2, G_2, C_2] = VF_iteration(Grid_2,V_0);

%% grid 3

N_3 = 10000;

Grid_3 = linspace(k_min,k_max,N_3);

V_0 = spline(Grid_2,V_2,Grid_3);

[V_3, G_3, C_3] = VF_iteration(Grid_3,V_0);

toc

%% comparison with no mult grid

tic
V_0 = zeros(1,N_3);
[V_s, G_s, C_s] = VF_iteration(Grid_3,V_0);
toc


%% plots

figure (1)
plot(Grid_3,G_3); hold on; grid on;

figure (2)
plot(Grid_s,G_s); hold on; grid on;

figure (3)
plot(Grid_1,C_1); hold on; grid on;
plot(Grid_2,C_2);
xlabel('capital')
ylabel('consumption')
legend('100 points', '1000 points')
set(gca,'FontSize',20)

figure (4)
plot(Grid_2,C_2); hold on; grid on;
plot(Grid_3,C_3);
xlabel('capital')
ylabel('consumption')
legend('1k points', '10k points')
set(gca,'FontSize',20)

%% dados

% sem acelerador

% it =
% 
%    461
% 
% 
% it =
% 
%    375
% 
% 
% it =
% 
%     82
% 
% Elapsed time is 143.852061 seconds.
% 
% it =
% 
%    454
% 
% Elapsed time is 782.596389 seconds.


