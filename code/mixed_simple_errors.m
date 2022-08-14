clear; close all; clc;

% produces the plots for the Taylor-Hood elements for the L2 error 
% and the H1 error on the displacement, the L2 error of the pressure and 
% the errors on the functionals Q1, Q2. The same code was used also for
% producing the plots of the MINI-element.

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% READ ERRORS

% simple form

h_vect_s = fftomatlab_vector('h_vect_simple.dat');
err_L2_s = fftomatlab_vector('errorL2_simple.dat');
err_H1_s = fftomatlab_vector('errorH1_simple.dat');
err_Q1_s = fftomatlab_vector('errorQ1_simple.dat');
err_Q2_s = fftomatlab_vector('errorQ2_simple.dat');

% mixed form

h_vect_m = fftomatlab_vector('h_vect_mixed.dat');
err_L2_m = fftomatlab_vector('errorL2_mixed.dat');
err_H1_m = fftomatlab_vector('errorH1_mixed.dat');
err_P_m = fftomatlab_vector('errorP_mixed.dat');
err_Q1_m = fftomatlab_vector('errorQ1_mixed.dat');
err_Q2_m = fftomatlab_vector('errorQ2_mixed.dat');

% h_vect_s and h_vect_m are the same.

%% PLOTS : errors and orders of convergence

% L2 and H1 errors of the displacement

figure
loglog(h_vect_s,8e-5*h_vect_s.^(0.7),'--','linewidth',1.5)
hold on
loglog(h_vect_s,err_H1_s,'-x','linewidth',1.5)
loglog(h_vect_m,err_H1_m,'-x','linewidth',1.5)
loglog(h_vect_s,err_L2_s,'-x','linewidth',1.5)
loglog(h_vect_m,err_L2_m,'-x','linewidth',1.5)
loglog(h_vect_s,0.09e-4*h_vect_s.^(1.4),'--','linewidth',1.5)
set(gca,'fontsize',18)
title('Taylor-Hood, $\lambda = 10^{17}$, $\nu \sim 0.5$')
xlabel('h')
ylabel('$\vert \vert u - u_{h} \vert \vert$')
legend('$h^{0.7}$'...,
    ,'simple form, $H^1$ error'...
    ,'mixed form, $H^1$ error'...
    ,'simple form, $L^2$ error'...
    ,'mixed form, $L^2$ error'...
    ,'$h^{1.4}$'...
    ,'location','best')

% L2 error on the pressure

figure
loglog(h_vect_s,err_P_m,'-x','linewidth',1.5)
hold on
loglog(h_vect_s,5e1*h_vect_s.^(0.6),'--','linewidth',1.5)
set(gca,'fontsize',18)
title('Taylor-Hood, $\lambda = 10^{17}$, $\nu \sim 0.5$')
xlabel('h')
ylabel('$\vert \vert p - p_{h} \vert \vert_{L^2(\Omega)}$')
legend('mixed form','$h^{0.6}$','location','best')

% errors on the functionals

figure
loglog(h_vect_s,err_Q2_s,'-x','linewidth',1.5)
hold on
loglog(h_vect_s,err_Q1_s,'-x','linewidth',1.5)
loglog(h_vect_m,err_Q2_m,'-x','linewidth',1.5)
loglog(h_vect_m,err_Q1_m,'-x','linewidth',1.5)
loglog(h_vect_s,0.01e-4*h_vect_s.^(1.5),'--','linewidth',1.5)
set(gca,'fontsize',18)
title('Taylor-Hood, $\lambda = 10^{17}$, $\nu \sim 0.5$')
xlabel('h')
ylabel('$\vert Q(u) - Q(u_h)\vert$')
legend('simple form $Q_2$'...
    ,'simple form $Q_1$'...
    ,'mixed form $Q_2$'...
    ,'mixed form $Q_1$'...
    ,'$h^{1.5}$'...
    ,'location','best')
