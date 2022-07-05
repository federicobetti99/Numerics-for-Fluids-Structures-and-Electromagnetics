%% importing the data from Free Fem ++
h = fftomatlab_vector('h_vect_mixed.dat');
errorP = fftomatlab_vector('errorP_mixed.dat');
errorL2_mixed = fftomatlab_vector('errorL2_mixed.dat');
errorH1_mixed = fftomatlab_vector('errorH1_mixed.dat');
errorL2_simple = fftomatlab_vector('errorL2_simple.dat');
errorH1_simple = fftomatlab_vector('errorH1_simple.dat');
errorQ1_mixed = fftomatlab_vector('errorQ1_mixed.dat');
errorQ2_mixed = fftomatlab_vector('errorQ2_mixed.dat');
errorQ1_simple = fftomatlab_vector('errorQ1_simple.dat');
errorQ2_simple = fftomatlab_vector('errorQ2_simple.dat');
%% defining possible convergence rates
h_to_one = 4 * h.^(1);
% plotting the results
loglog(h, errorP, '-x',h, h_to_one, '--',  'LineWidth', 1.5);
xlabel('$h$', 'interpreter', 'latex', 'Fontsize', 16);
ylabel('$\vert \vert p-p_h \vert \vert_{L^{2}(\Omega)}$', 'interpreter', 'latex','Fontsize', 16); 
% change to $\vert \vert u-u_h \vert \vert$ for functionals
title('$\lambda = 10^{5}$, $\mu = 8 * 10^5$, $\nu \sim 0.05$', 'interpreter', 'latex')
% change here \lambda properly
legend('mixed form','$h$','Location','Best', 'interpreter', 'latex');
figure
h_square = 0.0000002 * h.^(2);
h_to_one = 0.00004 * h.^(0.8);
h_middle = 0.000004 * h.^(1.5);
h_cube = 0.0000003 * h.^(3);
loglog(h, errorQ1_simple, '-x', h, errorQ2_simple, '-x', ...
    h, errorQ1_mixed, '-x', h, errorQ2_mixed, '-x', h, h_square, '--',  'LineWidth', 1.5);
xlabel('mesh size $h$', 'interpreter', 'latex', 'Fontsize', 16);
ylabel('$\vert Q(u)-Q(u_h) \vert$', 'interpreter', 'latex','Fontsize', 16); 
% change to $\vert Q(u)-Q(u_h) \vert \vert$ for functionals
title('$\lambda = 10^{5}$, $\mu = 8 * 10^5$, $\nu \sim 0.05$', 'interpreter', 'latex')
% change here \lambda properly
legend('simple form $Q_1$ error','simple form $Q_2$ error',...
      'mixed form $Q_1$ error','mixed form $Q_2$ error', ...
      '$h^{2}$','Location','Best', 'interpreter', 'latex');