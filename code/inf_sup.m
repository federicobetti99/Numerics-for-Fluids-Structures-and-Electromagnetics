clear; close all; clc;

% NUMERICAL INF-SUP CONDITION TEST : PROCESSING

% READ MATRICES

% Taylor-Hood

B1 = FFmatrix_fread('matB1.txt');
S1 = FFmatrix_fread('matS1.txt');
T1 = FFmatrix_fread('matT1.txt');

B2 = FFmatrix_fread('matB2.txt');
S2 = FFmatrix_fread('matS2.txt');
T2 = FFmatrix_fread('matT2.txt');

B25 = FFmatrix_fread('matB25.txt');
S25 = FFmatrix_fread('matS25.txt');
T25 = FFmatrix_fread('matT25.txt');

B3 = FFmatrix_fread('matB3.txt');
S3 = FFmatrix_fread('matS3.txt');
T3 = FFmatrix_fread('matT3.txt');

% P1c-P0 (unstable)

B1_P1 = FFmatrix_fread('matB1_P1.txt');
S1_P1 = FFmatrix_fread('matS1_P1.txt');
T1_P1 = FFmatrix_fread('matT1_P1.txt');

B2_P1 = FFmatrix_fread('matB2_P1.txt');
S2_P1 = FFmatrix_fread('matS2_P1.txt');
T2_P1 = FFmatrix_fread('matT2_P1.txt');

B25_P1 = FFmatrix_fread('matB25_P1.txt');
S25_P1 = FFmatrix_fread('matS25_P1.txt');
T25_P1 = FFmatrix_fread('matT25_P1.txt');

B3_P1 = FFmatrix_fread('matB3_P1.txt');
S3_P1 = FFmatrix_fread('matS3_P1.txt');
T3_P1 = FFmatrix_fread('matT3_P1.txt');

% MINI element

B1_b = FFmatrix_fread('matB1_b.txt');
S1_b = FFmatrix_fread('matS1_b.txt');
T1_b = FFmatrix_fread('matT1_b.txt');

B2_b = FFmatrix_fread('matB2_b.txt');
S2_b = FFmatrix_fread('matS2_b.txt');
T2_b = FFmatrix_fread('matT2_b.txt');

B25_b = FFmatrix_fread('matB25_b.txt');
S25_b = FFmatrix_fread('matS25_b.txt');
T25_b = FFmatrix_fread('matT25_b.txt');

B3_b = FFmatrix_fread('matB3_b.txt');
S3_b = FFmatrix_fread('matS3_b.txt');
T3_b = FFmatrix_fread('matT3_b.txt');

beta_h_TH = zeros(4,1); % beta of Taylor-Hood
beta_h_b = zeros(4,1); % MINI elements
beta_h_P1 = zeros(4,1); % P1-P0

%% solve the generalized eigenvalue problem
N = 2; % number of eigenvalues to compute if using eigs

% MINI-elements

% calculating the inf-sup constant for the first refinement of the MINI

[r1_b,c1_b] = find(B1_b); % filtering out the zeros
r1_b = unique(r1_b);
c1_b = unique(c1_b);
B1_b = B1_b(r1_b,c1_b);
% assemble GEP problem data
AM1_b = B1_b / T1_b * B1_b';
% solve GEP
[eigvecs1_b, spectrum1_b] = eigs (AM1_b, S1_b, N, 'sm');
spectrum1_b = diag(spectrum1_b);
beta_h_b(1) = sqrt(min(spectrum1_b));

% calculating the inf-sup constant for the second refinement of the MINI

[r2_b,c2_b] = find(B2_b); % filtering out the zeros
r2_b = unique(r2_b);
c2_b = unique(c2_b);
B2_b = B2_b(r2_b,c2_b);
% assemble GEP problem data
AM2_b = B2_b / T2_b * B2_b';
% solve GEP
[eigvecs2_b, spectrum2_b] = eigs (AM2_b, S2_b, N, 'sm');
spectrum2_b = diag(spectrum2_b);
beta_h_b(2) = sqrt(min(spectrum2_b));

% calculating the inf-sup constant for the third refinement of the MINI

[r25_b,c25_b] = find(B25_b); % filtering out the zeros
r25_b = unique(r25_b);
c25_b = unique(c25_b);
B25_b = B25_b(r25_b,c25_b);
% assemble GEP problem data
AM25_b = B25_b / T25_b * B25_b';
% solve GEP
[eigvecs25_b, spectrum25_b] = eigs (AM25_b, S25_b, N, 'sm');
spectrum25_b = diag(spectrum25_b);
beta_h_b(3) = sqrt(min(spectrum25_b));

% calculating the inf-sup constant for the fourth refinement of the MINI

[r3_b,c3_b] = find(B3_b); % filtering out the zeros
r3_b = unique(r3_b);
c3_b = unique(c3_b);
B3_b = B3_b(r3_b,c3_b);
% assemble GEP problem data
AM3_b = B3_b / T3_b * B3_b';
% solve GEP
[eigvecs3_b, spectrum3_b] = eigs (AM3_b, S3_b, N, 'sm');
spectrum3_b = diag(spectrum3_b);
beta_h_b(4) = sqrt(min(spectrum3_b));

% Taylor-Hood

% calculating the inf-sup constant for the first refinement of Taylor-Hood

[r1,c1] = find(B1); % filtering out the zeros
r1 = unique(r1);
c1 = unique(c1);
B1 = B1(r1,c1);
% assemble GEP problem data
AM1 = B1 / T1 * B1';
% solve GEP
[eigvecs1, spectrum1] = eigs (AM1, S1, N, 'sm');
spectrum1 = diag(spectrum1);
beta_h_TH(1) = sqrt(min(spectrum1));

% calculating the inf-sup constant for the second refinement of Taylor-Hood

[r2,c2] = find(B2); % filtering out the zeros
r2 = unique(r2);
c2 = unique(c2);
B2 = B2(r2,c2);
% assemble GEP problem data
AM2 = B2 / T2 * B2'; 
% solve GEP
[eigvecs2, spectrum2] = eigs (AM2, S2, N, 'sm');
spectrum2 = diag(spectrum2);
beta_h_TH(2) = sqrt(min(spectrum2));

% calculating the inf-sup constant for the third refinement of Taylor-Hood

[r25,c25] = find(B25); % filtering out the zeros
r25 = unique(r25);
c25 = unique(c25);
B25 = B25(r25,c25);
% assemble GEP problem data
AM25 = B25 / T25 * B25';
% solve GEP
[eigvecs25, spectrum25] = eigs (AM25, S25, N, 'sm');
spectrum25 = diag(spectrum25);
beta_h_TH(3) = sqrt(min(spectrum25));

% calculating the inf-sup constant for the fourth refinement of Taylor-Hood

[r3,c3] = find(B3); % filtering out the zeros
r3 = unique(r3);
c3 = unique(c3);
B3 = B3(r3,c3);
% assemble GEP problem data
AM3 = B3 / T3 * B3';
% solve GEP
[eigvecs3, spectrum3] = eigs (AM3, S3, N, 'sm');
spectrum3 = diag(spectrum3);
beta_h_TH(4) = sqrt(min(spectrum3));

% P1-P0 (unstable)

% calculating the inf-sup constant for the first refinement of P1-P0

[r1_P1,c1_P1] = find(B1_P1); % filtering out the zeros
r1_P1 = unique(r1_P1);
c1_P1 = unique(c1_P1);
B1_P1 = B1_P1(r1_P1,c1_P1);
% assemble GEP problem data
AM1_P1 = B1_P1 / T1_P1 * B1_P1';
% solve GEP
[eigvecs1_P1, spectrum1_P1] = eigs (AM1_P1, S1_P1, N, 'sm');
spectrum1_P1 = diag(spectrum1_P1);
beta_h_P1(1) = sqrt(min(spectrum1_P1));

% calculating the inf-sup constant for the second refinement of P1-P0

[r2_P1,c2_P1] = find(B2_P1); % filtering out the zeros
r2_P1 = unique(r2_P1);
c2_P1 = unique(c2_P1);
B2_P1 = B2_P1(r2_P1,c2_P1);
% assemble GEP problem data
AM2_P1 = B2_P1 / T2_P1 * B2_P1';
% solve GEP
[eigvecs2_P1, spectrum2_P1] = eigs (AM2_P1, S2_P1, N, 'sm');
spectrum2_P1 = diag(spectrum2_P1);
beta_h_P1(2) = sqrt(min(spectrum2_P1));

% calculating the inf-sup constant for the third refinement of P1-P0

[r25_P1,c25_P1] = find(B25_P1); % filtering out the zeros
r25_P1 = unique(r25_P1);
c25_P1 = unique(c25_P1);
B25_P1 = B25_P1(r25_P1,c25_P1);
% assemble GEP problem data
AM25_P1 = B25_P1 / T25_P1 * B25_P1';
% solve GEP
[eigvecs25_P1, spectrum25_P1] = eigs (AM25_P1, S25_P1, N, 'sm');
spectrum25_P1 = diag(spectrum25_P1);
beta_h_P1(3) = sqrt(min(spectrum25_P1));

% calculating the inf-sup constant for the fourth refinement of P1-P0

[r3_P1,c3_P1] = find(B3_P1); % filtering out the zeros
r3_P1 = unique(r3_P1);
c3_P1 = unique(c3_P1);
B3_P1 = B3_P1(r3_P1,c3_P1);
% assemble GEP problem data
AM3_P1 = B3_P1 / T3_P1 * B3_P1';
% solve GEP
[eigvecs3_P1, spectrum3_P1] = eigs (AM3_P1, S3_P1, N, 'sm');
spectrum3_P1 = diag(spectrum3_P1);
beta_h_P1(4) = sqrt(min(spectrum3_P1));

%% PLOT

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

figure
plot(0:3,beta_h_TH,'-x',0:3,beta_h_P1,'-x','linewidth',1.5)
hold on
plot(0:3,beta_h_b,'-x','linewidth',1.5)
xlabel('Refinements')
ylabel('$\beta_h$')
legend('$P_2^C-P_1^C$','$P_1^C-P_0$','$P_{1}^b-P_1^C$',...
    'location','northwest')
set(gca,'fontsize',18)
xlim ([-0.25 3.25])
