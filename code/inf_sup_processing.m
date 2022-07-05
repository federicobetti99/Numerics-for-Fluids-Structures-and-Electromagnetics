%% importing the matrices

B1 = FFmatrix_fread(':\Users\Federico Betti\Documents\Numerics for Fluids, Structures and Electromagnetism\matB1.txt');
S1 = FFmatrix_fread(':\Users\Federico Betti\Documents\Numerics for Fluids, Structures and Electromagnetism\matS1.txt');
T1 = FFmatrix_fread(':\Users\Federico Betti\Documents\Numerics for Fluids, Structures and Electromagnetism\matT1.txt');

%% calculating the beta_h value

N = 2; % number of eigenvalues to compute if using eigs
[r,c] = find(B1); % filtering out the zeros
r = unique(r);
c = unique(c);
B1 = B1(r,c);
AM = B1 / T1 * B1'; % operator / is matrix right division
[eigvecs,spectrum] = eigs(AM, S1, N, 'sm');
spectrum = diag(spectrum);
beta_h_1 = sqrt(min(spectrum));

%% plotting

beta_values = zeros(6,1);
h_values = zeros(6,1);
beta_values(1) = 0.2562;
beta_values(2) = 0.2454;
beta_values(3) = 0.2426;
beta_values(4) = 0.2412;
h_values(1) = 2.8512;
h_values(2) = 1.42565;
h_values(3) = 0.791485;
h_values(4) = 0.388588;

plot([1,2,3,4], beta_values)
label