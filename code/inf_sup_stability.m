B = FFmatrix_fread ("matB.txt");
T = FFmatrix_fread ("matT.txt");
S = FFmatrix_fread ("matS.txt");

N = 2;
[r,c] = find(B);
r = unique(r);
c = unique(c);
B = B(r,c);

AM = B / T * B'; % operator / is matrix right division
[ eigvecs , spectrum ] = eigs (AM , S, N, "sm");
spectrum = diag(spectrum);
beta_h = sqrt (min(spectrum));
