function v=fftomatlab_vector(file)

% v=fftomatlab_vector(file)
%
% reads a freefem vector and imports it in matlab 

fid=fopen(file);
% get the number of elements: FF puts at most 5 numbers on a line
el=str2num(fgetl(fid));
% the number of lines to be read
L=fix(el/5)+1;
l=1; % line counter
v=[];
while l <= L
    t=fgets(fid);
    v =[v str2num(t)];
    l=l+1;
end
