function [data_c, data_f, time, radius, cols_c, cols_f] = read_sol_1d_sacado(data_file)
% This is specifically designed to read the data
% computed from TDSolveSacado

% Get variable names
fid   = fopen(data_file,'r');
tline = fgetl(fid);
words = strsplit(tline);
newl  = 0;
var_c = 0;
var_f = 0;
name_c = {}; name_f = {}; indx_c = {}; indx_f = {};
for i = 1:length(words)
  is_newl = strcmp(words{i},'/newline');
  if is_newl ~= 1
    if newl == 1
      var_c = var_c + 1;
      name_c{var_c} = words{i};
      indx_c{var_c} = var_c;
    end
    if newl == 2
      var_f = var_f + 1;
      name_f{var_f} = words{i};
      indx_f{var_f} = var_f;
    end
  end
  newl = newl + is_newl;
end
cols_c = containers.Map(name_c,indx_c);
cols_f = containers.Map(name_f,indx_f);
Nd_c   = length(name_c);
Nd_f   = length(name_f);
% Get Nx
tline = fgetl(fid);
Nx = sscanf(tline,'%*f\t%*f\t%*f\t%d\t%d');
Nx_c = Nx(1);
Nx_f = Nx(2);
fclose(fid);

% Get the number of time stamps
extra_data_count = 1;
fid = fopen(data_file,'r');
n = 0;
while 1
  tline = fgetl(fid);
  n = n+1;
  if ischar(tline)==0
    break;
  end
end
Nt = round(n/(Nx_c+Nx_f+extra_data_count));
fclose(fid);

% Read the data
time   = zeros(1,Nt);
radius = zeros(1,Nt);
data_c = zeros(Nt,Nd_c,Nx_c);
data_f = zeros(Nt,Nd_f,Nx_f);
fid   = fopen(data_file,'r');
tline = fgetl(fid);
for t = 1:Nt
    tline     = fgetl(fid);
    scanvar   = sscanf(tline,'%f\t%*f\t%f\t%*d\t%*d');
    time(t)   = scanvar(1);
    radius(t) = scanvar(2);
    tmp_dat   = fscanf(fid,'%f',[Nd_c Nx_c]);
    data_c(t,:,:) = tmp_dat(:,:);
    tmp_dat = fscanf(fid,'%f',[Nd_f Nx_f]);
    tline   = fgetl(fid);
    data_f(t,:,:) = tmp_dat(:,:);
end
fclose(fid);

end
