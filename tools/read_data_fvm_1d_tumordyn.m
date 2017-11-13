function [data_c, data_f, time] = read_data_fvm_1d_tumordyn(prob_id, data_file)
% get Nx
fid   = fopen(data_file,'r');
tline = fgetl(fid);
tline = fgetl(fid);
Nx = sscanf(tline,'%*f\t%*f\t%d\t%d');
Nx_c = Nx(1);
Nx_f = Nx(2);
fclose(fid);

% Number of data pieces
Nd = 0;
if prob_id==0
  Nd = 5; % x_c G H N A
end
if prob_id==1
  Nd = 7; % x_c G H M A P1 P2
end

% get the number of time-stamps
extra_data_count = 1;
fid   = fopen(data_file,'r');
n = 0;
while 1
    tline = fgetl(fid);
    n     = n+1;
    if ischar(tline)==0
        break;
    end
end
Nt = round(n/(Nx_c+Nx_f+extra_data_count));
fclose(fid);

% read the data
time   = zeros(1,Nt);
data_c = zeros(Nt,Nd,Nx_c);
data_f = zeros(Nt,2,Nx_f);
fid   = fopen(data_file,'r');
tline = fgetl(fid);
for t = 1:Nt
    tline   = fgetl(fid);
    time(t) = sscanf(tline,'%f\t%*f\t%*d\t%*d');
    tmp_dat = fscanf(fid,'%f',[Nd Nx_c]);
    data_c(t,:,:) = tmp_dat(:,:);
    tmp_dat = fscanf(fid,'%f',[2 Nx_f]);
    tline   = fgetl(fid);
    data_f(t,:,:) = tmp_dat(:,:);
end
fclose(fid);

end
