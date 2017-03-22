function [data, time] = read_data_1d_tumordyn(data_file)
% get Nx
fid   = fopen(data_file,'r');
tline = fgetl(fid);
tline = fgetl(fid);
Nx    = sscanf(tline,'%*f\t%*f\t%d');
fclose(fid);

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
Nt = round(n/(Nx+extra_data_count));
fclose(fid);

% read the data
time  = zeros(1,Nt);
data  = zeros(Nt,6,Nx);
fid   = fopen(data_file,'r');
tline = fgetl(fid);
for t = 1:Nt
    tline   = fgetl(fid);
    time(t) = sscanf(tline,'%f\t%*f\t%*d');
    tmp_dat = fscanf(fid,'%f',[6 Nx]);
    tline   = fgetl(fid);
    data(t,:,:) = tmp_dat(:,:);
end
fclose(fid);

end
