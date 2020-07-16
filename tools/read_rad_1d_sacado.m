function [data, cols] = read_rad_1d_sacado(data_file)
% This is specifically designed to read the data
% computed from TDSolveSacado
% Format:
% tsol, data_r, data_v[Nx], data_u[Nx]

name = {'T', 'R', 'V', 'U'};
indx = {1, 2, 3, 4};

cols = containers.Map(name,indx);

data = load(data_file);

end
