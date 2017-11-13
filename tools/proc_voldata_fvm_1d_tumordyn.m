function data = proc_voldata_fvm_1d_tumordyn(prob_id,data_file)

[data_c, data_f, time] = read_data_fvm_1d_tumordyn(prob_id,data_file);
Nt = length(time);

Nx = size(data_c,3);

if prob_id==0 || prob_id==1
  data = zeros(Nt,8); % PDGF: time tot.G tot.H tot.N tot perc.G perc.H perc.N
                      % OPN:  time tot.G tot.H tot.M tot perc.G perc.H perc.M
  for t=1:Nt
    data(t,1) = time(t);
    tot_G = 0.0;
    tot_H = 0.0;
    tot_N = 0.0; % tot_M
    dr = data_f(t,1,2);
    for n=1:Nx
      loc_r = data_c(t,1,n);
      loc_g = data_c(t,2,n);
      loc_h = data_c(t,3,n);
      loc_n = data_c(t,4,n); % loc_m
  
      tot_G = tot_G + 4.0 * pi * dr * loc_r * loc_r * loc_g;
      tot_H = tot_H + 4.0 * pi * dr * loc_r * loc_r * loc_h;
      tot_N = tot_N + 4.0 * pi * dr * loc_r * loc_r * loc_n; % tot_M
    end
    data(t,2) = tot_G;
    data(t,3) = tot_H;
    data(t,4) = tot_N; % tot_M
    data(t,5) = tot_G + tot_H + tot_N; % tot_G + tot_H + tot_M
    data(t,6) = tot_G/data(t,5);
    data(t,7) = tot_H/data(t,5);
    data(t,8) = tot_N/data(t,5); % tot_M/data(t,5);
  end
end

end
