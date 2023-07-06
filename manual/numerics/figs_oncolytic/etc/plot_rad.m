plot_configs;

frad = @(R0,LAMBDA,KAPPA,M0,T)R0*(((LAMBDA+KAPPA)*M0*exp(LAMBDA*T)+LAMBDA*(1-M0)*exp(-KAPPA*M0*T))/(LAMBDA+KAPPA*M0))^(1/3);

r0 = 10; lambda = 0.01; m0 = 0.45; t = 120;

% Format:
% Row 1-6: dt = 0.0375, 0.075, 0.15, 0.3, 0.6, 1.2
% Col 1-5: nx = 05, 10, 20, 40, 80

knames = {'rad_k.txt', 'rad_nok.txt'};
kvals  = {'100.0','0.0'};

nxs = {'05', '10', '20', '40', '80'};
spc = {'r*-', 'y*-', 'm*-', 'c*-', 'b*-'};
dts = {'0.0375', '0.075', '0.15', '0.3', '0.6', '1.2'};

for nt = 1 : length(dts)
  dt(nt) = str2double(dts{nt});
end

for nk = 2 : 2 %1 : length(knames)
  kname = knames{nk};
  kval  = kvals{nk};
  k     = str2double(kval);
  figure('units','normalized','outerposition',[0 0 1 1]);
  hold on;
  data = load(kname);
  rad_ref = frad(r0,lambda,k,m0,t);
  plot( dt, ones(size(dt))*rad_ref, 'k', 'LineWidth', curve_lwdth );
  legs{1} = 'Reference';
  for nn = 1 : length(nxs)
    plot( dt, data(:,nn), spc{nn}, 'LineWidth', curve_lwdth, 'MarkerSize', 16 );
    legs{nn+1} = ['Num = ' nxs{nn}];
  end
  legend( legs, 'Location', 'NW', 'FontSize', legen_fsize );
  xlabel( 'Time step size', 'FontSize', label_fsize );
  ylabel( 'Radius at T=120h', 'FontSize', label_fsize );
  set( gca, 'FontSize', atick_fsize );
  hold off;
end
