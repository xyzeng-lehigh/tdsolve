function make_mov_sol_1d_sacado(data_file,mov_id_c,mov_id_f,frate,interval,varargin)
% This is specifically designed to read the data
% computed from TDSolveSacado
% make_mov_sol_1d_sacado(data_file,mov_id_c,mov_id_f,frate,interval,varargin)
% varargin{1} = opt
% opt{1} : ylim for _id_c
% opt{2} : ylim for _id_f
% opt{3} : xlabel
% opt{4} : ylabel for _id_c
% opt{5} : ylabel for _id_f
% opt{6} : Font size for title/label
% opt{7} : Font size for ticks

if nargin == 5
  opt = {};
else
  opt = varargin{1};
end

[data_c, data_f, time, radius, cols_c, cols_f] = read_sol_1d_sacado(data_file);
Nd_c = size(data_c,2);
Nx_c = size(data_c,3);
Nd_f = size(data_f,2);
Nx_f = size(data_f,3);
Nt   = length(time);
Nm_c = length(mov_id_c);
Nm_f = length(mov_id_f);

figure('units','normalized','outerposition',[0 0 1 1]);
%figure('units','normalized','outerposition',[0 0 .5 .5]);

Nt_plot = 1:interval:Nt;
if Nt_plot(end)~=Nt
  Nt_plot = [Nt_plot Nt];
end

for nm = 1:(Nm_c+Nm_f)
  if nm <= Nm_c
    % Cell-centered data
    sum_id = mov_id_c{nm};
    Nx     = Nx_c;
    Nd     = Nd_c;
    data   = data_c;
    idx_x  = cols_c('X_c');
    cols   = cols_c;
    names  = keys(cols);
    index  = values(cols);
    if isempty(opt)
      opt_ylim = {};
      opt_ylab = {};
    else
      opt_ylim = opt{1}{nm};
      opt_ylab = opt{4}{nm};
    end
  else
    % Node-centered data
    sum_id = mov_id_f{nm-Nm_c};
    Nx     = Nx_f;
    Nd     = Nd_f;
    data   = data_f;
    idx_x  = cols_f('X_f'); 
    cols   = cols_f;
    names  = keys(cols);
    index  = values(cols);
    if isempty(opt)
      opt_ylim = {};
      opt_ylab = {};
    else
      opt_ylim = opt{2}{nm-Nm_c};
      opt_ylab = opt{5}{nm-Nm_c};
    end
  end
  pdata = zeros(Nt,Nx);
  tosum = zeros(Nt,Nx);
  dname = '';
  for i = 1:length(sum_id)
    tosum(:,:) = data(:,sum_id(i),:);
    pdata = pdata + tosum;
    for j = 1:Nd
      if index{j} == sum_id(i)
        cname = names{j};
      end
    end
    dname = strcat(dname,cname);
    if i < length(sum_id)
      dname = strcat(dname,'+');
    end
  end
  
  ymax = max(max(pdata(:,:)));
  ymin = min(min(pdata(:,:)));
  yspn = ymax-ymin;
  ymin = ymin - max(0.1*yspn,0.05);
  ymax = ymax + max(0.1*yspn,0.05);
  
  xmin = min(min(data(:,idx_x,:)));
  xmax = max(max(data(:,idx_x,:)));
  
  D = zeros(1,Nx);
  X = zeros(1,Nx);

  mov_title = sprintf('movie_%d.mp4', nm);
  vObj = VideoWriter(mov_title, 'MPEG-4');
  vObj.FrameRate = frate;
  open(vObj);

  for nt=1:length(Nt_plot)
    t = Nt_plot(nt);
    D(1,:) = pdata(t,:);
    X(1,:) = data(t,idx_x,:);
  
    plot(X,D,'r-x','LineWidth',1.2,'MarkerSize',8);
    xlim([xmin xmax]);
    if isempty(opt_ylim)
      ylim([ymin ymax]);
    else
      ylim(opt_ylim);
    end
    if ~(isempty(opt) || isempty(opt{3}))
      if isempty(opt) || isempty(opt{6})
        xlabel(opt{3});
      else
        xlabel(opt{3},'FontSize',opt{6});
      end
    end
    if ~isempty(opt_ylab)
      if isempty(opt) || isempty(opt{6})
        ylabel(opt_ylab);
      else
        ylabel(opt_ylab,'FontSize',opt{6});
      end
    end
    if isempty(opt) || isempty(opt{6})
      title([dname sprintf(', time-stamp %f',time(t))]);
    else
      title([dname sprintf(', time-stamp %f',time(t))],'FontSize',opt{6});
    end
    if ~(isempty(opt) || isempty(opt{7}))
      set(gca, 'FontSize', opt{7});
    end

    frame = getframe(gcf);
    writeVideo(vObj,frame);
  end

  close(vObj);
end

end
