function make_mov_cell_1d_sacado(data_file,mov_id,frate,interval)
% This is specifically designed to read the data
% computed from TDSolveSacado

[data, data_f, time, radius, cols, cols_f] = read_sol_1d_sacado(data_file);
Nd = size(data,2);
Nx = size(data,3);
Nt = length(time);
Nm = length(mov_id);

idx_x  = cols('X_c');
names  = keys(cols);
index  = values(cols);
pdata  = zeros(Nm,Nt,Nx);

clr = {'r','b','m','k','g','c'};

for nm = 1:Nm
  sum_id = mov_id{nm};
  tosum  = zeros(Nt,Nx);
  titem  = zeros(Nt,Nx);
  dname  = '';
  for i = 1:length(sum_id)
    titem(:,:) = data(:,sum_id(i),:);
    tosum = tosum + titem;
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
  pdata(nm,:,:) = tosum(:,:);
  legs{nm} = dname;
end

ymax = max(max(max(pdata)));
ymin = min(min(min(pdata)));
yspn = ymax-ymin;
ymin = ymin - max(0.1*yspn,0.05);
ymax = ymax + max(0.1*yspn,0.40);

xmin = min(min(data(:,idx_x,:)));
xmax = max(max(data(:,idx_x,:)));

mov_title = 'movie_cell.avi';
vObj = VideoWriter(mov_title);
vObj.FrameRate = frate;
open(vObj);

D = zeros(1,Nx);
X = zeros(1,Nx);

Nt_plot = 1:interval:Nt;
if Nt_plot(end)~=Nt
  Nt_plot = [Nt_plot Nt];
end

for nt = 1:length(Nt_plot)
  t = Nt_plot(nt);
  fg = figure(1);
  set(fg,'units','normalized','outerposition',[0 0 1 1]);
  %set(fg,'units','normalized','outerposition',[0 0 .5 .5]);
  hold on;

  for nm = 1:Nm
    D(1,:) = pdata(nm,t,:);
    X(1,:) = data(t,idx_x,:);
    if length(mov_id{nm})>1
      plot(X,D,[clr{mod(nm-1,length(clr))+1} '-'],'LineWidth',2);
    else
      plot(X,D,[clr{mod(nm-1,length(clr))+1} '-x'],'LineWidth',2,'MarkerSize',6);
    end
  end

  legend(legs,'FontSize',24);
  xlim([xmin xmax]);
  ylim([ymin ymax]);
  title(sprintf('Cell numbers, time-stamp %f',time(t)));

  frame = getframe(gcf);
  writeVideo(vObj,frame);

  drawnow;
  clf(fg);

end

close(vObj);

end
