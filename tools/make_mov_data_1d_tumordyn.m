function make_mov_data_1d_tumordyn(data_file,mov_id,frate)

[data,time] = read_data_1d_tumordyn(data_file);

Nx = size(data,3);
Nt = length(time);

Nm = length(mov_id);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);

for nm = 1:Nm
  sum_id = mov_id{nm};
  pdata = zeros(Nt,Nx);
  tosum = zeros(Nt,Nx);
  for i=1:length(sum_id)
    tosum(:,:) = data(:,sum_id(i),:);
    pdata = pdata + tosum;
  end

  ymax = max(max(pdata(1:750,:)));
  ymin = min(min(pdata(1:750,:)));
  %ymax = max(max(pdata(:,:)));
  %ymin = min(min(pdata(:,:)));
  yspn = ymax-ymin;
  ymin = ymin - max(0.1*yspn,0.05);
  ymax = ymax + max(0.1*yspn,0.05);

  xmin = min(min(data(:,1,:)));
  xmax = max(max(data(:,1,:)));
  
  D = zeros(1,Nx);
  X = zeros(1,Nx);

  mov_title = sprintf('movie_%d.avi', nm);
  vObj = VideoWriter(mov_title);
  vObj.FrameRate = frate;
  open(vObj);
  
  for t=1:Nt
      D(1,:) = pdata(t,:);
      X(1,:) = data(t,1,:);
      plot(X,D,'r-x','LineWidth',1.2,'MarkerSize',8);
      xlim([xmin xmax]);
      ylim([ymin ymax]);
      title(sprintf('Data, time-stamp %f',time(t)));
      
      frame = getframe(gcf);
      writeVideo(vObj,frame);
  end

  close(vObj);
end


end
