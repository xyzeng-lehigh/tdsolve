function make_mov_data_fvm_1d_tumordyn(data_file,mov_id,frate)

[data_c,data_f,time] = read_data_fvm_1d_tumordyn(data_file);

Nx_c = size(data_c,3);
Nx_f = size(data_f,3);
Nt = length(time);

Nm = length(mov_id);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);

for nm = 1:Nm
  sum_id = mov_id{nm};

  if length(sum_id)>1 && max(sum_id) >= 6
    disp('Error: Do not sum up the velocity');
    continue;
  end
  
  if max(sum_id) == 6
    pdata = zeros(Nt,Nx_f);
    pdata(:,:) = data_f(:,2,:);
    Nx = Nx_f;
  else
    pdata = zeros(Nt,Nx_c);
    tosum = zeros(Nt,Nx_c);
    for i=1:length(sum_id)
      tosum(:,:) = data_c(:,sum_id(i),:);
      pdata = pdata + tosum;
    end
    Nx = Nx_c;
  end
  
  ymax = max(max(pdata(:,:)));
  ymin = min(min(pdata(:,:)));
  yspn = ymax-ymin;
  ymin = ymin - max(0.1*yspn,0.05);
  ymax = ymax + max(0.1*yspn,0.05);
  
  xmin = min(min(data_f(:,1,:)));
  xmax = max(max(data_f(:,1,:)));
  
  D = zeros(1,Nx);
  X = zeros(1,Nx);

  mov_title = sprintf('movie_%d.avi', nm);
  vObj = VideoWriter(mov_title);
  vObj.FrameRate = frate;
  open(vObj);
  
  for t=1:Nt
      D(1,:) = pdata(t,:);
      if max(sum_id) == 6
        X(1,:) = data_f(t,1,:);
      else
        X(1,:) = data_c(t,1,:);
      end
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
