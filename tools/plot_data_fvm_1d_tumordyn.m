function plot_data_fvm_1d_tumordyn(varargin)

data_file = varargin{1};
sum_id    = varargin{2};

if nargin>2
  Tterm = varargin{3};
else
  Tterm = -1;
end

if nargin>3
  xrange = varargin{4};
else
  xrange = [1 -1];
end

if nargin>4
  yrange = varargin{5};
else
  yrange = [1 -1];
end

[data_c,data_f,time] = read_data_fvm_1d_tumordyn(data_file);

Nx_c = size(data_c,3);
Nx_f = size(data_f,3);
Nt = length(time);

%pause_sec = 0.1;
pause_sec = 0.01;
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);

if length(sum_id)>1 && max(sum_id) >= 6
  disp('Error: Do not sum up the velocity');
  return;
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

if xrange(1) > xrange(2)
  xrange = [xmin xmax];
end
if yrange(1) > yrange(2)
  yrange = [ymin ymax];
end

D = zeros(1,Nx);
X = zeros(1,Nx);

for t=1:Nt
  if Tterm > 0 & time(t) > Tterm
    break;
  end
  D(1,:) = pdata(t,:);
  if max(sum_id) == 6
    X(1,:) = data_f(t,1,:);
  else
    X(1,:) = data_c(t,1,:);
  end
  plot(X,D,'r-x','LineWidth',1.2,'MarkerSize',8);
  xlim(xrange);
  ylim(yrange);
  title(sprintf('Data, time-stamp %f',time(t)));
  pause(pause_sec);
end

end
