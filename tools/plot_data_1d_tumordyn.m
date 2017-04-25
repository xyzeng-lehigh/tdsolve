function plot_data_1d_tumordyn(data_file,sum_id)

[data,time] = read_data_1d_tumordyn(data_file);

Nx = size(data,3);
Nt = length(time);

%pause_sec = 0.1;
pause_sec = 0.01;
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);

pdata = zeros(Nt,Nx);
tosum = zeros(Nt,Nx);
for i=1:length(sum_id)
  tosum(:,:) = data(:,sum_id(i),:);
  pdata = pdata + tosum;
end

ymax = max(max(pdata(:,:)));
ymin = min(min(pdata(:,:)));
yspn = ymax-ymin;
ymin = ymin - max(0.1*yspn,0.05);
ymax = ymax + max(0.1*yspn,0.05);

xmin = min(min(data(:,1,:)));
xmax = max(max(data(:,1,:)));

D = zeros(1,Nx);
X = zeros(1,Nx);

for t=1:Nt
    D(1,:) = pdata(t,:);
    X(1,:) = data(t,1,:);
    plot(X,D,'r-x','LineWidth',1.2,'MarkerSize',8);
    xlim([xmin xmax]);
    %ylim([ymin ymax]);
    ylim([-.1 1.1]);
    title(sprintf('Data, time-stamp %f',time(t)));
    pause(pause_sec);
end

end
