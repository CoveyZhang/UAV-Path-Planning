function createfigure1(xdata1, zdata1, X1, Y1, Z1)
%CREATEFIGURE1(XDATA1, ZDATA1, X1, Y1, Z1)
%  XDATA1:  surface xdata
%  ZDATA1:  surface zdata
%  X1:  stem3 x
%  Y1:  stem3 y
%  Z1:  stem3 z

%  由 MATLAB 于 13-May-2018 13:58:25 自动生成

% 创建 figure
figure;
colormap('jet');

% 创建 axes
axes1 = axes;
hold(axes1,'on');

% 创建 mesh
mesh(xdata1,xdata1,zdata1);

% 创建 stem3
stem3(X1,Y1,Z1,'MarkerFaceColor','auto','Color',[0 0 0]);

% 创建 plot3
plot3(X1,Y1,Z1,'Color',[0 0 0]);

% 创建 mesh
mesh(xdata1,xdata1,zdata1);

% 创建 plot3
plot3(X1,Y1,Z1,'Color',[0 0 0]);

% 创建 xlabel
xlabel('x/m');

% 创建 zlabel
zlabel('z/m');

% 创建 ylabel
ylabel('y/m');

% 取消以下行的注释以保留坐标轴的 X 范围
% xlim(axes1,[0 100]);
% 取消以下行的注释以保留坐标轴的 Y 范围
% ylim(axes1,[0 100]);
% 取消以下行的注释以保留坐标轴的 Z 范围
% zlim(axes1,[150 1100]);
view(axes1,[-37.5 30]);
