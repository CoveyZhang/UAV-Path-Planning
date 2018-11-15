function createfigure1(xdata1, zdata1, X1, Y1, Z1)
%CREATEFIGURE1(XDATA1, ZDATA1, X1, Y1, Z1)
%  XDATA1:  surface xdata
%  ZDATA1:  surface zdata
%  X1:  stem3 x
%  Y1:  stem3 y
%  Z1:  stem3 z

%  �� MATLAB �� 13-May-2018 13:58:25 �Զ�����

% ���� figure
figure;
colormap('jet');

% ���� axes
axes1 = axes;
hold(axes1,'on');

% ���� mesh
mesh(xdata1,xdata1,zdata1);

% ���� stem3
stem3(X1,Y1,Z1,'MarkerFaceColor','auto','Color',[0 0 0]);

% ���� plot3
plot3(X1,Y1,Z1,'Color',[0 0 0]);

% ���� mesh
mesh(xdata1,xdata1,zdata1);

% ���� plot3
plot3(X1,Y1,Z1,'Color',[0 0 0]);

% ���� xlabel
xlabel('x/m');

% ���� zlabel
zlabel('z/m');

% ���� ylabel
ylabel('y/m');

% ȡ�������е�ע���Ա���������� X ��Χ
% xlim(axes1,[0 100]);
% ȡ�������е�ע���Ա���������� Y ��Χ
% ylim(axes1,[0 100]);
% ȡ�������е�ע���Ա���������� Z ��Χ
% zlim(axes1,[150 1100]);
view(axes1,[-37.5 30]);
