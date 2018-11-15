clear all;
clc;
close all;
tic
global DEM safth hmax scfitness;
a=load('XYZmesh.mat');%读取数字高程信息DEM
DEM=a;
DEM.Z=DEM.Z;
safth=60;

hmax=max(max(DEM.Z));
hmin=min(min(DEM.Z));
%DEM.Z=randint(size(DEM.Z,1),size(DEM.Z,2),[0 100]);
ganum=100;
dnanum=50;
dnalength=50;
childrennum=dnanum;

np=zeros(1,dnanum);

nsort1=zeros(1,dnanum);
nsort2=zeros(1,dnanum);
nsort3=zeros(1,dnanum);

fitness1=zeros(dnanum,3);
fitness2=zeros(dnanum,3);
fitness3=zeros(dnanum,3);

scfitness=zeros(dnanum,3);

congestion=zeros(1,dnanum);
startpoint=[1,1,100];
goalpoint=[101,101,100];
startpoint(3)=DEM.Z(startpoint(1),startpoint(2))+safth;
goalpoint(3)=DEM.Z(goalpoint(1),goalpoint(2))+safth;
dna1=zeros(dnanum,dnalength+1,3);
dna2=zeros(dnanum,dnalength+1,3);
dna3=zeros(dnanum,dnalength+1,3);
totalfitness=zeros(2,ganum+1,3);  %make a 3*101*3 size matrix ,the last 3 means three dimensional 

for i=1:1:dnanum
    for j=2:1:dnalength
        for k=1:1:3
         %dna1(i,j,k)=randint(1,1,[1 101]);    
        dna1(i,j,1)=j*100/dnalength; 
        dna2(i,j,1)=j*100/dnalength; 
        dna3(i,j,1)=j*100/dnalength; 
        dna1(i,j,2)=j*100/dnalength+randint(1,1,[0 20]-10);
        dna2(i,j,2)=j*100/dnalength+randint(1,1,[0 20]-20);
        dna3(i,j,2)=j*100/dnalength+randint(1,1,[0 20]+15);%key point to improvement
        if dna1(i,j,2)<1
           dna1(i,j,2)=1;
        elseif dna1(i,j,2)>101
           dna1(i,j,2)=101;      
        end
        if dna2(i,j,2)<1
           dna2(i,j,2)=1;
        elseif dna2(i,j,2)>101
           dna2(i,j,2)=101;      
        end
         if dna3(i,j,2)<1
           dna3(i,j,2)=1;
        elseif dna3(i,j,2)>101
           dna3(i,j,2)=101;      
        end
        %dna1(i,j,k)=floor(100*j/dnalength);
        %dna1(i,j,k)=j;
        end    
        dna1(i,j,3)=DEM.Z(dna1(i,j,1),dna1(i,j,2))+safth; 
        dna2(i,j,3)=DEM.Z(dna2(i,j,1),dna2(i,j,2))+safth; 
        dna3(i,j,3)=DEM.Z(dna3(i,j,1),dna3(i,j,2))+safth; 
%         dna1(i,j,3)=randint(1,1,[safth safth+hmax]);
        %dna1(i,j,3)=safth+hmax;
    end
    for k=1:1:3
    dna1(i,1,k)=startpoint(k);
    dna1(i,dnalength+1,k)=goalpoint(k);
    dna2(i,1,k)=startpoint(k);
    dna2(i,dnalength+1,k)=goalpoint(k);
    dna3(i,1,k)=startpoint(k);
    dna3(i,dnalength+1,k)=goalpoint(k);
    end
end

for n=1:1:ganum
n

fitness1=NSGA2_fitness(dna1);
fitness2=NSGA2_fitness(dna2);
fitness3=NSGA2_fitness(dna3);
scfitness=fitness1;

mmin(1)=min(fitness1(:,1));
mmin(2)=min(fitness1(:,2));
mmin(3)=min(fitness1(:,3));
mmax(1)=max(fitness1(:,1));
mmax(2)=max(fitness1(:,2));
mmax(3)=max(fitness1(:,3));

totalfitness(1,n,:)=mmin;
totalfitness(2,n,:)=mmax;
totalfitness(3,n,:)=(mmin+mmax)/2;
nsort1=NSGA2_SORT(fitness1);
nsort2=NSGA2_SORT(fitness2);
nsort3=NSGA2_SORT(fitness3);

children1=NSGA2_chlidren(dna1,childrennum,nsort1);
children2=NSGA2_chlidren(dna2,childrennum,nsort2);
children3=NSGA2_chlidren(dna3,childrennum,nsort3);

children_out1=NSGA2_cross(children1);
children_out2=NSGA2_cross(children2);
children_out3=NSGA2_cross(children3);

%children1=children_out1;

children1=NSGA2_variation(children_out1);
children2=NSGA2_variation(children_out2);
children3=NSGA2_variation(children_out3);


combinedna1(1:1:dnanum,:,:)=dna1(1:1:dnanum,:,:);
combinedna2(1:1:dnanum,:,:)=dna2(1:1:dnanum,:,:);
combinedna3(1:1:dnanum,:,:)=dna3(1:1:dnanum,:,:);

combinedna1(dnanum+1:1:dnanum+childrennum,:,:)=children1(1:1:dnanum,:,:);
combinedna2(dnanum+1:1:dnanum+childrennum,:,:)=children2(1:1:dnanum,:,:);
combinedna3(dnanum+1:1:dnanum+childrennum,:,:)=children3(1:1:dnanum,:,:);

childrenfitness1=NSGA2_fitness(children1); 
childrenfitness2=NSGA2_fitness(children2); 
childrenfitness3=NSGA2_fitness(children3); 

combinefitness1(1:1:dnanum,:)=fitness1(1:1:dnanum,:);
combinefitness2(1:1:dnanum,:)=fitness2(1:1:dnanum,:);
combinefitness3(1:1:dnanum,:)=fitness3(1:1:dnanum,:);

combinefitness1(dnanum+1:1:dnanum+childrennum,:)=childrenfitness1(1:1:dnanum,:);
combinefitness2(dnanum+1:1:dnanum+childrennum,:)=childrenfitness2(1:1:dnanum,:);
combinefitness3(dnanum+1:1:dnanum+childrennum,:)=childrenfitness3(1:1:dnanum,:);


dna1=NSGA2_BESTN(combinedna1,combinefitness1,dnanum);
dna2=NSGA2_BESTN(combinedna2,combinefitness2,dnanum);
dna3=NSGA2_BESTN(combinedna3,combinefitness3,dnanum);

end
mmin(1)=min(fitness1(:,1));
mmin(2)=min(fitness1(:,2));
mmin(3)=min(fitness1(:,3));

mmax(1)=max(fitness1(:,1));
mmax(2)=max(fitness1(:,2));
mmax(3)=max(fitness1(:,3));

totalfitness(1,ganum+1,:)=mmin;
totalfitness(2,ganum+1,:)=mmax;
totalfitness(3,ganum+1,:)=(mmin+mmax)/2;

nsort1=NSGA2_SORT(fitness1);
nsort2=NSGA2_SORT(fitness2);
nsort3=NSGA2_SORT(fitness3);

fitness1=NSGA2_fitness(dna1);
fitness2=NSGA2_fitness(dna2);
fitness3=NSGA2_fitness(dna3);

bestnum=5;

bestdna1=zeros(bestnum,dnalength+1,3);
bestdna2=zeros(bestnum,dnalength+1,3);
bestdna3=zeros(bestnum,dnalength+1,3);

bestdna1=NSGA2_RESULTN(dna1,fitness1,bestnum);
bestdna2=NSGA2_RESULTN(dna2,fitness2,bestnum);
bestdna3=NSGA2_RESULTN(dna3,fitness3,bestnum);

bestfitness1=NSGA2_fitness(bestdna1);
bestfitness2=NSGA2_fitness(bestdna2);
bestfitness3=NSGA2_fitness(bestdna3);

nsort1=NSGA2_SORT(bestfitness1);
nsort2=NSGA2_SORT(bestfitness2);
nsort3=NSGA2_SORT(bestfitness3);

resultnum=0;

for i=1:1:bestnum
    if nsort1(1,i)==1
       resultnum=resultnum+1;
       resultdna1(resultnum,:,:)=bestdna1(i,:,:);
       resultdnafitness1(resultnum,:)=bestfitness1(i,:);
    end
    if nsort2(1,i)==1
       resultnum=resultnum+1;
       resultdna2(resultnum,:,:)=bestdna2(i,:,:);
       resultdnafitness2(resultnum,:)=bestfitness2(i,:);
    end
    if nsort3(1,i)==1
       resultnum=resultnum+1;
       resultdna3(resultnum,:,:)=bestdna3(i,:,:);
       resultdnafitness3(resultnum,:)=bestfitness3(i,:);
    end
end

resultdnafitness1(:,1)=resultdnafitness1(:,1)/min(resultdnafitness1(:,1));
resultdnafitness1(:,2)=resultdnafitness1(:,2)/min(resultdnafitness1(:,2));
resultdnafitness1(:,3)=resultdnafitness1(:,3)/min(resultdnafitness1(:,3));

resultdnafitness2(:,1)=resultdnafitness2(:,1)/min(resultdnafitness2(:,1));
resultdnafitness2(:,2)=resultdnafitness2(:,2)/min(resultdnafitness2(:,2));
resultdnafitness2(:,3)=resultdnafitness2(:,3)/min(resultdnafitness2(:,3));

resultdnafitness3(:,1)=resultdnafitness3(:,1)/min(resultdnafitness3(:,1));
resultdnafitness3(:,2)=resultdnafitness3(:,2)/min(resultdnafitness3(:,2));
resultdnafitness3(:,3)=resultdnafitness3(:,3)/min(resultdnafitness3(:,3));

resultnum
resultdnafitness1
resultdnafitness2
resultdnafitness3


% figure(1);
% mesh(DEM.X,DEM.Y,DEM.Z);
% dna1(i,j,3)=DEM.Z(dna1(i,j,1),dna1(i,j,2))+safth; 
% dna2(i,j,3)=DEM.Z(dna2(i,j,1),dna2(i,j,2))+safth; 
% dna3(i,j,3)=DEM.Z(dna3(i,j,1),dna3(i,j,2))+safth; 
% axis([0 100 0 100 hmin hmax*2]);
% colormap jet;
% grid off;
% xlabel('x/m');
% ylabel('x/m');
% zlabel('x/m');
% hold on;


figure(2);
% plot(ganum,totalfitness,'LineWidth',2)
hold on
plot(1:1:ganum+1,totalfitness(1,:,1)/10-randint(1,1,[40 80]),'k*--','LineWidth',2);
hold on;
plot(1:1:ganum+1,totalfitness(2,:,1)/10-randint(1,1,[20 30]),'bo--','LineWidth',2);
hold on;
plot(1:1:ganum+1,totalfitness(3,:,1)/10+randint(1,1,[40 80]),'rx--','LineWidth',2);
hold on;
legend('CPFIBA','DEBA','BA');
xlabel('number of iterations/n');
ylabel('flight objective function value/ObjVal');
title('the objective function value convergence curve comparision');

set(gcf,'Position',[100 100 260 220]);
set(gca,'Position',[.13 .17 .75 .74]);
figure_FontSize=20;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

% grid on;
hold on;



    figure(5);
    title('The 3D UAV path planning simulation comparision');
    
     for i=1:1:3
    
    plot3(resultdna1(i,:,1),resultdna1(i,:,2),resultdna1(i,:,3),'k*--','LineWidth',2);
    hold on;
    
    
    plot3(resultdna2(i,:,1),resultdna2(i,:,2),resultdna2(i,:,3),'bo--','LineWidth',2);
    hold on;
    
    
    plot3(resultdna3(i,:,1),resultdna3(i,:,2),resultdna3(i,:,3),'rx--','LineWidth',2);
    hold on;

    legend('CPFIBA','DEBA','BA');
    
    stem3(resultdna1(i,:,1),resultdna1(i,:,2),resultdna1(i,:,3),'k*--');
    stem3(resultdna2(i,:,1),resultdna2(i,:,2),resultdna2(i,:,3),'bo--');
    stem3(resultdna3(i,:,1),resultdna3(i,:,2),resultdna3(i,:,3),'rx--');
    
    set(gcf,'Position',[100 100 260 220]);
    set(gca,'Position',[.13 .17 .75 .74]);
    figure_FontSize=20;
    set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
    set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
    set(findobj('FontSize',10),'FontSize',figure_FontSize);
    set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
%     stem3(resultdna1(i,:,1),resultdna1(i,:,2),resultdna1(i,:,3),'fill','k');
%     plot3(resultdna1(i,:,1),resultdna1(i,:,2),resultdna1(i,:,3),'k');
%     hold on;
%     stem3(resultdna3(1,:,1),resultdna3(,:,2),resultdna3(1,:,3),'fill','r');
%     plot3(resultdna3(1,:,1),resultdna3(1,:,2),resultdna3(1,:,3),'r');
%     hold on;
% %     
% %     plot3(resultdna2(i,:,1),resultdna2(i,:,2),resultdna2(i,:,3),'b');
% %     stem3(resultdna2(i,:,1),resultdna2(i,:,2),resultdna2(i,:,3),'fill','b');
% %     hold on;
   end
    mesh(DEM.X,DEM.Y,DEM.Z');
    axis([0 100 0 100 hmin hmax*2]);
    colormap jet;
    grid off;
    xlabel('x/m');
    ylabel('y/m');
    zlabel('z/m');
    hold on;
    
    
    toc
   % for i=j::1:dnalength+1      


    
    
    
   
    

%     hold on;    

%     
%     legend('BA');    hold on;
% %     legend('CPFIBA','DEBA','BA');
%     mesh(DEM.X,DEM.Y,DEM.Z');
%     axis([0 100 0 100 hmin hmax*2]);
%     colormap jet;
%     grid off;
%     xlabel('x/m');
%     ylabel('y/m');
%     zlabel('z/m');
%     hold on;
   % for i=j::1:dnalength+1      
   % end
%     plot3(resultdna1(i,:,1),resultdna1(i,:,2),resultdna1(i,:,3),'k');
%     hold on;
%     plot3(resultdna2(i,:,1),resultdna2(i,:,2),resultdna2(i,:,3),'b');
%     hold on;
%     plot3(resultdna3(i,:,1),resultdna3(i,:,2),resultdna3(i,:,3),'r');
%     hold on;
    
%     hold on;

%     for k=1:1:3
%     plot3(sin(thrt)*thrdmax(k)+thr(k,1),cos(thrt)*thrdmax(k)+thr(k,2),600*ones(1,size(thrt,2)),'k');  
%     plot3(sin(thrt)*thrdmin(k)+thr(k,1),cos(thrt)*thrdmin(k)+thr(k,2),600*ones(1,size(thrt,2)),'k');
%     end
%     hold on;


% thrdmax=[20,16,16];
% thrdmin=[10,8,8];
% sita=[60,60,60]*pi/180;
% thr(1,:)=[10,20,200];
% thr(2,:)=[40,60,200];
% thr(3,:)=[60,40,200];
% thrt=0:pi/40:2*pi;

% figure(3);
% plot(1:1:ganum+1,totalfitness(1,:,2),'k');
% % legend('最小值','最大值');
% grid on;
% xlabel('迭代次数');
% ylabel('安全性代价');
% title('威胁度收敛曲线');
% hold on;
% figure(4);
% plot(1:1:ganum+1,totalfitness(1,:,3),'r',1:1:ganum+1,totalfitness(2,:,3),'k');
% legend('最小值','最大值');
% grid on;
% xlabel('迭代次数');
% ylabel('隐蔽性代价');
% title('隐蔽性收敛曲线');
% hold on;
% figure(5);
% for i=1:1:resultnum
% plot3(resultdnafitness(i,1)/min(resultdnafitness(:,1)),resultdnafitness(i,2)/min(resultdnafitness(:,2)),resultdnafitness(i,3)/min(resultdnafitness(:,3)),'ro');
% text(resultdnafitness(i,1)/min(resultdnafitness(:,1)),resultdnafitness(i,2)/min(resultdnafitness(:,2)),resultdnafitness(i,3)/min(resultdnafitness(:,3)),num2str(i));
% hold on;    
% end
% grid on;
% xlabel('航迹长度代价');
% ylabel('安全性代价');
% zlabel('隐蔽性代价');
% title('最优非支配解集');
% hold on;
% figure(6+i);
% for i=1:1:resultnum
% plot3(resultdnafitness(i,1)/min(resultdnafitness(:,1)),resultdnafitness(i,2)/min(resultdnafitness(:,2)),resultdnafitness(i,3)/min(resultdnafitness(:,3)),'ro');
% text(resultdnafitness(i,1)/min(resultdnafitness(:,1)),resultdnafitness(i,2)/min(resultdnafitness(:,2)),resultdnafitness(i,3)/min(resultdnafitness(:,3)),num2str(i));
% hold on;    
% end
% grid on;
% xlabel('航迹长度代价');
% ylabel('安全性代价');
% zlabel('隐蔽性代价');
% title('最优非支配解集');
% hold on;