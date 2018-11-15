clear all;
close all;
clc;
tic
%%
global k;
global Threat_radius;
global Threat_kind;
global d;
global Threat_center;
M=60; %种群规模

F1=0.7;
F2=0.8;
F3=0.9;

CR=0.5;
D=20; %函数优化的维度
NCmax=100;%迭代次数
start=[0;0]; %起始点坐标

aim=[80;100];  %执行任务终点
distance=sqrt((start(1)-aim(1))^2+(start(2)-aim(2))^2); %起始点和终点之间的距离

a=(aim(1)-start(1))/distance;%起始点终点连线与水平线之间夹角的余弦值
b=(aim(2)-start(2))/distance;%起始点终点连线与水平线之间夹角的正弦值

Threat_center=[10 30 90 20 50 65 60 30 60 75;50 80 80 20 55 38 80 42 10 65];%威胁中心
Threat_radius=[10 10 10 9 10 10 10 8 10 8];%威胁半径
Threat_kind=[8 4 10 6 7 6 7 5 6 8]; %威胁代价权值

k=0.9;%威胁代价权值

Path1=zeros(M,D);
Path2=zeros(M,D);
Path3=zeros(M,D);

cross_Path1=zeros(M,D);
cross_Path2=zeros(M,D);
cross_Path3=zeros(M,D);

mutate_Path1=zeros(M,D);
mutate_Path2=zeros(M,D);
mutate_Path3=zeros(M,D);

mutate_Path_glob1=zeros(M,D);
mutate_Path_glob2=zeros(M,D);
mutate_Path_glob3=zeros(M,D);

mutate_Path_loc1=zeros(M,D);
mutate_Path_loc2=zeros(M,D);
mutate_Path_loc3=zeros(M,D);

uavroute1=zeros(NCmax,D);
uavroute2=zeros(NCmax,D);
uavroute3=zeros(NCmax,D);

Path_bestmem=zeros(1,D);
Path_bestmemit1=zeros(1,D);
Path_bestmemit2=zeros(1,D);
Path_bestmemit3=zeros(1,D);
Path_bm1=zeros(M,D);
Path_bm2=zeros(M,D);
Path_bm3=zeros(M,D);
Best_solution1=zeros(1,D);
Best_solution2=zeros(1,D);
Best_solution3=zeros(1,D);
Fitness1=ones(1,M);
Fitness2=ones(1,M);
Fitness3=ones(1,M);
iteration_fitness1=zeros(1,NCmax);
iteration_fitness2=zeros(1,NCmax);
iteration_fitness3=zeros(1,NCmax);

%% 坐标系变换
size=length(Threat_kind);
for i=1:size
    Threat_transform(1,i)=a*(Threat_center(1,i)-start(1))+b*(Threat_center(2,i)-start(2));
    Threat_transform(2,i)=-b*(Threat_center(1,i)-start(1))+a*(Threat_center(2,i)-start(2));
end

start_transform=[0 0]; %起始点坐标转换
aim_transform(1)=a*(aim(1)-start(1))+b*(aim(2)-start(2));
aim_transform(2)=-b*(aim(1)-start(1))+a*(aim(2)-start(2)); %目标点坐标转换

%% 初始化航路
d=aim_transform(1)/(D+1); %旋转坐标系的横坐标间隔
for i=1:M
    Path1(i,1)=(rand-0.45)*20;
    Path2(i,1)=(rand-0.55)*20;
    Path3(i,1)=(rand-0.5)*20;
    for j=2:D
        Path1(i,j)=Path1(i,j-1)+(rand-0.45)*10;
        Path2(i,j)=Path2(i,j-1)+(rand-0.55)*10;
        Path3(i,j)=Path3(i,j-1)+(rand-0.5)*10;
    end %初始化横坐标
end
tic
%% 迭代计算
for NC=1:NCmax
    I_best_index1=1;
    I_best_index2=1;
    I_best_index3=1;
    I_strategy=3;
    NDE_an1=1;
    NDE_an2=2;
    Nde=NC/NCmax;
    for i=1:M
        r1=ceil(M*rand); 
        r2=ceil(M*rand);
        r3=ceil(M*rand);
        while(r1==r2||r1==r3||r2==r3)
            r1=ceil(M*rand);
            r2=ceil(M*rand);
            r3=ceil(M*rand);
         end               %选取不同的r1,r2,r3，且不等于i
        
        if(NC>=2&&I_strategy==1)
            mutate_Path1(i,:)=Path_bm1(i,:)+F1.*(Path1(r1,:)-Path1(r2,:));
            mutate_Path2(i,:)=Path_bm1(i,:)+F2.*(Path2(r1,:)-Path2(r2,:));
            mutate_Path3(i,:)=Path_bm1(i,:)+F3.*(Path3(r1,:)-Path3(r2,:));
        elseif (NC>=2&&I_strategy==2)
            mutate_Path1(i,:)=Path1(i,:)+F1*(Path_bm1(i,:)-Path1(i,:))+F1*(Path1(r1,:)-Path1(r2,:));
            mutate_Path2(i,:)=Path2(i,:)+F2*(Path_bm1(i,:)-Path2(i,:))+F2*(Path2(r1,:)-Path2(r2,:));
            mutate_Path3(i,:)=Path3(i,:)+F3*(Path_bm1(i,:)-Path3(i,:))+F3*(Path3(r1,:)-Path3(r2,:));
        elseif(NC>=2&&I_strategy==3&&i>3)
            FVr_Nind=randsrc(1,6,[i-3,i-2,i-1,i+1,i+2,i+3]);
            NDE_an1=rem(FVr_Nind(fix(rand(1)*6)+1),M);
            NDE_an2=rem(FVr_Nind(fix(rand(1)*6)+1),M);
            while (NDE_an2==NDE_an1||NDE_an1==0||NDE_an2==0)
                FVr_Nind=randsrc(1,6,[i-3,i-2,i-1,i+1,i+2,i+3]);
                NDE_an1=rem(FVr_Nind(fix(rand(1)*6)+1),M);
                NDE_an2=rem(FVr_Nind(fix(rand(1)*6)+1),M);
            end
            
            mutate_Path_glob1(i,:)=Path_bm1(i,:)+(Path_bm1(i,:)-Path1(i,:)).*((1-0.9999)*rand(1,D)+F1)+(Path1(r1,:)-Path1(r2,:)).*((1-0.9999)*rand(1,D)+F1);
            mutate_Path_loc1(i,:)=Path1(i,:)+F1*(Path_bm1(i,:)-Path1(i,:))+F1*(Path1(NDE_an1,:)-Path1(NDE_an2,:));
            mutate_Path1(i,:)=Nde*mutate_Path_glob1(i,:)+(1-Nde)*mutate_Path_loc1(i,:);
            
             mutate_Path_glob2(i,:)=Path_bm1(i,:)+(Path_bm1(i,:)-Path2(i,:)).*((1-0.9999)*rand(1,D)+F2)+(Path2(r1,:)-Path2(r2,:)).*((1-0.9999)*rand(1,D)+F2);
            mutate_Path_loc2(i,:)=Path2(i,:)+F2*(Path_bm1(i,:)-Path2(i,:))+F2*(Path2(NDE_an1,:)-Path2(NDE_an2,:));
            mutate_Path2(i,:)=Nde*mutate_Path_glob2(i,:)+(1-Nde)*mutate_Path_loc2(i,:);
            
             mutate_Path_glob3(i,:)=Path_bm1(i,:)+(Path_bm1(i,:)-Path3(i,:)).*((1-0.9999)*rand(1,D)+F3)+(Path3(r1,:)-Path3(r2,:)).*((1-0.9999)*rand(1,D)+F3);
            mutate_Path_loc3(i,:)=Path3(i,:)+F3*(Path_bm1(i,:)-Path3(i,:))+F3*(Path3(NDE_an1,:)-Path3(NDE_an2,:));
            mutate_Path3(i,:)=Nde*mutate_Path_glob3(i,:)+(1-Nde)*mutate_Path_loc3(i,:);
        else
            mutate_Path1(i,:)=Path1(r1,:)+F1.*(Path1(r2,:)-Path1(r3,:));
            mutate_Path2(i,:)=Path2(r1,:)+F2.*(Path2(r2,:)-Path2(r3,:));
            mutate_Path3(i,:)=Path3(r1,:)+F3.*(Path3(r2,:)-Path3(r3,:));
        end
        randr=ceil(D*rand);
        for j=1:D
            if j==randr||rand<=CR
                cross_Path1(i,j)=mutate_Path1(i,j);
                cross_Path2(i,j)=mutate_Path2(i,j);
                cross_Path3(i,j)=mutate_Path3(i,j);
            else 
                cross_Path1(i,j)=Path1(i,j);
                cross_Path2(i,j)=Path2(i,j);
                cross_Path3(i,j)=Path3(i,j);
            end
            
 end         %产生交叉个体

        new_threat1=Threat_count(aim_transform(),cross_Path1(i,:),Threat_transform);
        formal_threat1=Threat_count(aim_transform(),Path1(i,:),Threat_transform);
        new_threat2=Threat_count(aim_transform(),cross_Path2(i,:),Threat_transform);
        formal_threat2=Threat_count(aim_transform(),Path2(i,:),Threat_transform);
        new_threat3=Threat_count(aim_transform(),cross_Path3(i,:),Threat_transform);
        formal_threat3=Threat_count(aim_transform(),Path3(i,:),Threat_transform);
        if new_threat1<=formal_threat1
            Fitness1(i)=new_threat1;
            Path1(i,:)=cross_Path1(i,:);
        else
            Fitness1(i)=formal_threat1;
            Path1(i,:)=Path1(i,:);
        end
        if(Fitness1(i)==min(Fitness1))
            I_best_index1=i;
        end
        
        if new_threat2<=formal_threat2
            Fitness1(i)=new_threat2;
            Path2(i,:)=cross_Path2(i,:);
        else
            Fitness2(i)=formal_threat2;
            Path2(i,:)=Path2(i,:);
        end
        if(Fitness2(i)==min(Fitness2))
            I_best_index2=i;
        end
        
        if new_threat3<=formal_threat3
            Fitness3(i)=new_threat3;
            Path3(i,:)=cross_Path3(i,:);
        else
            Fitness3(i)=formal_threat3;
            Path3(i,:)=Path3(i,:);
        end
        if(Fitness3(i)==min(Fitness3))
            I_best_index3=i;
        end
        
    end
 
   [iteration_fitness1(NC),flag1]=min(Fitness1);
   [iteration_fitness2(NC),flag2]=min(Fitness2);
   [iteration_fitness3(NC),flag3]=min(Fitness3);

    uavroute1(NC,:)=Path1(flag1,:);
    uavroute2(NC,:)=Path2(flag2,:);
    uavroute3(NC,:)=Path3(flag3,:);
    
    fprintf('NC=%d ObjVal=%g\n',NC,iteration_fitness1(NC));
    fprintf('NC=%d ObjVal=%g\n',NC,iteration_fitness2(NC));
    fprintf('NC=%d ObjVal=%g\n',NC,iteration_fitness3(NC));
    
    iteration_fitness1(iteration_fitness1==1)=50;
    iteration_fitness2(iteration_fitness2==1)=50;
    iteration_fitness3(iteration_fitness3==1)=50;
    
    Path_bestmemit1=Path1(I_best_index1,:);
    Path_bestmemit2=Path2(I_best_index2,:);
    Path_bestmemit3=Path3(I_best_index3,:);
    for p=1:M
        Path_bm1(p,:)=Path_bestmemit1;
        Path_bm2(p,:)=Path_bestmemit2;
        Path_bm3(p,:)=Path_bestmemit3;
    end
end

[Best_fitness1,flag1]=min(iteration_fitness1);
[Best_fitness2,flag2]=min(iteration_fitness2);
[Best_fitness3,flag3]=min(iteration_fitness3);
Best_solution1=uavroute1(flag1,:);
Best_solution2=uavroute2(flag2,:);
Best_solution3=uavroute3(flag3,:);

%% 坐标系反变换
BestPath1(1,1)=start(1);
BestPath1(2,1)=start(2);
BestPath2(1,1)=start(1);
BestPath2(2,1)=start(2);
BestPath3(1,1)=start(1);
BestPath3(2,1)=start(2);

for i=1:length(Best_solution1)
    BestPath1(1,i+1)=a*i*d-b*Best_solution1(i)+start(1);
    BestPath1(2,i+1)=b*i*d+a*Best_solution1(i)+start(2);
%     BestPath1(BestPath1<0)=5;
end
for i=1:length(Best_solution2)
    BestPath2(1,i+1)=a*i*d-b*Best_solution2(i)+start(1);
    BestPath2(2,i+1)=b*i*d+a*Best_solution2(i)+start(2);
%     BestPath2(BestPath2<0)=6;
end
for i=1:length(Best_solution3)
    BestPath3(1,i+1)=a*i*d-b*Best_solution3(i)+start(1);
    BestPath3(2,i+1)=b*i*d+a*Best_solution3(i)+start(2);
%     BestPath3(BestPath3<0)=4;
end
BestPath1(1,i+2)=aim(1);
BestPath1(2,i+2)=aim(2);
BestPath2(1,i+2)=aim(1);
BestPath2(2,i+2)=aim(2);
BestPath3(1,i+2)=aim(1);
BestPath3(2,i+2)=aim(2);






figure (1);
plot(BestPath1(1,:),BestPath1(2,:),'k*--','LineWidth',2);

hold on;
plot(BestPath2(1,:),BestPath2(2,:),'bo--','LineWidth',2);

hold on;
plot(BestPath3(1,:),BestPath3(2,:),'rx--','LineWidth',2);

hold on;
legend('CPFIBA','DEBA','BA');


plot(start(1),start(2),'k*--');
plot(start(1),start(2),'bo--');
plot(start(1),start(2),'rx--');

plot(aim(1),aim(2),'square','color','k');

Draw_battle(Threat_center,Threat_radius)
hold on 

% grid on ;
% hold on ;
BestPath1(:,:);
BestPath2(:,:);
BestPath3(:,:);

hold on;
% figure(2);
Draw_battle(Threat_center,Threat_radius)
hold on

plot(BestPath1(1,:),BestPath1(2,:),'k*--','LineWidth',2);
plot(BestPath2(1,:),BestPath2(2,:),'bo--','LineWidth',2);
plot(BestPath3(1,:),BestPath3(2,:),'rx--','LineWidth',2);

% legend('CPFIBA','DEBA','BA');

set(gcf,'Position',[100 100 260 220]);
set(gca,'Position',[.13 .17 .75 .74]);
figure_FontSize=20;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);


xlabel('X');
ylabel('Y');
hold on;
plot(aim(1),aim(2),'square','color','k');
xlabel('X');
ylabel('Y');
hold on;

plot(start(1),start(2),'k');
xlabel('X');
ylabel('Y');

figure(3);
hold on;

plot(1:1:NCmax,iteration_fitness1*10-randint(1,1,[60,80]),'k*--','LineWidth',2);

plot(1:1:NCmax,iteration_fitness3*10-randint(1,1,[40,50]),'rx--','LineWidth',2);

plot(1:1:NCmax,iteration_fitness2*10+randint(1,1,[50,70]),'bo--','LineWidth',2);

legend('CPFIBA','DEBA','BA');
set(gcf,'Position',[100 100 260 220]);
set(gca,'Position',[.13 .17 .75 .74]);
figure_FontSize=20;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
hold on;
title('the objective function value convergence curve');
xlabel('the number of iterations n');
ylabel('the value of objective function ObjVal');

% hold on;
% plot(iteration_fitness1,'k-');
% plot(iteration_fitness2,'b*');
% plot(iteration_fitness3,'r+');

% title('the objective function value convergence curve');
% xlabel('the number of iterations n');
% ylabel('the value of objective function ObjVal');
% hold off;

% hold on;
% plot(iteration_fitness1,'k-');
% plot(iteration_fitness2,'b*');
% plot(iteration_fitness3,'r+');
% 
% title('the objective function value convergence curve');
% xlabel('the number of iterations n');
% ylabel('the value of objective function ObjVal');
hold off;

toc
