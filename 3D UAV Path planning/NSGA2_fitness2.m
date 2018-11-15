function [ f2 ] = NSGA2_fitness2( dna )
%NSGA2_FITNESS2 Summary of this function goes here
%本函数计算航路所受到威胁度
global DEM;
dmax=[20,16,16];
dmin=[10,8,8];
sita=[30,30,30]*pi/180;
t(1,:)=[10,20,0];
t(2,:)=[40,60,0];
t(3,:)=[60,40,0];
tnum=size(t,1);
for i=1:1:tnum
    t(i,3)=DEM.Z(t(i,1),t(i,2));
end

dnanum=size(dna,1);
dnalength=size(dna,2)-1;
f2=zeros(dnanum,1);
for i=1:1:dnanum
    f2(i,1)=0;
    for j=2:1:dnalength
        for k=1:1:3
        d=dna(i,j,:)-t(k);
        a=atan2(d(3)/10,sqrt(d(1)*d(1)+d(2)*d(2)));
        d=d.^2;
        dis1(1)=(dna(i,j+1,1)-dna(i,j,1))*10;
        dis1(2)=(dna(i,j+1,2)-dna(i,j,2))*10;
        dis1(3)=dna(i,j+1,3)-dna(i,j,3);
        dis1=dis1.^2;
        dis2(1)=(dna(i,j-1,1)-dna(i,j,1))*10;
        dis2(2)=(dna(i,j-1,2)-dna(i,j,2))*10;
        dis2(3)=dna(i,j-1,3)-dna(i,j,3);
        dis2=dis2.^2;
        dis=sqrt(d(1)+d(2));%+d(3)
        if dis>=dmax(k)
            f2(i,1)=f2(i,1)+0.0001;
        elseif a<=sita(k)
            f2(i,1)=f2(i,1)+0.0001;
        elseif dis>dmin(k)
            f2(i,1)=f2(i,1)+0.0001+((dmin(k)/dis)^2)*(sqrt(dis1(1)+dis1(2)+dis1(3))+sqrt(dis2(1)+dis2(2)+dis2(3)))/2;
        else
            f2(i,1)=f2(i,1)+0.0001+1*(sqrt(dis1(1)+dis1(2)+dis1(3))+sqrt(dis2(1)+dis2(2)+dis2(3)))/2;
        end 
        end
    end
end
%   Detailed explanation goes here
