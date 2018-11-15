function [ f3 ] = NSGA2_fitness3( dna )
%NSGA2_FITNESS3 Summary of this function goes here
%本函数计算航路隐蔽性
global DEM safth hmax;
dnanum=size(dna,1);
dnalength=size(dna,2)-1;

f3=zeros(dnanum,1);
for i=1:1:dnanum
    f3(i,1)=0;
    for j=1:1:dnalength
        dis(1)=(dna(i,j+1,1)-dna(i,j,1))*10;
        dis(2)=(dna(i,j+1,2)-dna(i,j,2))*10;
        dis=dis.^2;
        d=DEM.Z(dna(i,j,1),dna(i,j,2));
        if dna(i,j,3)-d<safth
            f3(i,1)=f3(i,1)+(2*(hmax+safth))*sqrt(dis(1)+dis(2));
        else
            f3(i,1)=f3(i,1)+(dna(i,j,3)-safth)*sqrt(dis(1)+dis(2));
        end
        %{
        if dna(i,j,3)-d<safth
            f3(i,1)=dnalength*3*hmax;
        elseif f3(i,1)<dnalength*3*hmax
            f3(i,1)=f3(i,1)+dna(i,j,3);
        end
        %}
    end
end
%   Detailed explanation goes here
