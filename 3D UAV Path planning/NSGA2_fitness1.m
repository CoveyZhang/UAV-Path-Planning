function [ f1 ] = NSGA2_fitness1( dna)
%NSGA2_FITNESS1 Summary of this function goes here
%本函数计算航路长度代价
dnanum=size(dna,1);
dnalength=size(dna,2)-1;
f1=zeros(dnanum,1);
for i=1:1:dnanum
    f1(i,1)=0;
    for j=1:1:dnalength
        d(1)=(dna(i,j+1,1)-dna(i,j,1))*10;
        d(2)=(dna(i,j+1,2)-dna(i,j,2))*10;
        d(3)=dna(i,j+1,3)-dna(i,j,3);
        d=d.^2;
        f1(i,1)=f1(i,1)+sqrt(d(1)+d(2)+d(3));
    end
    g=1;
end

%   Detailed explanation goes here
