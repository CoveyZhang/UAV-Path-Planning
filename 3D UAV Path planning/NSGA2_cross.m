function [ children_out ] = NSGA2_cross( children )
%NSGA2_CROSS Summary of this function goes here
%本函数执行交叉操作
global DEM safth;
dnanum=size(children,1);
dnalength=size(children,2)-1;
children_out=children;
pcross=1;
t=zeros(dnanum,dnalength+1,3);
for m=1:1:floor(dnanum/2)
    if rand(1,1)<=pcross
       i=randint(1,1,[1 dnanum]);
       j=randint(1,1,[1 dnanum]);
       
       k=randint(1,1,[2 dnalength]);
       t(i,:,:)=children(i,:,:);
       children_out(i,k:1:dnalength,:)=children_out(j,k:1:dnalength,:);
       children_out(j,k:1:dnalength,:)=t(i,k:1:dnalength,:);
       %{
       k=rand(1,1);
       children_out(i,:,:)=floor(k*children(i,:,:)+(1-k)*children(j,:,:));
       children_out(j,:,:)=floor(k*children(j,:,:)+(1-k)*children(i,:,:));
       %}
    end
end


%   Detailed explanation goes here
