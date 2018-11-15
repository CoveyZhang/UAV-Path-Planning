function [ msort ] = NSGA2_SORT(fitness_sort_in)
%NSGA2_SORT Summary of this function goes here
%本函数执行dna适应度计算，非支配分层，以后可能层内添加拥挤度排序
%fitness
dnanum_sort=size(fitness_sort_in,1);

np=zeros(2,dnanum_sort);
msort=zeros(1,dnanum_sort)+dnanum_sort+1;
sp=zeros(dnanum_sort,dnanum_sort);
for i=1:1:dnanum_sort
    for j=1:1:dnanum_sort
        if i~=j
            if fitness_sort_in(i,1)>=fitness_sort_in(j,1) & fitness_sort_in(i,2)>=fitness_sort_in(j,2) & fitness_sort_in(i,3)>=fitness_sort_in(j,3)
               if (fitness_sort_in(i,1)-fitness_sort_in(j,1)+fitness_sort_in(i,2)-fitness_sort_in(j,2)+fitness_sort_in(i,3)-fitness_sort_in(j,3))>0    
                   np(1,i)=np(1,i)+1;
                   sp(1,j)=sp(1,j)+1;
                   sp(sp(1,j)+1,j)=i;
               end
               %{
            elseif fitness_sort_in(i,1)>=2*min(fitness_sort_in(:,1)) &  fitness_sort_in(j,1)<2*min(fitness_sort_in(:,1))
                   np(1,i)=np(1,i)+1;
                   sp(1,j)=sp(1,j)+1;
                   sp(sp(1,j)+1,j)=i;
               %}
            end         
        end
    end
end
np(2,:)=np(1,:);
for n=1:1:dnanum_sort
    np(1,:)=np(2,:);
    for i=1:1:dnanum_sort
        if np(1,i)==0
           np(2,i)=dnanum_sort;
           msort(1,i)=n;
           for j=1:1:sp(1,i)
               np(2,sp(j+1,i))=np(2,sp(j+1,i))-1;
           end
        end
    end
end
%   Detailed explanation goes here
