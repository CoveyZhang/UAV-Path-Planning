function [ children_out ] = NSGA2_variation( children )
%NSGA2_VARIATION Summary of this function goes here
%本函数执行变异操作

global DEM safth hmax;
dnanum=size(children,1);
dnalength=size(children,2)-1;
children_out=children;
pvariation=0.5;
n=floor(pvariation*(dnalength-1));
for i=1:1:dnanum
    for m=1:1:n
        j=randint(1,1,[2 dnalength]);
        k(1)=floor((children_out(i,j-1,1)+children_out(i,j+1,1))/2);
        k(2)=floor((children_out(i,j-1,2)+children_out(i,j+1,2))/2);
        k(3)=children_out(i,j,3);
        %k(1)=k(1)+randint(1,1,[0 20])-10;
        %k(2)=k(2)+randint(1,1,[0 20])-10;
       
        if k(1)>101
           k(1)=101;
        elseif k(1)<1
           k(1)=1;
        end
        if k(2)>101
           k(2)=101;
        elseif k(2)<1
           k(2)=1;
        end
        k(3)=DEM.Z(k(1),k(2))+safth;
        %{
        if k(3)>hmax+safth
           k(3)=hmax+safth;
        elseif k(3)<safth
           k(3)=safth;
        end
        %}
        children_out(i,j,:)=k;         
    end
end
%{
for m=1:1:n
    i=randint(1,1,[1 dnanum]);
    j=randint(1,1,[2 dnalength]);
    k(1)=children_out(i,j,1);
    k(2)=children_out(i,j,2);
    k(3)=children_out(i,j,3);
    k(1)=k(1)+randint(1,1,[0 202])-101;
    k(2)=k(2)+randint(1,1,[0 202])-101;
    k(3)=k(3)+randint(1,1,[0 2*hmax])-hmax;
    if k(1)>101
       k(1)=k(1)-101;
    elseif k(1)<1
       k(1)=k(1)+101;
    end
    if k(2)>101
       k(2)=k(2)-101;
    elseif k(2)<1
       k(2)=k(2)+101;
    end
    if k(3)>hmax+safth
       k(3)=k(3)-hmax;
    elseif k(3)<safth
       k(3)=k(3)+hmax;
    end
    children_out(i,j,:)=k;
end
%}

%   Detailed explanation goes here
