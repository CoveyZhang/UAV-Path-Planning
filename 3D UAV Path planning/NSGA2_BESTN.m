function [ dna_out ] = NSGA2_BESTN( dna_in,fitness,N)
%NSGA2_BESTN Summary of this function goes here
global DEM safth hmax;
dna_out=dna_in(1:1:N,:,:);
dnalength=size(dna_in,2)-1;
dnanum=size(fitness,1);
np=zeros(2,dnanum);
msort=zeros(1,dnanum)+dnanum+1;
sp=zeros(dnanum,dnanum);
%s=zeros(dnanum,3)+max(fitness(:,1));
for i=1:1:dnanum
    for j=1:1:dnanum
        if i~=j 
            if fitness(i,1)>=fitness(j,1) & fitness(i,2)>=fitness(j,2) & fitness(i,3)>=fitness(j,3)
               if (fitness(i,1)-fitness(j,1)+fitness(i,2)-fitness(j,2)+fitness(i,3)-fitness(j,3))>0 
                   np(1,i)=np(1,i)+1;
                   sp(1,j)=sp(1,j)+1;
                   sp(sp(1,j)+1,j)=i;
               end 
            %{
            elseif fitness(i,1)>=2*min(fitness(:,1)) & fitness(j,1)<2*min(fitness(:,1))
                   np(1,i)=np(1,i)+1;
                   sp(1,j)=sp(1,j)+1;
                   sp(sp(1,j)+1,j)=i; 
             
            elseif fitness(i,3)>=2*min(fitness(:,3)) & fitness(j,3)<2*min(fitness(:,3))
                   np(1,i)=np(1,i)+1;
                   sp(1,j)=sp(1,j)+1;
                   sp(sp(1,j)+1,j)=i;   
               %}
            end         
        end
        %{
        if i~=j
            if fitness(i,1)>=fitness(j,1)
                if fitness(i,2)>=fitness(j,2)
                    if fitness(i,3)>=fitness(j,3)
                       if (fitness(i,1)-fitness(j,1)+fitness(i,2)-fitness(j,2)+fitness(i,3)-fitness(j,3))>0 
                          np(1,i)=np(1,i)+1;
                          sp(1,j)=sp(1,j)+1;
                          sp(sp(1,j)+1,j)=i;
                       end    
                    end
                end
            end         
        end
        %}
    end
end
x=0;
np(2,:)=np(1,:);
for n=1:1:dnanum
    y=0;
    np(1,:)=np(2,:);
    for i=1:1:dnanum
        if np(1,i)==0
           np(2,i)=dnanum;
           msort(1,i)=n;
           y=y+1;
           for j=1:1:sp(1,i)
               np(2,sp(j+1,i))=np(2,sp(j+1,i))-1;
           end
        end
    end
    if x+y>N
        break;
    else
        x=x+y;
    end
end

k=0;
ld=zeros(2,y);
for m=1:1:dnanum  
    if msort(1,m)==n
    k=k+1;
    s(k,:)=fitness(m,:);
    ld(1,k)=m;
    end
end

s1=s(:,1);
[t,order]=sort(s1);
ld(2,order(1))=ld(2,order(1))+1;
ld(2,order(y))=ld(2,order(y))+1;
if y>2
for m=2:1:y-1  
ld(2,order(m))=ld(2,order(m))+(s1(order(m+1))-s1(order(m-1)))/(max(s1)-min(s1));
end   
end
s2=s(:,2);
[t,order]=sort(s2);
ld(2,order(1))=ld(2,order(1))+1;
ld(2,order(y))=ld(2,order(y))+1;
if y>2
for m=2:1:y-1  
ld(2,order(m))=ld(2,order(m))+(s2(order(m+1))-s2(order(m-1)))/(max(s2)-min(s2));
end   
end
s3=s(:,3);
[t,order]=sort(s3);
ld(2,order(1))=ld(2,order(1))+1;
ld(2,order(y))=ld(2,order(y))+1;
if y>2
for m=2:1:y-1  
ld(2,order(m))=ld(2,order(m))+(s3(order(m+1))-s3(order(m-1)))/(max(s3)-min(s3));
end   
end
[t,order]=sort(ld(2,:));
k=0;

for i=1:1:dnanum
  if msort(1,i)<n
  k=k+1;
  dna_out(k,:,:)=dna_in(i,:,:);
  end
end

if k<N
   for m=x+y-N+1:y
       k=k+1;
       dna_out(k,:,:)=dna_in(ld(1,order(m)),:,:);
   end
end
%{
for k=1:1:N
   for i=2:1:dnalength
      dna_out(k,i,1)=floor((dna_out(k,i+1,1)+dna_out(k,i-1,1))/2);
      dna_out(k,i,2)=floor((dna_out(k,i+1,2)+dna_out(k,i-1,2))/2);
      dna_out(k,i,3)=floor((dna_out(k,i+1,3)+dna_out(k,i-1,3))/2);
      if dna_out(k,i,3)<DEM.Z(dna_out(k,i,1),dna_out(k,i,2))+safth
         dna_out(k,i,3)=DEM.Z(dna_out(k,i,1),dna_out(k,i,2))+safth;
      end
   end
end
%}




%   Detailed explanation goes here
