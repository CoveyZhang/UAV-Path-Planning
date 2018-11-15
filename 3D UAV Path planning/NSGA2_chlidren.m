function [ children ] = NSGA2_chlidren( dna,childrennum,msort )
%NSGA2_CHLIDREN Summary of this function goes here
%本函数执行选择操作，产生子代种群
global scfitness; 
dnalength=size(dna,2)-1;
dnanum=size(dna,1);
children=zeros(childrennum,dnalength+1,3);

for n=1:1:childrennum
    i=randint(1,1,[1 dnanum]);
    j=randint(1,1,[1 dnanum]);
    children(n,:,:)=dna(i,:,:);
    if msort(1,i)>msort(1,j)
        children(n,:,:)=dna(j,:,:);
    end
    if msort(1,i)==msort(1,j)
    if scfitness(i,1)>=2*min(scfitness(:,1)) & scfitness(j,1)<2*min(scfitness(:,1))
        children(n,:,:)=dna(j,:,:);
    end
    end
end
%{
for i=1:1:childrennum
    for j=1:1:dnalength+1
        for k=1:1:3
        children(i,j,k)=randint(1,1,[0 100]);   
        end     
    end
end
%}