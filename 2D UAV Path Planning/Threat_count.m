function Threat= Threat_count(aim,path,Threat_center)

%   计算威胁对无人机航路产生的威胁代价

global Threat_radius;
global Threat_kind;
global d;
global k;
n=length(path);
distance(1)=sqrt(d^2+path(1)^2);
threat=zeros(n+1,1);
for i=1:length(Threat_radius)
    if(((path(1)-Threat_center(2,i))^2+(d-Threat_center(1,i))^2)<Threat_radius(i)^2)
        threat(1)=10+threat(1)+distance(1)^5*Threat_kind(i)/5*...
            (1/dis(Threat_center(1,i),Threat_center(2,i),0.1*d,0.1*path(1))^4+1/dis(Threat_center(1,i),Threat_center(2,i),0.3*d,0.3*path(1))^4....
        +1/dis(Threat_center(1,i),Threat_center(2,i),0.5*d,0.5*path(1))^4+...
        1/dis(Threat_center(1,i),Threat_center(2,i),0.7*d,0.7*path(1))^4+1/dis(Threat_center(1,i),Threat_center(2,i),0.9*d,0.9*path(1))^4);
    end
end
for i=2:n
    distance(i)=sqrt(d^2+(path(i)-path(i-1))^2);
    for j=1:length(Threat_radius)
        if(((path(i)-Threat_center(2,j))^2+(i*d-Threat_center(1,j))^2)<Threat_radius(j)^2)
            threat(i)=10+threat(i)+distance(i)^5*Threat_kind(j)/5*(1/dis(Threat_center(1,j),Threat_center(2,j),((i-0.9)*d),(0.9*path(i-1)+0.1*path(i)))^4+...
              1/dis(Threat_center(1,j),Threat_center(2,j),((i-0.7)*d),(0.7*path(i-1)+0.3*path(i)))^4+...
              1/dis(Threat_center(1,j),Threat_center(2,j),((i-0.5)*d),(0.5*path(i-1)+0.5*path(i)))^4+...
              1/dis(Threat_center(1,j),Threat_center(2,j),((i-0.3)*d),(0.3*path(i-1)+0.7*path(i)))^4+...
              1/dis(Threat_center(1,j),Threat_center(2,j),((i-0.1)*d),(0.1*path(i-1)+0.9*path(i)))^4);
        end
    end
end
distance(n+1)=dis(aim(1),aim(2),n*d,path(n));
for j=1:length(Threat_radius)
    if(abs((aim(2)-path(n))*Threat_center(1,j)-d*Threat_center(2,j)+path(n)*(n+1)*d-n*d*aim(2))<Threat_radius(j)*distance(n+1))
        threat(n+1)=threat(n+1)+distance(n+1)^5*Threat_kind(j)/5*(1/dis(Threat_center(1,j),Threat_center(2,j),((n+0.1)*d),(0.9*path(n)+0.1*aim(1)))^4+...
            1/dis(Threat_center(1,j),Threat_center(2,j),((n+0.3)*d),(0.7*path(n)+0.3*aim(1)))^4+...
            1/dis(Threat_center(1,j),Threat_center(2,j),((n+0.5)*d),(0.5*path(n)+0.5*aim(1)))^4+...
            1/dis(Threat_center(1,j),Threat_center(2,j),((n+0.7)*d),(0.3*path(n)+0.7*aim(1)))^4+...
            1/dis(Threat_center(1,j),Threat_center(2,j),((n+0.9)*d),(0.1*path(n)+0.9*aim(1)))^4);
    end
end
T=sum(threat);
D=sum(distance);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
Threat=k*T+(1-k)*D;

end


