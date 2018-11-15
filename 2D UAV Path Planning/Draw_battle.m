function Draw= Draw_battle( Threat_center,Threat_radius)
for i=1:size(Threat_center,2)
    scatter(Threat_center(1,i),Threat_center(2,i),'ok');
    hold on
end
for i=1:size(Threat_center,2)
    cir_x=[ ];
    cir_y=[ ];
    for j=0:pi/500:2*pi
        aa_x=Threat_radius(i)*cos(j)+Threat_center(1,i);
        aa_y=Threat_radius(i)*sin(j)+Threat_center(2,i);
        cir_x=[cir_x,aa_x];
        cir_y=[cir_y,aa_y];
    end
    plot(cir_x,cir_y,'k--','LineWidth',2);
    hold on
end


end

