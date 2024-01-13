function [period errorx errory] = periodos(pos_s,pixtoeras,s,pendiente,M)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

xg=[];
yg=[];
errorx=[];
errory=[];
% figure(5);
for k=pixtoeras:length(pos_s)-pixtoeras
    %Se busca el primer valor para arrancar
%     plot(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),s(2,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),'b')
    if (k==pixtoeras)
        medianaX(k-pixtoeras+1)=median(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras));
        pos_xo_yo=find(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras)>medianaX(k-pixtoeras+1));
        xo=s(1,pos_xo_yo(1)+pos_s(k)+pixtoeras);yo=s(2,pos_xo_yo(1)+pos_s(k)+pixtoeras);
        polinomio=[pendiente yo-pendiente*xo];
        xg=[xg xo]; yg=[yg yo];
%         hold
    else
        pol_evaly=polyval(polinomio,s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),1);
        pos_x1_y1=find(s(2,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras)>mean(pol_evaly));
        pos_x2_y2=find(s(2,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras)<mean(pol_evaly));
%         plot(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),pol_evaly,'r');
%         plot(s(1,pos_s(k)+pixtoeras+pos_x1_y1(1)-1),s(2,pos_s(k)+pixtoeras+pos_x1_y1(1)-1),'g*')
%         plot(s(1,pos_s(k)+pixtoeras+pos_x2_y2(end)-1),s(2,pos_s(k)+pixtoeras+pos_x2_y2(end)-1),'k*')
        xg=[xg (s(1,pos_s(k)+pixtoeras+pos_x1_y1(1)-1)+s(1,pos_s(k)+pixtoeras+pos_x2_y2(end)-1))/2];
        errorx=[errorx abs(s(1,pos_s(k)+pixtoeras+pos_x1_y1(1)-1)-s(1,pos_s(k)+pixtoeras+pos_x2_y2(end)-1))];
        yg=[yg (s(2,pos_s(k)+pixtoeras+pos_x1_y1(1)-1)+s(2,pos_s(k)+pixtoeras+pos_x2_y2(end)-1))/2];
        errory=[errory abs(s(1,pos_s(k)+pixtoeras+pos_x1_y1(1)-1)-s(1,pos_s(k)+pixtoeras+pos_x2_y2(end)-1))];
%         plot(xg(end),yg(end),'y*')
    end
    
end
period=sqrt(diff(xg).^2+diff(yg).^2);
end

