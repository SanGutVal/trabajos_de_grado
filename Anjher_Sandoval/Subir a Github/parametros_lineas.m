function [orientacion,medianXg,puntoiniY,puntofinY,coefY1,mu1,sttd,mux] = parametros_lineas(pos_s,pixtoeras,s,Mask,titulo,figu)
%Función para calcular los patrones de todas las líneas de contorno
%asociadas a las franjas interferométricas

  

flag=0;
     puntoiniY=[];
    puntofinY=[];
    orientacion=[];
    medianXg=[];
    coefY1=[];
    mu1=[];
    sttd=[];
    kg=[];
    mux=[];
    for k=pixtoeras:length(pos_s)-pixtoeras

%       El n�mero pixtoeras es la cantidad de pixeles que se quitan de los bordes
%--------------------------------------
    %Se busca la coordenada en X para poder limitar la frontera de la
    %M�scara

    medianaX(k-pixtoeras+1)=round(median(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras)));
    %Encontrada la coordena en X, busco el perfil vertical en la m�scara y los
    %puntos que corresponden a la interferencia
    if((isnan(medianaX(k-pixtoeras+1))==0)&& medianaX(k-pixtoeras+1)>0)
        puntosY=find(Mask(:,medianaX(k-pixtoeras+1))>0);
        if (length(puntosY)~=0)
            medianXg=[medianXg median(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras))];
            %Se encuentran los puntos en la frontera de esa l�nea
            puntoiniY=[puntoiniY puntosY(1)];
            puntofinY=[puntofinY puntosY(end)];
            %Se limitan los puntos en Y para poder dibujar la m�scara
            tempCoory=find(s(2,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras)<puntosY(end));
            puntos_dib=find(s(2,pos_s(k)+pixtoeras+tempCoory)>puntosY(1));
%             figure(figu)
%             plot(s(1,pos_s(k)+pixtoeras+puntos_dib),s(2,pos_s(k)+pixtoeras+puntos_dib))

%             if(flag==0)
%                 figure(figu)
%                 grid;hold on;
%                 title(titulo)
%                 flag=1;
%            end                    

%             drawnow;
       %Orientacion Espacial de las franjas periodos

                [coefY,S,mu]=polyfit(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),s(2,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),1); 
                coefY1=[coefY1 coefY(1)];%pendientes de la lineas
                mu1=[mu1 mu(1)];% media del alg de interpolacion de lo datos reales
                sttd=[sttd mu(2)];%desvianci�n estandar de las  
                kg=[kg k];
                %Rectas interpoladas (Y)
                fy=polyval(coefY,s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),S,mu);
                orientacion=[orientacion atan(coefY(1))];
                %figure(4),plot(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),f,'r')
                %hold,
                %plot(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),s(2,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras),'b')
                %legend('Interpolada','Original')
                %Interpolaci�n lineal de puntos en X
                fx=linspace(min(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras)),max(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras)),length(s(1,pos_s(k)+pixtoeras:pos_s(k+1)-pixtoeras)));
                mux=[mux mean(fx)];
            end
        end
    end
end

