distance=2.2;
separation=0.07;
lambda=632e-9;
l=1;
%This algorithm tries to find global phase from Takeda algorithm. So, it
%takes the images, executes FFT2, find fundamental frequencies and separate
%them of spectrum, after it applies hamming windowing and return to direct
%space where is calculated phase.
M=1600; N=1200;
close all
% path='./Videos/amb';
path='./Videos/01_04_HUM';
spi=1;
spf=2031;
vidObj_ph = VideoWriter([path '/fulc_phase.avi'],'Motion JPEG AVI');
vidObj_ph.FrameRate=5;
open(vidObj_ph);
% V=VideoReader([path '/31-03_amb.avi']);
V=VideoReader([path '/hmd.avi']);

nFrames=V.NumberOfFrames;
az = 90;
el = 0;

tic;
for i=spi:nFrames
    I=read(V,i);
    D1=double(I(:,:,1));
    D=(hanning(N)*hanning(M)').*D1;
    I1=fftshift(fft2(double(D)));
    if (i==spi)
        [a1,ind]=max(abs(I1(600,805:end)));   %Se toma en la mitad del espectro (600)
        %para verificar d�nde aperece la fundamental, tambi�n se toman 4 pixeles
        %despu�s de la DC para encontrar el centro de masa.
        ind1=804+ind;
        %Aqu� se verifica la posici�n del centro de masa del lobulo.
        [a2,ind2]=max(abs(I1(:,ind1)));
        %A partir del centro de masa del l�bulo se sacan las posiciones
        %alternas y as� recortar el l�bulo
        posicini=find(abs(I1(:,ind1))>a2/100);
        %Aqu� se toma el primer valor separado de fo/3 para poder encontrar el
        %�rea a filtrar
        [a3,ind3]=max(abs(I1(601,:)));
        
        deltau=round((ind+4)/3)+1;
    end
    disp(i)
    %Buscamos extraer la imagen sin el fondo continuo
    %Se busca una copia para llevar s�lo los cambios de la fundamental al
    %espacio directo
    I1takeda=I1;
    I1takeda(1:ind2-deltau,:)=0;
    I1takeda(ind2+deltau:end,:)=0;
    I1takeda(ind2-deltau:ind2+deltau,1:ind1-deltau)=0;
    I1takeda(ind2-deltau:ind2+deltau,ind1+deltau:end)=0;
%     i1takeda=ifft2(I1takeda);     No se utiliza porque se observa la
%     modulación de la onda cuadrada
    %Aqu� se pasa un filtro hanning de suavizado
    I1takeda(ind2-deltau:ind2+deltau-1,ind1-deltau:ind1+deltau-1)=...
    I1takeda(ind2-deltau:ind2+deltau-1,ind1-deltau:ind1+deltau-1).*(hanning(deltau*2)*hanning(deltau*2)');
    i1tak=ifft2(fftshift(I1takeda));
    %Esta es la fase para ambas im�genes
%     phasei1takeda=angle(i1takeda);
    phasei1tak=angle(i1tak);
    %Aqu� est� la amplitud de esa fase
%     magni1takeda=abs(i1takeda);
    magni1tak=abs(i1tak);       %Aqu� es donde se ve el filtro hanning.
    %Esto simplemente es para pasar el filtro. En reconstrucci�n 3D lo
    %utilizan para s�lo buscar el objeto.
         
    Mask=((magni1tak-min(magni1tak(:)))/(max(magni1tak(:))-min(magni1tak(:))))>=1/exp(2);
    Mask2=((magni1tak-min(magni1tak(:)))/(max(magni1tak(:))-min(magni1tak(:))))>0;
  
    se=strel('disk',30);
    b11=imerode(Mask,se);
    Mask1=imdilate(b11,se);
    
    sizeIm=size(D);
    phasei1tak_u=phasei1tak;
    [PhaseC,MaskF,tiempo]=unwrap2DClasico([ind1 ind2],phasei1tak_u,Mask2);

     %Trazar l�neas de fase a partir del contorno con alta resoluci�n
     
%      tt=find(Mask1(600,1:end)>0); %(Dentro de la máscara)
   
     s1=contourc(PhaseC,round((PhaseC(round(sizeIm(1)/2),1):PhaseC(round(sizeIm(1)/2),end))*2*pi));
     
%%Determinaci�n de Plano de Fase interpolado y comparacion de las fluctuaciones de Fase
% 
% [fitplx, Sx, mdx]=polyfit(xx, PhaseC,1);
% [fitply, Sy mdy]=polyfit(yy, PhaseC,1);
% Xeval=polyval(fitplx,xx, Sx, mdx);
% Yeval=polyval(fitply,yy,Sy,mdy);
% pl=Xeval+Yeval-fitply(2); % Plano interpolado linealmente
% pendx(i-spi+1)=fitplx(1);
% pendy(i-spi+1)=fitply(1);

pl1=PolyVal2D_1(Polyfit2D_1(PhaseC,Mask1,1:M,1:N,1),1:M,1:N,1);
pl=pl1.*Mask1;
figure(1),mesh((pl-PhaseC).*Mask1), colorbar% Las fluctuaciones de fase NO son superiores a pi/4, es decir,
title(['Fluctuaciones de la fase. Muestra ' num2str(i) '/' num2str(nFrames)])
plane_timex(i,:)=(pl(600,:)-PhaseC(600,:)).*Mask1(600,:);
plane_timey(i,:)=(pl(:,800)-PhaseC(:,800)).*Mask1(:,800);
if (i==spi)
%     axis tight manual
    set(gca,'nextplot','replacechildren');

end
view(az, el);
currFrame_ph=getframe(gcf);
writeVideo(vidObj_ph,currFrame_ph);
%no son mas de un per�odo de las franjas, por tanto, se pueden considerar paralelas las franjas.

%Líneas de contorno para plano inclinado interpolado
s2=contourc(pl1,round((pl1(round(sizeIm(1)/2),1):pl1(round(sizeIm(1)/2),end))*2*pi));

%  figure(7); imagesc(D); colormap('gray')
%     title(['Interpolaci�n de las franjas. Imagen ' num2str(i) ' de ' num2str(nFrames)])
%     hold on, plot(s2(1,:),s2(2,:),'g.')
%     hold off

    pos_s1=[];
    
   %Se busca separar cada l�nea (el valor de 10 es por el periodo de las franjas)
    pos_s1=find(diff(s1(2,:))>10);
    pos_s2=find(diff(s2(2,:))>10);
    %Con este algoritmo se dibujan las l�neas
    pix_eras=10;

%     figure(3);
%     figure(4);
   [orientacion,medianXg,puntoiniY,puntofinY,coefY1,mu1,sttd,mux] =parametros_lineas(pos_s1,pix_eras,s1,Mask1,'Líneas curvas encontradas',3);
   [orientacion_i,medianXg_i,puntoiniY_i,puntofinY_i,coefY1_i,mu1_i,sttd_i,mux_i] =parametros_lineas(pos_s2,pix_eras,s2,Mask1,'Líneas del plano interpolado',4);

   [period_i errorx_i errory_i]=periodos(pos_s2,pix_eras,s2,-1/median(coefY1_i),M);
    
    
    %Colecciones para estimaci�n en el tiempo Curvas
    orient(i-spi+1,1:length(orientacion))=orientacion;
    medX(i-spi+1,1:length(medianXg))=medianXg;
    piniY(i-spi+1,1:length(puntoiniY))=puntoiniY;%Ver puntos iniciales de las lineas en la mascara
    pfinY(i-spi+1,1:length(puntofinY))=puntofinY;%Ver puntos finales de las lineas en la mascara
    pendient(i-spi+1,1:length(coefY1))=coefY1;
    meaOrigX(i-spi+1,1:length(mu1))=mu1;
    stdOrigX(i-spi+1,1:length(sttd))=sttd;
    meaIntX(i-spi+1,1:length(mux))=mux;
    periodos_int(i,1:length(period_i))=period_i;
    errorx_int(i,1:length(errorx_i))=errorx_i;
    errory_int(i,1:length(errory_i))=errory_i;
%     figure(5), plot(orientacion(5:end)-mean(orientacion(5:end)))
%     title('Fluctuaciones de la Orientacion Espacial');
%     grid
%Determinacion del Periodo Espacial

   pf1(i-spi+1,1:length(medianXg)-1)=diff(medianXg); %Original
   pr1(i-spi+1,1:length(medianXg)-1)=pf1(i-spi+1,1:length(medianXg)-1).*sin(mean(orientacion));
   pf2(i-spi+1,1:length(mux)-1)=diff(mux); %Interpolado
   pr2(i-spi+1,1:length(mux)-1)=pf2(i-spi+1,1:length(mux)-1).*sin(mean(orientacion));
   
    %Determinaci�n de Contraste filtrando la f. DC
%     P1=I1;
%     P1(1:ind2-deltau,:)=0;
%     P1(ind2+deltau:end,:)=0;
%     P1(ind2-deltau:ind2+deltau,1:ind3-deltau)=0;
%     P1(ind2-deltau:ind2+deltau,ind3+deltau:end)=0;
%     A=ifft2(fftshift(P1)); %Fondo cont�nuo
%     V= 2*magni1tak./abs(A); %Visibilidad 
    A=D-2*magni1tak.*cos(PhaseC);
    Vis= 2*magni1tak./A; %Visibilidad
    Vis_g(i,:)=Vis(600,:).*Mask1(600,:);
%     figure(6),imagesc(Vis.*Mask1);
%     title('Visibilidad')
%     drawnow
    r1(i)=toc;
end
vidObj_ph.close();
close all
save([path '/parameters.mat'])