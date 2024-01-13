clc
clear all
close all

p=1;
k=64;
l=1;
%path='.\frame_';
path='C:\Users\Alex\Desktop\frame_';
sufix='.png';
spi=20;
spf=20;
p=1;
k=64;

for m=spi:spf
  
  I1temp=imread([path num2str(m) sufix]);
  I1=I1temp(:,:,1);
end
% I1=m1(:,:,1);


I1=double(I1);
I2=255*I1/(max(max(I1)));
I1=I2>255*exp(-2);
[r,c]=centroide(I1);
  xp1=p:0.1:k;
  yp1=I1temp(r,:,1);
  wx=linspace(-pi,pi,length((yp1))-1);
  YP1=fft(yp1(2:end));
  [a1,ind]=max(abs(YP1(4:round(length(yp1)/2))));
  wxtemp=fftshift(wx);
  wo=wxtemp(ind+3);
  figure,plot(wx,abs(fftshift(YP1)))
grid
Nventana=1000;
% Nventana=Ventana;
Ventana=round(8*(2*pi)/wo);
wxx=linspace(-pi,pi,Nventana); %vector de frecuencias

wxxtemp=fftshift(wxx);
fven=floor(length(yp1)/Ventana)*Ventana;
for k=1:Ventana:fven
   yp2=yp1(k:k+Ventana-1);
   YP2=fft(yp2,Nventana);
   [cof(l) ind]=max(abs(YP2(50:round(Nventana/2))/Ventana));
   
   [cof2(l) ind2]=max(abs(YP2(round(Nventana/2):length(YP2)-50)/Ventana));
   
   cof(l)=cof(l)+cof2(l);
   %cof2(l)=2*cof2(l);
   
   phase(l)=(angle(YP2(ind+49))-angle(YP2(ind2+round(Nventana/2))))/2;
   wos(l)=(wxxtemp(ind+49)-wxxtemp(ind2+round(Nventana/2)))/2;
   fondo(l)=abs(YP2(1)/Ventana);
l=l+1;
%figure,plot(wxxtemp,abs(fftshift(YP2)))
end
figure,plot(wos),title('Frecuencias Fundamentales')
grid
xpp=0:1e-3:Ventana;
for j=1:l-1
    ypp1(j,:)=fondo(j)+cof(j)*cos(wos(j)*xpp+phase(j)-(2*pi/5.7));
    vis(j)=(max(ypp1(j,:))-min(ypp1(j,:)))/(max(ypp1(j,:))+min(ypp1(j,:)));
    
    figure(4),plot(xpp+Ventana*(j-1),ypp1(j,:))
    hold on
    plot(Ventana*(j-1)+1:(Ventana*j),yp1(Ventana*(j-1)+1:(Ventana*j)),'r'),title('Empalme de perfil a la frecuencia Fundamental')
    grid
    pause(0.01);
end

hold off
figure,plot(vis),title('Visibilidad')
grid

