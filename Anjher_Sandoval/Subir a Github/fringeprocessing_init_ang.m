
function fringeprocessing_init_ang(spi,spf,path_imag)
%path='D:\Trabajo tesis\algoritmos\Imagenes\frames_ambiente\frame_';
%path='C:\Users\Alex\Desktop\potencia con parametros\oscuro-CERCA\potencia cerca\frame_';
%path='D:\Trabajo tesis\algoritmos\Imagenes\1_cT_79\frame_';
%path='D:\Trabajo tesis\algoritmos\Imagenes\1_cT_40\frame_';
path=['./' path_imag '/frame_'];
%path='D:\Trabajo tesis\algoritmos\Imagenes\atm_congelada\frame_';
%%IMAGENES PARA ESTUDIO DE ATMOSFERA CONGELADA
sufix='.png';
% spi=1;
% spf=20;
%%%%%%%%Ir=THg(path,sufix,spi,spf);
mg=[];
h=[];
long_gamma=[];
%It=zeros(260,256);
It=zeros(1200,1600);
for m=spi:spf
  
  I1temp=imread([path num2str(m) sufix]);
  I1=I1temp(:,:,1);
  I1=double(I1);

  It=It+I1;
end

It=It/(spf-spi+1);

FIL2=fftshift(fft2(It));%Transformada de Fourier de f
%P=FIL2;
%P(590:604,524:538)=0;
%P(588:602,1066:1080)=0;
 If=abs(ifft2(fftshift((FIL2))));
 [yc,xc]=centroide(If);%centroide para recortar la imagen acorde a la energ�a
          %figure(1),colormap('gray'),imagesc(It); 
         
          
          
for m=spi:spf

  I1temp=imread([path num2str(m) sufix]);
  I1=I1temp(:,:,1);
  I1=double(I1);
FIL=fftshift(fft2(I1));%Transformada de Fourier de f

%P1=FIL;
 %P1(615:631,728:744)=0;
If1=abs(ifft2(fftshift((FIL))));
 
%                      figure(2), subplot(2,2,1)
%                      imagesc(((If1)))
%                     title('Imagen filtrada')
%                     
%centroide para recortar la imagen acorde a la energ�a
                    
%I=imcrop(I1,[400 287 672 672]); %[x y anchox anchoy] valores tomados del centroide

%                       subplot(2,2,2)
%                       imagesc(I)
%                       title('Imagen filtrada recortada')
%Normalizaci�n


N=uint8(255*double(I1)/double(max(max(I1))));
%                        subplot(2,2,2)
%                        colormap('gray'), imagesc(N)
% title('Imagen normalizada')

b=im2bw(N,graythresh(N));
% b=im2bw(N,1/exp(2));      %Para imagen d�nde el spot est� completo
%Deteccion de bordes
bn=edge(b,'canny');    
%                        subplot(2,2,3)
%                        imagesc(bn)
%                        title('Detecci�n de bordes')
%                        subplot(2,3,4)
%                        colormap('gray'),imagesc(I1);
s=35;
[H,T,R]=hough(bn,'Theta',-2:0.01:1.99);
peaks=houghpeaks(H,s);
lines1=houghlines(bn,T,R,peaks);

% hold on
max_len=0;
grados=[];
gamma=[];
distancia=zeros(1,length(lines1));
    xa=[];xb=[];ya=[];yb=[];
   for k=1:length(lines1)
      xy=[lines1(k).point1; lines1(k).point2];
%                              plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% dibuja el principo y el final de cada segmento
%                               plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%                              plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
     s=lines1(k);
     xa(k)=s.point1(1,1);
     ya(k)=s.point1(1,2);
     xb(k)=s.point2(1,1);
     yb(k)=s.point2(1,2);
     %d=sqrt((xa-xb)^2+(ya-yb)^2);
     gamma1=s.theta(1)*pi/180;
     gamma=[gamma;gamma1];
   
% dibuja el segmento

    len=norm(lines1(k).point1 - lines1(k).point2);
        if( len > max_len)
        max_len = len;
%         xy_long = xy;
        e=k;
        end
    
   end
   
  perfilVertical(m,:)=I1temp(:,1,1);
  %h{m}=gamma;
  h=[h;gamma];  
  [maxh(m), frame(m)]=max(mean(gamma));
  long_gamma=[long_gamma; length(gamma)];
  %long_gamma(1)=0;
  valor_orient(m)=mean(gamma);
  num_lines(m)=length(lines1);
  xapromedio=mean(xa);
  yapromedio=mean(ya);
  xbpromedio=mean(xb);
  ybpromedio=mean(yb);
  
end

%save('Transformada_de_Hough_1-200.mat');
 
 
 %%
 %load('Transformada_de_Hough,.mat');

%  sc=(mean((Ing./K).^2)-(mean(Ing./K).^2))./((mean(Ing./K).^2))

save(['./' path_imag '/proc2_vertical_' path_imag '_' num2str(spi) '_' num2str(spf) '.mat'])

