
function fringeprocessing(spi,spf,path_imag)
%path='D:\Trabajo tesis\algoritmos\Imagenes\frames_ambiente\frame_';
%path='C:\Users\Alex\Desktop\potencia con parametros\oscuro-CERCA\potencia cerca\frame_';
%path='D:\Trabajo tesis\algoritmos\Imagenes\1_cT_79\frame_';
%path='D:\Trabajo tesis\algoritmos\Imagenes\1_cT_40\frame_';
path=['.\' path_imag '\frame_'];
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
 [yc,xc]=centroide(If);%centroide para recortar la imagen acorde a la energía
          %figure(1),colormap('gray'),imagesc(It); 
         
          
          
for m=spi:spf

  I1temp=imread([path num2str(m) sufix]);
  I1=I1temp(:,:,1);
  I1=double(I1);
  I2(:,:,m-spi+1)=I1temp(:,:,1);
	FIL=fftshift(fft2(I1));%Transformada de Fourier de f

	%P1=FIL;
	 %P1(615:631,728:744)=0;
	If1=abs(ifft2(fftshift((FIL))));
	 
	%                      figure(2), subplot(2,2,1)
	%                      imagesc(((If1)))
	%                     title('Imagen filtrada')
	%                     
	%centroide para recortar la imagen acorde a la energía
						
	%I=imcrop(I1,[400 287 672 672]); %[x y anchox anchoy] valores tomados del centroide

	%                       subplot(2,2,2)
	%                       imagesc(I)
	%                       title('Imagen filtrada recortada')
	%Normalización


	N=uint8(255*double(I1)/double(max(max(I1))));
	%                        subplot(2,2,2)
	%                        colormap('gray'), imagesc(N)
	% title('Imagen normalizada')

	b=im2bw(N,graythresh(N));
	% b=im2bw(N,1/exp(2));      %Para imagen dónde el spot está completo
	%Deteccion de bordes
	bn=edge(b,'canny');    
	%                        subplot(2,2,3)
	%                        imagesc(bn)
	%                        title('Detección de bordes')
	%                        subplot(2,3,4)
	%                        colormap('gray'),imagesc(I1);
	s=35;
	[H,T,R]=hough(bn,'Theta',-2:0.01:1.99);
	peaks=houghpeaks(H,s);
	lines=houghlines(bn,T,R,peaks);
	if (m==spi)
		lines_new=lines;
	else
		long_lines_new=length(lines_new);
		lines2=lines;
		[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
		[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
		[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
		[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
	end
	% hold on
	max_len=0;
	grados=[];
	gamma=[];
	distancia=zeros(1,length(lines));
    
   for k=1:length(lines)
      xy=[lines(k).point1; lines(k).point2];
%                              plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% dibuja el principo y el final de cada segmento
%                               plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%                              plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
     s=lines(k);
     xa=s.point1(1,1);
     ya=s.point1(1,2);
     xb=s.point2(1,1);
     yb=s.point2(1,2);
     %d=sqrt((xa-xb)^2+(ya-yb)^2);
     gamma1=s.theta(1)*pi/180;
     gamma=[gamma;gamma1];
   
% dibuja el segmento

    len=norm(lines(k).point1 - lines(k).point2);
        if( len > max_len)
        max_len = len;
%         xy_long = xy;
        e=k;
        end
    
   end

  %h{m}=gamma;
  h=[h;gamma];  
  [maxh(m), frame(m)]=max(mean(gamma));
  long_gamma=[long_gamma; length(gamma)];
  %long_gamma(1)=0;
  
  Ir1=imrotate(I1,180/pi*mean(gamma),'bilinear');
  %Ir(:,:,m-spi+1)=imcrop(Ir1,[4 4 250 250]);
  Ir(:,:,m-spi+1)=imcrop(Ir1,[4 4 1550 1150]);
end

%save('Transformada_de_Hough_1-200.mat');
 
 
 %%
 %load('Transformada_de_Hough,.mat');
%%%%%%%%
[f c temp]=size(Ir);
angletFA=zeros(f,(spf-spi)+1);
for j=1:(spf-spi)+1
    A=Ir(1,:,j); 
    % Búsqueda de primer perfil referencia
    for i=2:f
        tfA(i,:)=fft(Ir(i,:,j)); % Transformada de Fourier de las filas de la imagen
        [maxv,index(i,j)]=max(abs(tfA(i,10:round(length(tfA(i,:))/2-10))));  %Búsqueda de la frecuencia fundamental
        index(i,j)=index(i,j)+9;
        angletFA(i,j)=angle(tfA(i,index(i,j)));%Se encuentra la fase inicial
        A=A+(Ir(i,:,j)); %vector acumulado de perfiles
    end
    angletFAu(:,j)=unwrap(angletFA(:,j));
    Am(j,:)=A/length(Ir(:,1,j)); %Promedio de todos los perfiles (Vector)
    tfAm(j,:)=fft(Am(j,:));
    [maxAm, indexAm(j)]=max(abs(tfAm(j,10:round(length(tfAm(j,:))/2))));    %Se divide sobre dos por la simetría del espectros
    indexAm(j)=indexAm(j)+9;
    angleAm(j)=angle(tfAm(j,indexAm(j)));  %Angulo promedio de todos los perfiles de cada imagen
    angleAm(j)=angleAm(j)+2*pi;

end


% Periodo=6.4; %Franjas ct_79, ct_40 y st, Separación de espejos=7.2cm, Z=2.04m
Periodo=6.54; %franjas de atm_congelada, Separación de espejos=7.6cm, Z=2.2m
paso_ang=2*pi/Periodo;
Im_alin=zeros(size(Ir));
sam=zeros(f,(spf-spi)+1); %matriz de ceros para llenar con valores de angulos promedio de cada imagen
for j=1:(spf-spi)+1
    for i=1:f
        angle_rec(i,j)=angleAm(j)-angletFA(i,j); %Desviación del angulo de cada imagen con respecto al promedio
        sam(i,j)=round(angle_rec(i,j)/paso_ang);
        if sam(i,j)==0 
        Im_alin(i,1:length(Ir(1,:,j))-sam(i,j),j)=Ir(i,sam(i,j)+1:end,j);
        else
        Im_alin(i,1:length(Ir(1,:,j))-sam(i,j)+1,j)=Ir(i,sam(i,j):end,j);
        end
    % plot(Im_alin(i,1:length(Ir(1,:))-sam(i)+1)),title('Todos los perfiles')
    end
    Am2(:,j)=mean(Im_alin(:,:,j));
    B1(:,:,j)=repmat(Am2(:,j)',f,1);    %Imagen ideal a partir del perfil promedio

    %figure(3),colormap('gray'),imagesc(B1(:,:,j)), title('Reconstrucción franjas a partir de perfil promedio')
    %drawnow;
    
    dif(:,1:f-3,j)=Im_alin(:,1:f-3,j)-B1(:,1:f-3,j);
    flI(:,j)=mean(dif(:,:,j).^2);
    figen(:,j)=flI(:,j)./Am2(1:f-3,j).^2;
    scin(:,j)=figen(:,j);
    %%%
    
    %Intensidad en nivel de gris de la imagen en la camara
        Ing(j)=sum(sum(I2(:,:,j)));
    
       %sc=mean(
end

%save('alineacion1-200.mat');
angle_rec=imcrop(Ir1,[2 2 200 256]);
fas=unwrap(angle_rec);
% figure, plot(fas(:,1)),title('desviación de la fase del promedio para la primera imagen ')% Dibuja la desviación de la fase del promedio para la primera imagen 
% figure, plot(fas(5,:)),title('valor de fase promedio para todas las imagenes con respecto al 5 perfil ')% Dibuja valor de fase promedio para todas las imagenes con respecto al 5 perfil 
% figure, mesh(fas),colorbar, title('Fase promedio de todos los perfiles de todas las imagenes') %Fase promedio de todos los perfiles de todas las imagenes


% %Intensidad real
%     fid = fopen('D:\Trabajo tesis\algoritmos\Imagenes\prueba_31_03\amb\PM1918000.dat');
%     tline = fgetl(fid);
%     c=0;
%     while tline > 0 
%        c=c+1;
% 
%        if c>13 && ~(strcmp(tline,'End of Data'))
%            PW(c-12)=str2num(tline);
%         end
%         tline = fgetl(fid);
%     end
%     fclose(fid);
%     Pre=mean(PW);
% 
%     Ire(j)=Pre*(pi*(11.3e-3/2)^2);
% 
%  
%  K=mean(Ing./Ire);
%  %scintillation
%  sc=(mean((Ing./K).^2)-(mean(Ing./K).^2))./((mean(Ing./K).^2))

clear I2
clear Im_alin
save(['proceso' path_imag '_' num2str(spi) '_' num2str(spf) '.mat'])

       
%      K(j)=Ire(j)/Ing(j);
 
% sc=mean((Ing./K)^2)-(mean(Ing./K)^2)./((mean(Ing./K)^2))

%save sam_ct_79_1901-2000



% % figure,plot(Im_alin(255,:))
% hold
% plot(Am2,'k'),
% plot(Im_alin(25,:),'r')
% plot(Im_alin(450,:),'g')
% grid
% title('Grafica de comparación de perfiles')

%     B1(:,:,j)=imcrop(B1(:,:,j),[1 1 1217 463]);
%Promedio en intensidad

%figure,plot(flI./Am2.^2)
%figure,plot(flI(1:1210)./Am2(1:1210).^2)
% scin=mean(flI(1:1210)./Am2(1:1210).^2)
%figure,plot(scin), title('Scintillation')


%plot(angletFA+2*pi)
%hold
%plot(angleAm*ones(1,length(angletFA)),'r')
%Fluctuación del ángulo

% difAn=angletFA-angleAm*ones(1,length(angletFA));
% 
% for k=1:f
%     flA(k)=mean(difAn(:,k).^2);
% end
% figure,plot(flA(1:464)/angleAm.^2), title('Fluctuación del angulo') 






% figure(2), plot(Mg,'r')
% title('Angulos Promedios en cada imagen')

% VISUALIZACION DE RANGO DE FRANJAS
% figure, plot(abs(If(100,:))),grid
% figure, plot(abs(If(255,:))),grid
% figure, plot(abs(If(400,:))),grid
%                                  fs_x=1/5.3e-6;    
%                                  faxis_x=linspace(-fs_x/2,fs_x/2,length(N(255,:)));
%                                  figure, plot(faxis_x,N(100,:)),grid  


% faxis_x=linspace(-fs_x/2,fs_x/2,length(P(255,:)));
% figure, plot(faxis_x,log(P(255,:))),grid
% axis tight

% fs_y=1/5.3e-6;
% faxis_y=linspace(-fs_y/2,fs_y/2,length(P(:,545)));
% figure, plot(faxis_y,log(P(:,545))),grid
% axis tight
            

