
function takelines(spi,spf,path_imag)
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
end

 for y=1:length(lines_new)
  orientacion(y)=atan((lines_new(y).point1(:,2)-lines_new(y).point2(:,2))/(lines_new(y).point1(:,1)-lines_new(y).point2(:,1)));
  y1(y)=lines_new(y).point1(:,2);
  x1(y)=lines_new(y).point1(:,1);
  y2(y)=lines_new(y).point2(:,2);
  x2(y)=lines_new(y).point2(:,1);
 end 

save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2_lines_new.mat'],'y1','y2','x1','x2','lines_new');