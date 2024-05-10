function [filtro, filtroa]=anillof(A,long_aro,umbral,imtam)
%Funcion para crear anillos cocentricos y filtros pasabanda 2D en frecuencia
%
%FILTRO=ANILLO(A,LON_ARO,UMBRAL,IMTAM)
%
%A hace referencia a cuantos anillos se quiere (entre más pequeño menos
%frecuencia y menos anillos). Valor por defecto para filtros 0.3
%
%long_aro es la longitud de paso del filtro, hace referencia al tamaño
%(cuadrado) de la imagen efectiva del filtro
%
%umbral  es el umbral que se desea para el diseño del filtro
%
%imtam es el tamaño final de la imagen a filtrar

tamsize=size(imtam);

if (tamsize(2) < 2)
    disp('Error: el tamaño de la imagen no debe ser unitario');
elseif((long_aro > imtam(1)) || (long_aro > imtam(2)))
    disp('Error: el tamaño de la imagen debe ser mayor al del umbral');
end

x = linspace(-pi, pi, long_aro);
[xx,yy] = meshgrid(x);
im3 = sin(A*(xx.^2 + yy.^2));
for k=1:long_aro
    for m=1:long_aro
        if im3(k,m)>umbral
            im1(k,m)=1;
        else
            im1(k,m)=0;
        end
    end
end
filtro=zeros(imtam);
filtroa=zeros(imtam);
X=ceil(imtam(1)-long_aro)/2;
Y=ceil(imtam(2)-long_aro)/2;
filtro(X+1:X+long_aro,Y+1:Y+long_aro)=im1;
filtroa(X+1:X+long_aro,Y+1:Y+long_aro)=im1.*im3;