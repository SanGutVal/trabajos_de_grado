close all
clear all
clc
% Convertir videos a frames (Entra un .avi y salen directorios con
% imagenes)

%Se crean los directorios para insertar las imagenes por cada prueba
mkdir('AMB'); 
% mkdir('T120V2');
mkdir('T90V');
mkdir('T60V');
mkdir('T120V');
mkdir('HMD85');
mkdir('HMD75');
mkdir('AIRE3_5');
mkdir('AIRE4_5');
mkdir('AIRE5_5');
mkdir('AIRE6_5');

ss=getenv('computername');
compname=str2num([ss(length(ss)-1) ss(length(ss))])

%Ambiente, pruebas sin activar funciones en el turbulador.
V=VideoReader('.\31-03_amb.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\AMB\frame_',num2str(i),'.png']);
end

% disp('120v')
% %Temperatura 120V2
% V=VideoReader('.\TEST_120V_2.avi');
% nFrames=V.NumberOfFrames;
% 
% for i=75*(compname-1)+1:75*compname
%     disp(i)
%     I=read(V,i);
%     D=I(:,:,1);
%     
%     imwrite(D,['.\T120V2\frame_',num2str(i),'.png']);
% end

disp('90V')
%Temperatura 90V
V=VideoReader('.\31_03_90V.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\T90V\frame_',num2str(i),'.png']);
end

disp('60V')
%Temperatura 60V
V=VideoReader('.\31_03_60V.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\T60V\frame_',num2str(i),'.png']);
end

disp('120Vm')
%Temperatura 120V
V=VideoReader('.\31_03_120V.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\T120V\frame_',num2str(i),'.png']);
end
disp('Hum75%')
V=VideoReader('.\75%_hum.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\HMD75\frame_',num2str(i),'.png']);
end
disp('Hum85%')
%Humedad 85%
V=VideoReader('.\hum_04_04_85%.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\HMD85\frame_',num2str(i),'.png']);
end
disp('Aire3.5%')
%Aire 3.5 m/S
V=VideoReader('.\aire3_5.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\AIRE3_5\frame_',num2str(i),'.png']);
end
disp('Aire4.5%')
%Aire 4.5 m/S
V=VideoReader('.\aire4_5.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\AIRE4_5\frame_',num2str(i),'.png']);
end
disp('Aire5.5%')
%Aire 5.5 m/S
V=VideoReader('.\aire5_5.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\AIRE5_5\frame_',num2str(i),'.png']);
end
disp('Aire6_5')
%Aire 6.5 m/S
V=VideoReader('.\aire6_5.avi');
nFrames=V.NumberOfFrames;

for i=75*(compname-1)+1:75*compname
    disp(i)
    I=read(V,i);
    D=I(:,:,1);
    
    imwrite(D,['.\AIRE6_5\frame_',num2str(i),'.png']);
end

