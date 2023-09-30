distance=2.2;
separation=0.07;
lambda=632e-9;
l=1;
%This algorithm tries to find global phase from Takeda algorithm. So, it
%takes the images, executes FFT2, find fundamental frequencies and separate
%them of spectrum, after it applies hamming windowing and return to direct
%space where is calculated phase.

%M=1280; N=1024; %%Anjher
%M = 640; N = 480;
M = 480; N = 272;
%M = 320; N = 240;
%M = 160; N = 120;
t1=1:M;
t2=1:N;
%path='D:\Trabajo tesis\ALEXANDER\Codigos_Finales\15PSI\frame_'; %%Anjher
%path='D:\Anjher\UIS\Trabajo_de_Grado\scripts\imagenes_scripts\prueba_simulada3\ps';
%sufix='.jpg';
spi=107;
spf=107;
for i=spi:spf
   %I=imread([path num2str(i) sufix]);
   I = rejillaDifraccion([M,N], 75);
   D = double (I);
    %D=double(I(:,:,1));
    I1=fftshift(fft2(double(D)));
    if (i==spi)
       % [a1,ind]=max(abs(I1(N/2,M/2 + 5:end)));
        [a1,ind]=max(abs(I1(N/2,M/2 + 5:end)));
        %%[a1,ind]=max(abs(I1(600,805:end)));   %Se toma en la mitad del espectro (512) %%Anjher
        %para verificar dónde aparece la fundamental, también se toman 4 pixeles
        %después de la DC para encontrar el centro de masa.
        
        ind1= M/2 + 4 + ind;
        %ind1=804+ind;              %Anjher
        %Aquí se verifica la posición del centro de masa del lobulo.
        [a2,ind2]=max(abs(I1(:,ind1)));

        %A partir del centro de masa del lobulo se sacan las posiciones
        %alternas y así recortar el lóbulo
        posicini=find(abs(I1(:,ind1))>a2/100);
    %     deltau=round((posicini(end)-posicini(1))/2)+1;
        %Aquí se toma el primer valor separado de fo/3 para poder encontrar el
        %área a filtrar
        deltau=round((ind+4)/3)+1;
    end
    disp(i)
    %Buscamos extraer la imagen sin el fondo continuo
   
    
    %Se busca una copia para llevar sólo los cambios de la fundamental al
    %espacio directo  % Espacio recíproco, no espacio directo --  Omar
    I1takeda=I1;     %Filtro digital en el dominio de fourier
    I1takeda(1:ind2-deltau,:)=0;    % Primera línea del filtro -- Omar
    I1takeda(ind2+deltau:end,:)=0;
    I1takeda(ind2-deltau:ind2+deltau,1:ind1-deltau)=0;
    I1takeda(ind2-deltau:ind2+deltau,ind1+deltau:end)=0;
    i1takeda=ifft2(I1takeda);
    %Aquí se pasa un filtro hanning de suavizado  %Apolización -- Yezid
    I1takeda(ind2-deltau:ind2+deltau-1,ind1-deltau:ind1+deltau-1)=...
    I1takeda(ind2-deltau:ind2+deltau-1,ind1-deltau:ind1+deltau-1).*(hanning(deltau*2)*hanning(deltau*2)');
    i1tak=ifft2(fftshift(I1takeda));   %%Cambio del dominio de fourier (ya enventadano) al espacio directo
    %Esta es la fase para ambas imágenes
    phasei1takeda=angle(i1takeda);
    phasei1tak=angle(i1tak);
    
    %Aquí está la amplitud de esa fase
    magni1takeda=abs(i1takeda);
    magni1tak=abs(i1tak);       %Aquí es donde se ve el filtro hanning.
    %Esto simplemente es para pasar el filtro. En reconstrucción 3D lo
    %utilizan para sólo buscar el objeto.

    
    
    %FIL=fftshift(fft2(D));
%     filtro=xor(anillof(0.3,10,0,[1200 1600]),0);
%     P=FIL.*filtro;
%     N1=(double(magni1tak)/double(max(max(magni1tak))));
%     b1=im2bw(N1,0.1);    
% se=strel('disk',30);
% b11=imerode(b1,se);
% b111=imdilate(b11,se);
%     Mask=b111;
%     MaskOne=ones(size(D));
    Mask=((magni1tak-min(magni1tak(:)))/(max(magni1tak(:))-min(magni1tak(:))))>=1/exp(2);
    Mask2=((magni1tak-min(magni1tak(:)))/(max(magni1tak(:))-min(magni1tak(:))))>0;
  
    se=strel('disk',30);  %% disco convolución
    b11=imerode(Mask,se);   %disco para erosionar la imagen (poner los bordes más delgados)
    Mask1=imdilate(b11,se); %Con la imagen erosionada, la voy a dilatar: adelgazamiento de los bordes, después los hago más anchos
    
    sizeIm=size(D);
    phasei1tak_u=phasei1tak;
    [PhaseC,MaskF,tiempo]=unwrap2DClasico([ind1 ind2],phasei1tak_u,Mask2);  % [ind1 ind2] esto es la frecuencia espacial fundamental, cambié por Mask1, tenía Mask2

     %Trazar líneas de fase a partir del contorno con alta resolución
     %tt=find(Mask1(600,1:end)>0);                  %Anjher
     tt=find(Mask1( N/2 ,1:end)>0);
     s1=contourc(PhaseC,round((PhaseC(round(sizeIm(1)/2),1):PhaseC(round(sizeIm(1)/2),end))*2*pi));
%    s1=contourc(PhaseC.*Mask1,round((PhaseC(round(sizeIm(1)/2),tt(1)):PhaseC(round(sizeIm(1)/2),tt(end)))*2*pi));
      
    figure(1); imagesc(cos(PhaseC)); colormap('gray') % Cambié imagesc(D) por imagesc(PhaceC)
    title(['Interpolación de las franjas. Imagen ' num2str(i) ' de 2031'])
    hold on, plot(s1(1,:),s1(2,:),'r.')
    hold off
    return ;   %OMar
    pos_s1=[];
    
   %Se busca separar cada línea
%     pos_s1=find((s1(2,:))==1200);
    pos_s1=find(diff(s1(2,:))>10);
    %Con este algoritmo se dibujan las líneas
    pix_eras=10;

    flag=0;
    puntoiniY=[];
    puntofinY=[];
    orientacion=[];
    medianXg=[];
 figure(8);
    for k=pix_eras:length(pos_s1)-pix_eras
%         linea=s1(:,pos_s1(1,k):pos_s1(1,k+1));
%          eval(['L',num2str(k),'=linea;']);
%          
%         plot(linea(1,:),linea(2,:),'r.');
%         hold on
%         clear linea
        %El número pix_eras es la cantidad de pixeles que se quitan de los bordes
%       figure (2)
        %Se busca la coordenada en X para poder limitar la frontera de la
        %Máscara
        medianaX(k-pix_eras+1)=round(median(s1(1,pos_s1(k)+pix_eras:pos_s1(k+1)-pix_eras)));
        %Encontrada la coordena en X, busco el perfil vertical en la máscara y los
        %puntos que corresponden a la interferencia
        if(isnan(medianaX(k-pix_eras+1))==0)
            puntosY=find(Mask1(:,medianaX(k-pix_eras+1))>0);
            if (length(puntosY)~=0)
                medianXg=[medianXg median(s1(1,pos_s1(k)+pix_eras:pos_s1(k+1)-pix_eras))];
                %Se encuentran los puntos en la frontera de esa línea
                puntoiniY=[puntoiniY puntosY(1)];
                puntofinY=[puntofinY puntosY(end)];
                %Se limitan los puntos en Y para poder dibujar la máscara
                tempCoory=find(s1(2,pos_s1(k)+pix_eras:pos_s1(k+1)-pix_eras)<puntosY(end));
                puntos_dib=find(s1(2,pos_s1(k)+pix_eras+tempCoory)>puntosY(1));
                plot(s1(1,pos_s1(k)+pix_eras+puntos_dib),s1(2,pos_s1(k)+pix_eras+puntos_dib))
        %        plot(s1(1,pos_s1(k)+pix_eras:pos_s1(k+1)-pix_eras),s1(2,pos_s1(k)+pix_eras:pos_s1(k+1)-pix_eras));
        %        plot(s1(1,pos_s1(k+1)-pix_eras:pos_s1(k+1)-1),s1(2,pos_s1(k+1)-pix_eras:pos_s1(k+1)-1));
        %        lineas(k,:)=s1(:,pos_s1(1,k));
%                orientacion(k-pix_eras+1)=atan((s1(2,pos_s1(k+1)-pix_eras)-s1(2,pos_s1(k)+pix_eras))/(s1(1,pos_s1(k+1)-pix_eras)-s1(1,pos_s1(k)+pix_eras)));
               

%Orientacion Espacial de las franjas 

               [coefY,S,mu]=polyfit(s1(1,pos_s1(k)+pix_eras:pos_s1(k+1)-pix_eras),s1(2,pos_s1(k)+pix_eras:pos_s1(k+1)-pix_eras),1); 
               orientacion=[orientacion atan(coefY(1))];
               f=polyval(coefY,s1(1,pos_s1(k)+pix_eras:pos_s1(k+1)-pix_eras),S,mu);
               figure(2), plot(orientacion)
               
               if(flag==0)
                    grid;hold on;
                    title('Líneas curvas encontradas')
                    disp('Líneas')
                    flag=1;
                end
                drawnow;
    %         pause;
            end
        end
    end
    
%     %Lineas corregidas verdes
%     figure,imagesc(D),colormap(gray),hold on;
%     for g=1:length(pos_s1)-1
%         eval(['linea=L',num2str(g),';']);
%         izq=find(linea(1,:)<mean(mean(linea(1,:)))-8);
%         eval(['L',num2str(g),'(:,izq(1,:))=[];']);
%         eval(['linea=L',num2str(g),';']);
%         plot(linea(1,:),linea(2,:),'g.');
%         
%         clear linea izq
%         
%     end
%     
%     hold off
%   axis tight
    %Aquí estaba intentando el plano de fase, pero realmente no fue
    %posible, si sabes como encontrarlo sería interesante.
%     pendiente1=(PhaseC(round(sizeIm(1)/2),end)-(PhaseC(round(sizeIm(1)/2),1)))/M;
%     curva1=pendiente1*t1+PhaseC(round(sizeIm(1)/2),1);
%     pendiente2=(PhaseC(end,round(sizeIm(2)/2))-(PhaseC(1,round(sizeIm(2)/2))))/N;
%     curva2=pendiente2*t2+PhaseC(1,round(sizeIm(2)/2));

%PERIODO ESPACIAL

for s=1:length(medianXg)-1
pf=diff(medianXg);
orientpro(s)=(mean(orientacion));
pr(s)=pf(s).*sin(orientpro(s));

end
%CONTRASTE

%  wp=fftshift(I1)/(M*N);
%     A(i)=wp(1,1); %Fon22do continuo
%     X=1:1600;
%     B=((D-A(i))./(2*magni1tak+real(i1tak))).*Mask1;
%     SSS=(magni1tak/A(101));
%     figure, imagesc(SSS)
% % for x=0:length(phasei1tak_u(:,1))
%     [peak_valuetf, peak_locationtf] = findpeaks(phasei1tak_u(x,:),'minpeakdistance',2,'minpeakheight',mean(real(phasei1tak_u(x,:))));
%     vt=peak_locationtf;
%     dvt=diff(vt);
%     temp=dvt(find(dvt>1));
%     long_per(x,i)=length(temp);
%     Per(x,1:length(temp))=temp;
%  end


end
%return %Anjher
% figure(3), imagesc(Per)
% title('Períodos para: ')
% ylabel('Perfíl');
% xlabel('Indice de Período');

%edges=1:0.3:15;
edges=1:0.3:15;
countp=histc(pr',edges);
figure(4),bar(edges,sum(countp,2)/sum(sum(countp)),'histc');
axis([0 5 0 0.99]);
title('Periódos para:')
ylabel('Porcentaje');
xlabel('Valor del Período');
grid
% end 
% %Fluctuación de Fase Espacial
% for t=1:length(phasei1takeda(:,1))
%      stdd(t)=std(phasei1takeda(t,:));
% end
%  aa=mean(phasei1takeda);
%  aa1=unwrap(aa);
%  figure(1),errorbar(aa1(2:end),stdd(2:end))
%  axis([0 2030 -1 4]);
%  hold
%  plot(aa1(2:end),'r')
%  grid
%  xlabel('Samples');
%  ylabel('Initial Phase [Rad]');
%  hleg1=legend('Spatial Standar deviation ','temporal variation');
%  set(hleg1,'Location','SouthEast');
%  hold off
% 
% % p = polyfit(s1(1,:),s1(2,:),1)
% [p,S] = polyfit(s1(1,:),s1(2,:),1);
return  %Anjher
