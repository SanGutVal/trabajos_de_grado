function matA=PolyVal2D_1(TT,x,y,N)
% x,y son vectores que definen las coordenadas de los puntos donde se va a evaluar el polinomio. 
%Esta funcion asume que se desea evaluar la zona definida por x y y para
%generar una matriz de valores. Generalmente se usa asuminedo x=1:Num Colu
%y=1:Num Filas.
%Function pour evaluer le polynome TT de degree N avec x,y 
%xN xN-1y xN-2y2.....xyN-1 yN
%.......
%x4 x3y x2y2 xy3 y4
%x3 x2y xy2 y3
%x2 xy y2
%x y
%1

NDA=(N+1)*(N+2)/2;

validateattributes(TT, {'double'}, {'real','vector','nonsparse', 'nonempty',...
    'finite'},mfilename, 'TT', 1);
validateattributes(x, {'double'}, {'real','vector','nonsparse', 'nonempty',...
    'finite'},mfilename, 'x', 2);
validateattributes(y, {'double'}, {'real','vector','nonsparse', 'nonempty',...
    'finite'},mfilename, 'y', 3);
validateattributes(N, {'double'}, {'real', 'scalar', 'integer', 'positive',...
    'finite'},mfilename, 'N', 4);


if (length(TT)~=NDA)
   ui=errordlg('Dimension de TT no valida','PolyVal2D_1');
   waitfor(ui);
   return;
end
[X,Y]=meshgrid(x,y);
if N>=7 && length(X(:))>=500000
    X=[];Y=[];
    matA=zeros(length(y),length(x));
    disp('    polyval2D_1.m:evaluacion camino largo');
    for ix=1:length(x)
        valx=x(ix);
        LY=y;
        LX=valx*ones(size(LY));
        A=zeros(length(LY),NDA);
        cont=1;
        for n=N:-1:0
            for l=n:-1:0
                a=(LX.^l).*(LY.^(n-l));
                A(:,cont)=a;
                cont=cont+1;
            end
        end
        LA=A*TT;
        matA(:,ix)=LA;
    end
else    
    
    [height,width]=size(X);
    LX=reshape(X,[width*height 1]);X=[];
    LY=reshape(Y,[width*height 1]);Y=[];
    
    A=zeros(length(LY),NDA);
    cont=1;
    for n=N:-1:0
        for l=n:-1:0
            a=(LX.^l).*(LY.^(n-l));
            A(:,cont)=a;
            cont=cont+1;
        end
    end
    matA=A*TT;
    matA=reshape(matA,[height,width]);
end