function TT=Polyfit2D_1(mat,Mask,x,y,N)
warning off
%ESTA FUNCION USA VECTORES O MATRICES x,y COMO POSICIONES y LA MATRIZ mat.
%TT=polyfit2D_1(mat,Mask,x,y,N);
%Function pour obtenir le polinome de N degree 2D selon les suivantes ordres
%
%xN xN-1y xN-2y2.....xyN-1 yN
%.......
%x4 x3y x2y2 xy3 y4
%x3 x2y xy2 y3
%x2 xy y2
%x y
%1
%La matrix es vecteurisée dans le sens verticalle. On utilise Mask pour valider les points 
%Le vecteur de sortie TT donne les coefficients:
%m=length(TT);
%TT(m)=const
%TT(m-1)=y;  TT(m-2)=x;
%TT(m-3)=y2;  TT(m-4)=xy; TT(m-5)=x2;
%......
width=size(mat,2);
height=size(mat,1);
TT=[];
if ischar(x) || ischar(y) || ischar(Mask) || ischar(mat) || ischar(N)
   errordlg('Datos de entrada no validos...'); 
   return;
end
if (size(x,1)*size(x,2) == 1) || (size(y,1)*size(y,2) == 1)%x,y escalar
   errordlg('dimension de x y y  no valida.');
   return;
end

if (size(N,1)*size(N,2) ~= 1 || N<=0)%N debe ser escalar positivo
    errordlg('dimension de N  no valida.');
   return;
end

if (size(mat)~=size(Mask))
    errordlg('dimension de matrices  no valida.');
   return;
end

if isvector(x) && isvector(y)
    if (length(x)~=width || length(y)~=height)
        errordlg('dimension de matrices  no valida.'); 
        return;
    end
    [X,Y]=meshgrid(x,y);
    x=[];y=[];
elseif ( size(x)==size(y) && size(x)==size(mat) ) %opcion x,y matrices 
    X=x;Y=y;
    x=[];y=[];
else
    errordlg('dimension x y  no valida.'); 
    return;
    
end

LM=reshape(Mask,[width*height 1]);
indVal=find(LM==1);clear LM;
clear Mask;
LX=reshape(X,[width*height 1]);X=[];
X=LX(indVal);clear('LX');
LY=reshape(Y,[width*height 1]);Y=[];
Y=LY(indVal);clear LY;
LI=reshape(mat,[width*height 1]);
I=LI(indVal);clear LI;clear mat;


if N>=7 && length(I)>=500000
    disp('    polyfit2D_1.m:reduciendo el numero de puntos al 40%');
    n=round(0.4*length(I));
   ind=ceil(length(I).*rand(n,1)); 
   X=X(ind);Y=Y(ind);I=I(ind);
end
NDA=(N+1)*(N+2)/2;
A=zeros(length(Y),NDA);
cont=1;clear indVal;
for n=N:-1:0
    for l=n:-1:0
        a=(X.^l).*(Y.^(n-l));
        A(:,cont)=a;
        cont=cont+1;
    end
end
clear X;clear Y;clear a;
val=A'*A;
b=A'*I;
% TT=inv(A'*A)*A'*I;
TT=val\b;


clear A;clear I;
warning on