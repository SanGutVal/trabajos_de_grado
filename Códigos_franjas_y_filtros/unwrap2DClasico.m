%Function pour calculer la phase continue à partir de la phase discontinue et la masque
%de validation. Como critere de traitement de bruit il utilise la premiere derive.
%PI est le point initial de correction
%function [PhaseC,MaskF,tiempo]=unwrap2DClasico(PI,phased,Mask)
%function [PhaseC,MaskF,tiempo]=unwrap2DClasico(PI,phased,Mask,fig1)
%fig1=Visualization de correction
function [PhaseC,MaskF,tiempo]=unwrap2DClasico(varargin)

[PI,phased,Mask,fig1]= parseInputs(varargin{:});

if Mask(PI(2),PI(1))==0
   disp('PI en dehors de la masque');
   errordlg('PI en dehors de la masque','ERROR:unwrap2DClasico');
   return;
end


W=size(phased,2);H=size(phased,1);
Wi=W;Hi=H;
MaskT=zeros(H+2,W+2);phasedT=zeros(H+2,W+2);
MaskT(2:1+H,2:1+W)=Mask;phasedT(2:1+H,2:1+W)=phased;
PI(1)=PI(1)+1;PI(2)=PI(2)+1;
W=size(phasedT,2);H=size(phasedT,1);
phasec=phasedT;

xo=[-1,-1,-1,0,0,1,1,1];
yo=[-1,0,1,-1,1,-1,0,1];
NP=sum(sum(MaskT));
XY=zeros(NP,2);
XY(1,1)=PI(1);XY(1,2)=PI(2);

MaskT(PI(2),PI(1))=0;
phasec(PI(2),PI(1))=phasedT(PI(2),PI(1));
PIN=PI;
cont=1;cont1=1;opct=1;
NV=round(NP/30);
tic,
while (opct)
   PCI=phasec(PIN(2),PIN(1));
   PDI=phasedT(PIN(2),PIN(1));
   for i=1:8
      px=xo(i)+PIN(1);
      py=yo(i)+PIN(2);
      if(MaskT(py,px)==1)
         PDC=phasedT(py,px);
         D=(PDC-PDI)/2/pi;
         phasec(py,px)=PCI+2*pi*(D-round(D));
         XY(cont1+1,1)=px;
         XY(cont1+1,2)=py;
         cont1=cont1+1;
         MaskT(py,px)=0;
      end
   end
   cont=cont+1;
   if(cont>cont1)
      opct=0;
   else
         PIN(1)=XY(cont,1);PIN(2)=XY(cont,2);
   end
   
   if ~isempty(fig1)
       if(rem(cont,NV)==0)
         figure(fig1);imagesc(phasec);colormap gray;
         drawnow;
       end
   end
   
end
tiempo=toc;
PhaseC=phasec(2:1+Hi,2:1+Wi);
Mask1=MaskT(2:1+Hi,2:1+Wi);
MaskF=(1-Mask1).*Mask;


function [PI,phased,Mask,fig1] = parseInputs(varargin)

    iptchecknargin(3,4,nargin,mfilename);
    PI=[];phased=[];Mask=[];fig1=[];
    PI = varargin{1};
    iptcheckinput(PI, {'double'}, {'real', 'vector', ...
        'finite','positive'}, mfilename, 'PI', 1);
    phased = varargin{2};
    iptcheckinput(phased, {'double'}, {'real', '2d', 'nonsparse', 'nonempty',...
        'finite'},mfilename, 'phased', 2);
    Mask = varargin{3};
    iptcheckinput(Mask, {'numeric','logical'}, {'real', '2d', 'nonsparse', 'nonempty',...
        'finite'},mfilename, 'Mask', 3);
    if nargin >3
        fig1 = varargin{4};
    end
