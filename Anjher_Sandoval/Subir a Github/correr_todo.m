clear all
disp('15PSI')
%load('alineacion1-200.mat');
tic;

path2='15PSI';
tipovar='PRE';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir(:,:,1:75);
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 aaaa1=unwrap(aaaa);
 figure(1),errorbar(aaaa1(2:end),stdd(2:end))
 hold
 plot(aaaa1(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 drawnow();
 pause(1);
 saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
 drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

toc
clear all
tic;
%% Second test
path2='20PSI';
disp(path2);
tipovar='PRE';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

toc
clear all
tic;
%% Third test
path2='25PSI';
disp(path2);
tipovar='PRE';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');

     drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

toc
clear all
tic;
%% Fourth test
path2='AIRE3_5';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
     drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %%
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

     

toc
clear all
tic;
%% Fifth test
path2='AIRE4_5';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 drawnow();
 pause(1);

 saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');

     drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);

saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

toc
clear all
tic;
%% Sixth test
path2='AIRE5_5';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

toc
clear all
tic;
%% Seventh test
path2='AIRE6_5';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

     
toc
clear all
tic;
%% Eigth test
path2='AMB';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

toc
clear all
tic;
%% Ninth test
path2='HMD';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

toc
clear all
tic;
%% Tenth test
path2='HMD58';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
drawnow();
 pause(1);
 saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');

toc
clear all
tic;
%% Eleventh test
path2='HMD75';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');



toc
clear all
tic;
% Twelveth test
path2='HMD85';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');




toc
clear all
tic;
%% Thirtenth test
path2='T60V';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');




toc
clear all
tic;
%% Fourtenth test
path2='T90V';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');



toc
clear all
tic;
%% Fiftenth test
path2='T120V';
disp(path2);
tipovar='';
genericpath=[path2 '/proceso' tipovar path2];
load([genericpath '_1_75.mat']);
angletFAu_new=angletFAu;
index_new=index;
lines_new=lines;
Ir_new=Ir;
for kp1=76:75:2000
    load([genericpath '_' num2str(kp1) '_' num2str(kp1+74) '.mat']);
    angletFAu_new=[angletFAu_new angletFAu];
    index_new=[index_new index];
    long_lines_new=length(lines_new);
    lines2=lines;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
    [lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
    Ir_new(:,:,kp1:kp1+75-1)=Ir;
end
load([genericpath '_2026_2030.mat']);
angletFAu_new=[angletFAu_new angletFAu];
index_new=[index_new index];
long_lines_new=length(lines_new);
lines2=lines;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point1]=lines2(:).point1;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).point2]=lines2(:).point2;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).theta]=lines2(:).theta;
[lines_new(long_lines_new+1:long_lines_new+length(lines2)).rho]=lines2(:).rho;
Ir_new(:,:,2026:2030)=Ir(:,:,1:5);
%Orientacion de las franjas
angletFAu=angletFAu_new;
index=index_new;
lines1=lines_new;
clear lines_new
Ir=Ir_new;
clear Ir_new
save([path2 '/' tipovar path2 '_vars_index_and_angletFAu.mat'],'index','angletFAu');

%Primera tarea, verificaci�n de los �ngulos de Hough por medio de la
%Transformada de Fourier

 for t=1:length(angletFAu(1,:))
     
     angletFAu(:,t)=unwrap(angletFAu(:,t));
     
     dipix(t)=(angletFAu(length(angletFAu(:,1)),t)-angletFAu(2,t)); %diferencia entre el �ltimo y primer fase inicial en radianes
     %SE toma la frecuencia del perfil 5 (index(5,t), para determinar el
     %periodo
     dipix(t)=dipix(t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1))));%diferencial entre el ultimo y primer fase inicial en metros
     stdd(t)=std(angletFAu(:,t)/(2*pi*index(5,t)*(1/(2.8e-6)/length(angletFAu(:,1)))));
     ang(t)=atan(2.8e-6*length(angletFAu(1,:))/dipix(t)); %angulo real en radianes
     
        %*ones(size(ang));
 end
 aaaa=mean(angletFAu);
 figure(1),errorbar(aaaa(2:end),stdd(2:end))
 hold
 plot(aaaa(2:end),'r')
 grid
 xlabel('Samples');
 ylabel('Initial Phase [Rad]');

 err=std(ang);
 legend('Spatial Standar deviation ','temporal variation')
 
drawnow();
 pause(1);
 saveas(gcf,[path2 '/' tipovar path2 '_timeandspatial.fig']);
 if mean(ang)>0
     
     figure(2), 
     plot(ang-pi/2), title('Determinaci�n de la fase inicial')
     grid;
     xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 else 
     figure(2), plot(ang+pi/2), title('Determinaci�n de la fase inicial')
 grid;
  xlabel('Samples');
     ylabel('Initial Phase [Rad]');
 end
 
 drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_initial_phase.fig']);

%  errorbar(ang,err);
%  
%
%Segunda tarea, validez de paralelismo

 for y=1:length(lines1)
  orientacion(y)=atan((lines1(y).point1(:,2)-lines1(y).point2(:,2))/(lines1(y).point1(:,1)-lines1(y).point2(:,1)));
  y1(y)=lines1(y).point1(:,2);
  x1(y)=lines1(y).point1(:,1);
  y2(y)=lines1(y).point2(:,2);
  x2(y)=lines1(y).point2(:,1);
 end 

%Gr�ficas para entender los �ngulos y el porqu� de sus diferencias
orientacion=orientacion*180/pi;
 figure(4),subplot(3,1,1), plot(orientacion),title('SLOPES')% 
 grid;
%      xlabel('Lines');
     ylabel('Orientacion');

     subplot(3,1,2), plot(y2-y1),title('Y-pixel differences (Y_2-Y_1)')% 
 grid;
%      xlabel('Lines');
     ylabel('Y-Differences');
 
      subplot(3,1,3), plot(x2-x1),title('X-pixel differences (X_2-X_1)')% 
 grid;
     xlabel('Lines detected (Hough transform)');
     ylabel('X-Differences');
     
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Orientation.fig']);
save([path2 '/' tipovar path2 '_vars_y1_y2_x1_x2.mat'],'y1','y2','x1','x2');

 %
 %Tercera tarea Periodo o Frecuencia espacial de las franjas
for j=1:2030
       
 for x=5:length(Ir(:,1,1))
vt=find(Ir(x,:,j)>mean(Ir(x,:,j)));
dvt=diff(vt);
tr=x-4;
temp=2*dvt(find(dvt>1));
long_per(tr,j)=length(temp);
Per(tr,1:length(temp),j)=temp;
 end
 
 
%plot(Per(tr,j,1:length(temp)))
%  hold on
end
figure(5),imagesc(Per(:,:,2030)),colorbar
xlabel('ZeroCrossing');
ylabel('Profiles');
drawnow();
 pause(1);
saveas(gcf,[path2 '/' tipovar path2 '_Periods.fig']);
save([path2 '/' tipovar path2 '_vars_long_per_and_period.mat'],'long_per','Per');


