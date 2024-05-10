clc
clear all
close all

Nx = 2^8;
Ny = 2^8;
DeltaX = 3e-6;
    Fsx = DeltaX;                        % Spatial Sampling frequency along X-axis

DeltaY = DeltaX;
    Fsy = DeltaY;                       % Spatial Sampling frequency along Y-axis

x  = -Nx/2*DeltaX:DeltaX:Nx/2*DeltaX;   % X- Axis
y  = -Ny/2*DeltaY:DeltaY:Ny/2*DeltaY;   % Y- Axis

[X Y] = meshgrid(x,y);

Tx = (2*DeltaX);                        % Fringes period along X-axis
Ty = (2*DeltaY);                        % Fringes period along Y-axis

Fx = ((1)/Tx);                          % #frinjes/longitude -- Adjustable Espatial Frequency X-Axis
Fy = ((1)/Ty);                          % #frinjes/longitude -- Adjustable Espatial Frequency Y-Axis

fxy = 2.*cos(0.5.*Fx.*X + 0.5.*Fy.*Y);                       % Objective Function - (Modificar pesos de Fx y Fy para rotar)
    
    figure
    subplot(121)
    imagesc(x,y,fxy)
    h1 = gca;
    h1.YDir = 'normal';
    h1.XLabel.String = 'x[m]';
    h1.YLabel.String = 'y[m]';
    h1.Title.String  = 'f(x,y)';
    axis square


Fjw = fftshift(fft2((fxy)));                                 % Fourier Transform Objective Function
   
dfx = Fsx/length(x);
freqx = -Fsx/2:dfx:Fsx/2;
dfy = dfx;
freqy = freqx;

absFjw = abs(Fjw);
angFjw = angle(Fjw);

    subplot(122)
    imagesc(freqx,freqy,absFjw);
    h0 = gca;

    h0.XAxisLocation = 'origin';
    h0.YAxisLocation = 'origin';
    h0.YDir = 'normal';
    h0.XLabel.String = 'u[m^{-1}]';
    h0.YLabel.String = 'v[m^{-1}]';
    h0.Title.String  = '|F(u,v)|';
    axis square
    
        
[M I] = max(absFjw(:));
[I_row, I_col] = ind2sub(size(absFjw),I);

hold on

line([freqy(I_col) + dfx/2, - (freqy(I_col) + dfx/2)], [freqx(I_row)+ dfx/2, -(freqx(I_row)+ dfx/2)], 'Color', 'g','LineWidth',1.5,'LineStyle','--');

P1 = [freqy(I_col) + dfy/2, freqx(I_row)+ dfx/2];
P2 = [- (freqy(I_col) + dfy/2), -(freqx(I_row)+ dfx/2)];

% Alternative way to address angle inclination problem
% P1 = [    (I_col - length(freqx)/2) + length(freqx)/2 ,   (I_row - length(freqy)/2) + length(freqy)/2]
% P2 = [    -(I_col - length(freqx)/2) + length(freqx)/2 ,  -(I_row - length(freqy)/2) + length(freqy)/2]

m = (P2(1)-P1(1))/(P2(2)-P1(2)) ;

AnglePointsFreq = atan(m)*(180/(pi))            % Measured form F(u,v)
MeanInclination = (90 - AnglePointsFreq) + 90   % Fringes Tilt Angle F(x,y)

