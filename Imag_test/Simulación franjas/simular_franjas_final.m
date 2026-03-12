%---------------- PARÁMETROS GENERALES ----------------%
clear; clc;

rng(1);  % reproducibilidad

nx   = 216;
ny   = 146;

A0   = 0.50;         % nivel medio base
B0   = 0.45;         % contraste base
fx   = 1/8;          % ciclos/pixel
phi0 = 0;

nFrames = 600;
fps     = 30;
t       = (0:nFrames-1)/fps;

%---------------- VIBRACIÓN PRINCIPAL ----------------%
f_vib   = 3;                 % Hz
A_phi   = pi/8;              % amplitud de fase
phi_vib = A_phi*sin(2*pi*f_vib*t);

%---------------- JITTER LENTO ESTACIONARIO ----------------%
rho_phi   = 0.98;            % AR(1), cercano a 1 = lento
sigma_phi = pi/180 * 0.15;
phi_jit   = zeros(1,nFrames);
for k = 2:nFrames
    phi_jit(k) = rho_phi*phi_jit(k-1) + sigma_phi*randn;
end

phi_total = phi0 + phi_vib + phi_jit;

%---------------- MODULACIÓN SUAVE DE CONTRASTE ----------------%
% Mantener pequeña para no volverla artificial
mB = 0.03;   % 3% de modulación
B_t = B0 * (1 + mB*sin(2*pi*f_vib*t + pi/6));

%---------------- DESPLAZAMIENTO ESPACIAL PEQUEÑO ----------------%
Ax = 0.10;   % pixeles
dx = Ax*sin(2*pi*f_vib*t + pi/8);

%---------------- RUIDOS ----------------%
sigmaAdd   = 0.115;   % bajar primero para validar PSD
varSpeckle = 0.003;   % bajar primero para validar PSD

%---------------- CAMPO DE ILUMINACIÓN FIJO ----------------%
[x, y] = meshgrid(0:nx-1, 0:ny-1);

illum = 1 + 0.02*imgaussfilt(randn(ny,nx), 12);  % no cambia en el tiempo
illum = illum / mean(illum(:));

outFolder = 'franjas_sim_temporal_refinada';
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

for k = 1:nFrames

    xk = x + dx(k);
    Ik = A0 + B_t(k) * cos(2*pi*fx*xk + phi_total(k));

    % campo fijo de iluminación
    Ik = Ik .* illum;

    % recorte suave antes del ruido
    Ik = min(max(Ik,0),1);

    % ruido aditivo
    Ik = Ik + sigmaAdd*randn(size(Ik));

    % recortar antes de imnoise, sin renormalizar por frame
    Ik = min(max(Ik,0),1);

    % ruido multiplicativo simple
    Ik = imnoise(Ik, 'speckle', varSpeckle);   % J = I + nI
    Ik = min(max(Ik,0),1);

    % sin gamma al inicio: evita armónicos artificiales
    I_uint8 = im2uint8(Ik);

    filename = fullfile(outFolder, sprintf('frame_%03d.png', k));
    imwrite(I_uint8, filename);
end
