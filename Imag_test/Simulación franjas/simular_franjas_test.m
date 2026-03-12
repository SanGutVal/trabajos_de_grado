clc
clear
clearvars

%---------------- PARÁMETROS GENERALES ----------------%
nx   = 216;          % ancho en píxeles
ny   = 146;          % alto en píxeles

A    = 0.5;          % nivel medio (0–1)
B    = 0.5;          % contraste (0–1)
fx   = 1/8;          % frecuencia espacial ciclo/píxel en x (ajusta a tu patrón)
phi0 = 0;            % fase base

nFrames = 600;
fps     = 30;        % frames por segundo
t       = (0:nFrames-1)/fps;   % vector de tiempo en segundos

%---------------- VIBRACIÓN (SEÑAL) ----------------%
f_vib      = 3;          % frecuencia de vibración [Hz] (ejemplo)
A_phi      = pi/6;       % amplitud de modulación de fase (rad)
phi_vib    = A_phi * sin(2*pi*f_vib*t);   % fase sinusoidal temporal

%---------------- FASE LENTA ALEATORIA (DRIFTS) ----------------%
% Caminata aleatoria suave: low‑freq noise en fase
sigma_rw   = pi/180 * 0.5;   % paso rms ~0.5 grados
phi_rw     = zeros(1,nFrames);
for k = 2:nFrames
    phi_rw(k) = phi_rw(k-1) + sigma_rw*randn;
end

% Fase total por frame
phi_total = phi0 + phi_vib + phi_rw;  % 1×nFrames

%---------------- RUIDOS DE INTENSIDAD ----------------%
sigmaAdd   = 0.03;   % desviación estándar ruido gaussiano aditivo
varSpeckle = 0.01;   % varianza ruido multiplicativo speckle

%---------------- GEOMETRÍA DEL PATRÓN ----------------%
[x, y] = meshgrid(0:nx-1, 0:ny-1);

outFolder = 'franjas_sim_temporal';
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

for k = 1:nFrames

    % 1) Fase del frame k (incluye vibración periódica + drift lento)
    phi_k = phi_total(k);

    % 2) Patrón base (sin vibración espacial adicional)
    I = A + B * cos(2*pi*fx*x + phi_k);

    % 3) Opcional: hacer barras más "cuadradas" (no estrictamente necesario)
    gamma = 0.9;                    % < 1 aumenta contraste en claros
    I = (I - min(I(:))) / (max(I(:)) - min(I(:))); % normalizar 0–1
    I = I.^gamma;

    % 4) Ruido aditivo gaussiano blanco
    I_add = I + sigmaAdd*randn(size(I));

    % 5) Ruido multiplicativo tipo speckle (independiente por frame)
    I_speckle = imnoise(mat2gray(I_add), 'speckle', varSpeckle);

    % 6) Recorte a [0,1] y conversión a uint8
    I_speckle = min(max(I_speckle,0),1);
    I_uint8   = im2uint8(I_speckle);

    % 7) Guardar imagen
    filename = fullfile(outFolder, sprintf('frame_%03d.png', k));
    imwrite(I_uint8, filename);
end
