clc
clear
clearvars
% Parámetros de la rejilla de franjas
nx   = 216;          % ancho en píxeles
ny   = 146;          % alto en píxeles
A    = 0.5;          % nivel medio (0–1)
B    = 0.5;          % contraste (0–1)
f    = 1/8;          % frecuencia espacial (ciclos por píxel en x)
phi0 = 0;            % fase inicial

% Parámetros del "video"
nFrames      = 150;
phaseDrift   = 2*pi*0.001;   % cambio de fase determinista por frame
maxShiftPix  = 0.2;          % vibración aleatoria máxima (subpíxel)

% Ruido
sigmaAdd     = 0.03;         % desviación estándar ruido gaussiano aditivo
varSpeckle   = 0.02;         % varianza ruido multiplicativo (speckle)

outFolder = 'franjas_sim';
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

% Mallado
[x, y] = meshgrid(0:nx-1, 0:ny-1); %#ok<ASGLU>

for k = 1:nFrames

    % 1) Fase del frame k (pequeño drift para que no sean idénticos)
    phi_k = phi0 + (k-1)*phaseDrift;

    % 2) Vibración: desplazamiento subpíxel aleatorio en x
    dx = maxShiftPix*(2*rand-1);   % uniforme en [-maxShiftPix, maxShiftPix]
    x_shifted = x + dx;

    % 3) Patrón de franjas ideal (normalizado 0–1)
    I = A + B * cos(2*pi*f*x_shifted + phi_k);  % rango teórico [-B+A, B+A]
    I = (I - min(I(:))) / (max(I(:)) - min(I(:))); % normalizar exacto 0–1

    % 4) Ruido aditivo gaussiano
    noise_add = sigmaAdd * randn(size(I));
    I_noisy = I + noise_add;

    % 5) Ruido multiplicativo tipo speckle
    I_noisy = imnoise(im2double(mat2gray(I_noisy)), 'speckle', varSpeckle);

    % 6) Limitar a [0,1] y convertir a uint8
    I_noisy = min(max(I_noisy,0),1);
    I_uint8 = im2uint8(I_noisy);

    % 7) Guardar imagen
    filename = fullfile(outFolder, sprintf('frame_%03d.png', k));
    imwrite(I_uint8, filename);

end
