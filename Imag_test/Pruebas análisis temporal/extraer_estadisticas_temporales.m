clear; clc;

%% ===================== CONFIGURACION =====================
carpeta   = 'C:\Users\sangu\Documents\Trabajo de Grado\trabajos_de_grado\Imag_test\Simulación franjas\franjas_sim_temporal_refinada';
prefijo   = 'frame_';
extension = '.png';

frame_ini = 100;
n_frames  = 450;

% ROI fija centrada
ancho_roi = 80;   % columnas
alto_roi  = 48;   % filas

% Muestreo temporal
fps = 30;

% Demodulacion espacial en la ROI
fx_nominal = 1/8;                    % ciclos/pixel
estimate_fx_from_first_valid = false;
use_hann_in_space = true;

% Fluctuacion de fase
detrend_order = 1;                   % 0: quitar media, 1: quitar tendencia lineal

% Guardado
nombre_mat = sprintf('estadisticas_roi_%d_a_%d.mat', ...
                     frame_ini, frame_ini + n_frames - 1);

%% ===================== PREASIGNACION INICIAL =====================
frames_req = frame_ini : frame_ini + n_frames - 1;
n_req = numel(frames_req);

mu = [];
sd = [];

frames_valid = zeros(n_req,1);
time_s       = zeros(n_req,1);

phi_wrapped  = zeros(n_req,1);
phi_unwrapped = zeros(n_req,1);
phi_fluct    = zeros(n_req,1);
carrier_amp  = zeros(n_req,1);

m_global     = zeros(n_req,1);
sd_global    = zeros(n_req,1);

roi_box = [];
fx_used = fx_nominal;

count_valid = 0;
Wroi = [];
Hroi = [];

%% ===================== PASO 1: LOCALIZAR PRIMER FRAME VALIDO =====================
first_valid_found = false;

for idx = frames_req
    nombre = sprintf('%s%d%s', prefijo, idx, extension);
    ruta   = fullfile(carpeta, nombre);

    if ~isfile(ruta)
        continue;
    end

    img = imread(ruta);
    if size(img,3) == 3
        img = rgb2gray(img);
    end
    img = im2double(img);

    [H,W] = size(img);

    cx = round(W/2);
    cy = round(H/2);

    half_w = floor(ancho_roi/2);
    half_h = floor(alto_roi/2);

    xmin = max(cx - half_w, 1);
    xmax = min(xmin + ancho_roi - 1, W);

    ymin = max(cy - half_h, 1);
    ymax = min(ymin + alto_roi - 1, H);

    if (xmax - xmin + 1) ~= ancho_roi
        xmin = max(W - ancho_roi + 1, 1);
        xmax = W;
    end
    if (ymax - ymin + 1) ~= alto_roi
        ymin = max(H - alto_roi + 1, 1);
        ymax = H;
    end

    roi = img(ymin:ymax, xmin:xmax);
    [Hroi, Wroi] = size(roi);

    mu = zeros(n_req, Wroi);
    sd = zeros(n_req, Wroi);

    roi_box = [xmin xmax ymin ymax];

    if estimate_fx_from_first_valid
        c0 = mean(roi,1).';
        c0 = c0 - mean(c0);

        if use_hann_in_space
            ws = hann(Wroi, 'periodic');
            c0 = c0 .* ws;
        end

        C0 = fft(c0);
        fgrid = (0:Wroi-1)'/Wroi;

        searchIdx = 2:floor(Wroi/2);
        [~, imax] = max(abs(C0(searchIdx)));
        fx_used = fgrid(searchIdx(imax));
    end

    first_valid_found = true;
    break;
end

if ~first_valid_found
    error('No se encontró ningún frame válido en el rango solicitado.');
end

%% ===================== PASO 2: RECORRIDO PRINCIPAL =====================
x = (0:Wroi-1).';

for k = 1:n_req
    idx = frames_req(k);
    nombre = sprintf('%s%d%s', prefijo, idx, extension);
    ruta   = fullfile(carpeta, nombre);

    if ~isfile(ruta)
        warning('No se encontró %s. Se omite.', ruta);
        continue;
    end

    img = imread(ruta);
    if size(img,3) == 3
        img = rgb2gray(img);
    end
    img = im2double(img);

    roi = img(roi_box(3):roi_box(4), roi_box(1):roi_box(2));

    count_valid = count_valid + 1;
    frames_valid(count_valid) = idx;
    time_s(count_valid) = (idx - frame_ini) / fps;

    % Estadisticas por columna
    for j = 1:Wroi
        col = roi(:,j);
        mu(count_valid,j) = mean(col);
        sd(count_valid,j) = std(col, 1);   % desviacion estandar poblacional
    end

    % Estadisticas globales de ROI
    m_global(count_valid) = mean(roi(:));
    sd_global(count_valid) = std(roi(:), 1);

    % Demodulacion de fase en ROI
    ck = mean(roi, 1).';
    ck0 = ck - mean(ck);

    if use_hann_in_space
        ws = hann(Wroi, 'periodic');
        ck0 = ck0 .* ws;
    end

    ref = exp(-1j * 2*pi * fx_used * x);
    Ak = sum(ck0 .* ref);

    phi_wrapped(count_valid) = angle(Ak);
    carrier_amp(count_valid) = abs(Ak);
end

%% ===================== COMPACTAR RESULTADOS =====================
frames_valid = frames_valid(1:count_valid);
time_s       = time_s(1:count_valid);

mu = mu(1:count_valid, :);
sd = sd(1:count_valid, :);

m_global = m_global(1:count_valid);
sd_global = sd_global(1:count_valid);

phi_wrapped = phi_wrapped(1:count_valid);
carrier_amp = carrier_amp(1:count_valid);

if count_valid < 4
    error('Muy pocos frames válidos para análisis temporal.');
end

phi_unwrapped = unwrap(phi_wrapped);

if detrend_order == 0
    phi_fluct = phi_unwrapped - mean(phi_unwrapped);
elseif detrend_order == 1
    phi_fluct = detrend(phi_unwrapped, 1);
else
    error('detrend_order debe ser 0 o 1.');
end

change_profile_L2 = [0; vecnorm(diff(mu,1,1), 2, 2)];
change_phase_abs  = [0; abs(diff(phi_fluct))];

%% ===================== ESTRUCTURA DE SALIDA =====================
results = struct();

results.info.carpeta = carpeta;
results.info.prefijo = prefijo;
results.info.extension = extension;
results.info.frame_ini = frame_ini;
results.info.n_frames_solicitados = n_frames;
results.info.n_frames_validos = count_valid;
results.info.fps = fps;

results.roi.ancho = Wroi;
results.roi.alto  = Hroi;
results.roi.xmin = roi_box(1);
results.roi.xmax = roi_box(2);
results.roi.ymin = roi_box(3);
results.roi.ymax = roi_box(4);

results.demod.fx_used = fx_used;
results.demod.use_hann_in_space = use_hann_in_space;
results.demod.detrend_order = detrend_order;

results.time.frames_valid = frames_valid;
results.time.time_s = time_s;

results.signals.mu = mu;                         % [nFrames x Wroi]
results.signals.sd = sd;                         % [nFrames x Wroi]
results.signals.m_global = m_global;             % [nFrames x 1]
results.signals.sd_global = sd_global;           % [nFrames x 1]
results.signals.phi_wrapped = phi_wrapped;       % [rad]
results.signals.phi_unwrapped = phi_unwrapped;   % [rad]
results.signals.phi_fluct = phi_fluct;           % [rad]
results.signals.carrier_amp = carrier_amp;       % amplitud carrier
results.signals.change_profile_L2 = change_profile_L2;
results.signals.change_phase_abs = change_phase_abs;

save(nombre_mat, 'results');

fprintf('Guardado %s\n', nombre_mat);
fprintf('Frames validos: %d de %d solicitados\n', count_valid, n_frames);
fprintf('ROI real: %d x %d px\n', Wroi, Hroi);
fprintf('fx usada: %.6f ciclos/pixel\n', fx_used);
