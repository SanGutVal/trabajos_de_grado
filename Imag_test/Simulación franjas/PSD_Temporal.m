%% =========================================================
%  ANALISIS TEMPORAL DE FASE EN SECUENCIA DE FRANJAS
%  Salidas principales:
%   1) fluctuaciones de fase [rad] vs tiempo [s]
%   2) score de cambio entre frames consecutivos
%   3) PSD temporal de la fase con Welch
% ==========================================================
clear; clc;
%close all;

%% ---------------- CONFIGURACION ----------------
frameFolder  = 'franjas_sim_temporal_refinada';   % carpeta de frames
framePattern = 'frame_*.png';
fps          = 30;

% Frecuencia espacial esperada de la franja
% Si se conoce de la simulacion, fijarla directamente:
fx_nominal = 1/8;     % ciclos/pixel

% para estimarla del primer frame, poner true
estimate_fx_from_first_frame = true;

% Perfil horizontal
use_hann_in_space = true;

% Regularizacion temporal por grafo (opcional pero recomendada)
use_graph_smoothing = true;
graphRadius = 3;
sigmaTime   = 1.5;
lambdaGraph = 4;

% PSD Welch
welchWindow  = 128;
welchOverlap = 0.5;
nfft         = 1024;

% Guardado
saveCSV = true;
saveMAT = true;

%% ---------------- CARGA DE FRAMES ----------------
files = dir(fullfile(frameFolder, framePattern));
assert(~isempty(files), 'No se encontraron frames en la carpeta indicada.');

[~, idxSort] = sort({files.name});
files = files(idxSort);

nFrames = numel(files);
t = (0:nFrames-1).' / fps;

I0 = imread(fullfile(frameFolder, files(1).name));
if size(I0,3) == 3
    I0 = rgb2gray(I0);
end
I0 = im2double(I0);
[ny, nx] = size(I0);

x = (0:nx-1).';

%% ---------------- ESTIMACION OPCIONAL DE fx ----------------
if estimate_fx_from_first_frame
    c0 = mean(I0,1).';
    c0 = c0 - mean(c0);

    if use_hann_in_space
        ws = hann(nx, 'periodic');
        c0w = c0 .* ws;
    else
        c0w = c0;
    end

    C0 = fft(c0w);
    fgrid = (0:nx-1)'/nx;   % ciclos/pixel

    % buscar pico positivo evitando DC
    searchIdx = 2:floor(nx/2);
    [~, imax] = max(abs(C0(searchIdx)));
    fx = fgrid(searchIdx(imax));
else
    fx = fx_nominal;
end

%% ---------------- EXTRACCION DE FASE POR FRAME ----------------
phi_wrapped = zeros(nFrames,1);
amp_carrier = zeros(nFrames,1);
dc_level    = zeros(nFrames,1);

for k = 1:nFrames
    Ik = imread(fullfile(frameFolder, files(k).name));
    if size(Ik,3) == 3
        Ik = rgb2gray(Ik);
    end
    Ik = im2double(Ik);

    % Perfil horizontal promedio por columna
    ck = mean(Ik,1).';
    dc_level(k) = mean(ck);

    % Remover DC
    ck0 = ck - mean(ck);

    % Ventana espacial para reducir leakage
    if use_hann_in_space
        ws = hann(nx, 'periodic');
        ck0 = ck0 .* ws;
    end

    % Demodulacion compleja en la frecuencia espacial de interes
    ref = exp(-1j*2*pi*fx*x);
    Ak = sum(ck0 .* ref);

    phi_wrapped(k) = angle(Ak);
    amp_carrier(k) = abs(Ak);
end

%% ---------------- DESENVOLVIMIENTO Y FLUCTUACION ----------------
phi_unwrapped = unwrap(phi_wrapped);

% Fluctuacion en radianes:
% quitar solo componente media + tendencia lineal lenta
phi_fluct_raw = detrend(phi_unwrapped, 1);

% Score de cambio entre frames consecutivos
phase_change_score = [0; abs(diff(phi_fluct_raw))];   % rad/frame

%% ---------------- SUAVIZADO TEMPORAL POR GRAFO ----------------
if use_graph_smoothing
    W = zeros(nFrames, nFrames);

    % similitud basada en tiempo y amplitud del carrier
    a = amp_carrier / max(amp_carrier + eps);
    for i = 1:nFrames
        j1 = max(1, i-graphRadius);
        j2 = min(nFrames, i+graphRadius);
        for j = j1:j2
            if i == j
                continue;
            end
            w_time = exp(-((i-j)^2)/(2*sigmaTime^2));
            w_amp  = exp(-((a(i)-a(j))^2)/(2*(0.15^2)));
            W(i,j) = w_time * w_amp;
        end
    end

    W = max(W, W.');
    D = diag(sum(W,2));
    L = D - W;

    phi_fluct = (eye(nFrames) + lambdaGraph*L) \ phi_fluct_raw;
else
    phi_fluct = phi_fluct_raw;
    W = [];
    L = [];
end

phase_change_score_smooth = [0; abs(diff(phi_fluct))];

%% ---------------- PSD TEMPORAL ----------------
xpsd = phi_fluct - mean(phi_fluct);

wlen = min(welchWindow, nFrames);
if mod(wlen,2) ~= 0
    wlen = wlen - 1;
end
noverlap = floor(welchOverlap * wlen);
win = hann(wlen, 'periodic');

[Pphi, f] = pwelch(xpsd, win, noverlap, nfft, fps, 'onesided');

% Ignorar DC al buscar pico dominante
valid = f >= 0.2;
[~, ipk] = max(Pphi(valid));
f_valid = f(valid);
f_peak = f_valid(ipk);

%% ---------------- TABLAS Y GUARDADO ----------------
resultsTable = table( ...
    (1:nFrames).', t, phi_wrapped, phi_unwrapped, phi_fluct_raw, phi_fluct, ...
    amp_carrier, dc_level, phase_change_score, phase_change_score_smooth, ...
    'VariableNames', {'frame','time_s','phi_wrapped_rad','phi_unwrapped_rad', ...
    'phi_fluct_raw_rad','phi_fluct_rad','carrier_amp','dc_level', ...
    'change_score_raw_rad','change_score_smooth_rad'});

psdTable = table(f, Pphi, 'VariableNames', {'freq_Hz','PSD_phase_rad2_per_Hz'});

if saveCSV
    writetable(resultsTable, fullfile(frameFolder, 'phase_time_series.csv'));
    writetable(psdTable, fullfile(frameFolder, 'phase_psd.csv'));
end

if saveMAT
    save(fullfile(frameFolder, 'phase_analysis_results.mat'), ...
        'resultsTable','psdTable','fx','f_peak','W','L');
end

%% ---------------- FIGURAS ----------------
figure('Color','w','Name','Analisis temporal de fase');

subplot(3,1,1);
plot(t, phi_fluct_raw, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0); hold on;
plot(t, phi_fluct, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Tiempo [s]');
ylabel('Fluctuacion [rad]');
title(sprintf('Fluctuaciones de fase vs tiempo (fx = %.4f ciclos/pixel)', fx));
legend('Estimación Directa','Suavizada','Location','northeast');

subplot(3,1,2);
plot(t, phase_change_score, 'Color', "#FFD966", 'LineWidth', 1.0); hold on;
plot(t, phase_change_score_smooth, 'Color', "#227D66", 'LineWidth', 1.5);
grid on;
xlabel('Tiempo [s]');
ylabel('Cambio [rad/frame]');
title('Score de cambio entre frames consecutivos');
legend('Estimación Directa','Suavizada','Location','northeast');

subplot(3,1,3);
plot(f, 10*log10(Pphi + eps), 'k', 'LineWidth', 1.5); hold on;
xline(f_peak, '--r', sprintf('Pico %.3f Hz', f_peak), 'LineWidth', 1.2);
grid on;
xlabel('Frecuencia [Hz]');
ylabel('PSD [dB/rad^2/Hz]');
title('PSD temporal de la fase');

fprintf('Frecuencia espacial usada fx = %.6f ciclos/pixel\n', fx);
fprintf('Pico dominante de la PSD = %.4f Hz\n', f_peak);
