
% - Score principal: divergencia entre histogramas de medias por columna
% - Score opcional fusionado: histogramas + similitud del perfil espacial
% - Regularización profesional: suavizado sobre grafo temporal
% - PSD final: Welch sobre score temporal

clear; clc;

%% ===================== CONFIGURACION =====================
frameFolder   = 'franjas_sim_temporal';   %
framePattern  = 'frame_*.png';
fps           = 30;

% Histograma del perfil horizontal (medias por columna)
nBins         = 64;
histRange     = [0 1];

% Score entre frames consecutivos
alpha_jsd     = 0.60;   % peso Jensen-Shannon
beta_w1       = 0.25;   % peso Wasserstein-1
gamma_prof    = 0.15;   % peso termino de perfil (1-rho)
useProfileTerm = true;  % false = score puramente por histogramas

% Grafo temporal
graphRadius   = 3;      % conecta cada frame con vecinos +/- graphRadius
sigmaTime     = 1.5;    % controla decaimiento por distancia temporal
lambdaGraph   = 8;      % regularizacion de Laplaciano

% PSD (Welch)
welchWindow   = 64;     % si N es menor, se ajusta automaticamente
welchOverlap  = 0.5;    % fraccion de traslape
fMinPeak      = 0.2;    % evitar DC y muy baja frecuencia al buscar pico

% Guardado
saveCSV       = true;
saveMAT       = true;

%% ===================== CARGA DE FRAMES =====================
files = dir(fullfile(frameFolder, framePattern));
assert(~isempty(files), 'No se encontraron frames en la carpeta indicada.');

% Orden natural por nombre
[~, idxSort] = sort({files.name});
files = files(idxSort);

nFrames = numel(files);

% Leer primera imagen para dimensiones
I0 = imread(fullfile(frameFolder, files(1).name));
if size(I0,3) == 3
    I0 = rgb2gray(I0);
end
I0 = im2double(I0);
[ny, nx] = size(I0);

t = (0:nFrames-1).' / fps;
edges = linspace(histRange(1), histRange(2), nBins+1);
binCenters = 0.5 * (edges(1:end-1) + edges(2:end));
binWidth = edges(2) - edges(1);

% Prealocacion
colProfiles = zeros(nx, nFrames);     % perfil horizontal por frame
H = zeros(nBins, nFrames);            % histograma del perfil por frame

%% ===================== EXTRACCION DE FEATURES =====================
for k = 1:nFrames
    Ik = imread(fullfile(frameFolder, files(k).name));
    if size(Ik,3) == 3
        Ik = rgb2gray(Ik);
    end
    Ik = im2double(Ik);

    % Perfil horizontal: media por columna
    ck = mean(Ik, 1).';               % nx x 1
    ck = min(max(ck, histRange(1)), histRange(2));

    % Histograma probabilistico del perfil
    hk = histcounts(ck, edges, 'Normalization', 'probability').';

    colProfiles(:, k) = ck;
    H(:, k) = hk;
end

%% ===================== SCORE CONSECUTIVO ENTRE FRAMES =====================
score_jsd  = zeros(nFrames,1);
score_w1   = zeros(nFrames,1);
score_prof = zeros(nFrames,1);
score_hist = zeros(nFrames,1);   % solo histogramas
score_fused = zeros(nFrames,1);  % histogramas + perfil

for k = 2:nFrames
    p = H(:,k-1);
    q = H(:,k);

    d_jsd = jensen_shannon_div(p, q);           % ~ [0,1]
    d_w1  = wasserstein1_1d(p, q, binWidth);    % en soporte [0,1], ~ [0,1]
    rho   = safe_corr(colProfiles(:,k-1), colProfiles(:,k));
    d_prof = 1 - rho;                           % 0 = identico, 2 = opuesto

    score_jsd(k)  = d_jsd;
    score_w1(k)   = d_w1;
    score_prof(k) = d_prof;

    score_hist(k) = alpha_jsd*d_jsd + beta_w1*d_w1;

    if useProfileTerm
        score_fused(k) = alpha_jsd*d_jsd + beta_w1*d_w1 + gamma_prof*d_prof;
    else
        score_fused(k) = score_hist(k);
    end
end

%% ===================== GRAFO TEMPORAL ENTRE FRAMES =====================
% Nodos: frames
% Aristas: vecinos temporales cercanos
% Peso: similitud basada en distancia de histogramas + opcional perfil

W = zeros(nFrames, nFrames);

% Escala robusta para similitud
baseDist = score_fused(2:end);
sigmaFeat = median(baseDist(baseDist > 0));
if isempty(sigmaFeat) || sigmaFeat == 0
    sigmaFeat = 1e-6;
end

for i = 1:nFrames
    j1 = max(1, i - graphRadius);
    j2 = min(nFrames, i + graphRadius);

    for j = j1:j2
        if i == j
            continue;
        end

        d_jsd_ij = jensen_shannon_div(H(:,i), H(:,j));
        d_w1_ij  = wasserstein1_1d(H(:,i), H(:,j), binWidth);
        rho_ij   = safe_corr(colProfiles(:,i), colProfiles(:,j));
        d_prof_ij = 1 - rho_ij;

        if useProfileTerm
            d_ij = alpha_jsd*d_jsd_ij + beta_w1*d_w1_ij + gamma_prof*d_prof_ij;
        else
            d_ij = alpha_jsd*d_jsd_ij + beta_w1*d_w1_ij;
        end

        w_feat = exp(-(d_ij^2) / (2*sigmaFeat^2));
        w_time = exp(-((i-j)^2) / (2*sigmaTime^2));

        W(i,j) = w_feat * w_time;
    end
end

W = max(W, W.');              % asegurar simetria
D = diag(sum(W,2));
L = D - W;                    % Laplaciano del grafo

%% ===================== SUAVIZADO GRAFICO DEL SCORE =====================
Ireg = eye(nFrames);

score_hist_graph  = (Ireg + lambdaGraph * L) \ score_hist;
score_fused_graph = (Ireg + lambdaGraph * L) \ score_fused;

%% ===================== PSD TEMPORAL =====================
% Omitimos el primer punto porque score(1)=0 por construccion
x_hist_raw   = detrend(score_hist(2:end), 0);
x_hist_graph = detrend(score_hist_graph(2:end), 0);

x_fused_raw   = detrend(score_fused(2:end), 0);
x_fused_graph = detrend(score_fused_graph(2:end), 0);

Npsd = numel(x_fused_graph);
wlen = min(welchWindow, Npsd);
if mod(wlen,2) ~= 0
    wlen = wlen - 1;
end
if wlen < 8
    error('Muy pocos frames para estimar PSD con Welch de forma razonable.');
end

win = hann(wlen, 'periodic');
noverlap = floor(welchOverlap * wlen);
nfft = max(256, 2^nextpow2(Npsd));

[P_hist_raw, f]   = pwelch(x_hist_raw,   win, noverlap, nfft, fps, 'onesided');
[P_hist_graph, ~] = pwelch(x_hist_graph, win, noverlap, nfft, fps, 'onesided');

[P_fused_raw, ~]   = pwelch(x_fused_raw,   win, noverlap, nfft, fps, 'onesided');
[P_fused_graph, ~] = pwelch(x_fused_graph, win, noverlap, nfft, fps, 'onesided');

% Pico dominante (evitar DC)
valid = f >= fMinPeak;
[~, i1] = max(P_hist_graph(valid));
[~, i2] = max(P_fused_graph(valid));

f_valid = f(valid);
peak_hist_graph  = f_valid(i1);
peak_fused_graph = f_valid(i2);

%% ===================== RESULTADOS Y GUARDADO =====================
resultsTable = table( ...
    (1:nFrames).', t, score_jsd, score_w1, score_prof, ...
    score_hist, score_hist_graph, score_fused, score_fused_graph, ...
    'VariableNames', {'frame','time_s','jsd','w1','profile_change', ...
                      'score_hist','score_hist_graph', ...
                      'score_fused','score_fused_graph'});

psdTable = table(f, P_hist_raw, P_hist_graph, P_fused_raw, P_fused_graph, ...
    'VariableNames', {'freq_Hz','PSD_hist_raw','PSD_hist_graph', ...
                      'PSD_fused_raw','PSD_fused_graph'});

if saveCSV
    writetable(resultsTable, fullfile(frameFolder, 'temporal_scores.csv'));
    writetable(psdTable,     fullfile(frameFolder, 'temporal_psd.csv'));
end

if saveMAT
    save(fullfile(frameFolder, 'temporal_psd_results.mat'), ...
        'H','colProfiles','resultsTable','psdTable','W','L', ...
        'peak_hist_graph','peak_fused_graph','fps','binCenters');
end

%% ===================== FIGURAS =====================
figure('Color','w','Name','Analisis temporal por histogramas');

subplot(3,1,1);
imagesc(t, binCenters, H);
axis xy;
xlabel('Tiempo [s]');
ylabel('Bin de intensidad');
title('Evolucion temporal de histogramas del perfil horizontal');
colorbar;

subplot(3,1,2);
plot(t, score_hist, '-', 'LineWidth', 1.0); hold on;
plot(t, score_hist_graph, '-', 'LineWidth', 1.5);
plot(t, score_fused, '--', 'LineWidth', 1.0);
plot(t, score_fused_graph, '-', 'LineWidth', 1.8);
grid on;
xlabel('Tiempo [s]');
ylabel('Score de cambio');
title('Cambio entre frames consecutivos');
legend('Hist raw','Hist graph','Fused raw','Fused graph', 'Location','best');

subplot(3,1,3);
plot(f, 10*log10(P_hist_graph + eps), 'LineWidth', 1.4); hold on;
plot(f, 10*log10(P_fused_graph + eps), 'LineWidth', 1.6);
xline(peak_hist_graph, '--');
xline(peak_fused_graph, '--');
grid on;
xlabel('Frecuencia [Hz]');
ylabel('PSD [dB/Hz]');
title('PSD temporal del score');
legend(sprintf('Hist graph (pico %.3f Hz)', peak_hist_graph), ...
       sprintf('Fused graph (pico %.3f Hz)', peak_fused_graph), ...
       'Location','best');

fprintf('Pico dominante PSD (hist+grafo):  %.4f Hz\n', peak_hist_graph);
fprintf('Pico dominante PSD (fused+grafo): %.4f Hz\n', peak_fused_graph);

%% ===================== FUNCIONES LOCALES =====================
function d = jensen_shannon_div(p, q)
    p = p(:); q = q(:);
    p = p / max(sum(p), eps);
    q = q / max(sum(q), eps);
    m = 0.5 * (p + q);

    d = 0.5 * kl_div(p, m) + 0.5 * kl_div(q, m);

    % Por estabilidad numérica
    d = max(d, 0);
    d = min(d, 1);
end

function d = kl_div(p, q)
    idx = (p > 0) & (q > 0);
    d = sum(p(idx) .* log2(p(idx) ./ q(idx)));
end

function d = wasserstein1_1d(p, q, dx)
    p = p(:); q = q(:);
    p = p / max(sum(p), eps);
    q = q / max(sum(q), eps);

    Fp = cumsum(p);
    Fq = cumsum(q);

    d = sum(abs(Fp - Fq)) * dx;
end

function rho = safe_corr(x, y)
    x = x(:); y = y(:);

    sx = std(x);
    sy = std(y);

    if sx < eps || sy < eps
        rho = 1;
        return;
    end

    C = corrcoef(x, y);
    rho = C(1,2);

    if isnan(rho)
        rho = 0;
    end
end