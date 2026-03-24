clear; clc; close all;

%% ===================== CARGA =====================
nombre_mat = 'estadisticas_roi_100_a_549.mat';
load(nombre_mat, 'results');

% Atajos
fps   = results.info.fps;
t     = results.time.time_s(:);
mu    = results.signals.mu;
sd    = results.signals.sd;

m_global        = results.signals.m_global(:);
sd_global       = results.signals.sd_global(:);
phi_wrapped     = results.signals.phi_wrapped(:);
phi_unwrapped   = results.signals.phi_unwrapped(:);
phi_fluct       = results.signals.phi_fluct(:);
carrier_amp     = results.signals.carrier_amp(:);
change_profile  = results.signals.change_profile_L2(:);
change_phase    = results.signals.change_phase_abs(:);

frames_valid = results.time.frames_valid(:);
fx_used      = results.demod.fx_used;

[nF, Wroi] = size(mu);

fprintf('Archivo cargado: %s\n', nombre_mat);
fprintf('Frames válidos: %d\n', nF);
fprintf('ROI: %d x %d px\n', results.roi.ancho, results.roi.alto);
fprintf('fx usada: %.6f ciclos/pixel\n', fx_used);

%% ===================== PARAMETROS PSD =====================
fmin_peak = 0.2;
fmax_peak = fps/2;

segmentLength = min(128, max(32, floor(nF/4)));
if mod(segmentLength,2) ~= 0
    segmentLength = segmentLength - 1;
end

%% ===================== PSD DE SENALES PRINCIPALES =====================
[f_phi, P_phi, meta_phi] = psd_temporal(phi_fluct, fps, ...
    'detrendOrder', 0, ...
    'segmentLength', segmentLength, ...
    'overlapFrac', 0.5);

[f_mg, P_mg, meta_mg] = psd_temporal(m_global, fps, ...
    'detrendOrder', 1, ...
    'segmentLength', segmentLength, ...
    'overlapFrac', 0.5);

[f_sg, P_sg, meta_sg] = psd_temporal(sd_global, fps, ...
    'detrendOrder', 1, ...
    'segmentLength', segmentLength, ...
    'overlapFrac', 0.5);

[f_cp, P_cp, meta_cp] = psd_temporal(change_profile, fps, ...
    'detrendOrder', 1, ...
    'segmentLength', segmentLength, ...
    'overlapFrac', 0.5);

[f_ch, P_ch, meta_ch] = psd_temporal(change_phase, fps, ...
    'detrendOrder', 1, ...
    'segmentLength', segmentLength, ...
    'overlapFrac', 0.5);

%% ===================== PSD PROMEDIO SOBRE COLUMNAS =====================
P_cols = [];
f_cols = [];

for j = 1:Wroi
    xj = mu(:,j);
    [fj, Pj] = psd_temporal(xj, fps, ...
        'detrendOrder', 1, ...
        'segmentLength', segmentLength, ...
        'overlapFrac', 0.5);

    if isempty(P_cols)
        P_cols = zeros(numel(Pj), Wroi);
        f_cols = fj(:);
    end

    P_cols(:,j) = Pj(:);
end

P_mean_cols = mean(P_cols, 2);

%% ===================== PICOS DOMINANTES =====================
[fpk_phi,  ppk_phi]  = dominant_peak(f_phi,  P_phi,  fmin_peak, fmax_peak);
[fpk_mg,   ppk_mg]   = dominant_peak(f_mg,   P_mg,   fmin_peak, fmax_peak);
[fpk_sg,   ppk_sg]   = dominant_peak(f_sg,   P_sg,   fmin_peak, fmax_peak);
[fpk_cp,   ppk_cp]   = dominant_peak(f_cp,   P_cp,   fmin_peak, fmax_peak);
[fpk_ch,   ppk_ch]   = dominant_peak(f_ch,   P_ch,   fmin_peak, fmax_peak);
[fpk_cols, ppk_cols] = dominant_peak(f_cols, P_mean_cols, fmin_peak, fmax_peak);

%% ===================== GRAFICAS TEMPORALES =====================
figure('Color','w','Name','Series temporales ROI');

subplot(3,1,1);
plot(t, phi_fluct, 'b', 'LineWidth', 1.4);
grid on;
xlabel('Tiempo [s]');
ylabel('Fluctuacion [rad]');
title(sprintf('Fluctuaciones de fase en la ROI (fx = %.4f ciclos/pixel)', fx_used));

subplot(3,1,2);
plot(t, m_global, 'k', 'LineWidth', 1.2); hold on;
plot(t, sd_global, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
grid on;
xlabel('Tiempo [s]');
ylabel('Intensidad');
title('Descriptores globales de intensidad en la ROI');
legend('Media global ROI', 'Desv. est. global ROI', 'Location', 'best');

subplot(3,1,3);
plot(t, carrier_amp, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.3);
grid on;
xlabel('Tiempo [s]');
ylabel('|A_k|');
title('Amplitud de la portadora en la ROI');

%% ===================== GRAFICAS DE CAMBIO =====================
figure('Color','w','Name','Cambios entre frames');

subplot(2,1,1);
plot(t, change_phase, 'Color', '#1791BA', 'LineWidth', 1.3);
grid on;
xlabel('Tiempo [s]');
ylabel('|Δfase| [rad]');
title('Cambio entre frames consecutivos basado en fase');

subplot(2,1,2);
plot(t, change_profile, 'Color', [0.85 0.5 0.1], 'LineWidth', 1.3);
grid on;
xlabel('Tiempo [s]');
ylabel('||\mu_k-\mu_{k-1}||_2');
title('Cambio entre frames consecutivos basado en perfiles');

%% ===================== PSD COMPARATIVAS =====================
figure('Color','w','Name','PSD comparativas ROI');

subplot(2,2,1);
plot(f_phi, 10*log10(P_phi + eps), 'b', 'LineWidth', 1.5); hold on;
xline(fpk_phi, '--r', sprintf('%.3f Hz', fpk_phi), 'LineWidth', 1.2);
grid on;
xlabel('Frecuencia [Hz]');
ylabel('PSD [dB/rad^2/Hz]');
title('PSD de la fase fluctuante');

subplot(2,2,2);
plot(f_mg, 10*log10(P_mg + eps), 'k', 'LineWidth', 1.3); hold on;
plot(f_sg, 10*log10(P_sg + eps), 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3);
xline(fpk_mg, '--k', sprintf('Media %.3f Hz', fpk_mg), 'LineWidth', 1.0);
xline(fpk_sg, '--', sprintf('SD %.3f Hz', fpk_sg), 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
grid on;
xlabel('Frecuencia [Hz]');
ylabel('PSD [dB]');
title('PSD de descriptores globales');
legend('Media global ROI', 'SD global ROI', 'Location', 'best');

subplot(2,2,3);
plot(f_cols, 10*log10(P_mean_cols + eps), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.4); hold on;
xline(fpk_cols, '--r', sprintf('%.3f Hz', fpk_cols), 'LineWidth', 1.2);
grid on;
xlabel('Frecuencia [Hz]');
ylabel('PSD [dB]');
title('PSD promedio de medias por columna');

subplot(2,2,4);
plot(f_cp, 10*log10(P_cp + eps), 'Color', "#FFD966", 'LineWidth', 1.3); hold on;
plot(f_ch, 10*log10(P_ch + eps), 'Color', "#227D66", 'LineWidth', 1.3);
xline(fpk_cp, '--', sprintf('Perfil %.3f Hz', fpk_cp), 'Color', "#FFD966", 'LineWidth', 1.0);
xline(fpk_ch, '--', sprintf('Fase %.3f Hz', fpk_ch), 'Color', "#227D66", 'LineWidth', 1.0);
grid on;
xlabel('Frecuencia [Hz]');
ylabel('PSD [dB]');
title('PSD de cambios entre frames');
legend('Cambio perfil', 'Cambio fase', 'Location', 'best');

%% ===================== RESUMEN NUMERICO =====================
summaryTable = table( ...
    ["phi_fluct"; "m_global"; "sd_global"; "mu_prom_columnas"; "change_profile_L2"; "change_phase_abs"], ...
    [fpk_phi; fpk_mg; fpk_sg; fpk_cols; fpk_cp; fpk_ch], ...
    [10*log10(ppk_phi+eps); 10*log10(ppk_mg+eps); 10*log10(ppk_sg+eps); ...
     10*log10(ppk_cols+eps); 10*log10(ppk_cp+eps); 10*log10(ppk_ch+eps)], ...
    'VariableNames', {'signal','peak_freq_Hz','peak_level_dB'});

disp(summaryTable);

save('resumen_psd_roi.mat', 'summaryTable', ...
     'f_phi','P_phi','f_mg','P_mg','f_sg','P_sg','f_cols','P_mean_cols', ...
     'f_cp','P_cp','f_ch','P_ch', 'meta_phi','meta_mg','meta_sg','meta_cp','meta_ch');

%% ===================== FUNCIONES LOCALES =====================
function [fpk, ppk] = dominant_peak(f, P, fmin, fmax)
    idx = (f >= fmin) & (f <= fmax);
    ff = f(idx);
    PP = P(idx);

    if isempty(ff)
        fpk = NaN;
        ppk = NaN;
        return;
    end

    [ppk, imax] = max(PP);
    fpk = ff(imax);
end
