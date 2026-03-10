% Cargar resultados de medias/desviaciones por columna
nombre_mat = 'estadisticas_temporales_101_a_115.mat';  % ajustar
load(nombre_mat, 'mu', 'sd', 'frames_vec', 'fps');

[nF, Wroi] = size(mu);

%% Serie temporal de la media global de la ROI
m_global = mean(mu, 2);    % media sobre columnas por frame

% Grafica tendencia temporal
t = (0:nF-1) / fps;        % tiempo en segundos

figure;
plot(t, m_global, '-o');
xlabel('Tiempo (s)');
ylabel('Intensidad media ROI');
title('Evolución temporal de la intensidad media en la ROI');
grid on;

% PSD temporal de la media global
[f_t, Pxx] = psd_temporal(m_global, fps);

figure;
plot(f_t, 10*log10(Pxx));
xlabel('Frecuencia temporal (Hz)');
ylabel('PSD media ROI (dB)');
title('PSD temporal de la intensidad media en la ROI');
grid on;

%% PSD temporal promedio sobre columnas
% PSD de cada columna y promedio
nfft = nF;
P_cols = zeros(nfft, Wroi);
for j = 1:Wroi
    xj = mu(:,j);
    [Pj, f_t2] = psd_temporal(xj, fps);
    % Asegurar mismo tamaño (por si hay NaN o similares)
    P_cols(:,j) = Pj(:);
end

P_mean_cols = mean(P_cols, 2);

figure;
plot(f_t2, 10*log10(P_mean_cols));
xlabel('Frecuencia temporal (Hz)');
ylabel('PSD promedio columnas (dB)');
title('PSD temporal promedio sobre columnas (medias por columna)');
grid on;

%% Medida de cambio frame a frame (distancia entre perfiles)
% Perfil de medias por frame: mu(k,:) ∈ R^{Wroi}
D = zeros(nF-1,1);
for k = 2:nF
    diff_k = mu(k,:) - mu(k-1,:);
    D(k-1) = norm(diff_k, 2);   % norma L2 entre perfiles consecutivos
end

t_D = (1:nF-1) / fps;   % tiempo asociado a diferencias

figure;
plot(t_D, D, '-o');
xlabel('Tiempo (s)');
ylabel('||mu_k - mu_{k-1}||_2');
title('Cambio entre frames consecutivos (norma L2 de medias)');
grid on;

% PSD temporal de la secuencia de cambios D
[f_D, P_D] = psd_temporal(D, fps);

figure;
plot(f_D, 10*log10(P_D));
xlabel('Frecuencia temporal (Hz)');
ylabel('PSD de cambios entre frames (dB)');
title('PSD temporal de la secuencia de cambios entre frames');
grid on;