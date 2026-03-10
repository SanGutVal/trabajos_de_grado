function [f_t, Pxx] = psd_temporal(x, fs)
% x: vector columna de longitud N (serie temporal)
% fs: frecuencia de muestreo temporal (fps)
% Salida:
%   f_t: vector de frecuencias (Hz)
%   Pxx: densidad espectral de potencia

    x = x(:);              % asegurar columna
    x = detrend(x);        % quitar tendencia lineal

    N = length(x);
    if N < 16
        warning('Muy pocas muestras (%d) para PSD robusta.', N);
    end

    % Ventana y parámetros de Welch (ajustables)
    win      = hanning(N);
    noverlap = floor(N/4);
    nfft     = N;

    [Pxx, f_t] = pwelch(x, win, noverlap, nfft, fs);
end