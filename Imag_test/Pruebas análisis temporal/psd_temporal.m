function [f_t, Pxx, meta] = psd_temporal(x, fs, varargin)
% PSD_TEMPORAL  Estima PSD temporal con Welch de forma robusta
%
% Uso:
%   [f_t, Pxx, meta] = psd_temporal(x, fs)
%   [f_t, Pxx, meta] = psd_temporal(x, fs, 'detrendOrder', 1, 'segmentLength', 128)
%
% Entradas:
%   x  : serie temporal
%   fs : frecuencia de muestreo [Hz]
%
% Salidas:
%   f_t  : frecuencias [Hz]
%   Pxx  : PSD
%   meta : estructura con parametros usados

    p = inputParser;
    addParameter(p, 'detrendOrder', 1);
    addParameter(p, 'segmentLength', []);
    addParameter(p, 'overlapFrac', 0.5);
    addParameter(p, 'nfft', []);
    parse(p, varargin{:});

    detrendOrder = p.Results.detrendOrder;
    segmentLength = p.Results.segmentLength;
    overlapFrac = p.Results.overlapFrac;
    nfft = p.Results.nfft;

    x = x(:);
    x = x(~isnan(x));

    N = numel(x);
    if N < 16
        warning('Muy pocas muestras (%d) para PSD robusta.', N);
    end

    if detrendOrder == 0
        x = x - mean(x);
    elseif detrendOrder == 1
        x = detrend(x, 1);
    else
        error('detrendOrder debe ser 0 o 1.');
    end

    if isempty(segmentLength)
        segmentLength = min(128, floor(N/4));
    end
    segmentLength = max(segmentLength, 16);
    segmentLength = min(segmentLength, N);

    if mod(segmentLength, 2) ~= 0
        segmentLength = segmentLength - 1;
    end

    if isempty(nfft)
        nfft = max(256, 2^nextpow2(segmentLength));
    end

    win = hann(segmentLength, 'periodic');
    noverlap = floor(overlapFrac * segmentLength);

    [Pxx, f_t] = pwelch(x, win, noverlap, nfft, fs);

    meta = struct();
    meta.N = N;
    meta.segmentLength = segmentLength;
    meta.noverlap = noverlap;
    meta.nfft = nfft;
    meta.fs = fs;
    meta.detrendOrder = detrendOrder;
end
