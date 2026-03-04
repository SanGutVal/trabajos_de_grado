function [freqs_top, pxx_top, snr_db] = analizar_roi_perfil(region, plotear)
    if nargin < 2
        plotear = true;
    end

    mask_valid = region > 0;
    if ~any(mask_valid(:))
        warning('ROI vacía o todos ceros; no se puede analizar.');
        freqs_top = []; pxx_top = []; snr_db = NaN;
        return;
    end

    [h_roi, w_roi] = size(region);
    perfil = zeros(w_roi,1);
    for x = 1:w_roi
        col = region(:,x);
        mcol = mask_valid(:,x);
        if any(mcol)
            perfil(x) = mean(col(mcol));
        else
            perfil(x) = 0;
        end
    end

    perfil = detrend(perfil);
    perfil = imgaussfilt(perfil',0.8)';    % suavizado ligero
    perfil = perfil .* hann(length(perfil));

    if plotear
        figure; plot(perfil);
        title('Perfil horizontal promedio');
        xlabel('Columna de píxel'); ylabel('Intensidad');
    end

    % PSD con pwelch (f en ciclos/píxel)
    N = length(perfil);
    win = hanning(N);
    noverlap = floor(N/4);
    nfft = N;
    [pxx,f] = pwelch(perfil, win, noverlap, nfft, 1);  % fs=1[web:45][web:51]

    if plotear
        figure; plot(f,10*log10(pxx)); grid on;
        title('PSD de la ROI (pwelch)');
        xlabel('Frecuencia (ciclos/píxel)'); ylabel('Potencia (dB)');
    end

    % Definir intervalos de frecuencia (ciclos/píxel)
    intervals = [0    0.1;
                 0.1  0.25;
                 0.25 0.4];

    freqs_top = nan(1,3);
    pxx_top   = nan(1,3);

    for i = 1:3
        f_min = intervals(i,1);
        f_max = intervals(i,2);

        idx = (f >= f_min) & (f < f_max);
        idx(1) = false;  % excluir DC explícitamente

        if ~any(idx)
            % No hay muestras en ese intervalo
            freqs_top(i) = NaN;
            pxx_top(i)   = NaN;
            continue;
        end

        % Buscar máximo local dentro del intervalo
        [pk, loc_rel] = max(pxx(idx));
        loc_abs = find(idx);
        loc_abs = loc_abs(loc_rel);

        freqs_top(i) = f(loc_abs);
        pxx_top(i)   = pk;
    end

    % SNR aproximado del intervalo principal
    % Se puede ajustar a otro intervalo si tiene más sentido físico.
    main_band = 2;   % Segundo intervalo como "principal"
    if ~isnan(pxx_top(main_band))
        idx_band = (f >= intervals(main_band,1)) & (f < intervals(main_band,2));
        idx_band(1) = false;
        pxx_band = pxx(idx_band);
        if ~isempty(pxx_band)
            ruido_med = mean(pxx_band);
            snr_lin = pxx_top(main_band) / max(ruido_med,eps);
            snr_db = 10*log10(snr_lin);
        else
            snr_db = NaN;
        end
    else
        snr_db = NaN;
    end

    if plotear
        fprintf('Máximos por intervalo:\n');
        for i = 1:3
            fprintf('  [%0.2f, %0.2f): f = %.4f c/px, P = %.3g\n', ...
                intervals(i,1), intervals(i,2), freqs_top(i), pxx_top(i));
        end
        fprintf('SNR (intervalo 0.1–0.25): %.1f dB\n', snr_db);
    end
end