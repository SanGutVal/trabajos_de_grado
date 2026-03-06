% ------------ CONFIGURACIÓN GENERAL ------------
archivo = 'frame_101.jpg';          % Nombre de la imagen
usar_roi_manual = false;            % true = Pedir ROI manual
ancho_roi_fix = 216;
alto_roi_fix  = 146;

% ------------ CARGA DE IMAGEN ------------
img = imread(archivo);
if size(img,3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end
img_d = double(img_gray);
[h,w] = size(img_gray);

figure; imshow(img_gray,[]); title('Imagen completa');
drawnow;

% ------------ ROI FIJA CENTRADA ------------
cx = round(w/2);
cy = round(h/2);
half_w = round(ancho_roi_fix/2);
half_h = round(alto_roi_fix/2);

xmin = max(cx-half_w+1,1);
xmax = min(cx+half_w,  w);
ymin = max(cy-half_h+1,1);
ymax = min(cy+half_h,  h);

roi_fix = img_d(ymin:ymax, xmin:xmax);

figure; imshow(roi_fix,[]); 
title(sprintf('ROI fija centrada (%dx%d)', size(roi_fix,2), size(roi_fix,1)));

% Analizar ROI fija
fprintf('----- ROI FIJA CENTRADA -----\n');
analizar_roi_perfil(roi_fix);

% ------------ ROI MANUAL (OPCIONAL) ------------
if usar_roi_manual
    figure; imshow(img_gray,[]); title('Dibuja ROI manual (polígono), doble clic para cerrar');
    mask = roipoly;                          % máscara binaria ROI manual [web:39]
    close;
    roi_manual = img_d .* mask;
    % recortar al bounding box para no procesar toda la imagen
    [r,c] = find(mask);
    ymin2 = min(r); ymax2 = max(r);
    xmin2 = min(c); xmax2 = max(c);
    roi_manual = roi_manual(ymin2:ymax2, xmin2:xmax2);

    figure; imshow(roi_manual,[]); title('ROI manual recortada');

    fprintf('----- ROI MANUAL -----\n');
    analizar_roi_perfil(roi_manual);
end

% ------------ FUNCIÓN DE ANÁLISIS ------------
function analizar_roi_perfil(region)
    % region: matriz double de la ROI

    % Máscara válida: ignorar ceros puros por si vienen de fuera de la ROI
    mask_valid = region > 0;
    if ~any(mask_valid(:))
        warning('ROI vacía o todos ceros; no se puede analizar.');
        return;
    end

    % Perfil horizontal promedio solo sobre píxeles válidos
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

    % Preprocesado
    perfil = detrend(perfil);              % quita DC/rampa
    perfil = imgaussfilt(perfil',0.8)';    % suavizado ligero
    perfil = perfil .* hann(length(perfil));

    figure; plot(perfil);
    title('Perfil horizontal promedio');
    xlabel('Columna de píxel'); ylabel('Intensidad');

    % PSD con pwelch (una sola ventana = equivalente a FFT con ventana) [web:45][web:51]
    N = length(perfil);
    win = hanning(N);
    noverlap = floor(N/4);
    nfft = N;
    [pxx,f] = pwelch(perfil, win, noverlap, nfft, 1);  % fs = 1 => f en ciclos/píxel

    figure; plot(f,10*log10(pxx)); grid on;
    title('PSD de la ROI (pwelch)');
    xlabel('Frecuencia (ciclos/píxel)'); ylabel('Potencia (dB)');

    % Picos (ignorando DC)
    pxx_no_dc = pxx(2:end);
    f_no_dc   = f(2:end);

    [pk, loc] = findpeaks(pxx_no_dc, 'MinPeakProminence',0.05*max(pxx_no_dc), ...
                                      'NPeaks',3, 'SortStr','descend');
    if isempty(loc)
        warning('No se detectaron picos significativos en la PSD.');
        return;
    end

    freqs = f_no_dc(loc);
    periodos = 1./freqs;

    % SNR aproximado usando primer pico
    ruido_med = mean(pxx_no_dc(1:max(loc(1)-1,1)));
    snr_lin = pk(1) / max(ruido_med,eps);
    snr_db = 10*log10(snr_lin);
    if snr_lin > 3
        conf_str = 'alta';
    else
        conf_str = 'baja';
    end

    fprintf('Picos principales (ignorando DC):\n');
    for i = 1:length(freqs)
        fprintf('  Pico %d: %.4f c/px (periodo %.2f px), potencia %.3g\n', ...
            i, freqs(i), periodos(i), pk(i));
    end
    fprintf('SNR estimado del primer pico: %.1f dB (confianza %s)\n', snr_db, conf_str);
end
