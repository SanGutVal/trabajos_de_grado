clc
clear
clearvars
%% -------- CONFIGURACIÓN --------
carpeta = 'C:\Users\sangu\Documents\Trabajo de Grado\trabajos_de_grado\Imag_test';
prefijo = 'frame_';            % e.g. 'frame_22'
extension = '.jpg';            % 

frame_ini = 101;                % primer índice
n_frames  = 10;                % cantidad de frames a procesar

usar_roi_manual = false;       % flase = ROI fija en el centro de la imagen
ancho_roi_fix = 216;
alto_roi_fix  = 146;

nombre_mat = sprintf('resultados_frames_%d_a_%d.mat', ...
                     frame_ini, frame_ini+n_frames-1);

%% -------- PRE-ASIGNAR RESULTADOS --------
resultado_base = struct( ...
    'frame',            [], ...
    'archivo',          '', ...
    'freqs_fix',        [], ...
    'potencias_fix',    [], ...
    'snr_db_fix',       NaN, ...
    'freqs_manual',     [], ...
    'potencias_manual', [], ...
    'snr_db_manual',    NaN);

resultados = repmat(resultado_base, n_frames, 1);
%% -------- BUCLE PRINCIPAL --------
for k = 1:n_frames
    idx = frame_ini + k - 1;
    nombre = sprintf('%s%d%s', prefijo, idx, extension);
    ruta   = fullfile(carpeta, nombre);

    if ~isfile(ruta)
        warning('No se encontró el archivo: %s. Se salta.', ruta);
        continue;
    end

    fprintf('\n===== Procesando %s =====\n', nombre);
    img = imread(ruta);
    if size(img,3) == 3
        img_gray = rgb2gray(img);
    else
        img_gray = img;
    end
    img_d = double(img_gray);
    [h,w] = size(img_gray);

    % --- ROI fija centrada ---
    cx = round(w/2); cy = round(h/2);
    half_w = round(ancho_roi_fix/2);
    half_h = round(alto_roi_fix/2);

    xmin = max(cx-half_w+1,1);
    xmax = min(cx+half_w,  w);
    ymin = max(cy-half_h+1,1);
    ymax = min(cy+half_h,  h);

    roi_fix = img_d(ymin:ymax, xmin:xmax);

    [freqs_fix, pows_fix, snr_fix] = analizar_roi_perfil(roi_fix, false);

    % --- ROI manual opcional (misma máscara para todos) ---
    freqs_manual   = [];
    pows_manual    = [];
    snr_manual     = NaN;

    if usar_roi_manual
        if k == 1
            figure; imshow(img_gray,[]);
            title('Dibujar ROI para todos los frames');
            mask_global = roipoly; close;
            [r,c] = find(mask_global);
            ymin2 = min(r); ymax2 = max(r);
            xmin2 = min(c); xmax2 = max(c);
        end
        roi_manual = img_d(ymin2:ymax2, xmin2:xmax2) .* ...
                     mask_global(ymin2:ymax2, xmin2:xmax2);

        [freqs_manual, pows_manual, snr_manual] = ...
            analizar_roi_perfil(roi_manual, false);
    end

    % Guardar en struct
    resultados(k).frame            = idx;
    resultados(k).archivo          = nombre;
    resultados(k).freqs_fix        = freqs_fix(:).';     % fila
    resultados(k).potencias_fix    = pows_fix(:).';
    resultados(k).snr_db_fix       = snr_fix;
    resultados(k).freqs_manual     = freqs_manual(:).';
    resultados(k).potencias_manual = pows_manual(:).';
    resultados(k).snr_db_manual    = snr_manual;
end

%% -------- FIGURAS RESUMEN (PICO PRINCIPAL ROI FIJA) --------
% Extraer frames válidos (por si faltó algún archivo)
frames = [resultados.frame];
mask_valid = ~arrayfun(@(r) isempty(r.freqs_fix), resultados);

frames_valid = frames(mask_valid);
freq1 = nan(size(frames_valid));
pow1  = nan(size(frames_valid));

j = 1;
for k = 1:numel(resultados)
    if isempty(resultados(k).freqs_fix)
        continue;
    end
    freq1(j) = resultados(k).freqs_fix(1);      % primer pico (mayor potencia)
    pow1(j)  = resultados(k).potencias_fix(1);  % potencia asociada
    j = j + 1;
end

% Frecuencia principal vs frame
figure;
plot(frames_valid, freq1, '-o');
xlabel('Número de frame');
ylabel('Frecuencia principal (ciclos/píxel)');
title('Evolución temporal de la frecuencia principal (ROI fija)');
grid on;

% Potencia del pico principal vs frame
figure;
plot(frames_valid, pow1, '-o');
xlabel('Número de frame');
ylabel('Potencia del pico principal (PSD)');
title('Evolución temporal de la potencia del pico principal (ROI fija)');
grid on;
fprintf('\nResultados guardados en %s\n', nombre_mat);

%% -------- GUARDAR A .MAT --------
save(nombre_mat, 'resultados');
fprintf('\nResultados guardados en %s\n', nombre_mat);