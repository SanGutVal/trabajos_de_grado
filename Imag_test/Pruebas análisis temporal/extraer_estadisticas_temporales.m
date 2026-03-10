% Configuración de archivos y ROI
carpeta   = '.';          % carpeta donde están los frames
prefijo   = 'frame_';
extension = '.jpg';       % o '.png'

frame_ini = 101;            % primer índice
n_frames  = 15;          % número de frames a procesar

ancho_roi = 54;           % columnas (ejemplo)
alto_roi  = 36;           % filas

% Frecuencia de muestreo temporal (video)
fps = 30;                 % frames por segundo

%% Preasignación
mu = [];      % medias por columna: [n_frames x Wroi]
sd = [];      % desviaciones estándar por columna
frames_vec = zeros(n_frames,1);

%% Bucle principal
for k = 1:n_frames
    idx = frame_ini + k - 1;
    nombre = sprintf('%s%d%s', prefijo, idx, extension);
    ruta   = fullfile(carpeta, nombre);

    if ~isfile(ruta)
        warning('No se encontró %s, se omite.', ruta);
        continue;
    end

    img = imread(ruta);
    if size(img,3) == 3
        img_gray = rgb2gray(img);
    else
        img_gray = img;
    end
    img_d = double(img_gray);
    [H,W] = size(img_gray);

    % ROI fija centrada
    cx = round(W/2); cy = round(H/2);
    half_w = round(ancho_roi/2);
    half_h = round(alto_roi/2);

    xmin = max(cx-half_w+1,1);
    xmax = min(cx+half_w,  W);
    ymin = max(cy-half_h+1,1);
    ymax = min(cy+half_h,  H);

    roi = img_d(ymin:ymax, xmin:xmax);  % tamaño aprox alto_roi x ancho_roi

    if k == 1
        [Hroi,Wroi] = size(roi);
        if Hroi ~= alto_roi || Wroi ~= ancho_roi
            warning('ROI real (%dx%d) difiere de lo esperado (%dx%d).', ...
                Hroi, Wroi, alto_roi, ancho_roi);
        end
        mu = zeros(n_frames, Wroi);
        sd = zeros(n_frames, Wroi);
    end

    frames_vec(k) = idx;

    % Medias y desviaciones por columna dentro de la ROI
    for j = 1:Wroi
        col = roi(:,j);
        mu(k,j) = mean(col);
        sd(k,j) = std(col, 0);  % desviación estándar poblacional
    end
end

%% Guardar resultados temporales
nombre_mat = sprintf('estadisticas_temporales_%d_a_%d.mat', ...
                     frame_ini, frame_ini + n_frames - 1);

save(nombre_mat, 'mu', 'sd', 'frames_vec', ...
     'ancho_roi', 'alto_roi', 'frame_ini', 'n_frames', 'fps');

fprintf('Guardado %s con mu, sd y frames_vec.\n', nombre_mat);