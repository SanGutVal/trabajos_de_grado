%% EXTRACCION + AUDITORIA DE FRAMES
clear; clc;

% ---------------- CONFIGURACION ----------------
videoFile = 'video_interferometria.avi';   % cambiar nombre
outDir    = 'frames_30s_audit';
tReq      = 30;                            % segundos a extraer
tolFracDt = 0.25;                          % tolerancia relativa para saltos temporales
guardarPNG = true;                         % mantener true para auditoria visible

% ---------------- VALIDACION INICIAL ----------------
if ~isfile(videoFile)
    error('No se encuentra el archivo: %s', videoFile);
end

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

v = VideoReader(videoFile);

fpsNom   = v.FrameRate;
durVideo = v.Duration;
tStop    = min(tReq, durVideo);
dtNom    = 1 / fpsNom;
tolDt    = tolFracDt * dtNom;

% Conteo nominal esperado en el intervalo [0, tStop)
expectedNominal = ceil(tStop * fpsNom - 1e-9);

% Preasignacion conservadora
nAlloc = max(expectedNominal + 20, 100);

frameIdx   = zeros(nAlloc,1);
tFrame     = nan(nAlloc,1);   % CurrentTime justo antes de leer el frame
tAfterRead = nan(nAlloc,1);   % CurrentTime justo despues de leer
meanI      = nan(nAlloc,1);
stdI       = nan(nAlloc,1);
madPrev    = nan(nAlloc,1);   % mean absolute difference vs frame anterior
dupPrev    = false(nAlloc,1);
fileNames  = strings(nAlloc,1);

inputClass = "";
sampleSavedBitDepth = NaN;
prevGray = [];
k = 0;

% ---------------- EXTRACCION SECUENCIAL ----------------
v.CurrentTime = 0;

while hasFrame(v) && v.CurrentTime < tStop
    tBefore = v.CurrentTime;
    frame = readFrame(v);
    tAfter  = v.CurrentTime;

    k = k + 1;

    % Convertir a gris preservando el tipo de dato
    if ndims(frame) == 3
        if isequal(frame(:,:,1), frame(:,:,2)) && isequal(frame(:,:,2), frame(:,:,3))
            frameGray = frame(:,:,1);
        else
            frameGray = rgb2gray(frame);
        end
    else
        frameGray = frame;
    end

    if k == 1
        inputClass = string(class(frameGray));
    end

    % Guardado lossless
    fname = sprintf('frame_%06d.png', k);
    if guardarPNG
        imwrite(frameGray, fullfile(outDir, fname));
        if k == 1
            infoPNG = imfinfo(fullfile(outDir, fname));
            sampleSavedBitDepth = infoPNG.BitDepth;
        end
    end

    % Metricas basicas por frame
    frameIdx(k)   = k;
    tFrame(k)     = tBefore;
    tAfterRead(k) = tAfter;
    meanI(k)      = mean(double(frameGray(:)));
    stdI(k)       = std(double(frameGray(:)));
    fileNames(k)  = fname;

    if ~isempty(prevGray)
        dupPrev(k) = isequal(frameGray, prevGray);
        madPrev(k) = mean(abs(double(frameGray(:)) - double(prevGray(:))));
    end

    prevGray = frameGray;
end

% Recorte
frameIdx   = frameIdx(1:k);
tFrame     = tFrame(1:k);
tAfterRead = tAfterRead(1:k);
meanI      = meanI(1:k);
stdI       = stdI(1:k);
madPrev    = madPrev(1:k);
dupPrev    = dupPrev(1:k);
fileNames  = fileNames(1:k);

% ---------------- AUDITORIA TEMPORAL ----------------
dtFrame = [NaN; diff(tFrame)];
irregularDt = false(k,1);
irregularDt(2:end) = abs(dtFrame(2:end) - dtNom) > tolDt;

frameGapEst = nan(k,1);
frameGapEst(2:end) = max(0, round(dtFrame(2:end) / dtNom) - 1);

% Deriva frente a una malla temporal nominal
tNominal = tFrame(1) + (0:k-1)' * dtNom;
timeDrift = tFrame - tNominal;

% Frames sospechosos
idxSuspicious = find(irregularDt | dupPrev | (frameGapEst > 0));

% Resumen
meta = struct();
meta.videoFile              = videoFile;
meta.outputFolder           = outDir;
meta.width                  = v.Width;
meta.height                 = v.Height;
meta.bitsPerPixelVideo      = v.BitsPerPixel;
meta.frameRateNominal       = fpsNom;
meta.dtNominal              = dtNom;
meta.videoDuration_s        = durVideo;
meta.requestedSeconds_s     = tReq;
meta.extractedWindow_s      = tStop;
meta.expectedNominalFrames  = expectedNominal;
meta.actualExtractedFrames  = k;
meta.missingNominalEstimate = max(0, expectedNominal - k);
meta.duplicateCount         = nnz(dupPrev);
meta.irregularDtCount       = nnz(irregularDt);
meta.estimatedGapEvents     = nnz(frameGapEst > 0);
meta.totalEstimatedMissingFromGaps = nansum(frameGapEst);
meta.firstFrameTime_s       = tFrame(1);
meta.lastFrameTime_s        = tFrame(end);
meta.maxAbsTimeDrift_s      = max(abs(timeDrift));
meta.inputClass             = inputClass;
meta.savedPNG_BitDepth      = sampleSavedBitDepth;

% Tabla de auditoria
auditTable = table( ...
    frameIdx, fileNames, tFrame, tAfterRead, dtFrame, timeDrift, ...
    meanI, stdI, madPrev, dupPrev, irregularDt, frameGapEst, ...
    'VariableNames', { ...
    'frame_idx','file_name','t_frame_s','t_after_read_s','dt_frame_s','time_drift_s', ...
    'mean_intensity','std_intensity','mad_prev','duplicate_prev','irregular_dt','estimated_missing_before_this' });

% Guardado
save(fullfile(outDir, 'audit_extraccion.mat'), 'meta', 'auditTable', 'idxSuspicious');
writetable(auditTable, fullfile(outDir, 'audit_extraccion.csv'));

% Reporte de texto
fid = fopen(fullfile(outDir, 'reporte_auditoria.txt'), 'w');
fprintf(fid, 'Archivo: %s\n', meta.videoFile);
fprintf(fid, 'Resolucion: %d x %d\n', meta.width, meta.height);
fprintf(fid, 'BitsPerPixel video: %g\n', meta.bitsPerPixelVideo);
fprintf(fid, 'Clase de dato de entrada: %s\n', meta.inputClass);
fprintf(fid, 'BitDepth PNG guardado: %g\n', meta.savedPNG_BitDepth);
fprintf(fid, 'FPS nominal: %.9f\n', meta.frameRateNominal);
fprintf(fid, 'Duracion video (s): %.9f\n', meta.videoDuration_s);
fprintf(fid, 'Ventana extraida (s): %.9f\n', meta.extractedWindow_s);
fprintf(fid, 'Frames nominales esperados: %d\n', meta.expectedNominalFrames);
fprintf(fid, 'Frames extraidos realmente: %d\n', meta.actualExtractedFrames);
fprintf(fid, 'Estimacion faltantes por conteo: %d\n', meta.missingNominalEstimate);
fprintf(fid, 'Duplicados consecutivos: %d\n', meta.duplicateCount);
fprintf(fid, 'Saltos temporales irregulares: %d\n', meta.irregularDtCount);
fprintf(fid, 'Eventos de gap estimados: %d\n', meta.estimatedGapEvents);
fprintf(fid, 'Frames faltantes estimados por gaps: %g\n', meta.totalEstimatedMissingFromGaps);
fprintf(fid, 'Max |drift temporal| (s): %.9g\n', meta.maxAbsTimeDrift_s);
fprintf(fid, 'Numero de frames sospechosos: %d\n', numel(idxSuspicious));
fclose(fid);

% Consola
fprintf('\n=== AUDITORIA COMPLETADA ===\n');
fprintf('Archivo: %s\n', meta.videoFile);
fprintf('Frames nominales esperados: %d\n', meta.expectedNominalFrames);
fprintf('Frames extraidos realmente: %d\n', meta.actualExtractedFrames);
fprintf('Estimacion faltantes por conteo: %d\n', meta.missingNominalEstimate);
fprintf('Duplicados consecutivos: %d\n', meta.duplicateCount);
fprintf('Saltos temporales irregulares: %d\n', meta.irregularDtCount);
fprintf('Frames faltantes estimados por gaps: %g\n', meta.totalEstimatedMissingFromGaps);
fprintf('BitDepth PNG guardado: %g\n', meta.savedPNG_BitDepth);
fprintf('Carpeta salida: %s\n', meta.outputFolder);
