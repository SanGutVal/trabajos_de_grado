%% FILTRADO DE FRAMES AUDITADOS PARA ANALISIS POSTERIOR
clear; clc;

% ---------------- CONFIGURACION ----------------
auditMatFile = fullfile('frames_30s_audit', 'audit_extraccion.mat');
sourceFramesDir = '';                 % '' usar meta.outputFolder del .mat
outDir = 'frames_30s_limpios';

% ---------------- CARGA ----------------
if ~isfile(auditMatFile)
    error('No se encuentra el archivo de auditoria: %s', auditMatFile);
end

S = load(auditMatFile);

if ~isfield(S, 'auditTable') || ~isfield(S, 'meta')
    error('El archivo .mat no contiene las variables esperadas: auditTable y meta.');
end

auditTable = S.auditTable;
meta = S.meta;

if isempty(sourceFramesDir)
    if isfield(meta, 'outputFolder') && isfolder(meta.outputFolder)
        sourceFramesDir = meta.outputFolder;
    else
        error('No se pudo inferir la carpeta de frames desde meta.outputFolder.');
    end
end

if ~isfolder(sourceFramesDir)
    error('No existe la carpeta de frames fuente: %s', sourceFramesDir);
end

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% ---------------- VALIDACIONES ----------------
requiredVars = ["frame_idx","file_name","t_frame_s","dt_frame_s", ...
                "time_drift_s","mad_prev","duplicate_prev", ...
                "irregular_dt","estimated_missing_before_this"];

for i = 1:numel(requiredVars)
    if ~ismember(requiredVars(i), string(auditTable.Properties.VariableNames))
        error('Falta la variable requerida en auditTable: %s', requiredVars(i));
    end
end

n = height(auditTable);

% Asegurar tipo string para nombres de archivo
if ~isa(auditTable.file_name, 'string')
    auditTable.file_name = string(auditTable.file_name);
end

% ---------------- CRITERIO DE DESCARTE ----------------
missingFile = false(n,1);
for i = 1:n
    missingFile(i) = ~isfile(fullfile(sourceFramesDir, auditTable.file_name(i)));
end

gapFlag = auditTable.estimated_missing_before_this > 0;
dupFlag = auditTable.duplicate_prev;
irrFlag = auditTable.irregular_dt;

% Regla principal: descartar cualquier frame sospechoso
discardMask = irrFlag | dupFlag | gapFlag | missingFile;
keepMask    = ~discardMask;

% ---------------- TABLAS DE SALIDA ----------------
keptTable = auditTable(keepMask, :);
droppedTable = auditTable(discardMask, :);

% Agregar razones de descarte
dropReason = strings(height(droppedTable),1);
for i = 1:height(droppedTable)
    reasons = strings(0,1);

    if droppedTable.irregular_dt(i)
        reasons(end+1) = "irregular_dt";
    end
    if droppedTable.duplicate_prev(i)
        reasons(end+1) = "duplicate_prev";
    end
    if droppedTable.estimated_missing_before_this(i) > 0
        reasons(end+1) = "gap_before";
    end
    if ~isfile(fullfile(sourceFramesDir, droppedTable.file_name(i)))
        reasons(end+1) = "missing_file";
    end

    if isempty(reasons)
        reasons = "unspecified";
    end

    dropReason(i) = strjoin(reasons, ';');
end

droppedTable.discard_reason = dropReason;

% ---------------- COPIA Y RENOMBRADO DE FRAMES LIMPIOS ----------------
nKeep = height(keptTable);
clean_file_name = strings(nKeep,1);
clean_idx = (1:nKeep).';

for i = 1:nKeep
    src = fullfile(sourceFramesDir, keptTable.file_name(i));
    dstName = sprintf('frame_clean_%06d.png', i);
    dst = fullfile(outDir, dstName);

    ok = copyfile(src, dst);
    if ~ok
        error('No se pudo copiar el archivo: %s', src);
    end

    clean_file_name(i) = dstName;
end

keptTable.clean_idx = clean_idx;
keptTable.clean_file_name = clean_file_name;

% Reordenar columnas para lectura mas clara
priorityVars = ["clean_idx","clean_file_name","frame_idx","file_name","t_frame_s", ...
                "dt_frame_s","time_drift_s","mad_prev","duplicate_prev", ...
                "irregular_dt","estimated_missing_before_this"];

existingPriorityVars = priorityVars(ismember(priorityVars, string(keptTable.Properties.VariableNames)));
otherVars = string(keptTable.Properties.VariableNames);
otherVars = otherVars(~ismember(otherVars, existingPriorityVars));
keptTable = keptTable(:, [cellstr(existingPriorityVars), cellstr(otherVars)]);

% ---------------- METADATOS DE FILTRADO ----------------
filterMeta = struct();
filterMeta.auditMatFile = auditMatFile;
filterMeta.sourceFramesDir = sourceFramesDir;
filterMeta.outputFolder = outDir;
filterMeta.originalFrameCount = n;
filterMeta.keptFrameCount = nKeep;
filterMeta.droppedFrameCount = height(droppedTable);
filterMeta.keptFraction = nKeep / max(n,1);
filterMeta.droppedFraction = height(droppedTable) / max(n,1);
filterMeta.count_irregular_dt = nnz(irrFlag);
filterMeta.count_duplicate_prev = nnz(dupFlag);
filterMeta.count_gap_before = nnz(gapFlag);
filterMeta.count_missing_file = nnz(missingFile);
filterMeta.frameRateNominal = meta.frameRateNominal;
filterMeta.dtNominal = meta.dtNominal;
filterMeta.videoFile = meta.videoFile;

% Vector temporal limpio para uso posterior
t_clean = keptTable.t_frame_s;
clean_idx = keptTable.clean_idx;

% ---------------- GUARDADO ----------------
save(fullfile(outDir, 'frames_limpios_index.mat'), ...
     'filterMeta', 'keptTable', 'droppedTable', 't_clean', 'clean_idx');

writetable(keptTable, fullfile(outDir, 'frames_limpios_index.csv'));
writetable(droppedTable, fullfile(outDir, 'frames_descartados.csv'));

% Reporte TXT
fid = fopen(fullfile(outDir, 'reporte_filtrado.txt'), 'w');
fprintf(fid, 'SCRIPT 2 - FILTRADO DE FRAMES AUDITADOS\n\n');
fprintf(fid, 'Archivo de auditoria: %s\n', filterMeta.auditMatFile);
fprintf(fid, 'Video original: %s\n', filterMeta.videoFile);
fprintf(fid, 'Carpeta fuente: %s\n', filterMeta.sourceFramesDir);
fprintf(fid, 'Carpeta salida: %s\n\n', filterMeta.outputFolder);

fprintf(fid, 'Frames originales: %d\n', filterMeta.originalFrameCount);
fprintf(fid, 'Frames conservados: %d\n', filterMeta.keptFrameCount);
fprintf(fid, 'Frames descartados: %d\n', filterMeta.droppedFrameCount);
fprintf(fid, 'Fraccion conservada: %.6f\n', filterMeta.keptFraction);
fprintf(fid, 'Fraccion descartada: %.6f\n\n', filterMeta.droppedFraction);

fprintf(fid, 'Causas de descarte detectadas en la auditoria:\n');
fprintf(fid, '  irregular_dt: %d\n', filterMeta.count_irregular_dt);
fprintf(fid, '  duplicate_prev: %d\n', filterMeta.count_duplicate_prev);
fprintf(fid, '  gap_before: %d\n', filterMeta.count_gap_before);
fprintf(fid, '  missing_file: %d\n', filterMeta.count_missing_file);
fclose(fid);

% ---------------- MENSAJES ----------------
fprintf('\n=== FILTRADO COMPLETADO ===\n');
fprintf('Frames originales: %d\n', filterMeta.originalFrameCount);
fprintf('Frames conservados: %d\n', filterMeta.keptFrameCount);
fprintf('Frames descartados: %d\n', filterMeta.droppedFrameCount);
fprintf('Carpeta de salida: %s\n', filterMeta.outputFolder);
fprintf('Indice limpio: %s\n', fullfile(outDir, 'frames_limpios_index.csv'));
fprintf('Descartados: %s\n', fullfile(outDir, 'frames_descartados.csv'));
