%% -------- CARGAR RESULTADOS --------
nombre_mat = 'resultados_frames_101_a_116.mat';   % ajustar nombre
load(nombre_mat, 'resultados');

% Extraer solo elementos con datos válidos en ROI fija
frames = [resultados.frame];
mask_valid = ~arrayfun(@(r) isempty(r.freqs_fix) || all(isnan(r.freqs_fix)), resultados);

frames_valid = frames(mask_valid);
R = resultados(mask_valid);
nF = numel(R);

freqs_fix = nan(nF,3);      % columnas = intervalos [0-0.1], [0.1-0.25], [0.25-0.4]
pows_fix  = nan(nF,3);

for k = 1:nF
    f = R(k).freqs_fix;
    p = R(k).potencias_fix;
    % asegura longitud 3 (por si en algún frame hubo menos)
    L = min(3, numel(f));
    freqs_fix(k,1:L) = f(1:L);
    pows_fix(k,1:L)  = p(1:L);
end

%% -------- GRÁFICA: FRECUENCIA vs FRAME --------
figure;
plot(frames_valid, freqs_fix(:,1), '-o', 'LineWidth',1.5); hold on;
plot(frames_valid, freqs_fix(:,2), '-s', 'LineWidth',1.5);
plot(frames_valid, freqs_fix(:,3), '-^', 'LineWidth',1.5);
hold off; grid on;
xlabel('Número de frame');
ylabel('Frecuencia (ciclos/píxel)');
title('Evolución temporal de los máximos locales');
legend({'[0, 0.1)', '[0.1, 0.25)', '[0.25, 0.4)'}, 'Location','best');

%% -------- GRÁFICA: POTENCIA vs FRAME --------
figure;
plot(frames_valid, pows_fix(:,1), '-o', 'LineWidth',1.5); hold on;
plot(frames_valid, pows_fix(:,2), '-s', 'LineWidth',1.5);
plot(frames_valid, pows_fix(:,3), '-^', 'LineWidth',1.5);
hold off; grid on;
xlabel('Número de frame');
ylabel('Potencia del máximo local (PSD)');
title('Evolución temporal de la potencia');
legend({'[0, 0.1)', '[0.1, 0.25)', '[0.25, 0.4)'}, 'Location','best');
