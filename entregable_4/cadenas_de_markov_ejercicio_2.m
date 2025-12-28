%% Ejercicio 2 - Clasificación y Análisis de Matrices de Transición
clear; clc;

% --- Definición de Datos ---
names = {'Pa', 'Pb', 'Pc', 'Pd', 'Pe', 'Pf', 'Pg', 'Ph', 'Pi', 'Pj', 'Pk', 'Pl'};
matrices = {
    [1 0; 0.2 0.8], [1 0; 0 1], [0.7 0.3; 0.4 0.6], ...
    [0.9 0 0.1; 0.6 0.2 0.2; 0 0.6 0.4], [0.8 0 0.2; 0 1 0; 0.9 0 0.1], ...
    [0.25 0.05 0.7; 0.5 0.3 0.2; 0 0 1], [0.8 0.2 0; 0.2 0.4 0.4; 0 0.3 0.7], ...
    [0.4 0 0.4 0.2; 0.1 0.8 0 0.1; 0 0 1 0; 0.2 0 0.7 0.1], ...
    [0.8 0.2 0 0; 0.3 0.7 0 0; 0 0 0.9 0.1; 0 0 0.4 0.6], ...
    [0.5 0.1 0.3 0.1; 1 0 0 0; 0.2 0.1 0.7 0; 0 0 1 0], ...
    [0.2 0.3 0.4 0.1; 0 1 0 0; 0.2 0.3 0.4 0.1; 0.1 0.4 0.2 0.3], ...
    [0.9 0 0.1 0; 0.1 0.6 0.1 0.2; 0.3 0 0.7 0; 0 0 0 1]
};

%% 1. Comprobación de Candidatas a Absorbentes
fprintf('\n=== 1. ANÁLISIS DE CANDIDATAS A ABSORBENTES ===\n');
for i = 1:length(matrices)
    if any(diag(matrices{i}) == 1)
        fprintf('%s es una matriz candidata (contiene 1 en la diagonal).\n', names{i});
    end
end

%% 2. Cálculo de Matriz Fundamental para Absorbentes Confirmadas
% Según informe: Pa, Pb, Pf, Ph, Pk son las seleccionadas
true_abs_idx = ismember(names, {'Pa', 'Pb', 'Pf', 'Ph', 'Pk'});
true_abs_names = names(true_abs_idx);
true_abs_mats = matrices(true_abs_idx);

fprintf('\n=== 2. MATRIZ FUNDAMENTAL (F = (I - Q)^-1) ===\n');
for i = 1:length(true_abs_mats)
    [F, Q] = calculate_fundamental_matrix(true_abs_mats{i});
    
    fprintf('\n--- Resultados para %s ---\n', true_abs_names{i});
    fprintf('Matriz Q (Estados transitorios):\n');
    disp(Q);
    fprintf('Matriz Fundamental F:\n');
    disp(F);
end

%% 3. Comprobación de Matrices Regulares
% Matrices restantes del análisis previo
rem_idx = ismember(names, {'Pc', 'Pd', 'Pe', 'Pg', 'Pi', 'Pj', 'Pl'});
rem_names = names(rem_idx);
rem_mats = matrices(rem_idx);

fprintf('\n=== 3. COMPROBACIÓN DE MATRICES REGULARES ===\n');
for i = 1:length(rem_mats)
    fprintf('Análisis de %s: ', rem_names{i});
    is_regular_matrix(rem_mats{i});
end

%% 4. Vectores Fijos (Matrices Regulares Confirmadas)
% Según lógica: Pc, Pd, Pg, Pj
reg_idx = ismember(names, {'Pc', 'Pd', 'Pg', 'Pj'});
reg_names = names(reg_idx);
reg_mats = matrices(reg_idx);

fprintf('\n=== 4. VECTORES FIJOS (ESTACIONARIOS) ===\n');
for i = 1:length(reg_mats)
    fprintf('\nVector fijo para %s:\n', reg_names{i});
    pi = stationary_vector(reg_mats{i});
    disp(pi);
end

%% 5. Resumen Final de Clasificación
fprintf('\n=== RESUMEN FINAL DE CLASIFICACIÓN ===\n');
fprintf('Absorbentes: %s\n', strjoin(true_abs_names, ', '));
fprintf('Regulares:   %s\n', strjoin(reg_names, ', '));
fprintf('Otras:       %s\n', strjoin(setdiff(names, [true_abs_names, reg_names]), ', '));

% -------------------------------------------------------------------------
% MÉTODOS AUXILIARES
% -------------------------------------------------------------------------

function [F, Q] = calculate_fundamental_matrix(M)
    abs_states = [];
    n = size(M,1);
    for i = 1:n
        % Identifica si el estado es puramente absorbente
        if M(i,i) == 1 && all(M(i,[1:i-1,i+1:n]) == 0)
            abs_states(end+1) = i;
        end
    end
    trans_states = setdiff(1:n, abs_states);
    % Si hay estados transitorios, calcula la matriz F
    if isempty(trans_states)
        Q = M(trans_states, trans_states);
        F = inv(eye(size(Q)) - Q);
    end
end

function is_regular_matrix(M)
    k_max = 50;
    Mk = M;
    regular = false;
    for k = 1:k_max
        if all(Mk(:) > 0)
            fprintf('SÍ (k=%d)\n', k);
            regular = true;
            break;
        end
        Mk = Mk * M;
    end
    if ~regular
        fprintf('NO (dentro del rango k=1:%d)\n', k_max);
    end
end

function pi = stationary_vector(P)
    n = size(P, 1);
    % Resolvemos (P' - I)pi' = 0 con la restricción suma(pi) = 1
    A = [P' - eye(n); ones(1, n)];
    b = [zeros(n, 1); 1];
    pi = (A \ b)';
end