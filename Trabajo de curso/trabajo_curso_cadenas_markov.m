% ================================================================
%  Trabajo de curso - Cadenas de Markov y Mantenimiento Industrial
%  Métodos Numéricos en Ingeniería Computacional
% ================================================================

clear; clc; close all;

% ================================================================
%  Definición del sistema
% ================================================================

% Estados del sistema
states = {'Optimo','Leve','Moderado','Severo','Fallo'};
transient_states = [1, 2, 3, 4];
absorbing_state  = 5;

% Matriz de transición base (sin mantenimiento)
P_base = [
    0.85  0.15  0     0     0               % Estado 1: Óptimo
    0     0.80  0.20  0     0               % Estado 2: Desgaste leve
    0     0     0.75  0.25  0               % Estado 3: Desgaste moderado
    0     0     0     0.70  0.30            % Estado 4: Desgaste severo
    0     0     0     0     1               % Estado 5: Fallo
];

% ================================================================
%  Diagrama de estados
% ================================================================

s = [1 1 2 2 3 3 4 4 5];
t = [1 2 2 3 3 4 4 5 5];
weights = [0.85 0.15 0.80 0.20 0.75 0.25 0.70 0.30 1.00];

G = digraph(s, t, weights, states);

figure;
plot(G, 'Layout','layered', 'EdgeLabel', G.Edges.Weight);
title('Diagrama de estados del sistema sin mantenimiento');

% ================================================================
%  Función auxiliar para análisis
% ================================================================

compute_markov_metrics = @(P) ...
    struct( ...
        'Q', P(transient_states, transient_states), ...
        'R', P(transient_states, absorbing_state) ...
    );

% ================================================================
%  Escenario 1: Sin mantenimiento
% ================================================================

disp('==================== SIN MANTENIMIENTO ====================');

results_base = compute_markov_metrics(P_base);
N_base = inv(eye(size(results_base.Q)) - results_base.Q);

T_base = N_base * ones(length(transient_states),1);
F_base = N_base * results_base.R;

disp('Matriz P'); disp(P_base);
disp('Matriz Q'); disp(results_base.Q);
disp('Matriz R'); disp(results_base.R);
disp('Matriz fundamental N'); disp(N_base);
disp('Tiempo medio hasta fallo'); disp(T_base);
disp('Probabilidad de absorción'); disp(F_base);

% ================================================================
%  Escenario 2: Mantenimiento imperfecto
% ================================================================

disp('================ MANTENIMIENTO IMPERFECTO =================');

P_imperfect = P_base;
P_imperfect(4,4) = 0.85;
P_imperfect(4,5) = 0.15;

results_imp = compute_markov_metrics(P_imperfect);
N_imp = inv(eye(size(results_imp.Q)) - results_imp.Q);

T_imp = N_imp * ones(length(transient_states),1);
F_imp = N_imp * results_imp.R;

disp('Matriz P'); disp(P_imperfect);
disp('Matriz Q'); disp(results_imp.Q);
disp('Matriz R'); disp(results_imp.R);
disp('Matriz fundamental N'); disp(N_imp);
disp('Tiempo medio hasta fallo'); disp(T_imp);
disp('Probabilidad de absorción'); disp(F_imp);

% ================================================================
%  Escenario 3: Mantenimiento perfecto
% ================================================================

disp('================= MANTENIMIENTO PERFECTO ==================');

P_perfect = P_base;
P_perfect(4,:) = [0 0.80 0 0 0.20];   % Retorno a estado leve

results_perf = compute_markov_metrics(P_perfect);
N_perf = inv(eye(size(results_perf.Q)) - results_perf.Q);

T_perf = N_perf * ones(length(transient_states),1);
F_perf = N_perf * results_perf.R;

disp('Matriz P'); disp(P_perfect);
disp('Matriz Q'); disp(results_perf.Q);
disp('Matriz R'); disp(results_perf.R);
disp('Matriz fundamental N'); disp(N_perf);
disp('Tiempo medio hasta fallo'); disp(T_perf);
disp('Probabilidad de absorción'); disp(F_perf);

% ================================================================
%  Análisis de sensibilidad
% ================================================================

disp('================= ANALISIS DE SENSIBILIDAD =================');

p_fail = [0.15 0.10 0.05];
T_sensitivity = zeros(length(p_fail),1);

for i = 1:length(p_fail)
    P_tmp = P_base;
    P_tmp(4,4) = 1 - p_fail(i);
    P_tmp(4,5) = p_fail(i);

    Q_tmp = P_tmp(transient_states, transient_states);
    N_tmp = inv(eye(size(Q_tmp)) - Q_tmp);

    T_tmp = N_tmp * ones(length(transient_states),1);
    T_sensitivity(i) = T_tmp(1);
end

sens_table = table(p_fail', T_sensitivity,'VariableNames', {'Prob_Fallo', 'Tiempo_Medio'});
disp(sens_table);

figure;
plot(p_fail, T_sensitivity, '-o','LineWidth',2);
xlabel('Probabilidad de fallo desde estado severo');
ylabel('Tiempo medio hasta fallo (estado óptimo)');
grid on;
