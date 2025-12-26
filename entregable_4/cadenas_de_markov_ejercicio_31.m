%% SOLUCIÓN: CADENAS DE MARKOV - LABERINTO DEL RATÓN
% Análisis de estados absorbentes y tiempos de tránsito
% -------------------------------------------------------------------------

% 1. Definición del Modelo
% Estados: 1=A, 2=B, 3=C, 4=D, 5=E, 6=F
% Matriz de transición global (M)
M = [1,   0,   0,   0,   0,   0;   % Nodo A (Absorbente)
     0.5, 0,   0.5, 0,   0,   0;   % Nodo B
     0,   1/3, 0,   1/3, 1/3, 0;   % Nodo C
     0,   0,   1,   0,   0,   0;   % Nodo D
     0,   0,   0.5, 0,   0,   0.5; % Nodo E
     0,   0,   0,   0,   0,   1];  % Nodo F (Absorbente)

% Segmentación de la matriz para el análisis de absorción
% Identificamos índices: Transitorios (idxT) y Absorbentes (idxA)
idxA = [1, 6]; 
idxT = [2, 3, 4, 5];

% Submatrices del sistema
Q_mat = M(idxT, idxT);
R_mat = M(idxT, idxA);

% 2. Cálculo de la Matriz Fundamental (N)
% N = (I - Q)^-1 representa el tiempo esperado en cada estado transitorio
N_fund = inv(eye(size(Q_mat)) - Q_mat);

% --- Resolución de Apartados ---

% (a) Visitas a C desde D
% En idxT: B=1, C=2, D=3, E=4. Queremos saber de D(3) a C(2)
fprintf('Resultados del análisis:\n');
visitas_D_C = N_fund(3, 2);
fprintf('a) Visitas esperadas a la habitación C iniciando en D: %.2f\n', visitas_D_C);

% (b) Tiempo de absorción desde E
% Sumamos la fila correspondiente a E (índice 4 en idxT)
tiempo_desde_E = sum(N_fund(4, :));
fprintf('b) Número de habitaciones que se espera visitar desde E antes de llegar a un estado de absorción (A o F): %.2f\n', tiempo_desde_E);

% (c) y (d) Matriz de probabilidades de absorción (B_prob = N * R)
B_prob = N_fund * R_mat;

% (c) Probabilidad de acabar en A (col 1) desde C (fila 2 de idxT)
prob_C_A = B_prob(2, 1);
fprintf('c) Probabilidad de terminar en A si se inicia en C: %.4f\n', prob_C_A);

% (d) Probabilidad de acabar en F (col 2) desde E (fila 4 de idxT)
prob_E_F = B_prob(4, 2);
fprintf('d) Probabilidad de terminar en F si se inicia en E: %.4f\n', prob_E_F);

% Mostrar matrices calculadas para verificación
disp('-----------------------------------------');
disp('Matriz Fundamental (N):'); disp(N_fund);
disp('Matriz de Absorción (B):'); disp(B_prob);
