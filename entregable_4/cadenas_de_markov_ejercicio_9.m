%% EJERCICIO 9: MODELO DE OCUPACIÓN EN BARBERÍA
% Análisis mediante Cadenas de Markov
% -------------------------------------------------------------------------

% --- (a) Construcción de la Matriz de Transición ---
% Definimos la matriz T donde T(i,j) es la probabilidad de pasar de i a j clientes.
% Estados: 0, 1, 2, y 3 (máximo permitido en el sistema)

T = [ 0.2, 0.5, 0.2, 0.1 ;   % Desde 0: Solo entran clientes
      0.2, 0.5, 0.2, 0.1 ;   % Desde 1: Sale uno y entran n
      0.0, 0.2, 0.5, 0.3 ;   % Desde 2: Sale uno, queda uno y entran n
      0.0, 0.0, 0.2, 0.8 ];  % Desde 3: Sale uno, quedan dos y entran n (saturación)

fprintf('--- Apartado (a): Matriz de Transición T ---\n');
disp(T);

% --- (b) Distribución tras 30 minutos (2 intervalos de 15 min) ---
% Estado inicial: Equiprobable entre 0 y 1 cliente.
p_ini = [0.5, 0.5, 0, 0]; 

% Evolución: p_final = p_ini * T^pasos
p_30min = p_ini * (T^2);

fprintf('--- Apartado (b): Distribución a los 30 minutos ---\n');
fprintf('P(X=0)=%.4f, P(X=1)=%.4f, P(X=2)=%.4f, P(X=3)=%.4f\n\n', p_30min);

% --- (c) Distribución Estacionaria (Largo Plazo) ---
% Resolvemos el sistema de ecuaciones para el estado de equilibrio.
% Se utiliza la transpuesta y se impone la condición de que la suma sea 1.

n = length(T);
MatrizEcuac = (T' - eye(n)); 
MatrizEcuac(n, :) = 1; % Condición de normalización (suma de probs = 1)

VectorB = [zeros(n-1, 1); 1];
dist_estacionaria = (MatrizEcuac \ VectorB)';

fprintf('--- Apartado (c): Distribución Estacionaria ---\n');
for i = 1:n
    fprintf('Estado %d clientes: %.4f\n', i-1, dist_estacionaria(i));
end