%% TEMPERATURPROJEKT - MATRISLÖSNING (Direkt metod)
% Löser Laplaces ekvation genom att ställa upp systemet A*x = b
% och lösa det direkt (utan iterationer).

clear; clc; close all;

%% 1. PARAMETRAR
% ---------------------------------------------------------
% A. Nätupplösning
N = 30;             % Antal inre punkter i x- och y-led

% B. Randvillkor
T_top    = 100;     % Temperatur övre kant
T_bottom = 0;       % Temperatur nedre kant
T_left   = 50;      % Temperatur vänster kant
T_right  = 50;      % Temperatur höger kant

% (Inga iterationsinställningar behövs här eftersom vi löser det exakt direkt)
% ---------------------------------------------------------

%% 2. BYGG MATRIS OCH HÖGERLED
fprintf('Bygger matris för N=%d (Totalt %d ekvationer)...\n', N, N^2);
tic; % Starta tidtagning

% Totalt antal okända punkter
n_unknowns = N * N;

% Vi använder "sparse" (gles matris) för att spara minne.
% Matrisen A kommer ha 5 diagonalband (mittenpunkten + 4 grannar).
A = spalloc(n_unknowns, n_unknowns, 5*n_unknowns);
b = zeros(n_unknowns, 1);

% Vi loopar igenom varje INRE punkt i rutnätet
% (i = radindex, j = kolumnindex)
for j = 1:N 
    for i = 1:N
        
        % Räkna ut vilket nummer (k) denna punkt har i ekvationssystemet (1 till N^2)
        k = i + (j-1) * N;
        
        % Själva ekvationen: 4*T(i,j) - (grannarna) = 0
        % Vi sätter koefficienten 4 på diagonalen
        A(k, k) = 4;
        
        % --- HANTERA GRANNAR ---
        
        % 1. UPPÅT (i-1)
        if i > 1
            % Om inre granne: Koefficient -1 i matrisen
            A(k, k-1) = -1;
        else
            % Om rand (vägg): Flytta värdet till högerledet (b)
            b(k) = b(k) + T_top;
        end
        
        % 2. NEDÅT (i+1)
        if i < N
            A(k, k+1) = -1;
        else
            b(k) = b(k) + T_bottom;
        end
        
        % 3. VÄNSTER (j-1)
        if j > 1
            A(k, k-N) = -1;
        else
            b(k) = b(k) + T_left;
        end
        
        % 4. HÖGER (j+1)
        if j < N
            A(k, k+N) = -1;
        else
            b(k) = b(k) + T_right;
        end
        
    end
end

%% 3. LÖS EKVATIONSSYSTEMET
% Här sker magin: MATLAB löser Ax = b
T_vec = A \ b;

time_taken = toc; % Stoppa tidtagning

%% 4. RESULTAT OCH VISUALISERING
% Omforma lösningsvektorn (1D) tillbaka till en platta (2D)
T_inner = reshape(T_vec, N, N);

% Skapa den fullständiga plattan med ränder för visualisering
T_plot = zeros(N+2, N+2);

% Fyll i ränderna
T_plot(1, :)     = T_top;
T_plot(end, :)   = T_bottom;
T_plot(:, 1)     = T_left;
T_plot(:, end)   = T_right;

% Fyll i mitten med vår lösning
T_plot(2:end-1, 2:end-1) = T_inner;

% Skriv ut info
fprintf('Lösning klar på %.4f sekunder.\n', time_taken);

% Rita
figure(1); clf;
imagesc(T_plot); 
colormap('jet'); 
colorbar;
axis equal; axis tight;
title(['Temperaturfördelning (Direkt lösning, N=', num2str(N), ')']);
xlabel('x'); ylabel('y');

% Konturlinjer
hold on;
contour(T_plot, 10, 'k', 'LineWidth', 1.5);
hold off;




