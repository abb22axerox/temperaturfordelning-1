
%% --- Globala Inställningar (Ändra här för dina experiment) ---
n = 50;                 % Nätupplösning (Fråga 1, 5, 6)
shape = 'L-shape';       % Välj: 'square' eller 'L-shape' (Fråga 3, 4)
method = 'direct';      % Välj: 'direct' eller 'monte_carlo' (Fråga 5, 6)

% Randvillkor (Fråga 2)
BC_top = 100; BC_bottom = 0; BC_left = 50; BC_right = 50;

%% --- 1. Geometri och Maskning ---
% Skapa en matris som definierar objektets form
mask = ones(n+1, n+1);
if strcmp(shape, 'L-shape')
    % Tar bort nedre högra kvadranten för att skapa ett L
    mask(round(n/2):end, round(n/2):end) = 0; 
end

% Hitta alla inre punkter (där vi ska beräkna temperatur)
% En punkt är "inre" om den är 1 i masken och inte ligger på ytterkanten
inner_mask = zeros(n+1, n+1);
inner_mask(2:n, 2:n) = mask(2:n, 2:n);
inner_indices = find(inner_mask == 1);
N = length(inner_indices);

% Skapa en karta för att gå från 2D-koordinat till ett index i vektorn u
node_map = zeros(n+1, n+1);
node_map(inner_indices) = 1:N;

%% --- 2. Numerisk Lösning (Direkt Metod) ---
if strcmp(method, 'direct')
    fprintf('Löser systemet med direkt metod (Gles matris)...\n');
    tic;
    
    A = sparse(N, N); % Skapar en gles matris för effektivitet
    b = zeros(N, 1);
    
    % Gå igenom alla inre noder och ställ upp ekvationerna
    for k = 1:N
        [r, c] = find(node_map == k);
        A(k, k) = 4;
        
        % Kolla grannar: [Rad, Kol, Randvärde]
        neighbors = [r, c-1, BC_left;   % Vänster
                     r, c+1, BC_right;  % Höger
                     r-1, c, BC_bottom; % Ner
                     r+1, c, BC_top];   % Upp
        
        for g = 1:4
            nr = neighbors(g,1); 
            nc = neighbors(g,2);
            val = neighbors(g,3);
            
            if node_map(nr, nc) > 0
                % Grannen är en inre punkt - lägg till i matrisen A
                A(k, node_map(nr, nc)) = -1;
            elseif mask(nr, nc) == 1
                % Grannen är en randpunkt - lägg till i vektor b
                b(k) = b(k) + val;
            end
        end
    end
    
    u_vec = A \ b;
    solve_time = toc;

%% --- 3. Monte Carlo-metoden ---
elseif strcmp(method, 'monte_carlo')
    fprintf('Löser systemet med Monte Carlo (Väntetid beror på n)...\n');
    tic;
    u_vec = zeros(N, 1);
    walks_per_node = 50; % Ju fler walks, desto högre noggrannhet
    
    for k = 1:N
        [start_r, start_c] = find(node_map == k);
        temp_sum = 0;
        
        for w = 1:walks_per_node
            r = start_r; c = start_c;
            while node_map(r, c) > 0 % Så länge vi är inne i plattan
                move = randi(4);
                if move == 1, r = r - 1;
                elseif move == 2, r = r + 1;
                elseif move == 3, c = c - 1;
                else, c = c + 1;
                end
            end
            % Vi har nått en rand, kolla vilket värde den har
            if r == 1, temp_sum = temp_sum + BC_bottom;
            elseif r == n+1, temp_sum = temp_sum + BC_top;
            elseif c == 1, temp_sum = temp_sum + BC_left;
            elseif c == n+1, temp_sum = temp_sum + BC_right;
            end
        end
        u_vec(k) = temp_sum / walks_per_node;
    end
    solve_time = toc;
end

%% --- 4. Visualisering ---
% Skapa en resultatmatris fylld med NaN (så att tomma ytor inte ritas)
U = nan(n+1, n+1);
U(inner_indices) = u_vec;

% Lägg till randvärden för visualisering
U(1, mask(1,:)==1) = BC_bottom;
U(n+1, mask(n+1,:)==1) = BC_top;
U(mask(:,1)==1, 1) = BC_left;
U(mask(:,n+1)==1, n+1) = BC_right;

% Rita grafen
figure;
surf(U);
shading interp; % Gör ytan mjukare
colorbar;
xlabel('x (Position)'); ylabel('y (Position)'); zlabel('Temperatur');
title(sprintf('Form: %s | Metod: %s | Tid: %.4f s', shape, method, solve_time));

fprintf('Klart! Beräkningstid: %.4f sekunder.\n', solve_time);