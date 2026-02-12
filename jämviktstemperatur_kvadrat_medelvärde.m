% --- 1. Inställningar ---
n = 40;                 % Antal delsteg (mesh size)
TL = 50; TR = 50;       % Temperatur Vänster (Left) och Höger (Right)
TT = 100; TB = 0;       % Temperatur Uppe (Top) och Nere (Bottom)

% --- 2. Förberedelser ---
m = n - 1;              % Antal inre punkter per rad/kolumn
N = m^2;                % Totalt antal inre punkter (storleken på vektorn u)
A = eye(N) * 4;         % Skapar huvuddiagonalen med fyror (4*u_ij)
b = zeros(N, 1);        % Högerledet (vektorn med kända randtemperaturer)

% --- 3. Bygg matrisen A och vektorn b ---
% Vi går igenom varje inre punkt (i,j) och ställer upp ekvationen:
% 4*u_ij - u_grannar = kända_temperaturer
for i = 1:m
    for j = 1:m
        k = i + (j-1)*m; % Omvandlar 2D-koordinat (i,j) till ett radnummer k
        
        % Kolla granne till Vänster (u_i-1,j)
        if i == 1, b(k) = b(k) + TL; else A(k, k-1) = -1; end
        
        % Kolla granne till Höger (u_i+1,j)
        if i == m, b(k) = b(k) + TR; else A(k, k+1) = -1; end
        
        % Kolla granne Nere (u_i,j-1)
        if j == 1, b(k) = b(k) + TB; else A(k, k-m) = -1; end
        
        % Kolla granne Uppe (u_i,j+1)
        if j == m, b(k) = b(k) + TT; else A(k, k+m) = -1; end
    end
end

% --- 4. Beräkna temperaturerna ---
u_vektor = A \ b;       % Löser ekvationssystemet direkt

% --- 5. Omvandla vektorn till en platta och rita ---
% Vi lägger in de beräknade värdena i mitten av en tom platta
U = zeros(n+1, n+1);
U(2:n, 2:n) = reshape(u_vektor, m, m);

% Lägg till randvärdena på kanterna för visualisering
U(:, 1) = TB; U(:, n+1) = TT; 
U(1, :) = TL; U(n+1, :) = TR;

% Visa resultatet
surf(U)
xlabel('x'); ylabel('y'); zlabel('Temperatur u(x,y)');
title('Diskret temperaturfördelning');