% Settings
n = 100;                   % Number of grid points per side
TL = 15; TR = 0;          % Left and Right edge
TT = 0; TB = 15;          % Top and Bottom edge
iterations_number = 2000;  % Number of iterations

%% Matrix solution
m = n - 2;              % Inner points for direction
N = m * m;              % Total unknown points

A = sparse(N,N);      % Wiring matrix
b = zeros(N,1);       % Edge temperature

% Build Matrix
for i = 1:m
    for j = 1:m
        
        k = i + (j - 1)*m;   % convert (i,j) → index
        
        A(k, k) = 4;
        
        % Left neighbour
        if i > 1
            A(k, k - 1) = -1;
        else
            b(k) = b(k) + TL;
        end
        
        % Right neighbour
        if i < m
            A(k, k + 1) = -1;
        else
            b(k) = b(k) + TR;
        end
        
        % Lower neighbour
        if j > 1
            A(k, k - m) = -1;
        else
            b(k) = b(k) + TB;
        end
        
        % Top neighbour
        if j < m
            A(k, k + m) = -1;
        else
            b(k) = b(k) + TT;
        end
        
    end
end

% Solve equation-system
tic
u = A\b;
matrix_solution_time = toc;

% Insert in mesh
U_1 = zeros(n,n);
U_1(2:end-1,2:end-1) = reshape(u,m,m);

% Edge values
U_1(:,1)=TL; U_1(:,end)=TR;
U_1(1,:)=TB; U_1(end,:)=TT;


%% Jacobi-iteration
% Create grid
U_2 = zeros(n, n);

% Set boundaries
U_2(:,1)   = TL; % Left
U_2(:,end) = TR; % Right
U_2(1,:)   = TB; % Bottom
U_2(end,:) = TT; % Top

% Temperature distribution (Jacobi-iterations)
tic
tolerance = 1e-6;
diff = Inf;
while diff > tolerance
    U_old = U_2;
    for i = 2:n-1
        for j = 2:n-1
            U_2(i,j) = 0.25*(U_2(i+1,j)+U_2(i-1,j)+U_2(i,j+1)+U_2(i,j-1));
        end
    end
    diff = max(max(abs(U_2-U_old)));
end
jacobi_iteration_time = toc;

%% Result
figure

% --- Subplot 1 ---
subplot(1,2,1)
pcolor(U_1); 
shading interp; 
colorbar; 
axis equal tight; 
title('Matrix solution')
xlabel(sprintf('Computing time = %.6f', matrix_solution_time))

% --- Subplot 2 ---
subplot(1,2,2)
pcolor(U_2); 
shading interp; 
colorbar; 
axis equal tight; 
title('Jacobi-iteration')
xlabel(sprintf('Computing time = %.6f', jacobi_iteration_time))
