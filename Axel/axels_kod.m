%% Settings
n_values = [20 40 60 80 100];  % Problemstorlekar att jämföra
TL = 15; TR = 0;               % Vänster och höger kant
TT = 0; TB = 15;               % Topp och botten
iterations_number = 50;        % Slumpvandringar per punkt (Monte Carlo)

% Förbered tidtagning
matrix_times = zeros(size(n_values));
jacobi_times  = zeros(size(n_values));
monte_times   = zeros(size(n_values));

% "Warm-up" för MATLABs interna bibliotek (förhindrar spikar i första mätningen)
dummy_A = sparse(eye(100)); dummy_b = ones(100,1); dummy_u = dummy_A\dummy_b;

%% Loop över problemstorlekar
for idx = 1:length(n_values)
    
    n = n_values(idx);
    m = n - 2;
    N = m * m;
    
    %% ------------------------ Matrix solution ------------------------
    A = sparse(N,N);
    b = zeros(N,1);
    
    for i = 1:m
        for j = 1:m
            k = i + (j-1)*m;
            A(k,k) = 4;
            if i>1,  A(k,k-1) = -1; else, b(k) = b(k) + TL; end
            if i<m,  A(k,k+1) = -1; else, b(k) = b(k) + TR; end
            if j>1,  A(k,k-m) = -1; else, b(k) = b(k) + TB; end
            if j<m,  A(k,k+m) = -1; else, b(k) = b(k) + TT; end
        end
    end
    
    tic
    u = A\b;
    matrix_times(idx) = toc;
    
    if idx == length(n_values)
        U_1 = zeros(n,n);
        U_1(2:end-1,2:end-1) = reshape(u,m,m);
        U_1(:,1) = TL; U_1(:,end) = TR;
        U_1(1,:) = TB; U_1(end,:) = TT;
    end
    
    %% ------------------------ Jacobi iteration ------------------------
    U_2 = zeros(n,n);
    U_2(:,1) = TL; U_2(:,end) = TR;
    U_2(1,:) = TB; U_2(end,:) = TT;
    
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
        diff = max(abs(U_2(:) - U_old(:)));
    end
    jacobi_times(idx) = toc;
    
    if idx == length(n_values)
        U_2_plot = U_2;
    end
    
    %% ------------------------ Monte Carlo ------------------------
    U_3 = zeros(n,n);
    U_3(:,1) = TL; U_3(:,end) = TR;
    U_3(1,:) = TB; U_3(end,:) = TT;
    
    tic
    for i = 2:n-1
        for j = 2:n-1
            temp_sum = 0;
            for walk = 1:iterations_number
                x = i; y = j;
                while x > 1 && x < n && y > 1 && y < n
                    r = rand();
                    if r < 0.25, x = x + 1;
                    elseif r < 0.50, x = x - 1;
                    elseif r < 0.75, y = y + 1;
                    else, y = y - 1;
                    end
                end
                
                % Check boundaries
                if x == 1, temp_sum = temp_sum + TB;
                elseif x == n, temp_sum = temp_sum + TT;
                elseif y == 1, temp_sum = temp_sum + TL;
                else, temp_sum = temp_sum + TR;
                end
            end
            U_3(i,j) = temp_sum / iterations_number;
        end
    end
    monte_times(idx) = toc;
    
    if idx == length(n_values)
        U_3_plot = U_3;
    end
end

%% ------------------------ Resultat ------------------------
figure('Name', 'Analys av Laplace-lösare', 'Color', 'w')

% Temperaturfördelning (Ren layout utan tids-labels)
subplot(2,3,1)
pcolor(U_1); shading interp; colorbar; axis equal tight; title('Matrislösning')

subplot(2,3,2)
pcolor(U_2_plot); shading interp; colorbar; axis equal tight; title('Jacobi')

subplot(2,3,3)
pcolor(U_3_plot); shading interp; colorbar; axis equal tight; title('Monte Carlo')

% Prestanda (Skalningsgrafer)
subplot(2,3,4)
plot(n_values, matrix_times,'-o','LineWidth',1.5); grid on
xlabel('Nätupplösning n'); ylabel('Tid (s)'); title('Skalning: Matris')

subplot(2,3,5)
plot(n_values, jacobi_times,'-o','LineWidth',1.5); grid on
xlabel('Nätupplösning n'); ylabel('Tid (s)'); title('Skalning: Jacobi')

subplot(2,3,6)
plot(n_values, monte_times,'-o','LineWidth',1.5); grid on
xlabel('Nätupplösning n'); ylabel('Tid (s)'); title('Skalning: Monte Carlo')