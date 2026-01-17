clear all; close all; clc;

%% Parameters
c0 = 343; rho0 = 1.21;
Lx = 1; Ly = 0.5;
dx = 0.002; dy = dx;
Nx = round(Lx/dx); Ny = round(Ly/dy);

% ABH Geometry
R = 0.047; Dmax = 0.07;
d = 0.0035; dc = d/dx; dtt = 0.004;
n_branches = 20;   
L = n_branches*(d+dtt); Lc = L/dx; 
x0 = 0.4; branch_x0 = x0/dx; branch_x1 = branch_x0+Lc; 
abh_start_x = x0; abh_end_x = x0 + L;

% Optimized parameters
m = 1.8; f0 = 1800;

duct_yc = round(Ny/2);
duct_y0 = duct_yc - floor(round(2*R/dy)/2);
duct_y1 = duct_yc + floor(round(2*R/dy)/2);

branch_x_positions = round(linspace(branch_x0, branch_x1, n_branches));
xi = linspace(0,1,n_branches);
branch_depths_top = round((d/dy) + (round(Dmax/dy) - d/dy) * xi.^m);
branch_depths_bot = round((d/dy) + (round(Dmax/dy) - d/dy) * xi.^m);

% Probe Position
probe_y_meters = (R - 0.01); 
probe_y = duct_yc + round(probe_y_meters/dy); 

Nprobes = 100;
probe_x = linspace(0.2, 0.78, Nprobes);
probe_ix = max(1, min(Nx, round(probe_x/dx)));

%% Simulation Loop
U_cases = [0, 20]; 
v_loc_results = zeros(Nprobes, length(U_cases)); 

for u_idx = 1:length(U_cases)
    U = U_cases(u_idx);
    Mach = U / c0;
    fprintf('Running U = %d m/s...\n', U);
    
    dt = 0.3 * min(dx,dy)/(c0 + abs(U) + 20); 
    Nt = ceil(0.015/dt); t_vec = (0:Nt-1)*dt;
    n_min = round(0.003 / dt);
    
    p = zeros(Nx, Ny); vx = zeros(Nx+1, Ny); vy = zeros(Nx, Ny+1);
    p_record = zeros(Nprobes, Nt);
    
    % Create ABH Geometry
    mask = zeros(Nx, Ny);
    for i = 1:Nx, for j = duct_y0:duct_y1, mask(i,j) = 1; end; end
    for b = 1:n_branches
        ix = branch_x_positions(b); ox0 = max(1, ix - floor(dc/2)); ox1 = min(Nx, ix + floor(dc/2));
        for i = ox0:ox1
            for j = duct_y1+1:min(Ny, duct_y1+branch_depths_top(b)), mask(i,j)=1; end
            for j = max(1, duct_y0-branch_depths_bot(b)):duct_y0-1, mask(i,j)=1; end
        end
    end
    
    mask_vx = zeros(Nx+1, Ny); mask_vy = zeros(Nx, Ny+1);
    for i = 1:Nx+1, for j = 1:Ny
        if (i-1>=1 && mask(i-1,j)==1 && i<=Nx && mask(i,j)==1), mask_vx(i,j)=1; end
    end; end
    for i = 1:Nx, for j = 1:Ny+1
        if (j-1>=1 && mask(i,j-1)==1 && j<=Ny && mask(i,j)==1), mask_vy(i,j)=1; end
    end; end
    
    % 1. PML Setup
    pml_size = 20;  % PML layer size in grid points
    sigma_max = 2000;  % Maximum absorption coefficient
    
    % PML absorption profile for pressure (Nx x Ny)
    sigma_x_p = zeros(Nx, Ny);
    sigma_y_p = zeros(Nx, Ny);
    
    for i = 1:pml_size
        % Left PML
        sigma_value = sigma_max * ((pml_size - i + 1) / pml_size)^2;
        sigma_x_p(i, :) = sigma_value;
        
        % Right PML
        sigma_value = sigma_max * ((pml_size - i + 1) / pml_size)^2;
        sigma_x_p(Nx-i+1, :) = sigma_value;
    end
    
    for j = 1:pml_size
        % Bottom PML
        sigma_value = sigma_max * ((pml_size - j + 1) / pml_size)^2;
        sigma_y_p(:, j) = max(sigma_y_p(:, j), sigma_value);
        
        % Top PML
        sigma_value = sigma_max * ((pml_size - j + 1) / pml_size)^2;
        sigma_y_p(:, Ny-j+1) = max(sigma_y_p(:, Ny-j+1), sigma_value);
    end
    
    % PML for vx (Nx+1 x Ny) - different size!
    sigma_x_vx = zeros(Nx+1, Ny);
    sigma_y_vx = zeros(Nx+1, Ny);
    
    for i = 1:pml_size
        % Left PML for vx
        sigma_value = sigma_max * ((pml_size - i + 1) / pml_size)^2;
        sigma_x_vx(i, :) = sigma_value;
        
        % Right PML for vx
        sigma_value = sigma_max * ((pml_size - i + 1) / pml_size)^2;
        sigma_x_vx(Nx+1-i+1, :) = sigma_value;
    end
    
    for j = 1:pml_size
        sigma_value = sigma_max * ((pml_size - j + 1) / pml_size)^2;
        sigma_y_vx(:, j) = max(sigma_y_vx(:, j), sigma_value);
        sigma_y_vx(:, Ny-j+1) = max(sigma_y_vx(:, Ny-j+1), sigma_value);
    end
    
    % PML for vy (Nx x Ny+1) - different size!
    sigma_x_vy = zeros(Nx, Ny+1);
    sigma_y_vy = zeros(Nx, Ny+1);
    
    for i = 1:pml_size
        sigma_value = sigma_max * ((pml_size - i + 1) / pml_size)^2;
        sigma_x_vy(i, :) = sigma_value;
        sigma_x_vy(Nx-i+1, :) = sigma_value;
    end
    
    for j = 1:pml_size
        sigma_value = sigma_max * ((pml_size - j + 1) / pml_size)^2;
        sigma_y_vy(:, j) = max(sigma_y_vy(:, j), sigma_value);
        sigma_y_vy(:, Ny+1-j+1) = max(sigma_y_vy(:, Ny+1-j+1), sigma_value);
    end
    
    % PML absorption factors
    absorb_x_p = exp(-sigma_x_p * dt);
    absorb_y_p = exp(-sigma_y_p * dt);
    absorb_x_vx = exp(-sigma_x_vx * dt);
    absorb_y_vx = exp(-sigma_y_vx * dt);
    absorb_x_vy = exp(-sigma_x_vy * dt);
    absorb_y_vy = exp(-sigma_y_vy * dt);
    
    % 2. Branch Damping Setup (simpler version)
    damping = ones(Nx, Ny);
    for i = 1:Nx
        for j = 1:Ny
            if mask(i,j) == 1 && (j < duct_y0 || j > duct_y1)
                damping(i,j) = 0.97;  % Constant damping in branches
            end
        end
    end
    
    % FDTD Update with PML and damping
    src_x = round(0.1/dx);
    
    for n = 1:Nt
        t = (n-1)*dt;
        
        % Velocity updates with PML
        dpdx = zeros(Nx+1, Ny);
        dpdx(2:Nx,:) = (p(2:Nx,:) - p(1:Nx-1,:))/dx;
        dvxdx = zeros(Nx+1, Ny);
        dvxdx(2:Nx,:) = (vx(2:Nx,:) - vx(1:Nx-1,:))/dx;
        
        vx = vx - (dt/rho0).*dpdx.*mask_vx - dt*U.*dvxdx.*mask_vx;
        
        % Apply PML absorption to vx (correct dimensions)
        vx = vx .* absorb_x_vx .* absorb_y_vx;
        
        dpdy = zeros(Nx, Ny+1);
        dpdy(:,2:Ny) = (p(:,2:Ny) - p(:,1:Ny-1))/dy;
        vy = vy - (dt/rho0).*dpdy.*mask_vy;
        
        % Apply PML absorption to vy (correct dimensions)
        vy = vy .* absorb_x_vy .* absorb_y_vy;
        
        % Pressure update with PML and damping
        divv = (vx(2:Nx+1,:) - vx(1:Nx,:))/dx + (vy(:,2:Ny+1) - vy(:,1:Ny))/dy;
        p_grad_x = zeros(Nx, Ny);
        p_grad_x(2:Nx,:) = (p(2:Nx,:) - p(1:Nx-1,:))/dx;
        
        p = p - (rho0*c0^2*dt).*divv - dt*U.*p_grad_x;
        
        % Apply PML absorption to pressure
        p = p .* absorb_x_p .* absorb_y_p;
        
        % Apply branch damping
        p = p .* damping;
        
        % Apply mask (set pressure to zero outside geometry)
        p = p .* mask;
        
        % Source
        if mask(src_x, duct_yc)==1
            src_value = sin(2*pi*f0*t);
            if t < 0.002
                src_value = (t/0.002) * src_value;
            end
            p(src_x, duct_yc) = p(src_x, duct_yc) + src_value;
        end
        
        % Record at probe line
        for ip = 1:Nprobes
            p_record(ip, n) = p(probe_ix(ip), probe_y);
        end
        
        % Live Visualization
        if mod(n, 50) == 0
            figure(1);
            subplot(2,2,[1 3]);
            imagesc((0:Nx-1)*dx, (0:Ny-1)*dy, p');
            hold on;
            [Y, X] = meshgrid(1:Ny, 1:Nx);
            scatter(X(mask==1)*dx, Y(mask==1)*dy, 1, 'k', 'filled', 'AlphaData', 0.3);
            xline(abh_start_x, 'r--', 'LineWidth', 2);
            xline(abh_end_x, 'r--', 'LineWidth', 2);
            plot(probe_x, probe_y_meters*ones(size(probe_x)), 'y-', 'LineWidth', 2);
            
            % Mark PML regions
            rectangle('Position', [0, 0, pml_size*dx, Ly], 'EdgeColor', 'g', 'LineWidth', 1, 'LineStyle', '--');
            rectangle('Position', [Lx-pml_size*dx, 0, pml_size*dx, Ly], 'EdgeColor', 'g', 'LineWidth', 1, 'LineStyle', '--');
            
            hold off;
            axis equal tight;
            xlabel('x (m)'); ylabel('y (m)');
            title(sprintf('Pressure Field (U=%d m/s, t=%.3f ms)', U, t*1000));
            colorbar;
            clim([-0.05 0.05]);
            
            subplot(2,2,2);
            plot((0:Nx-1)*dx, p(:, duct_yc), 'b-', 'LineWidth', 1.5);
            xline(abh_start_x, 'r--', 'LineWidth', 2);
            xline(abh_end_x, 'r--', 'LineWidth', 2);
            xlabel('x (m)'); ylabel('Pressure');
            title('Centerline Pressure');
            grid on;
            xlim([0 Lx]);
            ylim([-0.1 0.1]);
            
            subplot(2,2,4);
            plot(probe_x, p_record(:, n), 'r-', 'LineWidth', 1.5);
            xline(abh_start_x, 'r--', 'LineWidth', 2);
            xline(abh_end_x, 'r--', 'LineWidth', 2);
            xlabel('x (m)'); ylabel('Pressure');
            title(sprintf('Pressure at y=%.3f m', probe_y_meters));
            grid on;
            xlim([min(probe_x) max(probe_x)]);
            ylim([-0.1 0.1]);
            
            drawnow;
        end
    end
    
    % CORRECT Wave Decomposition
    p_steady = p_record(:, n_min:end);
    P_total = zeros(Nprobes, 1);
    exp_mat = exp(-1i * 2 * pi * f0 .* t_vec(n_min:end));
    
    for ip = 1:Nprobes
        P_total(ip) = sum(p_steady(ip, :) .* exp_mat) * dt;
    end
    
    P_inc_vec = zeros(Nprobes, 1);
    k_fwd = (2*pi*f0) / (c0 * (1 + Mach)); 
    k_bwd = (2*pi*f0) / (c0 * (1 - Mach)); 
    
    % FIXED: Use proper matrix method
    for ip = 2:Nprobes-1
        x1 = probe_x(ip-1);
        x2 = probe_x(ip);
        x3 = probe_x(ip+1);
        
        P1 = P_total(ip-1);
        P2 = P_total(ip);
        P3 = P_total(ip+1);
        
        % Solve for forward and backward waves
        A = [exp(-1i*k_fwd*x1), exp(1i*k_bwd*x1);
             exp(-1i*k_fwd*x3), exp(1i*k_bwd*x3)];
        
        b = [P1; P3];
        
        if cond(A) < 1e10
            X = A\b;
            A_fwd = X(1);
            P_inc_vec(ip) = A_fwd * exp(-1i*k_fwd*x2);
        else
            P_inc_vec(ip) = P2;
        end
    end
    P_inc_vec(1) = P_inc_vec(2);
    P_inc_vec(end) = P_inc_vec(end-1);
    
    % Phase velocity from phase gradient
    phiU = unwrap(angle(P_inc_vec));
    for ip = 1:Nprobes
        i0 = max(1, ip-3);
        i1 = min(Nprobes, ip+3);
        pfit = polyfit(probe_x(i0:i1), phiU(i0:i1), 1);
        if abs(pfit(1)) > 1e-8
            v_loc_results(ip, u_idx) = (2*pi*f0) / abs(pfit(1));
        end
    end
    v_loc_results(:, u_idx) = medfilt1(v_loc_results(:, u_idx), 5);
end

%% Plot Results
figure(2);
set(gcf, 'Position', [100, 100, 1200, 400]);

% Phase Velocity
subplot(1,2,1);
plot(probe_x, v_loc_results(:,1), 'b-', 'LineWidth', 2.5);
hold on;
plot(probe_x, v_loc_results(:,2), 'r-', 'LineWidth', 2.5);
yline(c0, 'k--', 'LineWidth', 2);
xline(abh_start_x, 'm--', 'LineWidth', 2.5);
xline(abh_end_x, 'm--', 'LineWidth', 2.5);
hold off;
xlabel('Position (m)');
ylabel('Phase Velocity (m/s)');
title(sprintf('Phase Velocity with PML & Damping (f=%d Hz, m=%.1f)', f0, m));
legend('U = 0 m/s', 'U = 20 m/s', 'c_0', 'ABH Start', 'ABH End');
grid on;
ylim([100 400]);
xlim([min(probe_x) max(probe_x)]);

% Normalized (β)
subplot(1,2,2);
beta = v_loc_results ./ c0;
plot(probe_x, beta(:,1), 'b-', 'LineWidth', 2.5);
hold on;
plot(probe_x, beta(:,2), 'r-', 'LineWidth', 2.5);
xline(abh_start_x, 'm--', 'LineWidth', 2.5);
xline(abh_end_x, 'm--', 'LineWidth', 2.5);
yline(1, 'k--', 'LineWidth', 2);
hold off;
xlabel('Position (m)');
ylabel('β = v_{phase}/c_0');
title('Normalized Phase Velocity');
legend('U = 0 m/s', 'U = 20 m/s', 'ABH Boundaries', 'β=1');
grid on;
ylim([0.3 1.5]);
xlim([min(probe_x) max(probe_x)]);

% Calculate β statistics in ABH region
abh_indices = (probe_x >= abh_start_x) & (probe_x <= abh_end_x);
min_beta_U0 = min(beta(abh_indices, 1));
min_beta_U20 = min(beta(abh_indices, 2));

fprintf('\n========================================\n');
fprintf('SIMULATION RESULTS WITH PML & DAMPING\n');
fprintf('========================================\n');
fprintf('PML Size: %d cells, Branch Damping: 0.97\n', pml_size);
fprintf('ABH Region: %.3f m to %.3f m\n', abh_start_x, abh_end_x);
fprintf('Frequency: %d Hz, Profile exponent: m = %.1f\n\n', f0, m);
fprintf('Minimum β in ABH region:\n');
fprintf('  U = 0 m/s:  β_min = %.3f (%.1f%% reduction)\n', min_beta_U0, (1-min_beta_U0)*100);
fprintf('  U = 20 m/s: β_min = %.3f (%.1f%% reduction)\n\n', min_beta_U20, (1-min_beta_U20)*100);

if min_beta_U0 < 1 && min_beta_U20 < 1
    fprintf('✅ SLOW SOUND ACHIEVED (β < 1) for both flow conditions!\n');
elseif min_beta_U0 < 1 || min_beta_U20 < 1
    fprintf('⚠ Partial slow sound achieved\n');
else
    fprintf('❌ Slow sound not achieved (β > 1)\n');
end
fprintf('========================================\n');