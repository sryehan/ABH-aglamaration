% =========================================================================
% SLOW SOUND ESTIMATION: Wave Decomposition Method (FIXED OSCILLATIONS)
% =========================================================================
% SOLUTION FOR POINT 1:
% 1. Uses Wave Decomposition to remove Reflections (Fixes Oscillations).
% 2. Calculates Phase Velocity from Incident Wave only.
% 3. Compares U=0 and U=20 m/s.
% 4. Cuts Transient Data (T < 3ms).
% =========================================================================
clear all; close all; clc;

%% -------------------- Parameters --------------------
c0 = 343;          
rho0 = 1.21;       
B0 = rho0*c0^2;    

Lx = 1;
Ly = 0.5;
dx = 0.001; dy = dx;
Nx = round(Lx/dx); Ny = round(Ly/dy);

% ABH Geometry
R = 0.047; Dmax = 0.07;
d = 0.0035; dc = d/dx; dtt = 0.004;
n_branches = 20;   
L = n_branches*(d+dtt); Lc = L/dx; 
x0 = 0.4; branch_x0 = x0/dx; branch_x1 = branch_x0+Lc; 

% === OPTION: Set to 0 for "Straight Duct Test" as per Professor ===
enable_ABH = 1;    
% ==================================================================

% Geometry Construction
duct_height = round(2*R/dy);
duct_yc = round(Ny/2);
duct_y0 = duct_yc - floor(duct_height/2);
duct_y1 = duct_yc + floor(duct_height/2);

branch_x_positions = round(linspace(branch_x0, branch_x1, n_branches));
branch_opening_width = dc;
max_branch_depth_top = round(Dmax/dy); max_branch_depth_bot = round(Dmax/dy);
m = 1/2.3;

xi = linspace(0,1,n_branches);
branch_depths_top = round((d/dy) + (max_branch_depth_top - d/dy) * xi.^m);
branch_depths_bot = round((d/dy) + (max_branch_depth_bot - d/dy) * xi.^m);

% Damping & PML
alpha_base = 5; 
alpha_branch = 200;
pml_size = round(0.06*Nx); %pml try from change spong to decrease sigma max 
sigma_max = 100000;

% Source
f0 = 1200;            
src_time =  @(t) (sin(2*pi*f0*t));
src_x = round(0.1/dx); 

src_y = round((duct_y0 + duct_y1)/2);

% Probes (100 Microphones)
Nprobes = 100;
probe_x = linspace(0.2, 0.8, Nprobes);  %0.7818 probes should be off here or increase the length L to 1
probe_ix = max(1, min(Nx, round(probe_x/dx))); 
probe_y = src_y; %max(1, min(Nx, round(probe_x/dx))); change here 0.01 distance from R

%% -------------------- Simulation Loop (0 & 20 m/s) --------------------
U_cases = [0, 20]; 
v_loc_results = zeros(Nprobes, length(U_cases)); 

fprintf('Starting Improved Slow Sound Analysis (Wave Decomposition)...\n');

for u_idx = 1:length(U_cases)
    U = U_cases(u_idx);
    fprintf('\n>>> Case %d/2: Running Simulation for U = %d m/s...\n', u_idx, U);
    
    % --- Stability & Time ---
    CFL = 0.3;
    dt = CFL * min(dx,dy)/(c0 + abs(U) + 20); 
    
    Tmax = 0.015;         
    Nt = ceil(Tmax/dt);
    t_vec = (0:Nt-1)*dt;
    
    Tmin = 0.003;         
    n_min = round(Tmin / dt);
    
    % Fields Init
    p = zeros(Nx, Ny); vx = zeros(Nx+1, Ny); vy = zeros(Nx, Ny+1);
    p_record = zeros(Nprobes, Nt);
    
    % Masks & Damping Init
    mask = zeros(Nx, Ny);
    for i = 1:Nx, for j = duct_y0:duct_y1, mask(i,j) = 1; end; end
    if enable_ABH
        for b = 1:n_branches
            ix = branch_x_positions(b);
            half_open = floor(branch_opening_width/2);
            ox0 = max(1, ix - half_open); ox1 = min(Nx, ix + half_open);
            depth_top = branch_depths_top(b);
            for i = ox0:ox1, for j = duct_y1+1 : min(Ny, duct_y1 + depth_top), mask(i,j) = 1; end; end
            depth_bot = branch_depths_bot(b);
            for i = ox0:ox1, for j = max(1, duct_y0 - depth_bot) : duct_y0-1, mask(i,j) = 1; end; end
        end
    end
    
    mask_vx = zeros(Nx+1, Ny); mask_vy = zeros(Nx, Ny+1);
    for i = 1:Nx+1, for j = 1:Ny, if (i-1>=1 && mask(i-1,j)==1 && i<=Nx && mask(i,j)==1), mask_vx(i,j)=1; end; end; end
    for i = 1:Nx, for j = 1:Ny+1, if (j-1>=1 && mask(i,j-1)==1 && j<=Ny && mask(i,j)==1), mask_vy(i,j)=1; end; end; end
    
    sigma_x_p = zeros(Nx,1); 
    for ii = 1:pml_size, s = sigma_max * ((pml_size - (ii-1)) / pml_size)^2; sigma_x_p(ii)=s; sigma_x_p(Nx-ii+1)=s; end
    SigmaP = repmat(sigma_x_p, 1, Ny); expP = exp(-SigmaP*dt);
    
    alpha_map = alpha_base * ones(Nx, Ny);
    if enable_ABH
        for i=1:Nx, for j=1:Ny, if mask(i,j)==1 && (j<duct_y0 || j>duct_y1), alpha_map(i,j)=alpha_branch; end; end; end
    end
    alpha_map(src_x, src_y) = alpha_base;
    alpha_vx = zeros(Nx+1,Ny); alpha_vx(2:Nx,:) = 0.5*(alpha_map(1:Nx-1,:)+alpha_map(2:Nx,:));
    alpha_vy = zeros(Nx,Ny+1); alpha_vy(:,2:Ny) = 0.5*(alpha_map(:,1:Ny-1)+alpha_map(:,2:Ny));
    
    invRho_dt = dt/rho0; Bdt = B0*dt;

    % --- Time Stepping ---
    for n = 1:Nt
        t = (n-1)*dt;
        
        % VX Update
        dpdx = zeros(Nx+1, Ny); dpdx(2:Nx, :) = (p(2:Nx,:) - p(1:Nx-1,:)) / dx;
        dvxdx_c = zeros(Nx+1, Ny); dvxdx_c(2:Nx, :) = (vx(2:Nx, :) - vx(1:Nx-1, :)) / dx;
        
        vx = vx - invRho_dt .* dpdx .* mask_vx ...
                - dt * U .* dvxdx_c .* mask_vx ...
                - dt .* (alpha_vx .* vx);
        
        % VY Update
        dpdy = zeros(Nx, Ny+1); dpdy(:, 2:Ny) = (p(:,2:Ny) - p(:,1:Ny-1)) / dy;
        vy = vy - invRho_dt .* dpdy .* mask_vy - dt .* (alpha_vy .* vy);
        
        % Pressure Update
        divv = (vx(2:Nx+1,:) - vx(1:Nx,:))/dx + (vy(:,2:Ny+1) - vy(:,1:Ny))/dy;
        dpdx_c = zeros(Nx, Ny); dpdx_c(2:Nx, :) = (p(2:Nx, :) - p(1:Nx-1, :)) / dx;
        
        p = p - Bdt .* divv - dt * U .* dpdx_c - dt .* (alpha_map .* p);
        if mask(src_x, src_y)==1, p(src_x, src_y) = p(src_x, src_y) + src_time(t); end
        
        % Record Probes
        for ip = 1:Nprobes
            p_record(ip, n) = p(probe_ix(ip), probe_y);
        end
        
        if mod(n, 2000) == 0, fprintf('  Progress: %.0f%%\n', (n/Nt)*100); end
    end
    
    % --- Post-Processing: WAVE DECOMPOSITION ---
    fprintf('  Processing with Wave Separation (Removing Reflections)...\n');
    t_steady = t_vec(n_min:end);
    p_steady = p_record(:, n_min:end);
    
    % 1. Get Total Complex Pressure (P_total) at 1200 Hz
    P_total = zeros(Nprobes, 1);
    exp_mat = exp(-1i * 2 * pi * f0 .* t_steady);
    for ip = 1:Nprobes
        P_total(ip) = sum(p_steady(ip, :) .* exp_mat) * dt;
    end
    
    % 2. Separate Incident Wave (P_inc) using Neighbor Probes
    P_inc_vec = zeros(Nprobes, 1);
    omega = 2 * pi * f0;
    k0 = omega / c0; 
    
    for ip = 2:Nprobes-1
        P_prev = P_total(ip-1); 
        P_next = P_total(ip+1); 
        dist = probe_x(ip+1) - probe_x(ip-1); 
        
        % Decomposition Formula
        H_fwd = exp(-1j * k0 * dist);
        H_bwd = exp(+1j * k0 * dist);
        Denom = H_fwd - H_bwd;
        
        if abs(Denom) > 1e-10
            A_inc_prev = (P_next - P_prev * H_bwd) / Denom;
            % Propagate to current probe
            dist_to_center = probe_x(ip) - probe_x(ip-1);
            P_inc_vec(ip) = A_inc_prev * exp(-1j * k0 * dist_to_center);
        else
            P_inc_vec(ip) = P_total(ip); 
        end
    end
    % Handle edges
    P_inc_vec(1) = P_inc_vec(2);
    P_inc_vec(end) = P_inc_vec(end-1);
    
    % 3. Calculate Phase Speed from INCIDENT Pressure
    phi = angle(P_inc_vec); 
    phiU = unwrap(phi); 
    W = 6; halfW = floor(W/2);
    
    for ip = 1:Nprobes
        i0 = max(1, ip-halfW); i1 = min(Nprobes, ip+halfW);
        xW = probe_x(i0:i1); yW = phiU(i0:i1);
        if length(xW) > 1
            pfit = polyfit(xW, yW, 1);
            k_loc = -pfit(1);
            if abs(k_loc) > 1e-8
                v_loc_results(ip, u_idx) = (2*pi*f0) / k_loc;
            end
        else
            v_loc_results(ip, u_idx) = NaN;
        end
    end
end

%% -------------------- Final Plotting --------------------
figure('Name', 'Slow Sound Comparison (Decomposed)', 'Color', 'w', 'Position', [100, 100, 900, 600]);

plot(probe_x, v_loc_results(:,1), 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'U = 0 m/s (No Flow)'); 
hold on;
plot(probe_x, v_loc_results(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'U = 20 m/s (With Flow)');

yline(c0, 'k--', 'Air Speed (343 m/s)', 'LineWidth', 1.5, 'DisplayName', 'Reference');
xline(0.4, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(0.4+L, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(0.45, 550, 'ABH REGION', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('Position along Duct (m)', 'FontSize', 12);
ylabel('Phase Velocity (m/s)', 'FontSize', 12);
title(sprintf('Local Phase Velocity (Incident Wave Only) at f = %d Hz', f0), 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 11);
grid on;
ylim([0 600]); 
xlim([0.2 0.8]);

fprintf('\nAnalysis Complete. Oscillations should be reduced.\n');