% =========================================================================
% OPTIMIZED SLOW SOUND: Consistent with Agglomeration Case
% =========================================================================
% KEY SETTINGS (From Paper's Best Case):
% 1. Frequency (f0): 1800 Hz [cite: 239]
% 2. Profile (m): 1.8 (Optimal rate of increase) [cite: 125]
% 3. Probe Position: Near wall (1cm from cavity) to catch strong signals.
% =========================================================================
clear all; close all; clc;

%% -------------------- Parameters --------------------
c0 = 343; rho0 = 1.21; B0 = rho0*c0^2;    
Lx = 1; Ly = 0.5;
dx = 0.001; dy = dx;
Nx = round(Lx/dx); Ny = round(Ly/dy);

% ABH Geometry Setup
R = 0.047; Dmax = 0.07;
d = 0.0035; dc = d/dx; dtt = 0.004;
n_branches = 20;   
L = n_branches*(d+dtt); Lc = L/dx; 
x0 = 0.4; branch_x0 = x0/dx; branch_x1 = branch_x0+Lc; 
abh_start_x = x0;  % ABH start position in meters
abh_end_x = x0 + L; % ABH end position in meters

% *** OPTIMIZED PARAMETERS ***
m = 1.8;   % Paper's optimal m [cite: 125]
f0 = 1800; % Frequency for best slow sound [cite: 239]
% ****************************

duct_yc = round(Ny/2);
duct_y0 = duct_yc - floor(round(2*R/dy)/2);
duct_y1 = duct_yc + floor(round(2*R/dy)/2);

branch_x_positions = round(linspace(branch_x0, branch_x1, n_branches));
xi = linspace(0,1,n_branches);
branch_depths_top = round((d/dy) + (round(Dmax/dy) - d/dy) * xi.^m);
branch_depths_bot = round((d/dy) + (round(Dmax/dy) - d/dy) * xi.^m);

% PML & Source
alpha_base = 5; alpha_branch = 200;
pml_size = round(0.06*Nx); sigma_max = 100000;
src_time = @(t) (min(t/0.002, 1) .* sin(2*pi*f0*t)); % Ramp-up
src_x = round(0.1/dx); src_y = duct_yc;

% Probe Position (1cm from top wall)
probe_y_meters = (R - 0.01); 
probe_y = duct_yc + round(probe_y_meters/dy); 

Nprobes = 100;
probe_x = linspace(0.2, 0.78, Nprobes); % Stop before PML
probe_ix = max(1, min(Nx, round(probe_x/dx))); 

% Create results folder with timestamp
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = sprintf('SlowSoundResults_%s', timestamp);
mkdir(results_folder);
fprintf('Results will be saved in folder: %s\n', results_folder);

%% -------------------- Simulation Loop (0 & 20 m/s) --------------------
U_cases = [0, 20]; 
v_loc_results = zeros(Nprobes, length(U_cases)); 
p_record_all = cell(length(U_cases), 1); % Store pressure records for all cases

for u_idx = 1:length(U_cases)
    U = U_cases(u_idx);
    Mach = U / c0;
    fprintf('\n>>> Running Optimized Slow Sound for U = %d m/s...\n', U);
    
    dt = 0.3 * min(dx,dy)/(c0 + abs(U) + 20); 
    Nt = ceil(0.015/dt); t_vec = (0:Nt-1)*dt;
    n_min = round(0.003 / dt);
    
    p = zeros(Nx, Ny); vx = zeros(Nx+1, Ny); vy = zeros(Nx, Ny+1);
    p_record = zeros(Nprobes, Nt);
    p_snapshot = zeros(Nx, Ny, min(100, Nt)); % Store snapshots for video
    
    % Mask Build
    mask = zeros(Nx, Ny);
    for i = 1:Nx, for j = duct_y0:duct_y1, mask(i,j) = 1; end; end
    for b = 1:n_branches
        ix = branch_x_positions(b); ox0 = max(1, ix - floor(dc/2)); ox1 = min(Nx, ix + floor(dc/2));
        for i = ox0:ox1, for j = duct_y1+1:min(Ny, duct_y1+branch_depths_top(b)), mask(i,j)=1; end; end
        for i = ox0:ox1, for j = max(1, duct_y0-branch_depths_bot(b)):duct_y0-1, mask(i,j)=1; end; end
    end
    
    mask_vx = zeros(Nx+1, Ny); mask_vy = zeros(Nx, Ny+1);
    for i = 1:Nx+1, for j = 1:Ny, if (i-1>=1 && mask(i-1,j)==1 && i<=Nx && mask(i,j)==1), mask_vx(i,j)=1; end; end; end
    for i = 1:Nx, for j = 1:Ny+1, if (j-1>=1 && mask(i,j-1)==1 && j<=Ny && mask(i,j)==1), mask_vy(i,j)=1; end; end; end
    
    sigma_x_p = zeros(Nx,1); 
    for ii = 1:pml_size, s = sigma_max * ((pml_size - (ii-1)) / pml_size)^2; sigma_x_p(ii)=s; sigma_x_p(Nx-ii+1)=s; end
    expP = exp(-repmat(sigma_x_p, 1, Ny)*dt);
    
    alpha_map = alpha_base * ones(Nx, Ny);
    for i=1:Nx, for j=1:Ny, if mask(i,j)==1 && (j<duct_y0 || j>duct_y1), alpha_map(i,j)=alpha_branch; end; end; end
    alpha_vx = zeros(Nx+1,Ny); alpha_vx(2:Nx,:) = 0.5*(alpha_map(1:Nx-1,:)+alpha_map(2:Nx,:));
    alpha_vy = zeros(Nx,Ny+1); alpha_vy(:,2:Ny) = 0.5*(alpha_map(:,1:Ny-1)+alpha_map(:,2:Ny));
    
    % FDTD Update with video snapshots
    snapshot_interval = ceil(Nt/100);
    snapshot_idx = 1;
    
    for n = 1:Nt
        t = (n-1)*dt;
        dpdx = zeros(Nx+1, Ny); dpdx(2:Nx,:) = (p(2:Nx,:) - p(1:Nx-1,:))/dx;
        dvxdx = zeros(Nx+1, Ny); dvxdx(2:Nx,:) = (vx(2:Nx,:) - vx(1:Nx-1,:))/dx;
        vx = vx - (dt/rho0).*dpdx.*mask_vx - dt*U.*dvxdx.*mask_vx - dt.*(alpha_vx.*vx);
        dpdy = zeros(Nx, Ny+1); dpdy(:,2:Ny) = (p(:,2:Ny) - p(:,1:Ny-1))/dy;
        vy = vy - (dt/rho0).*dpdy.*mask_vy - dt.*(alpha_vy.*vy);
        divv = (vx(2:Nx+1,:) - vx(1:Nx,:))/dx + (vy(:,2:Ny+1) - vy(:,1:Ny))/dy;
        p_grad_x = zeros(Nx, Ny); p_grad_x(2:Nx,:) = (p(2:Nx,:) - p(1:Nx-1,:))/dx;
        p = p - (B0*dt).*divv - dt*U.*p_grad_x - dt.*(alpha_map.*p);
        if mask(src_x, src_y)==1, p(src_x, src_y) = p(src_x, src_y) + src_time(t); end
        p = p .* expP; p(mask==0) = 0;
        
        % Record probe data
        for ip = 1:Nprobes, p_record(ip, n) = p(probe_ix(ip), probe_y); end
        
        % Store snapshots for video
        if mod(n, snapshot_interval) == 0 && snapshot_idx <= 100
            p_snapshot(:, :, snapshot_idx) = p;
            snapshot_idx = snapshot_idx + 1;
        end
    end
    
    % Store pressure records
    p_record_all{u_idx} = p_record;
    
    % Create propagation video
    fprintf('Creating propagation video for U = %d m/s...\n', U);
    createPropagationVideo(p_snapshot, mask, dx, dy, abh_start_x, abh_end_x, ...
                          results_folder, U, f0, m, dt, snapshot_interval);
    
    % WAVE DECOMPOSITION (Mach Corrected)
    p_steady = p_record(:, n_min:end);
    P_total = zeros(Nprobes, 1);
    exp_mat = exp(-1i * 2 * pi * f0 .* t_vec(n_min:end));
    for ip = 1:Nprobes, P_total(ip) = sum(p_steady(ip, :) .* exp_mat) * dt; end
    
    P_inc_vec = zeros(Nprobes, 1);
    k_fwd = (2*pi*f0) / (c0 * (1 + Mach)); 
    k_bwd = (2*pi*f0) / (c0 * (1 - Mach)); 
    
    for ip = 2:Nprobes-1
        P_prev = P_total(ip-1); P_next = P_total(ip+1); dist = probe_x(ip+1) - probe_x(ip-1); 
        Denom = exp(-1j*k_fwd*dist) - exp(+1j*k_bwd*dist);
        if abs(Denom) > 1e-8
            A_inc = (P_next - P_prev*exp(+1j*k_bwd*dist)) / Denom;
            P_inc_vec(ip) = A_inc * exp(-1j*k_fwd*(probe_x(ip)-probe_x(ip-1)));
        else
            P_inc_vec(ip) = P_total(ip); 
        end
    end
    P_inc_vec(1) = P_inc_vec(2); P_inc_vec(end) = P_inc_vec(end-1);
    
    % Final Velocity Calculation
    phiU = unwrap(angle(P_inc_vec)); 
    for ip = 1:Nprobes
        i0 = max(1, ip-3); i1 = min(Nprobes, ip+3);
        pfit = polyfit(probe_x(i0:i1), phiU(i0:i1), 1);
        if abs(pfit(1)) > 1e-8, v_loc_results(ip, u_idx) = (2*pi*f0) / -pfit(1); end
    end
    v_loc_results(:, u_idx) = medfilt1(v_loc_results(:, u_idx), 5); % Spike removal
end

%% -------------------- Create Comprehensive Plots --------------------
createAllPlots(probe_x, v_loc_results, p_record_all, t_vec, U_cases, ...
               abh_start_x, abh_end_x, f0, m, results_folder, c0);

%% -------------------- Display Summary --------------------
fprintf('\n========================================\n');
fprintf('SIMULATION COMPLETE\n');
fprintf('Results saved in: %s\n', results_folder);
fprintf('ABH Region: %.3f m to %.3f m\n', abh_start_x, abh_end_x);
fprintf('Optimized parameters: f0 = %d Hz, m = %.1f\n', f0, m);
fprintf('========================================\n');

%% -------------------- Video Creation Function --------------------
function createPropagationVideo(p_snapshot, mask, dx, dy, abh_start_x, abh_end_x, ...
                              folder, U, f0, m, dt, snapshot_interval)
    [Nx, Ny, n_frames] = size(p_snapshot);
    
    % Create video writer
    video_filename = fullfile(folder, sprintf('Propagation_U%d_f%d_m%.1f.mp4', U, f0, m));
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 20;
    v.Quality = 95;
    open(v);
    
    % Create figure for video
    fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 600], 'Color', 'w');
    
    % Pre-calculate global limits for consistent scaling
    p_all = p_snapshot(:);
    clim_value = max(abs(p_all));
    if clim_value < 1e-10
        clim_value = 1;
    end
    
    for frame = 1:n_frames
        p_frame = p_snapshot(:, :, frame);
        
        % Convert to physical coordinates
        x_phys = (0:Nx-1)*dx;
        y_phys = (0:Ny-1)*dy;
        
        % Create subplot layout
        subplot(2,2,[1 3]);
        
        % Safely calculate color limits
        current_max = max(abs(p_frame(:)));
        if current_max < 1e-10
            current_clim = [-1 1];
        else
            current_clim = [-1 1] * current_max;
        end
        
        imagesc(x_phys, y_phys, p_frame');
        hold on;
        
        % Overlay ABH structure (more efficient)
        [mask_y, mask_x] = find(mask' == 1);
        scatter(mask_x*dx, mask_y*dy, 1, 'k', 'filled', 'AlphaData', 0.3);
        
        % Add ABH region lines
        xline(abh_start_x, 'r--', 'LineWidth', 2, 'DisplayName', 'ABH Start');
        xline(abh_end_x, 'r--', 'LineWidth', 2, 'DisplayName', 'ABH End');
        
        hold off;
        axis equal tight;
        xlabel('x (m)', 'FontSize', 10);
        ylabel('y (m)', 'FontSize', 10);
        title(sprintf('Pressure Field (U=%d m/s, Frame %d/%d)', U, frame, n_frames), ...
              'FontSize', 11, 'FontWeight', 'bold');
        colorbar;
        
        % Use caxis instead of clim for compatibility
        caxis(current_clim);
        colormap(jet);
        
        % Middle horizontal line plot
        subplot(2,2,2);
        mid_y = round(Ny/2);
        plot(x_phys, p_frame(:, mid_y), 'b-', 'LineWidth', 1.5);
        hold on;
        xline(abh_start_x, 'r--', 'LineWidth', 2);
        xline(abh_end_x, 'r--', 'LineWidth', 2);
        hold off;
        xlabel('x (m)', 'FontSize', 10);
        ylabel('Pressure', 'FontSize', 10);
        title('Pressure along centerline', 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([min(x_phys) max(x_phys)]);
        
        % Probe line plot (near top wall)
        subplot(2,2,4);
        probe_y_idx = round(mid_y + (0.037/dy)); % 3.7cm from center (near top)
        plot(x_phys, p_frame(:, probe_y_idx), 'r-', 'LineWidth', 1.5);
        hold on;
        xline(abh_start_x, 'r--', 'LineWidth', 2);
        xline(abh_end_x, 'r--', 'LineWidth', 2);
        hold off;
        xlabel('x (m)', 'FontSize', 10);
        ylabel('Pressure', 'FontSize', 10);
        title(sprintf('Pressure at y=%.3f m', probe_y_idx*dy), ...
              'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([min(x_phys) max(x_phys)]);
        
        % Add overall title
        sgtitle(sprintf('Sound Propagation: U=%d m/s, f=%d Hz, m=%.1f | Time: %.3f ms', ...
                        U, f0, m, frame*snapshot_interval*dt*1000), ...
                'FontSize', 13, 'FontWeight', 'bold');
        
        % Capture frame
        frame_img = getframe(fig);
        writeVideo(v, frame_img);
    end
    
    close(v);
    close(fig);
    fprintf('  Video saved: %s\n', video_filename);
end

%% -------------------- All Plots Function --------------------
function createAllPlots(probe_x, v_loc_results, p_record_all, t_vec, U_cases, ...
                       abh_start_x, abh_end_x, f0, m, folder, c0)
    
    % Get number of probes
    Nprobes = length(probe_x);
    
    % ========== PLOT 1: Phase Velocity Comparison ==========
    fig1 = figure('Position', [100, 100, 1400, 600], 'Color', 'w');
    
    subplot(1,2,1);
    plot(probe_x, v_loc_results(:,1), 'b-o', 'LineWidth', 2, ...
         'DisplayName', 'U = 0 m/s', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); 
    hold on;
    plot(probe_x, v_loc_results(:,2), 'r-s', 'LineWidth', 2, ...
         'DisplayName', 'U = 20 m/s', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    yline(c0, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('c_0 = %.0f m/s', c0));
    xline(abh_start_x, 'm--', 'LineWidth', 2.5, 'DisplayName', 'ABH Start');
    xline(abh_end_x, 'm--', 'LineWidth', 2.5, 'DisplayName', 'ABH End');
    
    % Highlight slow sound region
    abh_indices = (probe_x >= abh_start_x) & (probe_x <= abh_end_x);
    if any(abh_indices)
        area_x = probe_x(abh_indices);
        area_y1 = v_loc_results(abh_indices, 1);
        area_y2 = v_loc_results(abh_indices, 2);
        
        % Create polygon for ABH region
        fill_x = [area_x, fliplr(area_x)];
        fill_y = [area_y1', fliplr(area_y2')];
        fill(fill_x, fill_y, 'g', 'FaceAlpha', 0.15, ...
             'EdgeColor', 'g', 'LineWidth', 1, 'DisplayName', 'ABH Region');
    end
    
    hold off;
    xlabel('Position (m)', 'FontSize', 12);
    ylabel('Phase Velocity (m/s)', 'FontSize', 12);
    title(sprintf('Phase Velocity: Optimized Slow Sound\nf=%d Hz, m=%.1f', f0, m), ...
          'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on; 
    ylim([0 500]); 
    xlim([min(probe_x) max(probe_x)]);
    
    % Add text annotations for slow sound effect
    min_vel_0 = min(v_loc_results(:,1));
    min_vel_20 = min(v_loc_results(:,2));
    reduction_0 = (1 - min_vel_0/c0) * 100;
    reduction_20 = (1 - min_vel_20/c0) * 100;
    
    text(0.5, 450, sprintf('U=0 m/s:\nMin: %.1f m/s\nReduction: %.1f%%', ...
         min_vel_0, reduction_0), ...
         'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'b', 'LineWidth', 1);
    text(0.5, 350, sprintf('U=20 m/s:\nMin: %.1f m/s\nReduction: %.1f%%', ...
         min_vel_20, reduction_20), ...
         'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'r', 'LineWidth', 1);
    
    % ========== PLOT 2: Velocity Reduction Ratio ==========
    subplot(1,2,2);
    velocity_ratio = v_loc_results ./ c0;
    plot(probe_x, velocity_ratio(:,1), 'b-', 'LineWidth', 2); 
    hold on;
    plot(probe_x, velocity_ratio(:,2), 'r-', 'LineWidth', 2);
    xline(abh_start_x, 'm--', 'LineWidth', 2.5);
    xline(abh_end_x, 'm--', 'LineWidth', 2.5);
    yline(1, 'k--', 'LineWidth', 1.5);
    
    % Highlight region where velocity < c0/2 (strong slow sound)
    slow_sound_threshold = 0.5;
    slow_indices = velocity_ratio(:,1) < slow_sound_threshold;
    if any(slow_indices)
        plot(probe_x(slow_indices), velocity_ratio(slow_indices,1), 'b-', ...
             'LineWidth', 3, 'DisplayName', 'v < c_0/2');
    end
    
    hold off;
    xlabel('Position (m)', 'FontSize', 12);
    ylabel('v_{phase} / c_0', 'FontSize', 12);
    title('Normalized Phase Velocity', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'U = 0 m/s', 'U = 20 m/s', 'ABH Boundaries', 'Reference (1.0)', 'Strong Slow Sound'}, ...
           'Location', 'best', 'FontSize', 9);
    grid on; 
    ylim([0 1.5]); 
    xlim([min(probe_x) max(probe_x)]);
    
    saveas(fig1, fullfile(folder, 'PhaseVelocity_Comparison.png'));
    saveas(fig1, fullfile(folder, 'PhaseVelocity_Comparison.fig'));
    
    % ========== PLOT 3: Time Signals at Different Positions ==========
    fig2 = figure('Position', [100, 100, 1200, 800], 'Color', 'w');
    
    probe_positions = [0.25, 0.45, 0.55, 0.75]; % Before, start, middle, after ABH
    [~, probe_indices] = min(abs(probe_x' - probe_positions), [], 1);
    
    for u_idx = 1:length(U_cases)
        subplot(2, length(U_cases), u_idx);
        colors = lines(length(probe_positions));
        for i = 1:length(probe_positions)
            idx = probe_indices(i);
            plot(t_vec*1000, p_record_all{u_idx}(idx,:), 'LineWidth', 1.5, ...
                 'Color', colors(i,:), 'DisplayName', sprintf('x=%.2f m', probe_x(idx)));
            hold on;
        end
        hold off;
        xlabel('Time (ms)', 'FontSize', 11);
        ylabel('Pressure', 'FontSize', 11);
        title(sprintf('Time Signals (U=%d m/s)', U_cases(u_idx)), ...
              'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        grid on;
        xlim([0 max(t_vec)*1000]);
        
        % Add ABH region indicator
        text(0.02, 0.95, sprintf('ABH: %.2f-%.2f m', abh_start_x, abh_end_x), ...
             'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', ...
             'EdgeColor', 'k');
    end
    
    % ========== PLOT 4: Spectrogram of signals ==========
    for u_idx = 1:length(U_cases)
        subplot(2, length(U_cases), u_idx+length(U_cases));
        idx_middle = find(probe_x >= (abh_start_x + abh_end_x)/2, 1);
        if isempty(idx_middle)
            idx_middle = round(length(probe_x)/2);
        end
        
        Fs = 1/(t_vec(2)-t_vec(1)); % Sampling frequency
        [S,F,T] = spectrogram(p_record_all{u_idx}(idx_middle,:), 256, 250, 256, Fs);
        imagesc(T*1000, F/1000, 20*log10(abs(S)));
        xlabel('Time (ms)', 'FontSize', 11);
        ylabel('Frequency (kHz)', 'FontSize', 11);
        title(sprintf('Spectrogram at x=%.2f m (U=%d m/s)', ...
              probe_x(idx_middle), U_cases(u_idx)), 'FontSize', 12, 'FontWeight', 'bold');
        colorbar;
        ylim([0 5]);
        colormap('jet');
        hold on;
        % Mark the source frequency
        plot(T*1000, f0/1000 * ones(size(T)), 'w--', 'LineWidth', 1.5);
        hold off;
    end
    
    sgtitle(sprintf('Time Domain Analysis - f=%d Hz, m=%.1f', f0, m), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    saveas(fig2, fullfile(folder, 'TimeSignals_Spectrograms.png'));
    saveas(fig2, fullfile(folder, 'TimeSignals_Spectrograms.fig'));
    
    % ========== PLOT 5: Spatial Amplitude Distribution ==========
    fig3 = figure('Position', [100, 100, 1000, 800], 'Color', 'w');
    
    % Calculate amplitude for each probe position
    amplitudes = zeros(Nprobes, length(U_cases));
    for u_idx = 1:length(U_cases)
        for ip = 1:Nprobes
            % Use last half of the signal for steady-state amplitude
            start_idx = round(size(p_record_all{u_idx}, 2)/2);
            amplitudes(ip, u_idx) = rms(p_record_all{u_idx}(ip, start_idx:end));
        end
    end
    
    subplot(2,2,1);
    plot(probe_x, amplitudes(:,1), 'b-', 'LineWidth', 2); 
    hold on;
    plot(probe_x, amplitudes(:,2), 'r-', 'LineWidth', 2);
    xline(abh_start_x, 'm--', 'LineWidth', 2.5);
    xline(abh_end_x, 'm--', 'LineWidth', 2.5);
    hold off;
    xlabel('Position (m)', 'FontSize', 12);
    ylabel('Amplitude (RMS)', 'FontSize', 12);
    title('Spatial Amplitude Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'U=0 m/s', 'U=20 m/s', 'ABH Boundaries'}, 'Location', 'best', 'FontSize', 10);
    grid on;
    xlim([min(probe_x) max(probe_x)]);
    
    % ========== PLOT 6: Phase Distribution ==========
    subplot(2,2,2);
    colors = {'b-', 'r-'};
    for u_idx = 1:length(U_cases)
        % Extract phase at source frequency
        start_idx = round(size(p_record_all{u_idx}, 2)/2);
        p_steady = p_record_all{u_idx}(:, start_idx:end);
        phase = zeros(Nprobes, 1);
        t_steady = t_vec(start_idx:end) - t_vec(start_idx);
        exp_mat = exp(-1i * 2 * pi * f0 .* t_steady);
        
        for ip = 1:Nprobes
            phase(ip) = angle(sum(p_steady(ip, :) .* exp_mat) * (t_vec(2)-t_vec(1)));
        end
        
        plot(probe_x, unwrap(phase), colors{u_idx}, 'LineWidth', 2); 
        hold on;
    end
    xline(abh_start_x, 'm--', 'LineWidth', 2.5);
    xline(abh_end_x, 'm--', 'LineWidth', 2.5);
    hold off;
    xlabel('Position (m)', 'FontSize', 12);
    ylabel('Phase (rad)', 'FontSize', 12);
    title('Phase Distribution along Probe Line', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'U=0 m/s', 'U=20 m/s', 'ABH Boundaries'}, 'Location', 'best', 'FontSize', 10);
    grid on;
    xlim([min(probe_x) max(probe_x)]);
    
    % ========== PLOT 7: Velocity Reduction Statistics ==========
    subplot(2,2,3);
    % Calculate reduction percentage
    velocity_reduction = (1 - v_loc_results/c0) * 100;
    
    % Create bar plot
    h1 = bar(1:Nprobes, velocity_reduction(:,1), 'FaceColor', 'b', 'FaceAlpha', 0.6);
    hold on;
    h2 = bar(1:Nprobes, velocity_reduction(:,2), 'FaceColor', 'r', 'FaceAlpha', 0.6);
    
    % Mark ABH region
    abh_start_idx = find(probe_x >= abh_start_x, 1);
    abh_end_idx = find(probe_x <= abh_end_x, 1, 'last');
    if ~isempty(abh_start_idx) && ~isempty(abh_end_idx)
        x_patch = [abh_start_idx, abh_end_idx, abh_end_idx, abh_start_idx];
        y_patch = [0, 0, 100, 100];
        patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'g', 'LineWidth', 1);
    end
    
    hold off;
    xlabel('Probe Index', 'FontSize', 12);
    ylabel('Velocity Reduction (%)', 'FontSize', 12);
    title('Velocity Reduction Percentage', 'FontSize', 14, 'FontWeight', 'bold');
    legend([h1, h2], {'U=0 m/s', 'U=20 m/s'}, 'Location', 'best', 'FontSize', 10);
    grid on;
    ylim([0 100]);
    xlim([1 Nprobes]);
    
    % ========== PLOT 8: Summary Statistics ==========
    subplot(2,2,4);
    stats_matrix = [
        min(v_loc_results(:,1)), mean(v_loc_results(:,1)), max(v_loc_results(:,1));
        min(v_loc_results(:,2)), mean(v_loc_results(:,2)), max(v_loc_results(:,2))
    ];
    
    hbar = bar(stats_matrix');
    set(hbar(1), 'FaceColor', 'b', 'FaceAlpha', 0.7);
    set(hbar(2), 'FaceColor', 'r', 'FaceAlpha', 0.7);
    
    ylabel('Velocity (m/s)', 'FontSize', 12);
    title('Velocity Statistics', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'U=0 m/s', 'U=20 m/s'}, 'Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'XTickLabel', {'Minimum', 'Mean', 'Maximum'});
    ylim([0 400]);
    
    % Add value labels on bars
    for i = 1:3
        for j = 1:2
            text(i + (j-1.5)*0.2, stats_matrix(j,i)*0.95, ...
                 sprintf('%.1f', stats_matrix(j,i)), ...
                 'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
                 'FontSize', 10, 'Color', 'w');
        end
    end
    
    sgtitle(sprintf('Comprehensive Slow Sound Analysis - f=%d Hz, m=%.1f', f0, m), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(fig3, fullfile(folder, 'Comprehensive_Analysis.png'));
    saveas(fig3, fullfile(folder, 'Comprehensive_Analysis.fig'));
    
    % ========== PLOT 9: Combined Final Summary ==========
    fig4 = figure('Position', [100, 100, 1200, 400], 'Color', 'w');
    
    subplot(1,3,1);
    imagesc(probe_x, 1:size(p_record_all{1},2), p_record_all{1}');
    xlabel('Position (m)', 'FontSize', 11);
    ylabel('Time Sample', 'FontSize', 11);
    title('Space-Time Plot (U=0 m/s)', 'FontSize', 12, 'FontWeight', 'bold');
    colorbar;
    hold on;
    xline(abh_start_x, 'w--', 'LineWidth', 2.5);
    xline(abh_end_x, 'w--', 'LineWidth', 2.5);
    hold off;
    colormap('jet');
    
    subplot(1,3,2);
    imagesc(probe_x, 1:size(p_record_all{2},2), p_record_all{2}');
    xlabel('Position (m)', 'FontSize', 11);
    ylabel('Time Sample', 'FontSize', 11);
    title('Space-Time Plot (U=20 m/s)', 'FontSize', 12, 'FontWeight', 'bold');
    colorbar;
    hold on;
    xline(abh_start_x, 'w--', 'LineWidth', 2.5);
    xline(abh_end_x, 'w--', 'LineWidth', 2.5);
    hold off;
    colormap('jet');
    
    subplot(1,3,3);
    plot(probe_x, v_loc_results(:,1), 'b-', 'LineWidth', 2); 
    hold on;
    plot(probe_x, v_loc_results(:,2), 'r-', 'LineWidth', 2);
    
    % Fill ABH region
    if any(abh_indices)
        fill_x = probe_x(abh_indices);
        fill([fill_x, fliplr(fill_x)], ...
             [zeros(1,sum(abh_indices)), 500*ones(1,sum(abh_indices))], ...
             'g', 'FaceAlpha', 0.1, 'EdgeColor', 'g', 'LineWidth', 1);
    end
    
    plot([abh_start_x, abh_end_x], [c0, c0], 'k--', 'LineWidth', 1.5);
    hold off;
    xlabel('Position (m)', 'FontSize', 11);
    ylabel('Phase Velocity (m/s)', 'FontSize', 11);
    title('Slow Sound Effect Summary', 'FontSize', 12, 'FontWeight', 'bold');
    legend({'U=0 m/s', 'U=20 m/s', 'ABH Region', 'c_0'}, 'Location', 'best', 'FontSize', 9);
    grid on;
    ylim([0 500]);
    xlim([min(probe_x) max(probe_x)]);
    
    sgtitle(sprintf('Optimized Slow Sound Analysis: f=%d Hz, m=%.1f', f0, m), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(fig4, fullfile(folder, 'Summary_Plot.png'));
    saveas(fig4, fullfile(folder, 'Summary_Plot.fig'));
    
    % ========== PLOT 10: ABH Geometry Visualization ==========
    fig5 = figure('Position', [100, 100, 800, 400], 'Color', 'w');
    
    % Create a simple visualization of the ABH geometry
    x_abh = linspace(abh_start_x, abh_end_x, 100);
    depth_profile = (Dmax - d) * ((x_abh - abh_start_x)/(abh_end_x - abh_start_x)).^m + d;
    
    subplot(1,2,1);
    plot(x_abh, depth_profile*1000, 'b-', 'LineWidth', 3);
    hold on;
    plot(x_abh, -depth_profile*1000, 'b-', 'LineWidth', 3);
    plot([abh_start_x, abh_end_x], [0, 0], 'k-', 'LineWidth', 2);
    fill([x_abh, fliplr(x_abh)], [depth_profile*1000, -fliplr(depth_profile)*1000], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
    xline(abh_start_x, 'r--', 'LineWidth', 2);
    xline(abh_end_x, 'r--', 'LineWidth', 2);
    hold off;
    xlabel('Position (m)', 'FontSize', 12);
    ylabel('Depth (mm)', 'FontSize', 12);
    title(sprintf('ABH Depth Profile (m=%.1f)', m), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    legend('Profile', '', 'Duct Center', 'ABH Region', 'Location', 'best');
    
    subplot(1,2,2);
    % Show velocity reduction in ABH region
    abh_vel_0 = v_loc_results(abh_indices, 1);
    abh_vel_20 = v_loc_results(abh_indices, 2);
    abh_x = probe_x(abh_indices);
    
    plot(abh_x, abh_vel_0, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    hold on;
    plot(abh_x, abh_vel_20, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    yline(c0, 'k--', 'LineWidth', 1.5);
    hold off;
    xlabel('Position in ABH (m)', 'FontSize', 12);
    ylabel('Phase Velocity (m/s)', 'FontSize', 12);
    title('Slow Sound in ABH Region', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'U=0 m/s', 'U=20 m/s', 'c_0'}, 'Location', 'best');
    grid on;
    ylim([0 400]);
    
    saveas(fig5, fullfile(folder, 'ABH_Geometry_Velocity.png'));
    saveas(fig5, fullfile(folder, 'ABH_Geometry_Velocity.fig'));
    
    % ========== Save data for further analysis ==========
    save(fullfile(folder, 'simulation_data.mat'), ...
         'probe_x', 'v_loc_results', 'p_record_all', 't_vec', 'U_cases', ...
         'abh_start_x', 'abh_end_x', 'f0', 'm', 'c0');
    
    % ========== Create a text summary file ==========
    summary_file = fullfile(folder, 'simulation_summary.txt');
    fid = fopen(summary_file, 'w');
    fprintf(fid, '========================================\n');
    fprintf(fid, 'SLOW SOUND SIMULATION SUMMARY\n');
    fprintf(fid, '========================================\n\n');
    fprintf(fid, 'Simulation Parameters:\n');
    fprintf(fid, '  Frequency (f0): %d Hz\n', f0);
    fprintf(fid, '  ABH profile exponent (m): %.1f\n', m);
    fprintf(fid, '  ABH Region: %.3f m to %.3f m\n', abh_start_x, abh_end_x);
    fprintf(fid, '  Sound speed (c0): %.1f m/s\n\n', c0);
    
    fprintf(fid, 'Velocity Results:\n');
    fprintf(fid, '  U = 0 m/s:\n');
    fprintf(fid, '    Minimum velocity: %.1f m/s (%.1f%% reduction)\n', ...
            min(v_loc_results(:,1)), (1-min(v_loc_results(:,1))/c0)*100);
    fprintf(fid, '    Mean velocity: %.1f m/s\n', mean(v_loc_results(:,1)));
    fprintf(fid, '    Maximum velocity: %.1f m/s\n\n', max(v_loc_results(:,1)));
    
    fprintf(fid, '  U = 20 m/s:\n');
    fprintf(fid, '    Minimum velocity: %.1f m/s (%.1f%% reduction)\n', ...
            min(v_loc_results(:,2)), (1-min(v_loc_results(:,2))/c0)*100);
    fprintf(fid, '    Mean velocity: %.1f m/s\n', mean(v_loc_results(:,2)));
    fprintf(fid, '    Maximum velocity: %.1f m/s\n\n', max(v_loc_results(:,2)));
    
    fprintf(fid, 'Files Generated:\n');
    fprintf(fid, '  1. Propagation videos (MP4 format)\n');
    fprintf(fid, '  2. Phase velocity comparison plots\n');
    fprintf(fid, '  3. Time signals and spectrograms\n');
    fprintf(fid, '  4. Comprehensive analysis plots\n');
    fprintf(fid, '  5. Summary plots\n');
    fprintf(fid, '  6. ABH geometry visualization\n');
    fprintf(fid, '  7. Data files (MAT format)\n');
    fprintf(fid, '  8. This summary file\n');
    fclose(fid);
    
    fprintf('\nAll plots saved in folder: %s\n', folder);
    fprintf('Generated %d plots and videos successfully!\n', 5); % 5 figures
end