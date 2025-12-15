% SCRIPT: NON-WSSUS ZONE DETECTOR (SPECTRAL DIVERGENCE)
% Method: Measures the "Information Distance" between channel snapshots.
% Goal: Automatically segment the path into Stationary Zones.

clear;
clc;
close all;

%% 1. Configuration
data_filename = 'Path1_newnew.xlsx'; 
c = 3e8; 

% --- Sensitivity Threshold ---
% Peaks above this value indicate a "Zone Boundary"
% 0.3 is a typical starting point for KL Divergence in wireless channels
divergence_threshold = 0.3; 

% --- Smoothing Window ---
% We smooth the score to remove tiny noise spikes and see the trend
smooth_window = 5; 

%% 2. Load Data
fprintf('Loading data from %s...\n', data_filename);
try
    T = readtable(data_filename);
catch
    error('File not found.');
end

col_names = T.Properties.VariableNames;
mobile_x_col = col_names{1}; freq_col = col_names{2};     
re_col = col_names{3};       im_col = col_names{4};       

mobile_x_list = sort(unique(T.(mobile_x_col)));
num_positions = length(mobile_x_list);
freq_list_hz = sort(unique(T.(freq_col))) * 1e9; 
num_freq_points = length(freq_list_hz);

% Pivot to Matrix
S21_matrix = zeros(num_positions, num_freq_points);
for i = 1:num_positions
    rows = T.(mobile_x_col) == mobile_x_list(i);
    T_subset = sortrows(T(rows, :), freq_col);
    if height(T_subset) == num_freq_points
        S21_matrix(i, :) = T_subset.(re_col) + 1j * T_subset.(im_col);
    end
end

% IFFT
window = hann(num_freq_points)';
impulse_response_matrix = ifft(S21_matrix .* window, [], 2);
num_valid = floor(num_freq_points/2) + 1;
h_matrix = impulse_response_matrix(:, 1:num_valid);

% Power Delay Profile (PDP)
PDP_matrix = abs(h_matrix).^2;

%% 3. Calculate Spectral Divergence (The "Change" Metric)
fprintf('Calculating Spectral Divergence...\n');

% 1. Normalize PDPs to be Probability Distributions (Sum = 1)
% This focuses the metric on the SHAPE of the delay profile, not the power.
PDP_norm = PDP_matrix ./ sum(PDP_matrix, 2);
eps_val = 1e-12; % Small value to prevent log(0)
PDP_norm = PDP_norm + eps_val; 

div_score = zeros(num_positions-1, 1);
pos_axis = mobile_x_list(2:end);

for i = 2:num_positions
    P_curr = PDP_norm(i, :);
    P_prev = PDP_norm(i-1, :);
    
    % Symmetric Kullback-Leibler (KL) Divergence
    % Formula: Sum( (P - Q) * log(P / Q) )
    kl = sum( (P_curr - P_prev) .* log(P_curr ./ P_prev) );
    
    div_score(i-1) = kl;
end

% Smooth the score
smooth_score = movmean(div_score, smooth_window);

%% 4. Detect Zones
[pks, locs] = findpeaks(smooth_score, 'MinPeakHeight', divergence_threshold);
boundary_pos = pos_axis(locs);

fprintf('\n--- DETECTED STATIONARY ZONES ---\n');
fprintf('Based on Divergence Threshold: %.2f\n', divergence_threshold);
fprintf('---------------------------------\n');

zone_start = mobile_x_list(1);
for k = 1:length(boundary_pos)
    zone_end = boundary_pos(k);
    fprintf('Zone %d: %.1fm to %.1fm\n', k, zone_start, zone_end);
    zone_start = zone_end;
end
fprintf('Zone %d: %.1fm to %.1fm (End)\n', length(boundary_pos)+1, zone_start, mobile_x_list(end));
fprintf('---------------------------------\n');

%% 5. Plotting
figure('Name', 'Non-WSSUS Zone Detection', 'Color', 'w', 'Position', [100 100 1000 800]);

% --- Plot 1: PDP Heatmap (Visual Reference) ---
subplot(2,1,1);
% Transpose for Yang Style (X=Position, Y=Delay)
delay_axis = linspace(0, 1/(freq_list_hz(2)-freq_list_hz(1))*1e9, num_freq_points);
delay_axis = delay_axis(1:num_valid);
imagesc(mobile_x_list, delay_axis, 10*log10(PDP_matrix)');
set(gca, 'YDir', 'normal');
colormap('jet'); colorbar;
caxis([-100 -80]); 
ylim([0 500]); 
title('Power Delay Profile (Visual Reference)');
ylabel('Delay (ns)');
% Draw detected boundaries
for k = 1:length(boundary_pos)
    xline(boundary_pos(k), 'w--', 'LineWidth', 2);
end

% --- Plot 2: Divergence Score (The Metric) ---
subplot(2,1,2);
area(pos_axis, smooth_score, 'FaceColor', [0.8 0.3 0.3], 'EdgeColor', 'r'); hold on;
yline(divergence_threshold, 'k--', 'Threshold');
plot(boundary_pos, pks, 'rv', 'MarkerFaceColor', 'k');

title('Spectral Divergence (Instability Index)');
xlabel('Position (m)'); ylabel('Divergence Score');
grid on; xlim([min(mobile_x_list) max(mobile_x_list)]);

% Add Labels
for k = 1:length(boundary_pos)
    text(boundary_pos(k), pks(k)+0.1, sprintf('Change\n@ %.1fm', boundary_pos(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end
