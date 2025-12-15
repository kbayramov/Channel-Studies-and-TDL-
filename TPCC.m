% SCRIPT TO CALCULATE STATIONARITY WITH SPATIAL AVERAGING
% Filters out Fast Fading to reveal true Environmental Stationarity.

clear;
clc;
close all;

%% 1. Configuration
data_filename = 'Path1_newnew.xlsx';
c = 3e8; 

% --- Parameters ---
corr_threshold = 0.8; 

% --- SPATIAL AVERAGING SETTINGS ---
% We average the PDPs over a small distance to remove fast fading.
% 40 wavelengths is a standard averaging window (approx 2 meters at 5.9GHz)
frequency_ghz = 1.815;
wavelength = c / (frequency_ghz * 1e9);
avg_window_meters = 20 * wavelength; % ~1.0 meter window

%% 2. Load Data & Generate Impulse Responses
fprintf('Loading and processing data...\n');
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

% Determine Spatial Step size
spatial_step = mean(diff(mobile_x_list));
window_samples = round(avg_window_meters / spatial_step);
% Ensure window is odd and at least 1
window_samples = max(1, window_samples + mod(window_samples+1, 2));

fprintf('Spatial Averaging Window: %.2f meters (%d samples)\n', ...
    avg_window_meters, window_samples);

% Create Matrix
S21_matrix = zeros(num_positions, num_freq_points);
for i = 1:num_positions
    rows = T.(mobile_x_col) == mobile_x_list(i);
    T_subset = sortrows(T(rows, :), freq_col);
    if height(T_subset) == num_freq_points
        S21_matrix(i, :) = T_subset.(re_col) + 1j * T_subset.(im_col);
    end
end

% IFFT to Time Domain
window = hann(num_freq_points)';
impulse_response_matrix = ifft(S21_matrix .* window, [], 2);
num_valid = floor(num_freq_points/2) + 1;
h_matrix = impulse_response_matrix(:, 1:num_valid);

%% 3. Calculate Power Delay Profiles (Instantaneous)
PDP_matrix = abs(h_matrix).^2;

%% 3.5. SPATIAL AVERAGING (Smoothing)
% We apply a moving average filter across the Position dimension (rows)
% This creates the "Average Power Delay Profile" (APDP) per location.

fprintf('Applying Spatial Averaging to remove fast fading...\n');
APDP_matrix = zeros(size(PDP_matrix));

% Use MATLAB's movmean function for sliding window average
% Dimensions: 1 = along columns (averaging across positions)
APDP_matrix = movmean(PDP_matrix, window_samples, 1);

%% 4. Calculate Temporal PDP Correlation Matrix (TPCC)
% NOTE: We now use APDP (Averaged) instead of raw PDP
fprintf('Calculating Correlation Matrix on Averaged Data...\n');

pdp_norms = sqrt(sum(APDP_matrix.^2, 2));
dot_products = APDP_matrix * APDP_matrix';
norm_cross_product = pdp_norms * pdp_norms';

TPCC_matrix = dot_products ./ norm_cross_product;
TPCC_matrix(logical(eye(size(TPCC_matrix)))) = 1;

%% 5. Calculate Stationarity Distance
stationary_dist_list = zeros(num_positions, 1);

for i = 1:num_positions
    current_corr_slice = TPCC_matrix(i, i:end);
    drop_idx = find(current_corr_slice < corr_threshold, 1);
    
    if isempty(drop_idx)
        dist = mobile_x_list(end) - mobile_x_list(i);
    else
        global_idx = i + drop_idx - 1;
        dist = mobile_x_list(global_idx) - mobile_x_list(i);
    end
    stationary_dist_list(i) = dist;
end

avg_stat_dist = mean(stationary_dist_list);
fprintf('New Average Stationarity Distance: %.2f meters\n', avg_stat_dist);

%% 6. Plotting
figure('Name', 'Spatially Averaged Stationarity', 'Color', 'w', 'Position', [100 100 1000 800]);

% --- Plot 1: Correlation Matrix (Smoothed) ---
subplot(2,2, [1 2]);
imagesc(mobile_x_list, mobile_x_list, TPCC_matrix);
colorbar;
colormap('jet');
caxis([0 1]);
axis xy;
title(['TPCC Matrix (Smoothed over ' num2str(avg_window_meters) 'm)']);
xlabel('Position (m)'); ylabel('Position (m)');
hold on; plot([min(mobile_x_list) max(mobile_x_list)], [min(mobile_x_list) max(mobile_x_list)], 'w--');

% --- Plot 2: Correlation Slice ---
subplot(2,2,3);
indices_to_plot = round(linspace(window_samples, num_positions-window_samples, 3));
hold on;
colors = ['r', 'g', 'b'];
for k = 1:3
    idx = indices_to_plot(k);
    pos = mobile_x_list(idx);
    corr_slice = TPCC_matrix(idx, :);
    plot(mobile_x_list, corr_slice, 'Color', colors(k), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Ref: %.1fm', pos));
end
yline(corr_threshold, 'k--', 'Threshold');
title('Correlation Decay (Smoothed)');
xlabel('Position (m)'); ylabel('Correlation');
legend; grid on; xlim([min(mobile_x_list) max(mobile_x_list)]); ylim([0 1]);

% --- Plot 3: Stationarity Distance ---
subplot(2,2,4);
plot(mobile_x_list, stationary_dist_list, 'k', 'LineWidth', 1.5);
yline(avg_stat_dist, 'r--', sprintf('Mean: %.2fm', avg_stat_dist));
title(['Stationarity Distance (Corr > ' num2str(corr_threshold) ')']);
xlabel('Position (m)'); ylabel('Valid Distance (m)');
grid on; xlim([min(mobile_x_list) max(mobile_x_list)]);
