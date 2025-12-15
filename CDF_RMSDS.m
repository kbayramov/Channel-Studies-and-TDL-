% SCRIPT TO GENERATE CDF OF RMS-DS (Comparing All Taps vs Active Taps)
% Based on methodology from Ansari et al., VTC2024-Fall
% Features: Loops all positions, calculates RMS-DS, plots CDF

clear;
clc;
close all;

%% 1. Configuration
data_filename = 'Path1_newnew.xlsx';
c = 3e8; 

% --- User Parameters for Selection ---
r_paths = 10;                % Max number of taps
min_spacing_ns = 25;         % Minimum time separation between taps
noise_threshold_db = 30;     % Dynamic threshold below peak for "All Taps" calculation 
                             % (Necessary to avoid noise floor blowing up RMS-DS)

%% 2. Load Data
fprintf('Loading data from %s...\n', data_filename);
try
    T = readtable(data_filename);
catch
    error('File not found. Please ensure filename is correct.');
end

col_names = T.Properties.VariableNames;
mobile_x_col = col_names{1}; freq_col = col_names{2};      
re_col = col_names{3};       im_col = col_names{4};        

mobile_x_list = sort(unique(T.(mobile_x_col)));
num_positions = length(mobile_x_list);
freq_list_hz = sort(unique(T.(freq_col))) * 1e9; 
num_freq_points = length(freq_list_hz);

%% 3. Process Data (IFFT)
fprintf('Processing IFFT for %d positions...\n', num_positions);
S21_matrix = zeros(num_positions, num_freq_points);
for i = 1:num_positions
    rows = T.(mobile_x_col) == mobile_x_list(i);
    T_subset = sortrows(T(rows, :), freq_col);
    if height(T_subset) == num_freq_points
        S21_matrix(i, :) = T_subset.(re_col) + 1j * T_subset.(im_col);
    end
end

window = hann(num_freq_points)';
impulse_response_matrix = ifft(S21_matrix .* window, [], 2);

%% 4. Axis Calculation
delta_f = freq_list_hz(2) - freq_list_hz(1);
max_range_m = c / delta_f;
max_delay_ns = (max_range_m / c) * 1e9;
time_axis_ns = linspace(0, max_delay_ns, num_freq_points);

%% 5. RMS-DS Calculation Loop
rms_ds_all = zeros(num_positions, 1);
rms_ds_active = zeros(num_positions, 1);

fprintf('Calculating RMS-DS for all snapshots...\n');

for i = 1:num_positions
    % Get PDP for current snapshot
    ir_abs = abs(impulse_response_matrix(i, :));
    pdp_linear = ir_abs.^2;
    pdp_db = 10 * log10(pdp_linear);
    [peak_power_db, ~] = max(pdp_db);
    
    % --- CASE 1: ALL TAPS (Noise Thresholded) ---
    % As per paper, we must zero out noise to get valid RMS-DS
    % Using a relative threshold (e.g. 25-30dB below peak) for "All Taps"
    valid_indices = pdp_db > (peak_power_db - noise_threshold_db);
    pwr_all = pdp_linear(valid_indices);
    tau_all = time_axis_ns(valid_indices);
    
    rms_ds_all(i) = calc_rms_ds(tau_all, pwr_all);
    
    % --- CASE 2: ACTIVE TAPS (Your Selection Logic) ---
    [pks, locs] = findpeaks(pdp_linear);
    
    if isempty(pks)
        rms_ds_active(i) = 0; % Flat line case
        continue;
    end
    
    % Sort Candidates
    [sorted_pks, sort_idx] = sort(pks, 'descend');
    sorted_locs = locs(sort_idx);
    
    % Iterative Selection (Min Spacing)
    selected_locs = [];
    selected_pks = [];
    
    for k = 1:length(sorted_locs)
        candidate_idx = sorted_locs(k);
        candidate_time = time_axis_ns(candidate_idx);
        
        is_far_enough = true;
        if ~isempty(selected_locs)
            existing_times = time_axis_ns(selected_locs);
            if any(abs(existing_times - candidate_time) < min_spacing_ns)
                is_far_enough = false;
            end
        end
        
        if is_far_enough
            selected_locs = [selected_locs; candidate_idx];
            selected_pks = [selected_pks; sorted_pks(k)];
        end
        
        if length(selected_locs) >= r_paths
            break;
        end
    end
    
    % Calculate RMS-DS for Active Taps
    tau_active = time_axis_ns(selected_locs);
    pwr_active = selected_pks';
    
    rms_ds_active(i) = calc_rms_ds(tau_active, pwr_active);
end

% Filter out NaNs or Zeros (if any snapshots failed)
valid_mask = rms_ds_all > 0 & rms_ds_active > 0;
rms_ds_all = rms_ds_all(valid_mask);
rms_ds_active = rms_ds_active(valid_mask);

%% 6. Generate CDF Plot (Paper Figure 4 style)
figure('Name', 'CDF of RMS-DS', 'Color', 'w');
hold on;

% Empirical CDF for "All Taps" (Blue in paper)
[f_all, x_all] = ecdf(rms_ds_all);
plot(x_all, f_all, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]); % Blue

% Empirical CDF for "Active Taps" (Red in paper)
[f_active, x_active] = ecdf(rms_ds_active);
plot(x_active, f_active, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); % Red

% Formatting
xlabel('RMS-DS [ns]');
ylabel('CDF');
title('CDF of RMS-DS: All Taps vs Active Taps');
legend('All Taps', 'Active Taps', 'Location', 'SouthEast');
grid on;
box on;

% Add Mean Difference annotation (optional, like paper text)
diff_mean = mean(rms_ds_all) - mean(rms_ds_active);
text(min(x_all), 0.8, sprintf('Mean Diff: %.2f ns', diff_mean), 'BackgroundColor', 'w');

fprintf('Processing complete.\n');

%% Helper Function: Calculate RMS Delay Spread
function sigma_tau = calc_rms_ds(tau, pwr)
    % Equations (5) and (6) from the paper
    % tau: delay vector (ns)
    % pwr: linear power vector (gamma^2)
    
    if isempty(pwr) || sum(pwr) == 0
        sigma_tau = 0;
        return;
    end
    
    % Normalize power
    total_pwr = sum(pwr);
    
    % Mean delay (First Moment)
    mean_delay = sum(tau .* pwr) / total_pwr;
    
    % Second Moment
    second_moment = sum((tau.^2) .* pwr) / total_pwr;
    
    % RMS Delay Spread (Square root of variance)
    val = second_moment - mean_delay^2;
    
    % numerical stability check
    if val < 0
        val = 0;
    end
    
    sigma_tau = sqrt(val);
end
