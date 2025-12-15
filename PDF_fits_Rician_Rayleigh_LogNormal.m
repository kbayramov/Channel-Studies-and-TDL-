% SCRIPT TO ANALYZE AMPLITUDE FADING & GENERATE HI TABLE (10 TAPS)
% Based on methodology from Ansari et al., VTC2024-Fall [Section IV-C]

clear;
clc;
close all;

%% 1. Configuration
data_filename = 'Path1_newnew.xlsx';

% --- User Parameters ---
num_taps_to_analyze = 10;    % Updated to analyze 10 taps
min_spacing_ns = 25;         % Minimum time separation
spatial_window_m = 1.5;      % "Minimal stationary length ds_min = 1.5m" [cite: 194]

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

% Calculate Spatial Resolution for Windowing
dx = mean(diff(mobile_x_list)); % Average step size in meters
if isnan(dx) || dx == 0, dx = 0.1; end % Fallback to avoid division by zero
window_size_samples = round(spatial_window_m / dx);
fprintf('Spatial Step: %.4f m. Window Size (w): %d samples.\n', dx, window_size_samples);

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

% Axis Calculation
c = 3e8;
delta_f = freq_list_hz(2) - freq_list_hz(1);
max_delay_ns = (1 / delta_f) * 1e9;
time_axis_ns = linspace(0, max_delay_ns, num_freq_points);

%% 4. Extract Amplitude History for Taps
% Matrix to store amplitude history: [Num_Snapshots x Num_Taps]
raw_tap_amplitudes_linear = nan(num_positions, num_taps_to_analyze);

fprintf('Extracting tap amplitudes across snapshots...\n');
for i = 1:num_positions
    ir_abs = abs(impulse_response_matrix(i, :));
    pdp_linear = ir_abs.^2;
    
    [pks, locs] = findpeaks(pdp_linear);
    if isempty(pks), continue; end
    
    % Sort by strength (to find the most significant paths first)
    [sorted_pks, sort_idx] = sort(pks, 'descend');
    sorted_locs = locs(sort_idx);
    
    % Select Taps (with Spacing Constraint)
    selected_locs = [];
    selected_pks = []; % Linear Power
    
    for k = 1:length(sorted_locs)
        curr_loc = sorted_locs(k);
        curr_time = time_axis_ns(curr_loc);
        
        is_far = true;
        if ~isempty(selected_locs)
            existing_times = time_axis_ns(selected_locs);
            if any(abs(existing_times - curr_time) < min_spacing_ns)
                is_far = false;
            end
        end
        
        if is_far
            selected_locs = [selected_locs; curr_loc];
            % Store AMPLITUDE (sqrt of power) because fitting is on Amplitude
            selected_pks = [selected_pks; sqrt(sorted_pks(k))]; 
        end
        
        if length(selected_locs) >= num_taps_to_analyze
            break;
        end
    end
    
    % Sort the collected taps by DELAY for consistent indexing
    % Tap 1 = First Arrival, Tap 2 = Second Arrival, etc.
    [~, delay_sort_idx] = sort(selected_locs, 'ascend');
    sorted_by_delay_amps = selected_pks(delay_sort_idx);
    
    % Store in matrix
    count = min(length(sorted_by_delay_amps), num_taps_to_analyze);
    raw_tap_amplitudes_linear(i, 1:count) = sorted_by_delay_amps(1:count);
end

%% 5. Fast/Slow Fading Separation & Fitting
% Prepare Storage for Table
table_data = cell(num_taps_to_analyze, 4); % Tap | Rician | Rayleigh | LogNormal

% Figure handle for Tap 1 Plot
fig_handle = figure('Name', 'Tap 1 Amplitude Fading', 'Color', 'w');

fprintf('Fitting distributions for %d taps...\n', num_taps_to_analyze);

for t = 1:num_taps_to_analyze
    % Get data for Tap 't'
    amp_data = raw_tap_amplitudes_linear(:, t);
    
    % Remove NaNs (snapshots where this tap wasn't found)
    amp_data = amp_data(~isnan(amp_data));
    
    % Skip if insufficient data
    if length(amp_data) < window_size_samples * 2
        fprintf('Tap %d: Not enough data points to analyze.\n', t);
        continue;
    end
    
    % --- A. Convert to dB ---
    R_t = 20 * log10(amp_data); 
    
    % --- B. Slow Fading (Moving Average) ---
    % Window size determined by spatial window (1.5m)
    R_s = movmean(R_t, window_size_samples);
    
    % --- C. Fast Fading (Subtraction) ---
    % Eq (8): R_t = R_s + R_f  => R_f = R_t - R_s
    R_f_db = R_t - R_s;
    
    % --- D. Convert back to Linear for Fitting ---
    % The paper fits distributions to the linear amplitude of the fast fading component
    r_f_linear = 10.^(R_f_db / 20);
    
    % --- E. Fit Distributions ---
    try
        pd_rice = fitdist(r_f_linear, 'Rician');
        pd_ray  = fitdist(r_f_linear, 'Rayleigh');
        pd_logn = fitdist(r_f_linear, 'Lognormal');
    catch
        warning('Fitting failed for Tap %d', t);
        continue;
    end
    
    % --- F. Calculate Histogram Intersection (HI) ---
    % Eq (11): HI = Sum( min(M_i, F_i) )
    num_bins = 50;
    [counts, edges] = histcounts(r_f_linear, num_bins, 'Normalization', 'pdf');
    bin_centers = edges(1:end-1) + diff(edges)/2;
    bin_width = edges(2) - edges(1);
    
    % Evaluate PDFs at bin centers
    y_rice = pdf(pd_rice, bin_centers);
    y_ray  = pdf(pd_ray,  bin_centers);
    y_logn = pdf(pd_logn, bin_centers);
    
    % Calculate HI
    HI_Rice = sum(min(counts, y_rice)) * bin_width;
    HI_Ray  = sum(min(counts, y_ray))  * bin_width;
    HI_Logn = sum(min(counts, y_logn)) * bin_width;
    
    % Store Data for Table
    table_data{t, 1} = t;
    table_data{t, 2} = HI_Rice;
    table_data{t, 3} = HI_Ray;
    table_data{t, 4} = HI_Logn;
    
    % --- G. Plotting (Only for Tap 1) ---
    if t == 1
        figure(fig_handle);
        hold on;
        
        % 1. Empirical Histogram
        h = histogram(r_f_linear, num_bins, 'Normalization', 'pdf', ...
            'FaceColor', [0.8 0.9 0.8], 'EdgeColor', 'none'); % Light Green
        
        % 2. Plot Fits
        x_grid = linspace(min(r_f_linear), max(r_f_linear), 1000);
        
        plot(x_grid, pdf(pd_rice, x_grid), 'r-', 'LineWidth', 2);
        plot(x_grid, pdf(pd_ray, x_grid), 'b-', 'LineWidth', 2);
        plot(x_grid, pdf(pd_logn, x_grid), 'k-', 'LineWidth', 2);
        
        title('Tap 1: Fast Fading Distribution (Fig. 5)');
        xlabel('Amplitude (linear)');
        ylabel('Probability Density');
        legend('Empirical', 'Rician', 'Rayleigh', 'Log normal');
        grid on;
        xlim([0, 2.5]); 
    end
end

%% 6. Display Table I (Extended to 10 Taps)
fprintf('\nTABLE I: Distance between measured data and PDF fits (HI)\n');
fprintf('----------------------------------------------------------\n');
fprintf('| Tap |   Rician   |  Rayleigh  | Log Normal |\n');
fprintf('----------------------------------------------------------\n');
for t = 1:num_taps_to_analyze
    if isempty(table_data{t,1})
        fprintf('| %-3d |    NaN     |    NaN     |    NaN     |\n', t);
    else
        fprintf('| %-3d |   %.3f    |   %.3f    |   %.3f    |\n', ...
            table_data{t,1}, table_data{t,2}, table_data{t,3}, table_data{t,4});
    end
end
fprintf('----------------------------------------------------------\n');
