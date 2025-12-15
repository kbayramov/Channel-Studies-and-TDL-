% SCRIPT TO GENERATE TDL PARAMETERS (TABLE II)
% Context: AGV (0.3 m/s, 1.8 GHz)
% Based on Ansari et al., VTC2024-Fall [Section IV-C]

clear;
clc;
close all;

%% 1. Configuration
data_filename = 'Path1_newnew.xlsx';

% --- AGV Context Parameters ---
velocity_m_s = 0.3;           % AGV Speed
carrier_freq_ghz = 1.8;       % Central Frequency
spatial_window_m = 1.5;       % Distance window to separate Slow/Fast fading 

% --- TDL Settings ---
num_taps_to_analyze = 10;     
min_spacing_ns = 25;          
multipath_threshold_db = 25;  % Dynamic threshold for Persistence (25dB below peak)

%% 2. Load Data & Process IFFT
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

% Calculate Spatial Step and Window Size
dx = mean(diff(mobile_x_list)); 
if isnan(dx) || dx == 0, dx = 0.1; end % Fallback
window_size_samples = round(spatial_window_m / dx);

fprintf('AGV Speed: %.1f m/s. Spatial Window: %.1fm (%d samples).\n', ...
    velocity_m_s, spatial_window_m, window_size_samples);

% IFFT Processing
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

% Time Axis
delta_f = freq_list_hz(2) - freq_list_hz(1);
max_delay_ns = (1 / delta_f) * 1e9;
time_axis_ns = linspace(0, max_delay_ns, num_freq_points);

%% 3. Identify Global Taps (Fixed Delay Grid)
% Use Average PDP to find the 10 dominant "Fixed" Taps
avg_pdp_linear = mean(abs(impulse_response_matrix).^2, 1);
[pks, locs] = findpeaks(avg_pdp_linear);
[sorted_pks, sort_idx] = sort(pks, 'descend');
sorted_locs = locs(sort_idx);

selected_locs = [];
for k = 1:length(sorted_locs)
    candidate_idx = sorted_locs(k);
    candidate_time = time_axis_ns(candidate_idx);
    
    is_far = true;
    if ~isempty(selected_locs)
        existing_times = time_axis_ns(selected_locs);
        if any(abs(existing_times - candidate_time) < min_spacing_ns)
            is_far = false;
        end
    end
    if is_far
        selected_locs = [selected_locs; candidate_idx];
    end
    if length(selected_locs) >= num_taps_to_analyze, break; end
end
[selected_locs, ~] = sort(selected_locs, 'ascend');
final_tap_delays_ns = time_axis_ns(selected_locs);

%% 4. Parameter Extraction
tap_energies = zeros(num_taps_to_analyze, 1);
P1_list = zeros(num_taps_to_analyze, 1);
P0_list = zeros(num_taps_to_analyze, 1);
s_rice_list = zeros(num_taps_to_analyze, 1);
sigma_rice_list = zeros(num_taps_to_analyze, 1);

% Pre-calc thresholds for persistence
snapshot_peaks_linear = max(abs(impulse_response_matrix).^2, [], 2);
thresholds_linear = snapshot_peaks_linear * 10^(-multipath_threshold_db/10); 

for t = 1:length(selected_locs)
    idx = selected_locs(t);
    tap_amps_raw = abs(impulse_response_matrix(:, idx)); % Raw Linear Amplitudes
    tap_powers = tap_amps_raw.^2;
    
    % --- A. Persistence (ON/OFF) [cite: 137] ---
    is_ON = tap_powers > thresholds_linear;
    count_ON = sum(is_ON);
    
    P1 = count_ON / num_positions;
    P0 = 1 - P1;
    P1_list(t) = P1;
    P0_list(t) = P0;
    
    % --- B. Energy & Fading Analysis ---
    if count_ON > window_size_samples % Need enough data for moving average
        % 1. Calculate Raw Energy (for Table "Energy" column)
        avg_on_energy = mean(tap_powers(is_ON));
        tap_energies(t) = avg_on_energy * P1; 
        
        % 2. Extract Fast Fading (Equation 7 & 8 in paper [cite: 184-188])
        % We operate on the ON samples (or fill OFF with interpolation)
        % For robustness, we analyze the continuous segments or just the ON points
        valid_amps = tap_amps_raw(is_ON);
        
        % Convert to dB
        R_total_db = 20*log10(valid_amps);
        
        % Isolate Slow Fading (Moving Average)
        R_slow_db = movmean(R_total_db, window_size_samples);
        
        % Isolate Fast Fading
        R_fast_db = R_total_db - R_slow_db;
        
        % Convert back to Linear
        r_fast_linear = 10.^(R_fast_db / 20);
        
        % 3. NORMALIZE to Unit Power (Crucial for correct s/sigma)
        % This ensures E[r^2] = 1, so s and sigma are standard parameters
        rms_val = sqrt(mean(r_fast_linear.^2));
        r_fast_norm = r_fast_linear / rms_val;
        
        % 4. Fit Rician Distribution
        try
            pd = fitdist(r_fast_norm, 'Rician');
            s_rice_list(t) = pd.s;      % Non-centrality parameter
            sigma_rice_list(t) = pd.sigma; % Scale parameter
        catch
            s_rice_list(t) = NaN;
            sigma_rice_list(t) = NaN;
        end
    else
        % Not enough ON samples
        tap_energies(t) = 0;
        s_rice_list(t) = NaN;
        sigma_rice_list(t) = NaN;
    end
end

% Normalize Energy Column (Sum = 1)
normalized_energies = tap_energies / sum(tap_energies);

%% 5. Display Table II
fprintf('\n=========================================================================\n');
fprintf('TABLE II: Parameters for TDL Channel based on Measurement Data\n');
fprintf('=========================================================================\n');
fprintf('| Tap | Delay(ns)| Energy(lin)|   P1    |   P0    |  SRice  |  Rice   |\n');
fprintf('|-----|----------|------------|---------|---------|---------|---------|\n');

for t = 1:num_taps_to_analyze
    s_val = s_rice_list(t);
    sig_val = sigma_rice_list(t);
    
    if isnan(s_val), s_str="NaN"; else, s_str=sprintf("%.3f",s_val); end
    if isnan(sig_val), sig_str="NaN"; else, sig_str=sprintf("%.3f",sig_val); end

    fprintf('| %-3d |  %6.1f  |   %.4f   |  %.3f  |  %.3f  |  %s  |  %s  |\n', ...
        t, final_tap_delays_ns(t), normalized_energies(t), ...
        P1_list(t), P0_list(t), s_str, sig_str); 
end
fprintf('-------------------------------------------------------------------------\n');
fprintf('Context: AGV v=%.1f m/s, f=%.1f GHz.\n', velocity_m_s, carrier_freq_ghz);
fprintf('Note: SRice/Rice parameters are fitted to Unit-Power Normalized Fast Fading.\n');
