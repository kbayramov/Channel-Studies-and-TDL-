% SCRIPT TO GENERATE TDL MODEL (NO THRESHOLD + MIN SPACING)
% Features: Full Range (250m), Spacing Constraint (25ns), No Noise Filter

clear;
clc;
close all;

%% 1. Configuration
data_filename = 'Path1_newnew.xlsx';
c = 3e8; 

% --- User Parameters ---
target_pos = 70;            % Position to analyze
r_paths = 10;               % Max number of taps
min_spacing_ns = 25;        % Minimum time separation between taps

% --- Plotting ---
plot_min_db = -100;

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

%% 3. Process Data (IFFT)
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

%% 4. Axis Calculation (Full Range)
delta_f = freq_list_hz(2) - freq_list_hz(1);
max_range_m = c / delta_f;
max_delay_ns = (max_range_m / c) * 1e9;

fprintf('Frequency Step: %.2f MHz\n', delta_f/1e6);
fprintf('Max Delay Limit: %.2f ns\n', max_delay_ns);

time_axis_ns = linspace(0, max_delay_ns, num_freq_points);

ir_valid = impulse_response_matrix; % Full buffer
ir_abs = abs(ir_valid);
ir_db = 20 * log10(ir_abs);

%% 5. Extract TDL Model (Spacing Constrained, No Threshold)
fprintf('\n--- Processing Position %.1f m ---\n', target_pos);
idx = find(mobile_x_list == target_pos, 1);

if ~isempty(idx)
    pdp_linear = ir_abs(idx, :).^2;
    pdp_db = 10 * log10(pdp_linear); 
    
    [peak_power_db, ~] = max(pdp_db);
    
    % 1. Find ALL Peaks (No MinPeakHeight threshold)
    [pks, locs] = findpeaks(pdp_linear);
    
    if isempty(pks)
        warning('No peaks found (Flat signal).');
    else
        % 2. Sort Candidates by Strength (Strongest First)
        [sorted_pks, sort_idx] = sort(pks, 'descend');
        sorted_locs = locs(sort_idx);
        
        % 3. Iterative Selection with Spacing Constraint
        selected_locs = [];
        selected_pks = [];
        
        for k = 1:length(sorted_locs)
            candidate_idx = sorted_locs(k);
            candidate_time = time_axis_ns(candidate_idx);
            
            % Check spacing against ALREADY selected taps
            is_far_enough = true;
            if ~isempty(selected_locs)
                existing_times = time_axis_ns(selected_locs);
                time_diffs = abs(existing_times - candidate_time);
                
                % If candidate is too close to ANY existing tap, reject it
                if any(time_diffs < min_spacing_ns)
                    is_far_enough = false;
                end
            end
            
            % If valid, add it
            if is_far_enough
                selected_locs = [selected_locs; candidate_idx];
                selected_pks = [selected_pks; sorted_pks(k)];
            end
            
            % Stop if we have enough taps
            if length(selected_locs) >= r_paths
                break;
            end
        end
        
        fprintf('Selected %d taps (Constraint: >%.1f ns spacing).\n', ...
            length(selected_locs), min_spacing_ns);
        
        % 4. Convert to TDL Format
        % Reference Delay = First Arriving Tap among SELECTED ones
        [min_loc_idx, ~] = min(selected_locs); 
        min_time_ns = time_axis_ns(min_loc_idx);
        
        final_times_ns = time_axis_ns(selected_locs);
        rel_delays_ns = final_times_ns - min_time_ns;
        
        % Relative Power
        final_pks_db = 10 * log10(selected_pks);
        rel_power_db = final_pks_db - max(final_pks_db);
        
        % Fix dimensions for table
        tap_indices = (1:length(selected_locs))';
        
        TDL_Table = table(tap_indices, ...
                          final_times_ns(:), ...
                          rel_delays_ns(:), ...
                          rel_power_db(:), ...
            'VariableNames', {'Tap', 'Abs_Delay_ns', 'Rel_Delay_ns', 'Rel_Power_dB'});
            
        TDL_Table = sortrows(TDL_Table, 'Rel_Delay_ns');
        disp(TDL_Table);
        
        % 5. Plotting
        figure('Name', 'TDL Model (No Threshold)', 'Color', 'w');
        
        % Plot Raw PDP
        plot(time_axis_ns, pdp_db - peak_power_db, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
        
        % Plot Selected Taps
        stem(TDL_Table.Abs_Delay_ns, TDL_Table.Rel_Power_dB, 'r', 'filled', 'LineWidth', 2, 'BaseValue', -100);
        
        title(sprintf('TDL @ %.1fm (Min Spacing: %.0fns, No Threshold)', target_pos, min_spacing_ns));
        xlabel('Excess Delay (ns)'); ylabel('Relative Power (dB)');
        legend('Raw PDP', 'Selected Taps');
        grid on; 
        
        xlim([0, max_delay_ns]); 
        
        % Dynamic Y-Limits
        if isempty(TDL_Table)
             ylim([-100, 5]);
        else
             min_pwr = min(TDL_Table.Rel_Power_dB);
             ylim([min_pwr - 10, 5]);
        end
    end
else
    fprintf('Position not found.\n');
end
