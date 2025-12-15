% SCRIPT: INTER-TAP CORRELATION COEFFICIENTS
% Objective: Analyze if taps fade together or independently.
% Method: Pearson Correlation of Normalized Power History.

clear;
clc;
close all;

%% 1. Configuration
data_filename = 'Path1_newnew.xlsx';
c = 3e8; 

% --- Analysis Parameters ---
num_taps_to_track = 10;     % We will correlate the top 10 taps
min_spacing_ns = 25;        % Spacing constraint
physics_threshold_db = 25;  % Power threshold

%% 2. Load Data & Prepare
fprintf('Loading data...\n');
try
    T = readtable(data_filename);
catch
    error('File not found.');
end

col_names = T.Properties.VariableNames;
mobile_x_col = col_names{1}; freq_col = col_names{2};     
re_col = col_names{3};       im_col = col_names{4};       

mobile_x_list = sort(unique(T.(mobile_x_col)));
freq_list_hz = sort(unique(T.(freq_col))) * 1e9; 
num_freq_points = length(freq_list_hz);
num_positions = length(mobile_x_list);

% IFFT Processing
fprintf('Processing IFFT...\n');
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

% Axis Setup
delta_f = freq_list_hz(2) - freq_list_hz(1);
time_axis_ns = linspace(0, 1/delta_f*1e9, num_freq_points);

%% 3. Extract Tap History (The "Tracking" Matrix)
% Matrix Size: [Num_Positions x Num_Taps_Tracked]
Tap_Power_History = nan(num_positions, num_taps_to_track);

fprintf('Extracting top %d taps for correlation analysis...\n', num_taps_to_track);

for i = 1:num_positions
    % Get PDP
    h_inst = impulse_response_matrix(i, :);
    pdp_linear = abs(h_inst).^2;
    
    % --- CRITICAL: NORMALIZE TO REMOVE PATH LOSS ---
    % We want correlation of FADING, not correlation of DISTANCE.
    total_energy = sum(pdp_linear);
    pdp_norm_db = 10*log10(pdp_linear / total_energy); 
    
    % Find Peaks (Physics Based)
    [peak_pwr, ~] = max(pdp_norm_db);
    cutoff = peak_pwr - physics_threshold_db;
    
    [pks, locs] = findpeaks(10.^(pdp_norm_db/10), 'MinPeakHeight', 10^(cutoff/10));
    pks_db = 10*log10(pks);
    
    % Sort & Spacing Filter
    [sorted_pks_db, sort_idx] = sort(pks_db, 'descend');
    sorted_locs = locs(sort_idx);
    
    selected_pks_db = [];
    selected_locs = [];
    
    for k = 1:length(sorted_locs)
        candidate_time = time_axis_ns(sorted_locs(k));
        is_distinct = true;
        if ~isempty(selected_locs)
            existing_times = time_axis_ns(selected_locs);
            if any(abs(existing_times - candidate_time) < min_spacing_ns)
                is_distinct = false;
            end
        end
        if is_distinct
            selected_locs = [selected_locs; sorted_locs(k)];
            selected_pks_db = [selected_pks_db; sorted_pks_db(k)];
        end
        if length(selected_locs) >= num_taps_to_track
            break;
        end
    end
    
    % Fill Matrix (Pad with noise floor if fewer taps found)
    num_found = length(selected_pks_db);
    Tap_Power_History(i, 1:num_found) = selected_pks_db;
    if num_found < num_taps_to_track
        % Fill remaining with a low floor (e.g. -100dB) to represent "missing"
        Tap_Power_History(i, num_found+1:end) = -100; 
    end
end

%% 4. Calculate Correlation Matrix
% We compute the correlation of the COLUMNS (Taps) over the ROWS (Positions)
Corr_Matrix = corrcoef(Tap_Power_History);

% Display in Command Window
disp('--- Inter-Tap Correlation Matrix (Top 5 Taps) ---');
disp(array2table(Corr_Matrix(1:5, 1:5), 'VariableNames', {'Tap1','Tap2','Tap3','Tap4','Tap5'}, ...
    'RowNames', {'Tap1','Tap2','Tap3','Tap4','Tap5'}));

%% 5. Visualization (Correlation Heatmap)
figure('Name', 'Tap Correlation Matrix', 'Color', 'w', 'Position', [100 100 900 700]);

% Using imagesc for the heatmap
imagesc(Corr_Matrix);
colormap(jet);
colorbar;
caxis([-1 1]); % Correlation ranges from -1 to 1

% Styling
title('Inter-Tap Correlation Coefficients (Normalized Fading)');
xlabel('Tap Index (Sorted by Power)');
ylabel('Tap Index (Sorted by Power)');
xticks(1:num_taps_to_track);
yticks(1:num_taps_to_track);
axis square;

% Annotate values on the grid
for r = 1:num_taps_to_track
    for c = 1:num_taps_to_track
        text(c, r, sprintf('%.2f', Corr_Matrix(r,c)), ...
            'HorizontalAlignment', 'center', ...
            'Color', 'k', 'FontWeight', 'bold', 'FontSize', 9);
    end
end

fprintf('Plot Generated. Values near 0 indicate independent fading.\n');
