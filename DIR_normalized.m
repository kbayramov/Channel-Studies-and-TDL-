% SCRIPT: DIR HEATMAP - ALIGNED & NORMALIZED
% Visualization: Position (Y) vs. Excess Range (X)
% Features:
%   1. Removes Path Loss (Peak = 0 dB)
%   2. Removes Time-of-Flight (Peak = 0 meters)

clear;
clc;
close all;

%% 1. Configuration
data_filename = 'S21_Re_Im_Path1.xlsx';
c = 3e8; 

% --- Plotting Limits ---
plot_max_db = 0;    % Peak is always 0
plot_min_db = -40;  % Show 40dB dynamic range

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

%% 4. Alignment & Normalization
delta_f = freq_list_hz(2) - freq_list_hz(1);
max_range = c / delta_f; 
range_axis_m = linspace(0, max_range, num_freq_points);

% Matrix to store the aligned (shifted) data
dir_db_aligned = zeros(size(impulse_response_matrix));

fprintf('Aligning snapshots to remove Time-of-Flight...\n');

for i = 1:num_positions
    % Get raw amplitude for this snapshot
    raw_response = abs(impulse_response_matrix(i, :));
    
    % 1. Find the Peak (assumed to be LOS)
    [max_val, max_idx] = max(raw_response);
    
    % 2. Normalize Magnitude (Peak -> 0 dB)
    % Convert to dB and shift so max is 0
    row_db = 20 * log10(raw_response) - 20 * log10(max_val);
    
    % 3. Align Delay (Shift Peak to Index 1)
    % circshift moves the data so 'max_idx' becomes index 1
    shift_amount = -(max_idx - 1);
    aligned_row = circshift(row_db, shift_amount);
    
    % (Optional) Mask the "wrap-around" part if desired. 
    % IFFT is circular, so wrapping is technically correct, 
    % but visually we often ignore the wrapped tail.
    % For now, we keep it to preserve all energy.
    
    dir_db_aligned(i, :) = aligned_row;
end

%% 5. Generate Heatmap
figure('Name', 'DIR Heatmap (Excess Range)', 'Color', 'w', 'Position', [100 100 1200 800]);

imagesc(range_axis_m, mobile_x_list, dir_db_aligned);

% --- Styling ---
colormap('jet');
colorbar;
caxis([plot_min_db, plot_max_db]); 
axis xy; 

% Limits
% We zoom in on the first 100m because everything is now compressed to the left
xlim([0, 100]); 
ylim([min(mobile_x_list), max(mobile_x_list)]);

% Labels
title('Impulse Response (Aligned to First Peak)');
xlabel('Excess Range (m) - Relative to LOS');
ylabel('AGV Position (m)');
zlabel('Relative Power (dB)');

% Add Reference Line at 0 (Main Path)
xline(0, 'w-', 'Main Path (Aligned)');

fprintf('Plot Generated. All positions aligned to x=0.\n');
