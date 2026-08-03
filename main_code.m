clc;
clear;
close all;

%% ===== 1. Select ABF file =====
[file, path] = uigetfile('*.abf', 'Select ABF file');

if isequal(file,0)
    error('No ABF file selected.');
end

abf_file = fullfile(path, file);

%% ===== 2. Read ABF file =====
[data, si, header] = abfload(abf_file);

if size(data,1) < size(data,2)
    data = data';
end

%% ===== 3. Read Channel 1 =====
signal = data(:,1);

fprintf('Using Channel 1 (IN 0)\n');

%% ===== 4. Time axis =====
sample_rate = 1/(si*1e-6);

N = length(signal);

time = (0:N-1)/sample_rate;

%% ===== 5. Smooth signal (ONLY for detection) =====
smooth_signal = smoothdata(signal,'movmean',10);

% nanopore downward events -> invert signal
sig = -smooth_signal;

%% ===== 6. Detection parameters =====
minProm = 50;
minWidth = 0.0001;
minDist = 0.0005;

%% ===== 7. Detect events =====
[~, locs] = findpeaks(sig, time, ...
    'MinPeakProminence', minProm, ...
    'MinPeakWidth', minWidth, ...
    'MinPeakDistance', minDist);

idx = round(locs * sample_rate) + 1;

%% ===== 8. Find real peaks in RAW signal =====
search_window = 20;

true_idx = zeros(size(idx));

for i = 1:length(idx)

    left = max(idx(i)-search_window,1);
    right = min(idx(i)+search_window,N);

    [~,local_min] = min(signal(left:right));

    true_idx(i) = left + local_min - 1;

end

true_peaks = signal(true_idx);

true_time = time(true_idx);

%% ===== 9. Global baseline =====
baseline = mean(signal(1:round(0.1*N)));

%% ===== 10. Dwell time (simplified FWHM) =====
fwhm = zeros(length(true_idx),1);

t1_all = zeros(length(true_idx),1);
t2_all = zeros(length(true_idx),1);

h_all = zeros(length(true_idx),1);

alpha = 0.5;

for i = 1:length(true_idx)

    peak_i = true_idx(i);

    peak_val = signal(peak_i);

    % half level
    half_level = baseline - alpha*(baseline - peak_val);

    left = max(peak_i-500,1);
    right = min(peak_i+500,N);

    %% left crossing
    left_idx = find(signal(left:peak_i) > half_level, ...
        1, 'last');

    %% right crossing
    right_idx = find(signal(peak_i:right) > half_level, ...
        1, 'first');

    if ~isempty(left_idx) && ~isempty(right_idx)

        t1 = time(left + left_idx -1);

        t2 = time(peak_i + right_idx -1);

        fwhm(i) = t2 - t1;

        t1_all(i) = t1;
        t2_all(i) = t2;

        h_all(i) = half_level;

    end

end

%% ===== 11. Remove invalid =====
valid = fwhm > 0 & fwhm < 1e-3;

fwhm = fwhm(valid);

t1_all = t1_all(valid);
t2_all = t2_all(valid);

h_all = h_all(valid);

fprintf('Valid events: %d / %d\n', ...
    length(fwhm), length(true_idx));

%% ===== 12. Mean dwell time =====
mean_dwell = mean(fwhm);

fprintf('Mean dwell time: %.6f s (%.2f us)\n', ...
    mean_dwell, mean_dwell*1e6);

%% ===== 13. Event depth =====
depth = zeros(length(true_idx),1);

for i = 1:length(true_idx)

    peak_val = signal(true_idx(i));

    depth(i) = baseline - peak_val;

end

depth = depth(valid);

%% ===== 14 Detect pre-event upward peaks =====

pre_idx = zeros(length(true_idx),1);
pre_peaks = zeros(length(true_idx),1);

% IR for every event
IR_all = zeros(length(true_idx),1);

% search window before event
pre_window_us = 500;

pre_window = round(pre_window_us * 1e-6 * sample_rate);

% upward peak must be at least 15 pA above baseline
peak_threshold = 15;

for i = 1:length(true_idx)

    peak_i = true_idx(i);

    left = max(peak_i - pre_window, 1);
    right = peak_i - 1;

    if right > left

        [local_max, local_idx] = max(signal(left:right));

        % Only keep peaks above baseline + threshold
        if local_max > (baseline + peak_threshold)

            pre_idx(i) = left + local_idx - 1;

            pre_peaks(i) = local_max;

            % IR = peak height above baseline
            IR_all(i) = local_max - baseline;

        end

    end

end

%% Keep only valid nanopore events

pre_idx   = pre_idx(valid);
pre_peaks = pre_peaks(valid);

IR = IR_all(valid);

%% Events that actually have an upward peak

has_peak = pre_idx > 0;

pre_time_plot  = time(pre_idx(has_peak));
pre_peaks_plot = pre_peaks(has_peak);

%% ===== 15 Calculate IT =====

IC = depth;

IT = (IC + IR) ./ IC;

fprintf('Mean IT = %.4f\n', mean(IT));
fprintf('Median IT = %.4f\n', median(IT));

%% ===== 16. Plot signal + FWHM =====
figure

plot(time, signal, 'b')

hold on

%% detected nanopore events
plot(true_time(valid), ...
    true_peaks(valid), ...
    'ro', ...
    'MarkerSize',6)

%% pre-event upward peaks
plot(pre_time_plot, ...
    pre_peaks_plot, ...
    'mo', ...
    'MarkerSize',5, ...
    'LineWidth',1.2)

%% FWHM lines
for i = 1:length(fwhm)

    plot([t1_all(i), t2_all(i)], ...
        [h_all(i), h_all(i)], ...
        'g', ...
        'LineWidth',1.5)

end

xlabel('Time (s)')

ylabel('Current (pA)')

title('Nanopore Signal')

legend('Raw Signal', ...
       'Event Peaks', ...
       'Pre-event Peaks', ...
       'FWHM')

grid on

%% ===== 17. Dwell Time Histogram =====
figure

histogram(fwhm,50)

xlabel('Dwell Time (s)')

ylabel('Frequency')

title('Dwell Time Distribution')

grid on

%% ===== 18. Depth Histogram =====
figure

histogram(depth,50)

xlabel('Event Depth (pA)')

ylabel('Frequency')

title('Nanopore Event Depth Distribution')

grid on

%% ===== 19 IT Histogram =====

figure

histogram(IT,50)

xlabel('IT')

ylabel('Frequency')

title('IT Distribution')

grid on

%% ===== 20. Scatter Plot =====
amplitude = baseline - true_peaks(valid);

figure

scatter(fwhm*1e6, amplitude, ...
    8, ...
    'b', ...
    'filled', ...
    'MarkerFaceAlpha',0.2)

xlabel('Dwell Time (\mus)')

ylabel('Peak Amplitude (pA)')

title('Scatter Plot: Peak Amplitude vs Dwell Time')

grid on

