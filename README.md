% Specifying the folder path and filename for Excel file
save_path = 'C:\Users\syedm\Desktop\New folder'; 
excel_filename = 'ECG_features1.xlsx';
excel_file = fullfile(save_path, excel_filename);

%Initialize matrix to store all the features
all_features = [];

for i = 1:70
    % Command for loading all the files in MATLAB
    filename = sprintf('a%d.dat', i);

    % Open the .dat file for reading the data in binary mode
    fid = fopen(filename, 'rb');


    % Read the binary data from the file
    binary_data = fread(fid, 'int16');

    % Close the file
    fclose(fid);

% Initializing sampling freq & time
fs = 210;
t = (0:length(binary_data)-1)/fs;

% filter data
notch_frequency = 60; % Hz
bandwidth = 2; % Hz
fs = 210; % Sampling frequency 

% Normalizing the notch frequency and bandwidth to the Nyquist rate
Wo = notch_frequency / (fs / 2);
BW = bandwidth / (fs / 2);

% Designing the 2nd order IIR notch filter
[b, a] = iirnotch(Wo, BW);

filtered_data = filter(b, a, binary_data);

% Decrease the amplitude of the filtered signal
amplitude_scale_factor = 0.5; % AMplitude factor
filtered_data_scaled = amplitude_scale_factor * filtered_data;

% Applying FT
n = length(filtered_data_scaled);
f = fft(filtered_data_scaled); 

f_axis = (0:n-1)*(fs/n); % Frequency axis

    % time vector
    t1 = 1:length(filtered_data_scaled);

    % Perform MODWT on the filtered signal
    wt = modwt(filtered_data_scaled, 4, 'sym4'); %using sym4 bcz it is very good at being symmetric and efficient

    % Initialize wtrec to be of the same size as wt
    wtrec = zeros(size(wt));

    % Retain only specific levels (e.g., levels 3 and 4)
    wtrec(3:4, :) = wt(3:4, :);

    % Performign inverse MODWT to reconstruct the signal
    y = imodwt(wtrec, 'sym4');

    % Square the absolute value of the reconstructed signal
    y1 = abs(y).^2;

    % Calculate the mean value of y
    avg = mean(y1);

    % Find peaks in y with specific criteria
    [Rpeaks, locs] = findpeaks(y1, t1, 'MinPeakHeight', 6 * avg, 'MinPeakDistance', 60);

    % location contains the locations of the detected R-peaks

    % Calculate RR intervals
    RR_intervals = diff(locs); % RR intervals in samples
    HR = 60 ./ RR_intervals; % Heart rate in beats per minute

    % Number of RR intervals
    num_intervals = length(RR_intervals);

    % Extract features from RR intervals
    mean_RR = mean(RR_intervals);
    std_RR = std(RR_intervals);
    median_RR = median(RR_intervals);
    min_RR = min(RR_intervals);
    max_RR = max(RR_intervals);
    range_RR = max_RR - min_RR;
    meanRpeaks = mean(Rpeaks);

    % Heart rate variability features
    Rmssd = sqrt(mean(diff(RR_intervals).^2)); % Root mean square of successive differences
    Nn50 = sum(abs(diff(RR_intervals)) > 0.05); % Number of pairs of successive NNs that differ by more than 50 ms
    pnn50 = 100 * Nn50 / num_intervals; % Proportion of NN50 divided by the total number of NNs
    Stdnn = std(RR_intervals); % Standard deviation of NN interval
    SD_HR = std(HR); % Standard deviation of heart rate

    % Combine features into a row vector
    features = horzcat(mean_RR, std_RR, median_RR, min_RR, max_RR, range_RR, meanRpeaks, Rmssd, Stdnn, Nn50, pnn50,SD_HR);

    % Store the features for this file
    all_features = [all_features; features];
end

% Write the features matrix to the Excel file
xlswrite(excel_file, all_features);

% Plot the filtered ECG signal with detected R peaks
subplot(3, 1, 1);
t_5sec = t(t <= 5);
plot(t_5sec, filtered_data_scaled(1:length(t_5sec)));
hold on;
t_5sec = t(t <= 5);
plot(locs, Rpeaks, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
title('ECG Signal with Detected R Peaks');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Filtered ECG', 'R Peaks');
grid on;
hold off;

% Plot the magnitude of the FFT
subplot(3, 1, 2);
plot(f_axis, abs(f)); % Plot the magnitude of the FFT
title('ECG Signal in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Plot the filtered data
subplot(3, 1, 3);
plot(t_5sec, filtered_data_scaled(1:length(t_5sec)));
title('Filtered Plot');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot the original binary data
figure;
subplot(2, 1, 1);
plot(t_5sec, binary_data(1:length(t_5sec)));
title('Original Binary Data');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot the R-peaks
subplot(2, 1, 2);
plot(t_5sec, y1(1:length(t_5sec)));
title('R Peaks');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
