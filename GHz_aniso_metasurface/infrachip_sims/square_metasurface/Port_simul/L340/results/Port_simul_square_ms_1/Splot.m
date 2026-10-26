% Define the file name 
filename1 = 'jones_parameters_1.txt';
filename2 = 'jones_parameters_2.txt';

% Import the data 
% Skip the first line because it contains headers
data1 = readmatrix(filename1, 'NumHeaderLines', 1);
data2 = readmatrix(filename2, 'NumHeaderLines', 1);

% Extract the frequency (x-axis) and the real/imaginary parts
freq1 = data1(:, 1);               % Frequency in THz
real1_S11 = data1(:, 2);           % Real part of S11
imag1_S11 = data1(:, 3);           % Imaginary part of S11
real1_S21 = data1(:, 4);           % Real part of S21
imag1_S21 = data1(:, 5);           % Imaginary part of S21
real1_S31 = data1(:, 6);           % Real part of S31
imag1_S31 = data1(:, 7);           % Imaginary part of S31

% Create the figure
figure(1)

% Plot all the variables on the same plot
hold on;
plot(freq1, real1_S11, 'DisplayName', 'Real(S11)', 'LineWidth', 2);
plot(freq1, imag1_S11, 'DisplayName', 'Imag(S11)', 'LineWidth', 2);
plot(freq1, real1_S21, 'DisplayName', 'Real(S21)', 'LineWidth', 2);
plot(freq1, imag1_S21, 'DisplayName', 'Imag(S21)', 'LineWidth', 2);
plot(freq1, real1_S31, 'DisplayName', 'Real(S31)', 'LineWidth', 2);
plot(freq1, imag1_S31, 'DisplayName', 'Imag(S31)', 'LineWidth', 2);

% Add labels and title
xlabel('Frequency (THz)');
ylabel('Magnitude');
title('S-Parameter 1');

% Display legend
legend;

% Hold off to stop adding to the current plot
hold off;

% Extract the frequency (x-axis) and the real/imaginary parts
freq2 = data2(:, 1);               % Frequency in THz
real2_S11 = data2(:, 2);           % Real part of S11
imag2_S11 = data2(:, 3);           % Imaginary part of S11
real2_S21 = data2(:, 4);           % Real part of S21
imag2_S21 = data2(:, 5);           % Imaginary part of S21
real2_S31 = data2(:, 6);           % Real part of S31
imag2_S31 = data2(:, 7);           % Imaginary part of S31

% Create the figure
figure(2)

% Plot all the variables on the same plot
hold on;
plot(freq2, real2_S11, 'DisplayName', 'Real(S11)', 'LineWidth', 2);
plot(freq2, imag2_S11, 'DisplayName', 'Imag(S11)', 'LineWidth', 2);
plot(freq2, real2_S21, 'DisplayName', 'Real(S21)', 'LineWidth', 2);
plot(freq2, imag2_S21, 'DisplayName', 'Imag(S21)', 'LineWidth', 2);
plot(freq2, real2_S31, 'DisplayName', 'Real(S31)', 'LineWidth', 2);
plot(freq2, imag2_S31, 'DisplayName', 'Imag(S31)', 'LineWidth', 2);

% Add labels and title
xlabel('Frequency (THz)');
ylabel('Magnitude');
title('S-Parameter 2');

% Display legend
legend;

% Hold off to stop adding to the current plot
hold off;
