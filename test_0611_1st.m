clear; clc;

%% Search
%file_path = 'E:/png2019/10/';
file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
%file_path = '/media/zerog/EA168F94168F6085/png2020/';

file_list = dir([file_path, '*.png']);

%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(160);
surf_center = 870 + 600/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 600/2;

%% Input parameter
Nx = 200;
Ny = 200;
Nt = 128;

dx = 3;
dy = 3;
dt = 1.43;

Lx = 600;
Ly = 600;
Lt = Nt * dt;

g = 9.81;
h = 30;

%% Windowing
window_xy = hann(Nx) * hann(Ny)';      % Nx x Ny
window_t  = hann(Nt);                  % Nt
window    = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% Save
nFile = length(file_list);

r_Date  = NaT(nFile,1);

r_LandE   = zeros(nFile,1);

r_surf_SNR   = zeros(nFile,1);
r_wave_SNR   = zeros(nFile,1);
r_surf_Singal = zeros(nFile,1);
r_wave_Signal = zeros(nFile,1);
r_surf_Noise  = zeros(nFile,1);
r_wave_Noise  = zeros(nFile,1);

r_surf_SNR1   = zeros(nFile,1);
r_wave_SNR1   = zeros(nFile,1);
r_surf_Signal1 = zeros(nFile,1);
r_wave_Signal1 = zeros(nFile,1);
r_surf_Noise1  = zeros(nFile,1);
r_wave_Noise1  = zeros(nFile,1);

%% Loop
tic
for i = 1 : nFile
    % ---------- Read png ---------- %
    png_path = fullfile(file_list(i).folder, file_list(i).name);
    dateStr  = file_list(i).name(5:end-4);
    Date_Val  = datetime(dateStr, 'InputFormat', 'yyyyMMdd_HHmm');

    png_long = imread(png_path);
    png_long = png_long(1 : 512*1080*128);
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    [png_surf, png_wave] = f_extract_zone( ...
        png_long, modifiy_theta, ...
        surf_theta1, surf_theta2, surf_center, ...
        wave_theta1, wave_theta2, wave_center);

    % ---------- Land energy ---------- %
    LandE_Val = sum(png_long(1:512, 355:385, :), 'all') + sum(png_long(1:512, 875:935, :), 'all');

    % ---------- Edit image ---------- %
    % surf
    img_surf = single(png_surf);
    img_surf = rot90(flip(img_surf));
    % wave
    img_wave = single(png_wave);
    img_wave = rot90(flip(img_wave));

    % ---------- Zerocrossing ---------- %
    % surf
    mu_surf = mean(img_surf, 3);
    mu_surf = repmat(mu_surf, [1, 1, 128]);
    img_surf = img_surf - mu_surf;

    % wave
    mu_wave = mean(img_wave, 3);
    mu_wave = repmat(mu_wave, [1, 1, 128]);
    img_wave = img_wave - mu_wave;

    % ---------- Window ---------- %
    % surf
    img_surf_win = img_surf .* window;
    % wave
    img_wave_win = img_wave .* window;

    % ---------- 3D FFT ---------- %
    % surf
    spectrum_surf_normal_win = fftshift(fftn(img_surf_win));
    spectrum_surf_normal_non = fftshift(fftn(img_surf));
    % wave
    spectrum_wave_normal_win = fftshift(fftn(img_wave_win));
    spectrum_wave_normal_non = fftshift(fftn(img_wave));

    % ---------- Spectrum ---------- %
    % surf
    spectrum_surf_power_win = abs(spectrum_surf_normal_win).^2 / Nx^2 / Ny^2 / Nt^2;
    spectrum_surf_power_non = abs(spectrum_surf_normal_non).^2 / Nx^2 / Ny^2 / Nt^2;
    % wave
    spectrum_wave_power_win = abs(spectrum_wave_normal_win).^2 / Nx^2 / Ny^2 / Nt^2;
    spectrum_wave_power_non = abs(spectrum_wave_normal_non).^2 / Nx^2 / Ny^2 / Nt^2;

    % ---------- BPF ---------- %
    dw = abs(w(2)-w(1));
    % surf
    dispersion_surf = sqrt(g .* K .* tanh(K .* h));% + Kx.*ux_surf + Ky.*uy_surf;
    bpMask_surf = ( abs( dispersion_surf - abs(W) ) < 16 * dw );
    spectrum_surf_power_win_bp = spectrum_surf_power_win .* bpMask_surf;
    spectrum_surf_power_non_bp = spectrum_surf_power_non .* bpMask_surf;
    % wave
    dispersion_wave = sqrt(g .* K .* tanh(K .* h));% + Kx.*ux_wave + Ky.*uy_wave;
    bpMask_wave = ( abs( dispersion_wave - abs(W) ) < 16 * dw );
    spectrum_wave_power_win_bp = spectrum_wave_power_win .* bpMask_wave;
    spectrum_wave_power_non_bp = spectrum_wave_power_non .* bpMask_wave;

    % ---------- SNR ---------- %
    % surf
    surf_signal = sum(spectrum_surf_power_win_bp(:));
    surf_noise = sum(spectrum_surf_power_win(:)) - surf_signal;
    surf_SNR_Val = surf_signal / surf_noise;
    % wave
    wave_signal = sum(spectrum_wave_power_win_bp(:));
    wave_noise = sum(spectrum_wave_power_win(:)) - wave_signal;
    wave_SNR_Val = wave_signal / wave_noise;

    % ---------- SNR ---------- %
    % surf
    surf_signal1 = sum(spectrum_surf_power_non_bp(:));
    surf_noise1 = sum(spectrum_surf_power_non(:)) - surf_signal1;
    surf_SNR_Val1 = surf_signal1 / surf_noise1;
    % wave
    wave_signal1 = sum(spectrum_wave_power_non_bp(:));
    wave_noise1 = sum(spectrum_wave_power_non(:)) - wave_signal1;
    wave_SNR_Val1 = wave_signal1 / wave_noise1;

    % ---------- Save ---------- %
    r_Date(i)  = Date_Val;

    r_LandE(i) = LandE_Val;

    r_surf_SNR(i)   = surf_SNR_Val;
    r_wave_SNR(i)   = wave_SNR_Val;
    r_surf_Singal(i) = surf_signal;
    r_wave_Signal(i) = wave_signal;
    r_surf_Noise(i) = surf_noise;
    r_wave_Noise(i) = wave_noise;

    r_surf_SNR1(i)   = surf_SNR_Val1;
    r_wave_SNR1(i)   = wave_SNR_Val1;
    r_surf_Signal1(i) = surf_signal1;
    r_wave_Signal1(i) = wave_signal1;
    r_surf_Noise1(i) = surf_noise1;
    r_wave_Noise1(i) = wave_noise1;

    disp([nFile, i]);
end
toc;