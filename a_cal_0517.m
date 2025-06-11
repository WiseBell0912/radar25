clear; close all; clc;

%% Search
%file_path = 'E:/png2019/10/';
%file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
file_path = '/media/zerog/EA168F94168F6085/png2020/';

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
Nx = 201;
Ny = 201;
Nt = 128;

dx = 3;
dy = 3;
dt = 1.43;

Lx = 600;
Ly = 600;
Lt = Nt * dt;

g = 9.81;
h = 30;

%% Frequency
dkx = 2 * pi / (256 * dx);
dky = 2 * pi / (256 * dy);
dw  = 2 * pi / (Nt * dt);

kx = (-floor(256/2) : ceil(256/2)-1) * dkx;
ky = (-floor(256/2) : ceil(256/2)-1) * dky;
w  = (-floor(Nt/2) : ceil(Nt/2)-1) * dw;

[Kx, Ky, W] = meshgrid(ky, kx, w);  % (ky, kx, w) 순서 주의
Kx = single(Kx);
Ky = single(Ky);
W  = single(W);
K  = sqrt(Kx.^2 + Ky.^2);

%% Windowing
window_xy = hann(Nx) * hann(Ny)';      % Nx x Ny
window_t  = hann(Nt);                  % Nt
window    = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% HPF
hpK = (2*pi/15 <= K) & (K <= 2*pi/5);
hpW = (2*pi/30 <= abs(W));
hpMask = hpK & hpW;

%% Save
nFile = length(file_list);

r_Date  = NaT(nFile,1);

r_surf_SNR   = zeros(nFile,1);
r_wave_SNR   = zeros(nFile,1);
r_surf_Current_num = zeros(nFile,1);
r_wave_Current_num = zeros(nFile,1);

%% Loop
tic
parfor i = 1 : nFile
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

    % ---------- Edit image ---------- %
    % surf
    img_surf = single(png_surf);
    img_surf = rot90(flip(img_surf));
    img_surf_win = img_surf .* window;
    % wave
    img_wave = single(png_wave);
    img_wave = rot90(flip(img_wave));
    img_wave_win = img_wave .* window;

    % ---------- FFT ---------- %
    % surf
    spectrum_surf = fftn(img_surf_win, [2^nextpow2(Nx) 2^nextpow2(Ny) 2^nextpow2(Nt)]);
    spectrum_surf = fftshift(spectrum_surf);
    spectrum_surf = abs(spectrum_surf).^2;
    spectrum_surf = spectrum_surf .* hpMask;
    % wave
    spectrum_wave = fftn(img_wave_win, [2^nextpow2(Nx) 2^nextpow2(Ny) 2^nextpow2(Nt)]);
    spectrum_wave = fftshift(spectrum_wave);
    spectrum_wave = abs(spectrum_wave).^2;
    spectrum_wave = spectrum_wave .* hpMask;

    % ---------- Spectrum ---------- %
    % surf
    [~, idx_max_surf] = max(spectrum_surf(:));
    [ky_idx_surf, kx_idx_surf, w_idx_surf] = ind2sub(size(spectrum_surf), idx_max_surf);
    % wave
    [~, idx_max_wave] = max(spectrum_wave(:));
    [ky_idx_wave, kx_idx_wave, w_idx_wave] = ind2sub(size(spectrum_wave), idx_max_wave);
    % view
    % for j = 1 : 1
    %     tiledlayout(2, 2);
    %     nexttile;
    %     imagesc(flip(spectrum_surf(:, :, w_idx_surf), 2));
    %     view(0, 90); axis tight; axis equal;
    %     nexttile;
    %     imagesc(img_surf(:, :, j));
    %     view(0, 90); axis tight; axis equal;
    %     nexttile;
    %     imagesc(flip(spectrum_wave(:, :, w_idx_wave), 2));
    %     view(0, 90); axis tight; axis equal;
    %     nexttile;
    %     imagesc(img_wave(:, :, j));
    %     view(0, 90); axis tight; axis equal;
    %     pause(1.43);
    % end

    % ---------- Current estimation ---------- %
    % surf
    threshold_surf        = 0.2 * max(spectrum_surf, [], 'all');
    filter_surf           = (spectrum_surf >= threshold_surf);

    current_num_surf = sum(filter_surf(:));

    current_Kx    = Kx(filter_surf);
    current_Ky    = Ky(filter_surf);
    current_K     = K(filter_surf);
    current_W     = W(filter_surf);
    current_sigma = sqrt(g .* current_K .* tanh(current_K .* h));

    a11 = sum( current_Kx.^2 );
    a12 = sum( current_Kx .* current_Ky );
    a21 = sum( current_Kx .* current_Ky );
    a22 = sum( current_Ky.^2 );

    b1 = sum( (current_W - current_sigma) .* current_Kx );
    b2 = sum( (current_W - current_sigma) .* current_Ky );

    a = [a11, a12; a21, a22];
    b = [b1; b2];

    U = a \ b;

    ux_surf = -U(2);
    uy_surf = -U(1);

    % wave
    threshold_wave        = 0.2 * max(spectrum_wave, [], 'all');
    filter_wave           = (spectrum_wave >= threshold_wave);

    current_num_wave = sum(filter_wave(:));

    current_Kx    = Kx(filter_wave);
    current_Ky    = Ky(filter_wave);
    current_K     = K(filter_wave);
    current_W     = W(filter_wave);
    current_sigma = sqrt(g .* current_K .* tanh(current_K .* h));

    a11 = sum( current_Kx.^2 );
    a12 = sum( current_Kx .* current_Ky );
    a21 = sum( current_Kx .* current_Ky );
    a22 = sum( current_Ky.^2 );

    b1 = sum( (current_W - current_sigma) .* current_Kx );
    b2 = sum( (current_W - current_sigma) .* current_Ky );

    a = [a11, a12; a21, a22];
    b = [b1; b2];

    U = a \ b;

    ux_wave = -U(2);
    uy_wave = -U(1);

    % ---------- BPF ---------- %
    % surf
    dispersion_surf = sqrt(g .* K .* tanh(K .* h)) + Kx .* ux_surf + Ky .* uy_surf;
    safe_surf = 4 * abs(w(1) - w(2));
    bpMask_surf = (abs(dispersion_surf - W) <= safe_surf);

    spectrum_surf_bp = spectrum_surf .* bpMask_surf;

    % wave
    dispersion_wave = sqrt(g .* K .* tanh(K .* h)) + Kx .* ux_wave + Ky .* uy_wave;
    safe_wave = 4 * abs(w(1) - w(2));
    bpMask_wave = (abs(dispersion_wave - W) <= safe_wave);

    spectrum_wave_bp = spectrum_wave .* bpMask_wave;

    % ---------- SNR ---------- %
    % surf
    signal_surf = sum(spectrum_surf_bp(:));
    noise_surf = sum(spectrum_surf(:)) - signal_surf;
    snr_surf = signal_surf / noise_surf;
    % surf
    signal_wave = sum(spectrum_wave_bp(:));
    noise_wave = sum(spectrum_wave(:)) - signal_wave;
    snr_wave = signal_wave / noise_wave;

    % ---------- Save ---------- %
    r_Date(i)  = Date_Val;

    r_surf_SNR(i)   = snr_surf;
    r_wave_SNR(i)   = snr_wave;

    r_surf_Current_num(i) = current_num_surf;
    r_wave_Current_num(i) = current_num_wave;

    disp([nFile, i]);

end
toc;

save("snr_y2020_0517.mat", 'r_Date', 'r_surf_SNR', 'r_wave_SNR', 'r_surf_Current_num', 'r_wave_Current_num');