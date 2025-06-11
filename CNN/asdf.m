clear; clc;

%% Search
%file_path = 'E:/png2019/10/';
%file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
file_path = '/media/zerog/EA168F94168F6085/png2019/';

%file_list = dir([file_path, '*.png']);
file_list = dir(fullfile(file_path, '*_201910*.png'));

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
kx = -pi/dx : 2*pi/Lx : pi/dx;
ky = -pi/dy : 2*pi/Ly : pi/dy;
w  = -pi/dt : 2*pi/Lt : pi/dt;
w(65) = [];

[Kx, Ky, W] = meshgrid(ky, kx, w);
Kx = single(Kx);
Ky = single(Ky);
W  = single(W);
K  = sqrt(Kx.^2 + Ky.^2);

%% Windowing
window_xy = hann(Nx) * hann(Ny)';      % Nx x Ny
window_t  = hann(Nt);                  % Nt
window    = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);

%% Save
nFile = length(file_list);

r_Date  = NaT(nFile,1);

r_LandE   = zeros(nFile,1);

r_surf_SNR   = zeros(nFile,1);
r_surf_Ux    = zeros(nFile,1);
r_surf_Uy    = zeros(nFile,1);
r_surf_Tp    = zeros(nFile,1);
r_surf_Pdir  = zeros(nFile,1);

r_wave_SNR   = zeros(nFile,1);
r_wave_Ux    = zeros(nFile,1);
r_wave_Uy    = zeros(nFile,1);
r_wave_Tp    = zeros(nFile,1);
r_wave_Pdir  = zeros(nFile,1);

r_surf_Singal = zeros(nFile,1);
r_wave_Signal = zeros(nFile,1);
r_surf_Noise  = zeros(nFile,1);
r_wave_Noise  = zeros(nFile,1);
r_surf_Std   = zeros(nFile,1);
r_wave_Std   = zeros(nFile,1);

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

    % ---------- Land energy ---------- %
    LandE_Val = sum(png_long(1:512, 355:385, :), 'all') + sum(png_long(1:512, 875:935, :), 'all');

    % ---------- Edit image ---------- %
    % surf
    img_surf = single(png_surf);
    img_surf = rot90(flip(img_surf));
    % wave
    img_wave = single(png_wave);
    img_wave = rot90(flip(img_wave));

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
    spectrum_surf_power_win = abs(spectrum_surf_normal_win).^2;
    spectrum_surf_power_non = abs(spectrum_surf_normal_non).^2;
    % wave
    spectrum_wave_power_win = abs(spectrum_wave_normal_win).^2;
    spectrum_wave_power_non = abs(spectrum_wave_normal_non).^2;

    % ---------- HPF ---------- %
    hpK = (2*pi/150 <= K) & (K <= 2*pi/5);
    hpW = (2*pi/25 <= abs(W));
    hpMask = hpK & hpW;
    % surf
    spectrum_surf_power_win_hp = spectrum_surf_power_win .* hpMask;
    spectrum_surf_power_non_hp = spectrum_surf_power_non .* hpMask;
    % wave
    spectrum_wave_power_win_hp = spectrum_wave_power_win .* hpMask;
    spectrum_wave_power_non_hp = spectrum_wave_power_non .* hpMask;

    % ---------- Current estimate ---------- %
    % % % mtf = K.^(-1.2) .* W.^(-0.6);
    % % % sigma = sqrt(g .* K .* tanh(K .* h));
    % % % % surf
    % % % threshold_surf = 0.2 * max(spectrum_surf_power_win_hp(:));
    % % % filter_surf = (spectrum_surf_power_win_hp > threshold_surf);
    % % % surf11 = sum( filter_surf .* spectrum_surf_power_win_hp .* mtf .* Kx.^2, "all", 'omitnan' );
    % % % surf12 = sum( filter_surf .* spectrum_surf_power_win_hp .* mtf .* Kx .* Ky, "all", 'omitnan' );
    % % % surf21 = sum( filter_surf .* spectrum_surf_power_win_hp .* mtf .* Ky .* Kx, "all", 'omitnan' );
    % % % surf22 = sum( filter_surf .* spectrum_surf_power_win_hp .* mtf .* Ky.^2, "all", 'omitnan' );
    % % % surf1  = sum( filter_surf .* spectrum_surf_power_win_hp .* mtf .* (W - sigma) .* Kx, "all", 'omitnan' );
    % % % surf2  = sum( filter_surf .* spectrum_surf_power_win_hp .* mtf .* (W - sigma) .* Ky, "all", 'omitnan' );
    % % % A = [surf11, surf12; surf21, surf22];
    % % % b = [surf1; surf2];
    % % % U = A \ b;
    % % % ux_surf = U(1);
    % % % uy_surf = U(2);
    % % % % wave
    % % % threshold_wave = 0.2 * max(spectrum_wave_power_win_hp(:));
    % % % filter_wave = (spectrum_wave_power_win_hp > threshold_wave);
    % % % wave11 = sum( filter_wave .* spectrum_wave_power_win_hp .* mtf .* Kx.^2, "all", 'omitnan' );
    % % % wave12 = sum( filter_wave .* spectrum_wave_power_win_hp .* mtf .* Kx .* Ky, "all", 'omitnan' );
    % % % wave21 = sum( filter_wave .* spectrum_wave_power_win_hp .* mtf .* Ky .* Kx, "all", 'omitnan' );
    % % % wave22 = sum( filter_wave .* spectrum_wave_power_win_hp .* mtf .* Ky.^2, "all", 'omitnan' );
    % % % wave1  = sum( filter_wave .* spectrum_wave_power_win_hp .* mtf .* (W - sigma) .* Kx, "all", 'omitnan' );
    % % % wave2  = sum( filter_wave .* spectrum_wave_power_win_hp .* mtf .* (W - sigma) .* Ky, "all", 'omitnan' );
    % % % A = [wave11, wave12; wave21, wave22];
    % % % b = [wave1; wave2];
    % % % U = A \ b;
    % % % ux_wave = U(1);
    % % % uy_wave = U(2);

    % ---------- Dispersion shell ---------- %
    % figure(2);
    % surf(Kx(:, :, 1), Ky(:, :, 1), sqrt(K(:, :, 1) .* g .* tanh(K(:, :, 1) .* h)) + Kx(:, :, 1).*ux_surf + Ky(:, :, 1).*uy_surf);
    % xlabel('Kx');
    % ylabel('Ky');
    % zlim([min(w), max(w)]);

    % ---------- BPF ---------- %
    dw = abs(w(2)-w(1));
    % surf
    dispersion_surf = sqrt(g .* K .* tanh(K .* h));% + Kx.*ux_surf + Ky.*uy_surf;
    bpMask_surf = ( abs( dispersion_surf - abs(W) ) < 4 * dw );
    spectrum_surf_power_win_bp = spectrum_surf_power_win_hp .* bpMask_surf;
    spectrum_surf_power_non_bp = spectrum_surf_power_non_hp .* bpMask_surf;
    % wave
    dispersion_wave = sqrt(g .* K .* tanh(K .* h));% + Kx.*ux_wave + Ky.*uy_wave;
    bpMask_wave = ( abs( dispersion_wave - abs(W) ) < 4 * dw );
    spectrum_wave_power_win_bp = spectrum_wave_power_win_hp .* bpMask_wave;
    spectrum_wave_power_non_bp = spectrum_wave_power_non_hp .* bpMask_wave;

    % ---------- Find max ---------- %
    % surf
    [~, idx_max_surf] = max(spectrum_surf_power_win_bp(:));
    [ky_idx_surf, kx_idx_surf, w_idx_surf] = ind2sub(size(spectrum_surf_power_win_bp), idx_max_surf);
    % wave
    [~, idx_max_wave] = max(spectrum_wave_power_win_bp(:));
    [ky_idx_wave, kx_idx_wave, w_idx_wave] = ind2sub(size(spectrum_wave_power_win_bp), idx_max_wave);

    % ---------- SNR ---------- %
    % surf
    surf_signal = sum(spectrum_surf_power_win_bp(:));
    surf_noise = sum(spectrum_surf_power_win_hp(:)) - surf_signal;
    surf_SNR_Val = surf_signal / surf_noise;
    % wave
    wave_signal = sum(spectrum_wave_power_win_bp(:));
    wave_noise = sum(spectrum_wave_power_win_hp(:)) - wave_signal;
    wave_SNR_Val = wave_signal / wave_noise;

    % ---------- Pdir ---------- %
    Kx2D = Kx(:, :, 1);
    Ky2D = Ky(:, :, 1);
    dir_deg = mod(mod(rad2deg(atan2(Ky2D, Kx2D)) - 180, 360) .* hpK(:, :, 1), 180);
    % surf
    spectrum_surf_Pdir = spectrum_surf_power_win_bp;
    spectrum_surf_Pdir = sum(spectrum_surf_Pdir, 3);
    [~, idx_sort_surf_Pdir] = sort(spectrum_surf_Pdir(:), 'descend');
    surf_Pdir_Val = mean( dir_deg(idx_sort_surf_Pdir(1:6)) , 'omitnan');
    % wave
    spectrum_wave_Pdir = spectrum_wave_power_win_bp;
    spectrum_wave_Pdir = sum(spectrum_wave_Pdir, 3);
    [~, idx_sort_wave_Pdir] = sort(spectrum_wave_Pdir(:), 'descend');
    wave_Pdir_Val = mean( dir_deg(idx_sort_wave_Pdir(1:6)) , 'omitnan');

    % ---------- Tp ---------- %
    % surf
    spectrum_surf_Tp = spectrum_surf_power_win_bp;
    spectrum_surf_Tp = squeeze(sum(sum(spectrum_surf_Tp, 1), 2));
    [~, idx_sort_surf_Tp] = sort(spectrum_surf_Tp, 'descend');
    surf_Tp_Val = mean(abs(2*pi ./ w(idx_sort_surf_Tp(1:6))), 'omitnan');
    % wave
    spectrum_wave_Tp = spectrum_wave_power_win_bp;
    spectrum_wave_Tp = squeeze(sum(sum(spectrum_wave_Tp, 1), 2));
    [~, idx_sort_wave_Tp] = sort(spectrum_wave_Tp, 'descend');
    wave_Tp_Val = mean(abs(2*pi ./ w(idx_sort_wave_Tp(1:6))), 'omitnan');

    % ---------- Check ---------- %
    % % % for j = 1 : 20
    % % %     figure(1);
    % % %     tiledlayout(3, 2);
    % % %     nexttile;
    % % %     imagesc(0:3:600, 0:3:600, img_surf(:, :, j));
    % % %     hold on;
    % % %     quiver(300, 300, 100*cos(deg2rad(surf_Pdir_Val)), 100*sin(deg2rad(surf_Pdir_Val)), 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    % % %     %quiver(300, 300, -ux_surf, -uy_surf, 10, 'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 2);
    % % %     hold off;
    % % %     axis equal; axis tight;
    % % %     nexttile;
    % % %     imagesc(0:3:600, 0:3:600, img_wave(:, :, j));
    % % %     hold on;
    % % %     quiver(300, 300, 100*cos(deg2rad(wave_Pdir_Val)), 100*sin(deg2rad(wave_Pdir_Val)), 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    % % %     %quiver(300, 300, -ux_wave, -uy_wave, 10, 'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 2);
    % % %     hold off;
    % % %     axis equal; axis tight;
    % % %     nexttile;
    % % %     imagesc(kx, ky, sum(spectrum_surf_power_win_hp, 3));
    % % %     title([num2str(2*pi/K(idx_max_surf)), ' [m] ', num2str(surf_Pdir_Val)]);
    % % %     hold on;
    % % %     %quiver(0, 0, ux_surf, uy_surf, 0.05, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    % % %     hold off;
    % % %     axis equal; axis tight;
    % % %     nexttile;
    % % %     imagesc(kx, ky, sum(spectrum_wave_power_win_hp, 3));
    % % %     title([num2str(2*pi/K(idx_max_wave)), ' [m] ', num2str(wave_Pdir_Val)]);
    % % %     hold on;
    % % %     %quiver(0, 0, ux_wave, uy_wave, 0.05, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    % % %     hold off;
    % % %     axis equal; axis tight;
    % % %     nexttile;
    % % %     imagesc(kx, ky, sum(spectrum_surf_power_win_bp, 3));
    % % %     title([num2str(2*pi/K(idx_max_surf)), ' [m] ', num2str(surf_Pdir_Val)]);
    % % %     axis equal; axis tight;
    % % %     nexttile;
    % % %     imagesc(kx, ky, sum(spectrum_wave_power_win_bp, 3));
    % % %     title([num2str(2*pi/K(idx_max_wave)), ' [m] ', num2str(wave_Pdir_Val)]);
    % % %     axis equal; axis tight;
    % % %     pause(0.1);
    % % % end

    % ---------- Save ---------- %
    r_Date(i)  = Date_Val;

    r_LandE(i) = LandE_Val;

    r_surf_SNR(i)   = surf_SNR_Val;
    %r_surf_Ux(i)    = ux_surf;
    %r_surf_Uy(i)    = uy_surf;
    r_surf_Tp(i)    = surf_Tp_Val;
    r_surf_Pdir(i)  = surf_Pdir_Val;

    r_wave_SNR(i)   = wave_SNR_Val;
    %r_wave_Ux(i)    = ux_wave;
    %r_wave_Uy(i)    = uy_wave;
    r_wave_Tp(i)    = wave_Tp_Val;
    r_wave_Pdir(i)  = wave_Pdir_Val;

    r_surf_Singal(i) = surf_signal;
    r_wave_Signal(i) = wave_signal;
    r_surf_Noise(i) = surf_noise;
    r_wave_Noise(i) = wave_noise;

    disp([nFile, i]);
end
toc;

%% Save
save( ...
    "snr_y202010_0518.mat", ...
    'r_Date', ...
    'r_LandE', ...
    'r_surf_Pdir', ...
    'r_surf_SNR', ...
    'r_surf_Tp', ...
    'r_surf_Ux', ...
    'r_surf_Uy', ...
    'r_wave_Pdir', ...
    'r_wave_SNR', ...
    'r_wave_Tp', ...
    'r_wave_Ux', ...
    'r_wave_Uy', ...
    'r_surf_Singal', ...
    'r_wave_Signal', ...
    'r_surf_Noise', ...
    'r_wave_Noise'...
    );