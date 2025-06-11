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
h = 50;

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
sc        = Nt * Nx * Ny / sum(window(:).^2);

%% HPF
hpK    = (K > 0.03);
hpW    = (abs(W) > 0.35);
hpMask = hpK & hpW;

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

r_surf_Entrophy = zeros(nFile,1);
r_wave_Entrophy = zeros(nFile,1);
r_surf_Flat  = zeros(nFile,1);
r_wave_Flat  = zeros(nFile,1);
r_surf_Std   = zeros(nFile,1);
r_wave_Std   = zeros(nFile,1);

%% 미리 1D화된 변수 준비(현재추정용)
Kx_1D    = reshape(Kx, 1, []);
Ky_1D    = reshape(Ky, 1, []);
K_1D     = sqrt(Kx_1D.^2 + Ky_1D.^2);
sigma1d  = sqrt(g * K_1D .* tanh(K_1D * h));

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
    img_surf = rot90(flip(flip(img_surf), 2));
    % wave
    img_wave = single(png_wave);
    img_wave = rot90(flip(flip(img_wave), 2));

    % ---------- Window ---------- %
    % surf
    img_surf_win = img_surf .* window;
    % wave
    img_wave_win = img_wave .* window;

    % ---------- FFT ---------- %
    % surf
    spectrum_surf_normal_win = fftshift(fftn(img_surf_win));
    spectrum_surf_normal_non = fftshift(fftn(img_surf));
    % wave
    spectrum_wave_normal_win = fftshift(fftn(img_wave_win));
    spectrum_wave_normal_non = fftshift(fftn(img_wave));

    % ---------- Threshold ---------- %
    % surf
    tmpPower_surf   = sc * abs(spectrum_surf_normal_win).^2 / (Nt^2 * Nx^2 * Ny^2);
    globalMax_surf  = max(tmpPower_surf(:));
    % wave
    tmpPower_wave   = sc * abs(spectrum_wave_normal_win).^2 / (Nt^2 * Nx^2 * Ny^2);
    globalMax_wave  = max(tmpPower_wave(:));

    % ---------- HPF ---------- %
    % surf
    spectrum_surf_normal_win_hp = spectrum_surf_normal_win .* hpMask;
    spectrum_surf_normal_non_hp = spectrum_surf_normal_non .* hpMask;
    spectrum_surf_power_win_hp = abs(spectrum_surf_normal_win_hp).^2 / (Nt^2* Nx^2* Ny^2);
    spectrum_surf_power_non_hp = abs(spectrum_surf_normal_non_hp).^2 / (Nt^2* Nx^2* Ny^2);
    % wave
    spectrum_wave_normal_win_hp = spectrum_wave_normal_win .* hpMask;
    spectrum_wave_normal_non_hp = spectrum_wave_normal_non .* hpMask;
    spectrum_wave_power_win_hp = abs(spectrum_wave_normal_win_hp).^2 / (Nt^2* Nx^2* Ny^2);
    spectrum_wave_power_non_hp = abs(spectrum_wave_normal_non_hp).^2 / (Nt^2* Nx^2* Ny^2);

    % ---------- Current estimation ---------- %
    % surf
    var   = zeros(1,6,'single');
    idxVar = 0;

    threshold = 1e-8 * globalMax_surf;

    for ii_ = (Nt/2) : Nt
        pt_1D = reshape(spectrum_surf_power_win_hp(:,:,ii_), 1, []);
        idx   = (pt_1D > threshold);

        Npt = sum(idx);
        if Npt>0
            oldEnd = idxVar;
            newEnd = idxVar + Npt;
            % 필요하면 동적 할당 대비
            if newEnd > size(var,1)
                var(end + 2*Npt, :) = 0;
            end

            var(oldEnd+1:newEnd,1) = w(ii_);
            var(oldEnd+1:newEnd,2) = sigma1d(idx);
            var(oldEnd+1:newEnd,3) = Kx_1D(idx);
            var(oldEnd+1:newEnd,4) = Ky_1D(idx);
            var(oldEnd+1:newEnd,5) = K_1D(idx);

            MTF1d = sqrt(K_1D(idx).^(-1.2) .* abs(w(ii_)).^(-0.6));
            var(oldEnd+1:newEnd,6) = pt_1D(idx) .* MTF1d;

            idxVar = newEnd;
        end
    end

    var = var(1:idxVar, :);

    % w<0인 경우 부호 반전
    idxNeg = (var(:,1) < 0);
    var(idxNeg,2) = -var(idxNeg,2);

    a1 = sum(var(:,6).*var(:,3).^2);
    b1 = sum(var(:,6).*var(:,3).*var(:,4));
    c1 = sum(var(:,6).*(var(:,1) - var(:,2)).*var(:,3));
    a2 = b1;
    b2 = sum(var(:,6).*var(:,4).^2);
    c2 = sum(var(:,6).*(var(:,1) - var(:,2)).*var(:,4));

    denom = (a1*b2 - a2*b1);
    if abs(denom) > 1e-12
        ux_cal = (c1*b2 - c2*b1) / denom;
        uy_cal = (c1*a2 - c2*a1) / (-denom);
    else
        ux_cal = 0;
        uy_cal = 0;
        disp('error or 1d wave');
    end

    ux_surf = ux_cal;
    uy_surf = uy_cal;

    % wave
    var   = zeros(1,6,'single');
    idxVar = 0;

    threshold = 1e-8 * globalMax_wave;

    for ii_ = (Nt/2) : Nt
        pt_1D = reshape(spectrum_wave_power_win_hp(:,:,ii_), 1, []);
        idx   = (pt_1D > threshold);

        Npt = sum(idx);
        if Npt>0
            oldEnd = idxVar;
            newEnd = idxVar + Npt;
            % 필요하면 동적 할당 대비
            if newEnd > size(var,1)
                var(end + 2*Npt, :) = 0;
            end

            var(oldEnd+1:newEnd,1) = w(ii_);
            var(oldEnd+1:newEnd,2) = sigma1d(idx);
            var(oldEnd+1:newEnd,3) = Kx_1D(idx);
            var(oldEnd+1:newEnd,4) = Ky_1D(idx);
            var(oldEnd+1:newEnd,5) = K_1D(idx);

            MTF1d = sqrt(K_1D(idx).^(-1.2) .* abs(w(ii_)).^(-0.6));
            var(oldEnd+1:newEnd,6) = pt_1D(idx) .* MTF1d;

            idxVar = newEnd;
        end
    end

    var = var(1:idxVar, :);

    % w<0인 경우 부호 반전
    idxNeg = (var(:,1) < 0);
    var(idxNeg,2) = -var(idxNeg,2);

    a1 = sum(var(:,6).*var(:,3).^2);
    b1 = sum(var(:,6).*var(:,3).*var(:,4));
    c1 = sum(var(:,6).*(var(:,1) - var(:,2)).*var(:,3));
    a2 = b1;
    b2 = sum(var(:,6).*var(:,4).^2);
    c2 = sum(var(:,6).*(var(:,1) - var(:,2)).*var(:,4));

    denom = (a1*b2 - a2*b1);
    if abs(denom) > 1e-12
        ux_cal = (c1*b2 - c2*b1) / denom;
        uy_cal = (c1*a2 - c2*a1) / (-denom);
    else
        ux_cal = 0;
        uy_cal = 0;
        disp('error or 1d wave');
    end

    ux_wave = ux_cal;
    uy_wave = uy_cal;

    % ---------- BPF ---------- %
    % surf
    wk_pos =  sqrt( g .* K .* tanh(K .* h) ) + Kx .* ux_surf + Ky .* uy_surf;
    wk_neg = -sqrt( g .* K .* tanh(K .* h) ) + Kx .* ux_surf + Ky .* uy_surf;
    bpv = 2 * (Nt/32) * 2*pi/Lt;
    bpMask = ( ( ((wk_pos - bpv) < W) & (W < (wk_pos + bpv)) ) | ( ((wk_neg - bpv) < W) & (W < (wk_neg + bpv)) ) );
    spectrum_surf_normal_win_bp = spectrum_surf_normal_win_hp .* bpMask;
    spectrum_surf_normal_non_bp = spectrum_surf_normal_non_hp .* bpMask;
    spectrum_surf_power_win_bp = spectrum_surf_power_win_hp .* bpMask;
    spectrum_surf_power_non_bp = spectrum_surf_power_non_hp .* bpMask;
    % wave
    wk_pos =  sqrt( g .* K .* tanh(K .* h) ) + Kx .* ux_wave + Ky .* uy_wave;
    wk_neg = -sqrt( g .* K .* tanh(K .* h) ) + Kx .* ux_wave + Ky .* uy_wave;
    bpv = 2 * (Nt/32) * 2*pi/Lt;
    bpMask = ( ( ((wk_pos - bpv) < W) & (W < (wk_pos + bpv)) ) | ( ((wk_neg - bpv) < W) & (W < (wk_neg + bpv)) ) );
    spectrum_wave_normal_win_bp = spectrum_wave_normal_win_hp .* bpMask;
    spectrum_wave_normal_non_bp = spectrum_wave_normal_non_hp .* bpMask;
    spectrum_wave_power_win_bp = spectrum_wave_power_win_hp .* bpMask;
    spectrum_wave_power_non_bp = spectrum_wave_power_non_hp .* bpMask;

    % ---------- SNR ---------- %
    % surf
    signal = sum(spectrum_surf_power_non_bp(:));
    noise  = sum(spectrum_surf_power_non_hp(:)) - signal;
    surf_SNR_Val = signal / max(noise,1e-12);
    % wave
    signal = sum(spectrum_wave_power_non_bp(:));
    noise  = sum(spectrum_wave_power_non_hp(:)) - signal;
    wave_SNR_Val = signal / max(noise,1e-12);

    % ---------- Tp ---------- %
    % surf
    spectrum_surf_w_non = 2 .* squeeze(sum(sum(spectrum_surf_power_non_bp, 1), 2));
    spectrum_surf_w_non = spectrum_surf_w_non .* (w > 0);
    [~, sorted_idx_surf] = sort(spectrum_surf_w_non, 'descend');
    surf_Tp_Val = mean(2 * pi ./ w(sorted_idx_surf(1:2)));
    % wave
    spectrum_wave_w_non = 2 .* squeeze(sum(sum(spectrum_wave_power_non_bp, 1), 2));
    spectrum_wave_w_non = spectrum_wave_w_non .* (w > 0);
    [~, sorted_idx_wave] = sort(spectrum_wave_w_non, 'descend');
    wave_Tp_Val = mean(2 * pi ./ w(sorted_idx_wave(1:2)));

    % ---------- Pdir ---------- %
    % surf
    spcetrum_surf_k_non = sum(spectrum_surf_power_non_bp(:, :, end/2:end), 3);
    [~, sorted_idx_surf] = sort(spcetrum_surf_k_non(:), 'descend');
    surf_Pdir_Val = mod(mean(rad2deg(atan2(Ky(sorted_idx_surf(1:5)), Kx(sorted_idx_surf(1:5))))), 360);
    % wave
    spcetrum_wave_k_non = sum(spectrum_wave_power_non_bp(:, :, end/2:end), 3);
    [~, sorted_idx_wave] = sort(spcetrum_wave_k_non(:), 'descend');
    wave_Pdir_Val = mod(mean(rad2deg(atan2(Ky(sorted_idx_wave(1:5)), Kx(sorted_idx_wave(1:5))))), 360);

    % ---------- Entropy ---------- %
    % surf
    spectrum_surf_power_win_hp_norm = spectrum_surf_power_win_hp / sum(spectrum_surf_power_win_hp(:));
    surf_Entrophy_Val = - sum(spectrum_surf_power_win_hp_norm(:) .* log(spectrum_surf_power_win_hp_norm(:) + eps));
    % wave
    spectrum_wave_power_win_hp_norm = spectrum_wave_power_win_hp / sum(spectrum_wave_power_win_hp(:));
    wave_Entrophy_Val = - sum(spectrum_wave_power_win_hp_norm(:) .* log(spectrum_wave_power_win_hp_norm(:) + eps));

    % ---------- Flat rate ---------- %
    % surf
    spectrum_surf_power_win_hp_norm = spectrum_surf_power_win_hp / sum(spectrum_surf_power_win_hp(:));
    surf_max = max(spectrum_surf_power_win_hp_norm(:));
    surf_mean = mean(spectrum_surf_power_win_hp_norm(:));
    surf_Flat_Val = surf_max / surf_mean;
    % wave
    spectrum_wave_power_win_hp_norm = spectrum_wave_power_win_hp / sum(spectrum_wave_power_win_hp(:));
    wave_max = max(spectrum_wave_power_win_hp_norm(:));
    wave_mean = mean(spectrum_wave_power_win_hp_norm(:));
    wave_Flat_Val = wave_max / wave_mean;

    % ---------- Std ---------- %
    % surf
    spectrum_surf_power_win_hp_norm = spectrum_surf_power_win_hp / sum(spectrum_surf_power_win_hp(:));
    surf_Std_Val = std(spectrum_surf_power_win_hp_norm(:));
    % wave
    spectrum_wave_power_win_hp_norm = spectrum_wave_power_win_hp / sum(spectrum_wave_power_win_hp(:));
    wave_Std_Val = std(spectrum_wave_power_win_hp_norm(:));

    % ---------- Save ---------- %
    r_Date(i)  = Date_Val;

    r_LandE(i) = LandE_Val;

    r_surf_SNR(i)   = surf_SNR_Val;
    r_surf_Ux(i)    = ux_surf;
    r_surf_Uy(i)    = uy_surf;
    r_surf_Tp(i)    = surf_Tp_Val;
    r_surf_Pdir(i)  = surf_Pdir_Val;

    r_wave_SNR(i)   = wave_SNR_Val;
    r_wave_Ux(i)    = ux_wave;
    r_wave_Uy(i)    = uy_wave;
    r_wave_Tp(i)    = wave_Tp_Val;
    r_wave_Pdir(i)  = wave_Pdir_Val;

    r_surf_Entrophy(i) = surf_Entrophy_Val;
    r_wave_Entrophy(i) = wave_Entrophy_Val;
    r_surf_Flat(i) = surf_Flat_Val;
    r_wave_Flat(i) = wave_Flat_Val;
    r_surf_Std(i) = surf_Std_Val;
    r_wave_Std(i) = wave_Std_Val;

    disp([nFile, i]);
end
toc;

%% Save
save( ...
    "snr_y2020_0513.mat", ...
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
    'r_surf_Entrophy', ...
    'r_wave_Entrophy', ...
    'r_surf_Flat', ...
    'r_wave_Flat', ...
    'r_surf_Std', ...
    'r_wave_Std' ...
    );