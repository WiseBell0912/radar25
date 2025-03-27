clear; close all; clc;

%% Search
file_path = 'E:/png2019/10/';
file_list = dir([file_path, '*.png']);

%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(160);
surf_center = 870 + 600/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 600/2;

%% input parameter
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

%% 하이패스 필터(HPF)용 마스크 미리 계산(반복문 밖)
hpK = (K > 0.025);
hpW = (abs(W) > 2*pi*0.03);
hpMask = hpK & hpW;

%% 결과 저장용
nFile = length(file_list);

for i = 1 : nFile
    %% 파일 읽기 및 Zone 추출
    png_path = fullfile(file_list(i).folder, file_list(i).name);
    dateStr  = file_list(i).name(5:end-4);
    dateVal  = datetime(dateStr, 'InputFormat', 'yyyyMMdd_HHmm');

    png_long = imread(png_path);
    png_long = png_long(1 : 512*1080*128);
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    [png_surf, png_wave] = f_extract_zone( ...
        png_long, modifiy_theta, ...
        surf_theta1, surf_theta2, surf_center, ...
        wave_theta1, wave_theta2, wave_center);

    img_surf = single(png_surf);
    img_wave = single(png_wave);

    png_long = [];
    png_surf = [];
    png_wave = [];

    %% FFT: 윈도우 적용 후 FFT 및 HPF 적용
    img_surf_windowed = img_surf .* window;
    image_surf_spectrum = fftshift(fftn(img_surf_windowed));
    image_surf_spectrum = abs(image_surf_spectrum).^2 / Nx / Ny / Nt;

    img_wave_windowed = img_wave .* window;
    image_wave_spectrum = fftshift(fftn(img_wave_windowed));
    image_wave_spectrum = abs(image_wave_spectrum).^2 / Nx / Ny / Nt;

    %% HPF
    image_surf_spectrum = image_surf_spectrum .* hpMask;
    image_wave_spectrum = image_wave_spectrum .* hpMask;

    %% Current velocity estimation
    sigma = sqrt(g .* K .* tanh(K .* h));

    a1 = sum( image_surf_spectrum .* Kx.^2 , 'all');
    a2 = sum( image_surf_spectrum .* Ky.^2 , 'all');
    a3 = sum( image_surf_spectrum .* Kx .* Ky , 'all');
    a4 = sum( image_surf_spectrum .* (W - sigma) .* Kx.^2 , 'all');
    a5 = sum( image_surf_spectrum .* (W - sigma) .* Ky.^2 , 'all');

    b1 = sum( image_wave_spectrum .* Kx.^2 , 'all');
    b2 = sum( image_wave_spectrum .* Ky.^2 , 'all');
    b3 = sum( image_wave_spectrum .* Kx .* Ky , 'all');
    b4 = sum( image_wave_spectrum .* (W - sigma) .* Kx.^2 , 'all');
    b5 = sum( image_wave_spectrum .* (W - sigma) .* Ky.^2 , 'all');

    surf_ux = ( a2 * a4 - a3 * a5 ) / ( a1 * a2 - a3^2 );
    surf_uy = ( a1 * a5 - a3 * a4 ) / ( a1 * a2 - a3^2 );

    wave_ux = ( b2 * b4 - b3 * b5 ) / ( b1 * b2 - b3^2 );
    wave_uy = ( b1 * b5 - b3 * b4 ) / ( b1 * b2 - b3^2 );

    %% BPF
    bpv = 64 * (2*pi/Lt);

    bpMask_surf = (sqrt( g .* K .* tanh( K .* h ) ) + (Kx .* surf_ux) + (Ky .* surf_uy) - bpv <= abs(W)) & (abs(W) <= sqrt( g .* K .* tanh( K .* h ) ) + (Kx .* surf_ux) + (Ky .* surf_uy) + bpv);
    bpMask_wave = (sqrt( g .* K .* tanh( K .* h ) ) + (Kx .* wave_ux) + (Ky .* wave_uy) - bpv <= abs(W)) & (abs(W) <= sqrt( g .* K .* tanh( K .* h ) ) + (Kx .* wave_ux) + (Ky .* wave_uy) + bpv);

    image_surf_spectrum = image_surf_spectrum .* bpMask_surf;
    image_wave_spectrum = image_wave_spectrum .* bpMask_wave;

    %% 피규어 설정: 2x3 타일로 두 영역의 3 스펙트럼씩 표시
    figure(1);
    tiledlayout(3,3)
    set(gcf, 'Position', [0, 0, 1250, 1250]);

    %% Wave 영역 스펙트럼 분석
    % 1. 주파수 스펙트럼 S(ω) (Wave)
    S_omega_wave = squeeze(sum(sum(image_wave_spectrum, 1), 2));
    nexttile;
    plot(w, S_omega_wave);
    xlabel('Angular Frequency \omega (rad/s)');
    ylabel('Energy');
    title(['Wave: Frequency Spectrum for file ', num2str(i)]);

    % 2. 파수 스펙트럼 S(kx, ky) (Wave)
    S_k_wave = sum(image_wave_spectrum, 3);
    nexttile;
    imagesc(kx, ky, S_k_wave);
    xlabel('k_x (rad/m)');
    ylabel('k_y (rad/m)');
    title(['Wave: Wavenumber Spectrum for file ', num2str(i)]);
    colorbar;

    % 3. 주파수-방향 스펙트럼 S(ω, θ) (Wave)
    theta2D = atan2(Ky(:,:,1), Kx(:,:,1));  % 첫 슬라이스에서 방위각 계산 (rad)
    theta_bins = -pi : deg2rad(5) : pi;
    nThetaBins = length(theta_bins) - 1;
    S_omega_theta_wave = zeros(nThetaBins, length(w));
    for n = 1:length(w)
        S_slice = image_wave_spectrum(:,:,n);
        for j = 1:nThetaBins
            idx = (theta2D >= theta_bins(j)) & (theta2D < theta_bins(j+1));
            S_omega_theta_wave(j, n) = sum(S_slice(idx));
        end
    end
    nexttile;
    imagesc(theta_bins(1:end-1)*180/pi, w, S_omega_theta_wave');
    xlabel('Theta (deg)');
    ylabel('Angular Frequency \omega (rad/s)');
    title(['Wave: Frequency-Direction Spectrum for file ', num2str(i)]);
    colorbar;

    %% Surf 영역 스펙트럼 분석
    % 1. 주파수 스펙트럼 S(ω) (Surf)
    S_omega_surf = squeeze(sum(sum(image_surf_spectrum, 1), 2));
    nexttile;
    plot(w, S_omega_surf);
    xlabel('Angular Frequency \omega (rad/s)');
    ylabel('Energy');
    title(['Surf: Frequency Spectrum for file ', num2str(i)]);

    % 2. 파수 스펙트럼 S(kx, ky) (Surf)
    S_k_surf = sum(image_surf_spectrum, 3);
    nexttile;
    imagesc(kx, ky, S_k_surf);
    xlabel('k_x (rad/m)');
    ylabel('k_y (rad/m)');
    title(['Surf: Wavenumber Spectrum for file ', num2str(i)]);
    colorbar;

    % 3. 주파수-방향 스펙트럼 S(ω, θ) (Surf)
    theta2D = atan2(Ky(:,:,1), Kx(:,:,1));  % 동일하게 계산
    theta_bins = -pi : deg2rad(2) : pi;
    nThetaBins = length(theta_bins) - 1;
    S_omega_theta_surf = zeros(nThetaBins, length(w));
    for n = 1:length(w)
        S_slice = image_surf_spectrum(:,:,n);
        for j = 1:nThetaBins
            idx = (theta2D >= theta_bins(j)) & (theta2D < theta_bins(j+1));
            S_omega_theta_surf(j, n) = sum(S_slice(idx));
        end
    end
    nexttile;
    imagesc(theta_bins(1:end-1)*180/pi, w, S_omega_theta_surf');
    xlabel('Theta (deg)');
    ylabel('Angular Frequency \omega (rad/s)');
    title(['Surf: Frequency-Direction Spectrum for file ', num2str(i)]);
    colorbar;

    %%
    nexttile([1 3]);
    imshow(['C:\Users\Hyeonjong Im\Documents\새 폴더\image/Image_', dateStr, '.png']);

    %%


end