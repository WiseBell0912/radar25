clear; clc;
for main = 10 : 12
clearvars -except main


%% Search
%file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
file_path = ['E:/png2019/',num2str(main), '/'];
%file_path = 'E:/png2019/10/';
file_list = dir([file_path, '*.png']);

%% Information of Zone
modifiy_theta = pi * 5/3;

surf_theta1 = deg2rad(90);
surf_theta2 = deg2rad(160);
surf_center = 870 + 600/2;

wave_theta1 = deg2rad(90);
wave_theta2 = deg2rad(145);
wave_center = 1150 + 600/2;

energy_theta1 = deg2rad(100);
energy_theta2 = deg2rad(185);
energy_center = 1750 + 500/2;

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
surf_h = 15;
wave_h = 30;

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
hpK = (K > 0.03);
hpW = (abs(W) > 2*pi*0.03);
hpMask = hpK & hpW;

%% 결과 저장용
nFile = length(file_list);

r_Date  = NaT(nFile,1);
r_LandE = zeros(nFile,1);
r_SurfE = zeros(nFile,1);
r_WaveE = zeros(nFile,1);

r_surf_K_max = zeros(nFile,1);
r_surf_K_mean = zeros(nFile,1);
r_surf_W_max = zeros(nFile,1);
r_surf_W_mean = zeros(nFile,1);

r_wave_K_max = zeros(nFile,1);
r_wave_K_mean = zeros(nFile,1);
r_wave_W_max = zeros(nFile,1);
r_wave_W_mean = zeros(nFile,1);

r_surf_SNR   = zeros(nFile,1);
r_surf_Ux    = zeros(nFile,1);
r_surf_Uy    = zeros(nFile,1);
r_surf_Tp   = zeros(nFile,1);
r_surf_Pdir = zeros(nFile,1);

r_wave_SNR   = zeros(nFile,1);
r_wave_Ux    = zeros(nFile,1);
r_wave_Uy    = zeros(nFile,1);
r_wave_Tp   = zeros(nFile,1);
r_wave_Pdir = zeros(nFile,1);

tic

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

    [png_surf, png_wave, png_energy] = f_extract_zone_new( ...
        png_long, modifiy_theta, ...
        surf_theta1, surf_theta2, surf_center, ...
        wave_theta1, wave_theta2, wave_center, ...
        energy_theta1, energy_theta2, energy_center);

    img_surf = single(png_surf);
    img_wave = single(png_wave);
    img_energy = single(png_energy);

    img_surf = rot90(flip(flip(img_surf), 2));
    img_wave = rot90(flip(flip(img_wave), 2));
    img_energy = rot90(flip(flip(img_energy), 2));

    %% Energy
    LandEVal = sum(img_energy, 'all') / (size(img_energy, 1) * size(img_energy, 2) * size(img_energy, 3));
    SurfEVal = sum(img_surf, 'all') / (size(img_surf, 1) * size(img_surf, 2) * size(img_surf, 3));
    WaveEVal = sum(img_wave, 'all') / (size(img_wave, 1) * size(img_wave, 2) * size(img_wave, 3));

    %% FFT
    img_surf_no_windo   = img_surf;
    img_surf_windowed   = img_surf .* window;

    image_surf_spectrum          = fftshift(fftn(img_surf_windowed));
    image_surf_spectrum          = abs(image_surf_spectrum).^2 / Nx^2 / Ny^2 / Nt^2;
    image_surf_spectrum_no_windo = fftshift(fftn(img_surf_no_windo));
    image_surf_spectrum_no_windo = abs(image_surf_spectrum_no_windo).^2 / Nx^2 / Ny^2 / Nt^2;


    img_wave_no_windo   = img_wave;
    img_wave_windowed   = img_wave .* window;

    image_wave_spectrum          = fftshift(fftn(img_wave_windowed));
    image_wave_spectrum          = abs(image_wave_spectrum).^2 / Nx^2 / Ny^2 / Nt^2;
    image_wave_spectrum_no_windo = fftshift(fftn(img_wave_no_windo));
    image_wave_spectrum_no_windo = abs(image_wave_spectrum_no_windo).^2 / Nx^2 / Ny^2 / Nt^2;

    %% HPF
    image_surf_spectrum_hp = image_surf_spectrum .* hpMask;
    image_wave_spectrum_hp = image_wave_spectrum .* hpMask;

    image_surf_spectrum_hp_no_windo = image_surf_spectrum_no_windo .* hpMask;
    image_wave_spectrum_hp_no_windo = image_wave_spectrum_no_windo .* hpMask;

    %% Smooth factor
    dummy1 = sum(image_surf_spectrum_hp, 3);
    surf_K_max = max(dummy1, [], 'all');
    surf_K_mean = mean(dummy1, 'all');
    dummy2 = sum(sum(image_surf_spectrum_hp, 2), 3);
    surf_W_max = max(dummy2);
    surf_W_mean = mean(dummy2);

    dummy1 = sum(image_wave_spectrum_hp, 3);
    wave_K_max = max(dummy1, [], 'all');
    wave_K_mean = mean(dummy1, 'all');
    dummy2 = sum(sum(image_wave_spectrum_hp, 2), 3);
    wave_W_max = max(dummy2);
    wave_W_mean = mean(dummy2);

    %% Current velocity estimation
    sigma_surf = sqrt(g .* K .* tanh(K .* surf_h));
    sigma_wave = sqrt(g .* K .* tanh(K .* wave_h));

    MTF_surf = (K.^(-1.2)) .* W.^(-0.6) .* (W > 0) .* (image_surf_spectrum_hp > 1e-6 * max(image_surf_spectrum(:)));
    MTF_surf(~isfinite(MTF_surf)) = 0;
    MTF_wave = (K.^(-1.2)) .* W.^(-0.6) .* (W > 0) .* (image_wave_spectrum_hp > 1e-6 * max(image_wave_spectrum(:)));
    MTF_wave(~isfinite(MTF_wave)) = 0;

    a1 = sum( image_surf_spectrum_hp .* MTF_surf .* Kx.^2 , 'all');
    a2 = sum( image_surf_spectrum_hp .* MTF_surf .* Ky.^2 , 'all');
    a3 = sum( image_surf_spectrum_hp .* MTF_surf .* Kx .* Ky , 'all');
    a4 = sum( image_surf_spectrum_hp .* MTF_surf .* (W - sigma_surf) .* Kx , 'all');
    a5 = sum( image_surf_spectrum_hp .* MTF_surf .* (W - sigma_surf) .* Ky , 'all');

    ux_surf = ( a2 * a4 - a3 * a5 ) / ( a1 * a2 - a3^2 );
    uy_surf = ( a1 * a5 - a3 * a4 ) / ( a1 * a2 - a3^2 );

    b1 = sum( image_wave_spectrum_hp .* MTF_wave .* Kx.^2 , 'all');
    b2 = sum( image_wave_spectrum_hp .* MTF_wave .* Ky.^2 , 'all');
    b3 = sum( image_wave_spectrum_hp .* MTF_wave .* Kx .* Ky , 'all');
    b4 = sum( image_wave_spectrum_hp .* MTF_wave .* (W - sigma_wave) .* Kx , 'all');
    b5 = sum( image_wave_spectrum_hp .* MTF_wave .* (W - sigma_wave) .* Ky , 'all');

    ux_wave = ( b2 * b4 - b3 * b5 ) / ( b1 * b2 - b3^2 );
    uy_wave = ( b1 * b5 - b3 * b4 ) / ( b1 * b2 - b3^2 );

    % % %% Check
    % % figure(1);
    % % tiledlayout(1, 2);
    % % nexttile;
    % % hold on;
    % % surf(img_surf(:, :, 1), 'EdgeAlpha', 0);
    % % quiver3(101,101,max(img_surf(:)),ux_surf,uy_surf, 0, 10, 'r');
    % % hold off;
    % % view(0, 90);
    % % axis equal;
    % % xlim([1, 201]); ylim([1, 201]);
    % % nexttile;
    % % hold on;
    % % surf(img_wave(:, :, 1), 'EdgeAlpha', 0);
    % % quiver3(101,101,max(img_wave(:)),ux_wave,uy_wave, 0, 10, 'r');
    % % hold off;
    % % view(0, 90);
    % % axis equal;
    % % xlim([1, 201]); ylim([1, 201]);

    %% BPF
    bpv = 1 * (2*pi/Lt);

    bpMask_surf = (sqrt( g .* K .* tanh( K .* surf_h ) ) + (Kx .* ux_surf) + (Ky .* uy_surf) - bpv <= abs(W)) & (abs(W) <= sqrt( g .* K .* tanh( K .* surf_h ) ) + (Kx .* ux_surf) + (Ky .* uy_surf) + bpv);
    bpMask_wave = (sqrt( g .* K .* tanh( K .* wave_h ) ) + (Kx .* ux_wave) + (Ky .* uy_wave) - bpv <= abs(W)) & (abs(W) <= sqrt( g .* K .* tanh( K .* wave_h ) ) + (Kx .* ux_wave) + (Ky .* uy_wave) + bpv);

    image_surf_spectrum_bp = image_surf_spectrum_hp .* bpMask_surf;
    image_wave_spectrum_bp = image_wave_spectrum_hp .* bpMask_wave;

    image_surf_spectrum_bp_no_windo = image_surf_spectrum_hp_no_windo .* bpMask_surf;
    image_wave_spectrum_bp_no_windo = image_wave_spectrum_hp_no_windo .* bpMask_wave;

    %% SNR
    signal_surf = sum(image_surf_spectrum_bp_no_windo, 'all');
    noise_surf = sum(image_surf_spectrum_hp_no_windo, 'all') - signal_surf;
    SNR_surfVal = signal_surf / noise_surf;

    signal_wave = sum(image_wave_spectrum_bp_no_windo, 'all');
    noise_wave = sum(image_wave_spectrum_hp_no_windo, 'all') - signal_wave;
    SNR_waveVal = signal_wave / noise_wave;

    %% Tp
    image_surf_W_spectrum = squeeze(sum(sum(image_surf_spectrum_bp .* (W > 0), 1), 2));
    % [~, maxidx_surf_Tp] = max(image_surf_W_spectrum);
    % Tp_surfVal = 2 * pi / w(maxidx_surf_Tp);
    [~, sorted_idx_surf] = sort(image_surf_W_spectrum, 'descend');
    Tp_surfVal = mean(2 * pi ./ w(sorted_idx_surf(1:2)));

    image_wave_W_spectrum = squeeze(sum(sum(image_wave_spectrum_bp .* (W > 0), 1), 2));
    % [~, maxidx_wave_Tp] = max(image_wave_W_spectrum);
    % Tp_waveVal = 2 * pi / w(maxidx_wave_Tp);
    [~, sorted_idx_wave] = sort(image_wave_W_spectrum, 'descend');
    Tp_waveVal = mean(2 * pi ./ w(sorted_idx_wave(1:2)));

    %% Pdir
    image_surf_K_spectrum = sum(image_surf_spectrum_bp, 3) .* (Ky(:, :, 1) < 0);
    % [~, maxidx_surf_Pdir] = max(image_surf_K_spectrum, [], "all");
    % Pdir_surfVal = mod(90 - rad2deg(atan2(Ky(maxidx_surf_Pdir), Kx(maxidx_surf_Pdir))) + 60, 360);
    [~, sorted_idx_surf] = sort(image_surf_K_spectrum(:), 'descend');
    Pdir_surfVal = mod(mean(rad2deg(atan2(Ky(sorted_idx_surf(1:5)), Kx(sorted_idx_surf(1:5))))), 360);

    image_wave_K_spectrum = sum(image_wave_spectrum_bp, 3) .* (Ky(:, :, 1) < 0);
    % [~, maxidx_wave_Pdir] = max(image_wave_K_spectrum, [], "all");
    % Pdir_waveVal = mod(90 - rad2deg(atan2(Ky(maxidx_wave_Pdir), Kx(maxidx_wave_Pdir))) + 60, 360);
    [~, sorted_idx_wave] = sort(image_wave_K_spectrum(:), 'descend');
    Pdir_waveVal = mod(mean(rad2deg(atan2(Ky(sorted_idx_wave(1:5)), Kx(sorted_idx_wave(1:5))))), 360);

    % % %% Check
    % % figure(1);
    % % tiledlayout(1, 2);
    % % nexttile;
    % % hold on;
    % % surf(sum(image_surf_spectrum_hp, 3), 'EdgeAlpha', 0);
    % % quiver3(101,101,max(img_surf(:)),ux_surf,uy_surf, 0, 10, 'r');
    % % hold off;
    % % view(0, 90);
    % % axis equal;
    % % xlim([1, 201]); ylim([1, 201]);
    % % nexttile;
    % % hold on;
    % % surf(sum(image_wave_spectrum_hp, 3), 'EdgeAlpha', 0);
    % % quiver3(101,101,max(img_wave(:)),ux_wave,uy_wave, 0, 10, 'r');
    % % hold off;
    % % view(0, 90);
    % % axis equal;
    % % xlim([1, 201]); ylim([1, 201]);

    %% Save
    r_Date(i)  = dateVal;
    r_LandE(i) = LandEVal;
    r_SurfE(i) = SurfEVal;
    r_WaveE(i) = WaveEVal;

    r_surf_K_max(i) = surf_K_max;
    r_surf_K_mean(i) = surf_K_mean;
    r_surf_W_max(i) = surf_W_max;
    r_surf_W_mean(i) = surf_W_mean;

    r_wave_K_max(i) = wave_K_max;
    r_wave_K_mean(i) = wave_K_mean;
    r_wave_W_max(i) = wave_W_max;
    r_wave_W_mean(i) = wave_W_mean;

    r_surf_SNR(i)   = SNR_surfVal;
    r_surf_Ux(i)    = ux_surf;
    r_surf_Uy(i)    = uy_surf;
    r_surf_Tp(i)   = Tp_surfVal;
    r_surf_Pdir(i) = Pdir_surfVal;

    r_wave_SNR(i)   = SNR_waveVal;
    r_wave_Ux(i)    = ux_wave;
    r_wave_Uy(i)    = uy_wave;
    r_wave_Tp(i)   = Tp_surfVal;
    r_wave_Pdir(i) = Pdir_surfVal;

    %% Check
    %disp([num2str(i) , '/', num2str(nFile)]);
    disp([num2str(main), '/', num2str(i) , '/', num2str(nFile)]);
    disp('surf zone');
    disp([ux_surf, uy_surf, SNR_surfVal]);
    disp('wave zone');
    disp([ux_wave, uy_wave, SNR_waveVal]);
end

toc


%save('snr_y1910_0406.mat', 'r_Date', 'r_LandE', 'r_SurfE', 'r_WaveE', 'r_surf_Pdir', 'r_surf_SNR', 'r_surf_Tp', 'r_surf_Ux', 'r_surf_Uy', 'r_wave_Pdir', 'r_wave_SNR', 'r_wave_Tp', 'r_wave_Ux', 'r_wave_Uy', 'r_surf_W_max', 'r_surf_W_mean', 'r_wave_K_max', 'r_wave_K_mean', 'r_wave_W_max', 'r_wave_W_mean', 'r_surf_K_max', 'r_surf_K_mean');
save(['snr_y19', num2str(main), '_0422.mat'], 'r_Date', 'r_LandE', 'r_SurfE', 'r_WaveE', 'r_surf_Pdir', 'r_surf_SNR', 'r_surf_Tp', 'r_surf_Ux', 'r_surf_Uy', 'r_wave_Pdir', 'r_wave_SNR', 'r_wave_Tp', 'r_wave_Ux', 'r_wave_Uy', 'r_surf_W_max', 'r_surf_W_mean', 'r_wave_K_max', 'r_wave_K_mean', 'r_wave_W_max', 'r_wave_W_mean', 'r_surf_K_max', 'r_surf_K_mean');


end

system('shutdown -s -t 600');