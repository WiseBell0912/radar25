clear; clc;

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
h = 50;

%% Frequency 계산 및 GPU 배열로 변환
kx = -pi/dx : 2*pi/Lx : pi/dx;
ky = -pi/dy : 2*pi/Ly : pi/dy;
w  = -pi/dt : 2*pi/Lt : pi/dt;
w(65) = [];

[Kx, Ky, W] = meshgrid(ky, kx, w);
Kx = single(Kx);  Ky = single(Ky);  W = single(W);
K  = sqrt(Kx.^2 + Ky.^2);

% GPU 배열로 변환
gpuKx = gpuArray(Kx);
gpuKy = gpuArray(Ky);
gpuW  = gpuArray(W);
gpuK  = gpuArray(K);

%% Windowing 계산 및 GPU 배열로 변환
window_xy = hann(Nx) * hann(Ny)';      % Nx x Ny
window_t  = hann(Nt);                  % Nt
window    = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);
gpu_window = gpuArray(single(window));  % 창을 single형으로 GPU로 이동

%% HPF용 마스크 미리 계산(반복문 밖, GPU로 변환)
hpK = (K > 0.025);
hpW = (abs(W) > 2*pi*0.03);
hpMask = hpK & hpW;
gpu_hpMask = gpuArray(hpMask);

%% 결과 저장용 (CPU 메모리)
nFile = length(file_list);
r_Date  = NaT(nFile,1);
r_LandE = zeros(nFile,1);

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
    tic
    %% 파일 읽기 및 Zone 추출 (CPU에서 수행)
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
    
    % GPU로 데이터 전송 (계산에 사용하기 전에)
    img_surf   = gpuArray(single(png_surf));
    img_wave   = gpuArray(single(png_wave));
    img_energy = gpuArray(single(png_energy));
    
    %% Land energy 계산 (GPU에서 수행 후 gather)
    LandEVal = gather( sum(img_energy, 'all') );
    
    %% FFT 연산 (GPU 연산)
    % surface zone
    img_surf_no_windo = rot90(flip(img_surf));
    img_surf_windowed = rot90(flip(img_surf)) .* gpu_window;
    
    image_surf_spectrum = fftshift(fftn(img_surf_windowed));
    image_surf_spectrum = abs(image_surf_spectrum).^2 / Nx / Ny / Nt;
    image_surf_spectrum_no_windo = fftshift(fftn(img_surf_no_windo));
    image_surf_spectrum_no_windo = abs(image_surf_spectrum_no_windo).^2 / Nx / Ny / Nt;
    
    % wave zone
    img_wave_no_windo = rot90(flip(img_wave));
    img_wave_windowed = rot90(flip(img_wave)) .* gpu_window;
    
    image_wave_spectrum = fftshift(fftn(img_wave_windowed));
    image_wave_spectrum = abs(image_wave_spectrum).^2 / Nx / Ny / Nt;
    image_wave_spectrum_no_windo = fftshift(fftn(img_wave_no_windo));
    image_wave_spectrum_no_windo = abs(image_wave_spectrum_no_windo).^2 / Nx / Ny / Nt;
    
    %% HPF (GPU 상에서)
    image_surf_spectrum_hp = image_surf_spectrum .* gpu_hpMask;
    image_wave_spectrum_hp = image_wave_spectrum .* gpu_hpMask;
    
    image_surf_spectrum_hp_no_windo = image_surf_spectrum_no_windo .* gpu_hpMask;
    image_wave_spectrum_hp_no_windo = image_wave_spectrum_no_windo .* gpu_hpMask;
    
    %% Current velocity estimation (GPU 상 계산)
    sigma_surf = sqrt(g .* gpuK .* tanh(gpuK .* h));
    sigma_wave = sqrt(g .* gpuK .* tanh(gpuK .* h));
    
    MTF_surf = (gpuK.^(-1.2)) .* (gpuW > 0) .* (image_surf_spectrum_hp > 1e-6 * max(image_surf_spectrum(:)));
    MTF_surf(~isfinite(MTF_surf)) = 0;
    MTF_wave = (gpuK.^(-1.2)) .* (gpuW > 0) .* (image_wave_spectrum_hp > 1e-6 * max(image_wave_spectrum(:)));
    MTF_wave(~isfinite(MTF_wave)) = 0;
    
    a1 = sum( image_surf_spectrum_hp .* MTF_surf .* (gpuKx.^2), 'all');
    a2 = sum( image_surf_spectrum_hp .* MTF_surf .* (gpuKy.^2), 'all');
    a3 = sum( image_surf_spectrum_hp .* MTF_surf .* (gpuKx .* gpuKy), 'all');
    a4 = sum( image_surf_spectrum_hp .* MTF_surf .* (gpuW - sigma_surf) .* (gpuKx.^2), 'all');
    a5 = sum( image_surf_spectrum_hp .* MTF_surf .* (gpuW - sigma_surf) .* (gpuKy.^2), 'all');
    
    ux_surf = gather( ( a2 * a4 - a3 * a5 ) / ( a1 * a2 - a3^2 ) );
    uy_surf = gather( ( a1 * a5 - a3 * a4 ) / ( a1 * a2 - a3^2 ) );
    
    b1 = sum( image_wave_spectrum_hp .* MTF_wave .* (gpuKx.^2), 'all');
    b2 = sum( image_wave_spectrum_hp .* MTF_wave .* (gpuKy.^2), 'all');
    b3 = sum( image_wave_spectrum_hp .* MTF_wave .* (gpuKx .* gpuKy), 'all');
    b4 = sum( image_wave_spectrum_hp .* MTF_wave .* (gpuW - sigma_wave) .* (gpuKx.^2), 'all');
    b5 = sum( image_wave_spectrum_hp .* MTF_wave .* (gpuW - sigma_wave) .* (gpuKy.^2), 'all');
    
    ux_wave = gather( ( b2 * b4 - b3 * b5 ) / ( b1 * b2 - b3^2 ) );
    uy_wave = gather( ( b1 * b5 - b3 * b4 ) / ( b1 * b2 - b3^2 ) );
    
    %% BPF (GPU 연산)
    bpv = 8 * (2*pi/Lt);
    
    bpMask_surf = (sqrt( g .* gpuK .* tanh( gpuK .* h ) ) + (gpuKx .* ux_surf) + (gpuKy .* uy_surf) - bpv <= abs(gpuW)) & ...
                  (abs(gpuW) <= sqrt( g .* gpuK .* tanh( gpuK .* h ) ) + (gpuKx .* ux_surf) + (gpuKy .* uy_surf) + bpv);
    bpMask_wave = (sqrt( g .* gpuK .* tanh( gpuK .* h ) ) + (gpuKx .* ux_wave) + (gpuKy .* uy_wave) - bpv <= abs(gpuW)) & ...
                  (abs(gpuW) <= sqrt( g .* gpuK .* tanh( gpuK .* h ) ) + (gpuKx .* ux_wave) + (gpuKy .* uy_wave) + bpv);
    
    image_surf_spectrum_bp = image_surf_spectrum_hp .* bpMask_surf;
    image_wave_spectrum_bp = image_wave_spectrum_hp .* bpMask_wave;
    
    image_surf_spectrum_bp_no_windo = image_surf_spectrum_hp_no_windo .* bpMask_surf;
    image_wave_spectrum_bp_no_windo = image_wave_spectrum_hp_no_windo .* bpMask_wave;
    
    %% SNR 계산 (GPU 상에서)
    MTF = (gpuK.^(-1.2)) .* (gpuW > 0);
    MTF(~isfinite(MTF)) = 0;
    
    signal_surf = sum( 2 .* image_surf_spectrum_bp_no_windo .* MTF, 'all');
    noise_surf = sum(image_surf_spectrum_no_windo, 'all') - sum(2 .* image_surf_spectrum_bp_no_windo .* (gpuW > 0), 'all');
    SNR_surfVal = gather( signal_surf / noise_surf );
    
    signal_wave = sum( 2 .* image_wave_spectrum_bp_no_windo .* MTF, 'all');
    noise_wave = sum(image_wave_spectrum_no_windo, 'all') - sum(2 .* image_wave_spectrum_bp_no_windo .* (gpuW > 0), 'all');
    SNR_waveVal = gather( signal_wave / noise_wave );
    
    %% Tp 계산 (GPU 상에서)
    image_surf_W_spectrum = squeeze(sum(sum(image_surf_spectrum_bp .* (gpuW > 0), 1), 2));
    [~, sorted_idx_surf] = sort(gather(image_surf_W_spectrum), 'descend');
    Tp_surfVal = mean(2 * pi ./ w(sorted_idx_surf(1:2)));
    
    image_wave_W_spectrum = squeeze(sum(sum(image_wave_spectrum_bp .* (gpuW > 0), 1), 2));
    [~, sorted_idx_wave] = sort(gather(image_wave_W_spectrum), 'descend');
    Tp_waveVal = mean(2 * pi ./ w(sorted_idx_wave(1:2)));
    
    %% Pdir 계산 (GPU 상에서)
    image_surf_K_spectrum = sum(image_surf_spectrum_bp, 3) .* (gpuKy(:, :, 1) < 0);
    [~, sorted_idx_surf] = sort(gather(image_surf_K_spectrum(:)), 'descend');
    Pdir_surfVal = mod(mean(rad2deg(atan2(gather(gpuKy(sorted_idx_surf(1:5))), gather(gpuKx(sorted_idx_surf(1:5)))))), 360);
    
    image_wave_K_spectrum = sum(image_wave_spectrum_bp, 3) .* (gpuKy(:, :, 1) < 0);
    [~, sorted_idx_wave] = sort(gather(image_wave_K_spectrum(:)), 'descend');
    Pdir_waveVal = mod(mean(rad2deg(atan2(gather(gpuKy(sorted_idx_wave(1:5))), gather(gpuKx(sorted_idx_wave(1:5)))))), 360);
    
    %% 결과 저장 (CPU 배열에 gather하여 저장)
    r_Date(i)  = dateVal;
    r_LandE(i) = LandEVal;
    
    r_surf_SNR(i)   = SNR_surfVal;
    r_surf_Ux(i)    = ux_surf;
    r_surf_Uy(i)    = uy_surf;
    r_surf_Tp(i)    = Tp_surfVal;
    r_surf_Pdir(i)  = Pdir_surfVal;
    
    r_wave_SNR(i)   = SNR_waveVal;
    r_wave_Ux(i)    = ux_wave;
    r_wave_Uy(i)    = uy_wave;
    r_wave_Tp(i)    = Tp_waveVal;
    r_wave_Pdir(i)  = Pdir_waveVal;
    
    %% 진행 상황 출력
    disp([num2str(i) , '/', num2str(nFile)]);
    disp('surf zone');
    disp([ux_surf, uy_surf, SNR_surfVal, Tp_surfVal, Pdir_surfVal]);
    disp('wave zone');
    disp([ux_wave, uy_wave, SNR_waveVal, Tp_waveVal, Pdir_waveVal]);
    toc
end

toc

save("snr_y1910_0406.mat", 'r_Date', 'r_LandE', 'r_surf_Pdir', 'r_surf_SNR', 'r_surf_Tp', 'r_surf_Ux', 'r_surf_Uy', 'r_wave_Pdir', 'r_wave_SNR', 'r_wave_Tp', 'r_wave_Ux', 'r_wave_Uy');
