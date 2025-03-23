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
w = -pi/dt : 2*pi/Lt : pi/dt;
w(65) = [];


[Kx, Ky, W] = meshgrid(ky, kx, w);
Kx = single(Kx);
Ky = single(Ky);
W = single(W);
K = sqrt(Kx.^2 + Ky.^2);

%% Windowing
window_xy = hann(Nx) * hann(Ny)';
window_t = hann(Nt);
window = repmat(window_xy, 1, 1, Nt) .* reshape(window_t, 1, 1, Nt);
sc = Nt * Nx * Ny / sum(sum(sum(window.^2)));

%% Save
Date = NaT(length(file_list), 1);
SNR = zeros(length(file_list), 1);
Ux = zeros(length(file_list), 1);
Uy = zeros(length(file_list), 1);

%% Image
tic
parfor i = 1 : length(file_list)

    png_path = [file_list(i).folder '/' file_list(i).name];
    date = file_list(i).name(5:end-4);
    date = datetime(date, 'InputFormat', 'yyyyMMdd_HHmm');

    png_long = imread(png_path);
    png_long = png_long(1:512*1080*128);
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    [png_surf, png_wave] = f_extract_zone(png_long, modifiy_theta, surf_theta1, surf_theta2, surf_center, wave_theta1, wave_theta2, wave_center);

    img = png_wave;
    %% FFT
    image_spectrum = fftn(img .* window);
    image_spectrum_no_window = fftn(img);

    pt = sc*abs(image_spectrum).^2 / Nt^2 / Nx^2 / Ny^2;

    %% High pass filter
    hpK = K > 0.03;
    hpK = fftshift(hpK);

    hpW = abs(W) > 0.35;
    hpW = fftshift(hpW);

    image_spectrum = image_spectrum .* hpK .* hpW;
    image_spectrum_no_window = image_spectrum_no_window .* hpK .* hpW;

    image_spectrum_hp = abs(image_spectrum).^2 / Nt^2 / Nx^2 / Ny^2;
    pt_hp_ = abs(image_spectrum_no_window).^2 / Nt^2 / Nx^2 / Ny^2;
    %% current estimation
    %%%%%%%%%%%%% full frequency algorithm
    Kx_1D = reshape(Kx, 1, []);
    Ky_1D = reshape(Ky, 1, []);
    K_1D = sqrt(Kx_1D.^2 + Ky_1D.^2);
    sigma1d = sqrt(g * K_1D .* tanh(K_1D * h));

    image_spectrum_hp = fftshift(image_spectrum_hp);

    var = []; % 초기화를 빈 배열로 진행

    for ii = Nt/2:Nt
        pt_1D = reshape(image_spectrum_hp(:,:,ii), 1, []);
        idx = pt_1D > 0.00000001 * max(pt, [], 'all');

        Npt = sum(idx);
        Nend = size(var, 1);

        if Npt ~= 0
            var(Nend+1:Nend+Npt,1) = w(ii);
            var(Nend+1:Nend+Npt,2) = sigma1d(idx);
            var(Nend+1:Nend+Npt,3) = Kx_1D(idx);
            var(Nend+1:Nend+Npt,4) = Ky_1D(idx);
            var(Nend+1:Nend+Npt,5) = K_1D(idx);

            MTF1d = sqrt(K_1D(idx).^(-1.2) .* abs(w(ii))^(-0.6));

            var(Nend+1:Nend+Npt,6) = pt_1D(idx) .* MTF1d;
        end
    end

    var = var(2:end,:);



    %%%%%%%%%%%%% full frequency algorithm


    idx = var(:,1)<0;
    var(idx,2) = -var(idx,2);

    a1 = sum(var(:,6).*var(:,3).^2);%kx^2
    b1 = sum(var(:,6).*var(:,3).*var(:,4));%kx*ky
    c1 = sum(var(:,6).*(var(:,1) - var(:,2)) .* var(:,3));%(w-sigma)kx
    a2 = b1;%kx*ky
    b2 = sum(var(:,6).*var(:,4).^2);%ky^2
    c2 = sum(var(:,6).*(var(:,1) - var(:,2)) .* var(:,4));%(w-sigma)ky

    if (a1*b2 - a2*b1) ~= 0
        ux_cal = (c1*b2 - c2*b1) / (a1*b2 - a2*b1);
        uy_cal = (c1*a2 - c2*a1) / (a1*b2 - a2*b1);
    else
        ux_cal = a1/c1;
        uy_cal = a2/c2;
        disp('error or 1d wave')
        ux
        uy
    end

    %% band pass filter
    wk1 = sqrt(g*K.*tanh(K*h)) + Kx*ux_cal + Ky*uy_cal;
    wk2 = -sqrt(g*K.*tanh(K*h)) + Kx*ux_cal + Ky*uy_cal;

    bpv = 2 * (Nt/32) * 2*pi/Lt; %Band pass filtering 하는 정도 / a * (dw)
    idx = ~((W > wk1 - bpv & W < wk1 + bpv) | (W > wk2 - bpv & W < wk2 + bpv));
    idx = fftshift(idx);

    image_spectrum(idx) = 0;
    image_spectrum_no_window(idx) = 0;
    pt_bp = abs(image_spectrum).^2 / Nt^2 / Nx^2 / Ny^2;
    pt_bp_ = abs(image_spectrum_no_window).^2 / Nt^2 / Nx^2 / Ny^2; %windowing 미적용 for SNR(jj)

    %% spectrum analysis
    % pt = fftshift(pt);
    pt_bp_ = fftshift(pt_bp_);

    Sk2d_bp = (sum(pt_bp_(:,:,round(end/2):end),3));

    pw1d = 2*squeeze(sum(sum(pt_bp_,1),2));
    pw1d(1:end/2) = 0;

    %% calculate SNR
    signal = sum(sum(sum(pt_bp_)));
    noise = sum(sum(sum(pt_hp_))) - signal;
    snr = (signal / noise);

    % result
    disp(['snr = ', num2str(snr)]);

    %%
    SNR(i) = snr;
    Date(i) = date;
    Ux(i) = ux_cal;
    Uy(i) = uy_cal;

    %% Checker
    disp([length(file_list), i]);

end
toc

save("snr_y1910_wave_new.mat", 'Date', 'SNR', 'Ux', 'Uy');