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

%% 임의 상수(추가)
beta = deg2rad(10);  % 예시: 10도
mdeg = 5;            % 예시: 5도

%% Frequency
kx = -pi/dx : 2*pi/Lx : pi/dx;  
ky = -pi/dy : 2*pi/Ly : pi/dy;  
w  = -pi/dt : 2*pi/Lt : pi/dt;  
% 원래 코드대로 w(65)를 제거
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

%% 하이패스 필터(HPF)용 마스크 미리 계산(반복문 밖)
hpK = fftshift(K > 0.03);
hpW = fftshift(abs(W) > 0.35);
hpMask = hpK & hpW;

%% 결과 저장용
nFile = length(file_list);
Date  = NaT(nFile,1);
SNR   = zeros(nFile,1);
Ux    = zeros(nFile,1);
Uy    = zeros(nFile,1);

% 새로 추가된 결과(Tp, Pdir)도 배열로 선언
Tp   = zeros(nFile,1);
Pdir = zeros(nFile,1);

%% 미리 1D화된 변수 준비(현재추정용)
Kx_1D    = reshape(Kx, 1, []);
Ky_1D    = reshape(Ky, 1, []);
K_1D     = sqrt(Kx_1D.^2 + Ky_1D.^2);
sigma1d  = sqrt(g * K_1D .* tanh(K_1D * h));

%% w_ft 정의(예시: 전체 w 사용)
w_ft = w;  
% 만약 양수 주파수만 고려한다면
% w_ft = w(w>0); 
% 그리고 pw1d의 양수영역과 인덱스를 맞춰야 함.

tic;

parfor i = 1 : nFile
    %% 파일 읽기
    png_path = fullfile(file_list(i).folder, file_list(i).name);
    dateStr  = file_list(i).name(5:end-4);
    dateVal  = datetime(dateStr, 'InputFormat', 'yyyyMMdd_HHmm');

    png_long = imread(png_path);
    png_long = png_long(1 : 512*1080*128);  % 파일 크기에 맞춰 잘라내는 로직(원 코드 동일)
    png_long = reshape(png_long, 512, 1080, 128);
    png_long = flip(png_long, 2);
    png_long = flip(png_long, 3);

    [png_surf, png_wave] = f_extract_zone( ...
        png_long, modifiy_theta, ...
        surf_theta1, surf_theta2, surf_center, ...
        wave_theta1, wave_theta2, wave_center);

    png_long = [];  % 큰 배열 해제

    img = single(png_surf);
    png_surf = [];
    png_wave = [];

    %% FFT
    img_windowed          = img .* window;
    image_spectrum        = fftn(img_windowed);
    image_spectrum_no_win = fftn(img);

    % 임시 스펙트럼 파워 (threshold 계산)
    tmpPower   = sc * abs(image_spectrum).^2 / (Nt^2 * Nx^2 * Ny^2);
    globalMax  = max(tmpPower(:));

    %% High pass filter
    image_spectrum        = image_spectrum        .* hpMask;
    image_spectrum_no_win = image_spectrum_no_win .* hpMask;

    image_spectrum_hp = abs(image_spectrum).^2 / (Nt^2* Nx^2* Ny^2);
    pt_hp_            = abs(image_spectrum_no_win).^2 / (Nt^2* Nx^2* Ny^2);

    %% current estimation
    image_spectrum_hp_shift = fftshift(image_spectrum_hp);

    var   = zeros(1,6,'single');
    idxVar = 0;

    threshold = 1e-8 * globalMax;

    for ii_ = (Nt/2) : Nt
        pt_1D = reshape(image_spectrum_hp_shift(:,:,ii_), 1, []);
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

    %% band pass filter
    wk1 = sqrt(g*K.*tanh(K*h)) + Kx*ux_cal + Ky*uy_cal;
    wk2 = -sqrt(g*K.*tanh(K*h)) + Kx*ux_cal + Ky*uy_cal;

    bpv = 2 * (Nt/32) * 2*pi/Lt;

    bpMask = ~ ( (W > (wk1 - bpv) & W < (wk1 + bpv)) | ...
                 (W > (wk2 - bpv) & W < (wk2 + bpv)) );
    bpMask = fftshift(bpMask);

    image_spectrum(bpMask)        = 0;
    image_spectrum_no_win(bpMask) = 0;

    pt_bp  = abs(image_spectrum).^2        / (Nt^2* Nx^2* Ny^2);
    pt_bp_ = abs(image_spectrum_no_win).^2 / (Nt^2* Nx^2* Ny^2);

    pt_bp_ = fftshift(pt_bp_);

    %% calculate SNR
    signal = sum(pt_bp_(:));
    noise  = sum(pt_hp_(:)) - signal;
    snrVal = signal / max(noise,1e-12);

    %% calculate Tp
    pw1d = 2*squeeze(sum(sum(pt_bp_,1),2)); 
    % 음의 주파수영역 제거
    pw1d(1:floor(end/2)) = 0; 
    [~, idxMax] = max(pw1d);

    % w_ft(idxMax)가 0이 아닌지 확인 필요
    if idxMax > 0 && idxMax <= numel(w_ft) && w_ft(idxMax) ~= 0
        TpVal = 2*pi / w_ft(idxMax);
    else
        TpVal = 0; % 혹은 NaN
    end

    %% calculate Pdir
    Sk2d_bp = sum(pt_bp_(:,:,round(end/2):end),3);
    idxDir  = (Sk2d_bp == max(Sk2d_bp(:)));

    % idxDir가 여러 개일 수 있음 -> 그 중 첫 번째만 사용하거나, 평균 처리 등
    if nnz(idxDir)==1
        kxDir = Kx(idxDir);
        kyDir = Ky(idxDir);
    else
        % 여러 점이면 첫 점 사용 (혹은 평균)
        [rowDir, colDir] = find(idxDir);
        kxDir = Kx(rowDir(1), colDir(1), 1);
        kyDir = Ky(rowDir(1), colDir(1), 1);
    end

    PdirVal = mod((90 - rad2deg(atan2(kyDir, kxDir))) - rad2deg(beta) - mdeg, 360);

    %% 결과 저장
    SNR(i)  = snrVal;
    Date(i) = dateVal;
    Ux(i)   = ux_cal;
    Uy(i)   = uy_cal;
    Tp(i)   = TpVal;
    Pdir(i) = PdirVal;

    disp([nFile, i, snrVal, TpVal, PdirVal]);
end

toc;

save("snr_y1910_surf_0327.mat", 'Date', 'SNR', 'Ux', 'Uy', 'Tp', 'Pdir');
