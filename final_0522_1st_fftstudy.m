clear; clc;

%% Search
%file_path = 'E:/png2019/10/';
file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
%file_path = '/media/zerog/EA168F94168F6085/png2019/';

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
Nx = 200;
Ny = 200;
Nt = 128;

dx = 3;
dy = 3;
dt = 1.43;

Lx = Nx * dx;
Ly = Ny * dy;
Lt = Nt * dt;

g = 9.81;
h = 30;

%% Frequency
function f = make_frequency_axis(N, d)
if mod(N, 2) == 0
    f = (-N/2 : N/2 - 1) / (N * d);    % 중심 기준
else
    f = (-(N-1)/2 : (N-1)/2) / (N * d);
end
f = f * 2 * pi;   % rad/m 또는 rad/s
end

kx = make_frequency_axis(Nx, dx);
ky = make_frequency_axis(Ny, dy);
w  = make_frequency_axis(Nt, dt);

[Ky, Kx, W] = meshgrid(kx, ky, w);
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

%% Loop
tic
for i = 696
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
    %img_surf = img_surf - mean(img_surf, 3);
    % wave
    img_wave = single(png_wave);
    img_wave = rot90(flip(img_wave));
    %img_wave = img_wave - mean(img_wave, 3);

    % ---------- 3D FFT ---------- %
    % surf
    F_surf = fftshift(fftn(img_surf));
    % wave
    F_wave = fftshift(fftn(img_wave));

    % ---------- Spectrum ---------- %
    % surf
    P_surf = abs(F_surf).^2 / Nx / Ny / Nt;
    P_surf = log(P_surf);
    % wave
    P_wave = abs(F_wave).^2 / Nx / Ny / Nt;
    P_wave = log(P_wave);

    % ---------- HPF ---------- %
    hpK = (0.0231642599 <= K) .* (K <= 0.1609926865);
    hpW = (2*pi/17 <= abs(W)) .* (abs(W) <= 2*pi/5); % 29번 ~ 55번 & 76번 ~ 101번
    hpMask = hpK .* hpW;
    % surf
    P_surf_hp = P_surf;% .* hpMask;
    % wave
    P_wave_hp = P_wave;% .* hpMask;


    % true 위치만 추출
    [x, y, z] = ind2sub(size(hpW), find(hpW));

    dkx = abs(kx(2) - kx(1));
    dky = abs(ky(2) - ky(1));
    dw = abs(w(2) - w(1));

    % 시각화
    figure(3);
    tiledlayout(1, 3);
    nexttile;
    [x, y, z] = ind2sub(size(hpW), find(hpW));
    scatter3((x*dkx-kx(end)), (y*dky-ky(end)), (z*dw-w(end)), 1, 'filled');  % 점 크기 10, 채움
    xlabel('Kx'); ylabel('Ky'); zlabel('W');
    title('Filter W : 5s ~ 17s');
    axis equal; grid off; view(0,0);
    xlim([min(kx), max(kx)]);
    ylim([min(ky), max(ky)]);
    zlim([min(w), max(w)]);
    nexttile;
    [x, y, z] = ind2sub(size(hpK), find(hpK));
    scatter3((x*dkx-kx(end)), (y*dky-ky(end)), (z*dw-w(end)), 1, 'filled');  % 점 크기 10, 채움
    xlabel('Kx'); ylabel('Ky'); zlabel('W');
    title('Filter K : w(5s) ~ w(17s)');
    axis equal; grid off; view(0,0);
    xlim([min(kx), max(kx)]);
    ylim([min(ky), max(ky)]);
    zlim([min(w), max(w)]);
    nexttile;
    [x, y, z] = ind2sub(size(hpMask), find(hpMask));
    scatter3((x*dkx-kx(end)), (y*dky-ky(end)), (z*dw-w(end)), 1, 'filled');  % 점 크기 10, 채움
    xlabel('Kx'); ylabel('Ky'); zlabel('W');
    title('Filter K & W');
    axis equal; grid off; view(0,0);
    xlim([min(kx), max(kx)]);
    ylim([min(ky), max(ky)]);
    zlim([min(w), max(w)]);

    % ---------- Find max ---------- %
    % surf
    [~, linear_idx_surf] = max(P_surf_hp(:));
    [ix_surf, iy_surf, it_surf] = ind2sub(size(P_surf), linear_idx_surf);
    % ix_surf = 101;
    % iy_surf = 101;
    % it_surf = 65;
    % wave
    [~, linear_idx_wave] = max(P_wave_hp(:));
    [ix_wave, iy_wave, it_wave] = ind2sub(size(P_wave), linear_idx_wave);
    % ix_wave = 101;
    % iy_wave = 101;
    % it_wave = 65;

    % ---------- Check 2 ---------- %
    figure(2);
    tiledlayout(1, 3);
    nexttile;
    imagesc(P_surf(:, :, 64));
    nexttile;
    imagesc(P_surf(:, :, 64) - P_surf(:, :, 65));
    nexttile;
    imagesc(P_surf(:, :, 66));

    % ---------- Check 3 ---------- %
    P_surf(101-3, 101-3, 65-5) - P_surf(101+3, 101+3, 65+5)

    % % ---------- Check 1 ---------- %
    % figure(1);
    % tiledlayout(2, 4);
    % % surf
    % nexttile;
    % imagesc(img_surf(:, :, 1));
    % axis equal; axis tight;
    % % surf w
    % nexttile;
    % hold on;
    % plot(w, squeeze(P_surf(ix_surf, iy_surf, :)));
    % plot(w, squeeze(P_surf_hp(ix_surf, iy_surf, :)));
    % hold off;
    % xlim([min(w), max(w)]); ylim([0, max(P_surf_hp(:))]);
    % title('surf w');
    % % surf kx
    % nexttile;
    % hold on;
    % plot(kx, squeeze(P_surf(:, iy_surf, it_surf)));
    % plot(kx, squeeze(P_surf_hp(:, iy_surf, it_surf)));
    % hold off;
    % xlim([min(kx), max(kx)]); ylim([0, max(P_surf_hp(:))]);
    % title('surf kx');
    % % surf ky
    % nexttile;
    % hold on;
    % plot(ky, squeeze(P_surf(ix_surf, :, it_surf)));
    % plot(ky, squeeze(P_surf_hp(ix_surf, :, it_surf)));
    % hold off;
    % xlim([min(ky), max(ky)]); ylim([0, max(P_surf_hp(:))]);
    % title('surf ky');
    % % wave img
    % nexttile;
    % imagesc(img_wave(:, :, 1));
    % axis equal; axis tight;
    % % wave w
    % nexttile;
    % hold on
    % plot(w, squeeze(P_wave(ix_wave, iy_wave, :)));
    % plot(w, squeeze(P_wave_hp(ix_wave, iy_wave, :)));
    % hold off
    % xlim([min(w), max(w)]); ylim([0, max(P_wave_hp(:))]);
    % title('wave w');
    % % surf kx
    % nexttile;
    % hold on
    % plot(kx, squeeze(P_wave(:, iy_wave, it_wave)));
    % plot(kx, squeeze(P_wave_hp(:, iy_wave, it_wave)));
    % hold off
    % xlim([min(kx), max(kx)]); ylim([0, max(P_wave_hp(:))]);
    % title('wave kx');
    % % surf ky
    % nexttile;
    % hold on;
    % plot(ky, squeeze(P_wave(ix_wave, :, it_wave)));
    % plot(ky, squeeze(P_wave_hp(ix_wave, :, it_wave)));
    % hold off;
    % xlim([min(ky), max(ky)]); ylim([0, max(P_wave_hp(:))]);
    % title('wave ky');

    disp([nFile, i]);
end
toc;