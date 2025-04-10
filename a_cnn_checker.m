clear; close all;

%% 1) 불러오기
load('a_cnn_result.mat');
load('data_grid.mat');
load("snr_y1910_0406.mat");
r_LandE = LandE;
r_SurfE = SurfE;
r_WaveE = WaveE;
r_Date = Date;
load("snr_y1911_0406.mat");
r_LandE = [r_LandE; LandE];
r_SurfE = [r_SurfE; SurfE];
r_WaveE = [r_WaveE; WaveE];
r_Date = [r_Date; Date];
load("snr_y1912_0406.mat");
r_LandE = [r_LandE; LandE];
r_SurfE = [r_SurfE; SurfE];
r_WaveE = [r_WaveE; WaveE];
r_Date = [r_Date; Date];

clear LandE SurfE WaveE Date wave_Uy wave_Ux wave_Tp wave_SNR wave_Pdir surf_Uy surf_Ux surf_Tp surf_SNR surf_Pdir

r_LandD = r_LandE / (167*67*128);
isRain_tf = (r_LandD > mean(r_LandD) + 2 * std(r_LandD));
isRain_categori = categorical(isRain_tf, [0 1], {'norain', 'rain'});

%% 2) 이미지 추출
file_path = 'E:/png2019/10/';
file_list = dir([file_path, '*.png']);
nFile = 1000;

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

for i = 500 : nFile
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

    figure(1);
    set(gcf, 'Position', [50, 980-700, 1920-100, 500]);
    tiledlayout(1,3);

    nexttile;
    modifiy_theta = pi * 5/3;
    r = linspace(800, 2333, 512);
    theta = linspace(0 - modifiy_theta, 2*pi - modifiy_theta, 1080);
    x = r' * cos(theta);
    y = r' * sin(theta);
    zone_lenght_x = 600;
    zone_length_y = 600;
    zx = -zone_lenght_x/2 : 3 : zone_lenght_x/2;
    zy = -zone_length_y/2 : 3 : zone_length_y/2;
    [Zx_sub, Zy_sub] = meshgrid(zx, zy);
    surf_Zx_center = Zx_sub * cos(surf_theta1) - Zy_sub * sin(surf_theta1);
    surf_Zy_center = Zx_sub * sin(surf_theta1) + Zy_sub * cos(surf_theta1);
    surf_Zx_real = surf_Zx_center + surf_center * cos(surf_theta2);
    surf_Zy_real = surf_Zy_center + surf_center * sin(surf_theta2);
    hold on;
    surf(x, y, png_long(:, :, 1), 'EdgeAlpha', 0);
    surf_xx = [surf_Zx_real(1, 1), surf_Zx_real(1, 201); surf_Zx_real(201, 1), surf_Zx_real(201, 201)];
    surf_yy = [surf_Zy_real(1, 1), surf_Zy_real(1, 201); surf_Zy_real(201, 1), surf_Zy_real(201, 201)];
    surf(surf_xx, surf_yy, [201, 201; 201, 201], 'FaceAlpha', 0, 'EdgeColor', 'red', 'LineWidth', 2);
    hold off
    title([dateStr(1:4), '/', dateStr(5:6), '/', dateStr(7:8), ' ', dateStr(10:11), ':', dateStr(12:13)], FontSize=15);
    axis equal;
    axis([-2500 2500 -2500 2500]);
    view(0, 90);

    nexttile;
    surf(data_grid_X, data_grid_Y, rot90(flip(single(png_surf(:, :, 1)))), 'EdgeAlpha', 0);
    view(0, 90);
    set(gca,'xtick',[],'ytick',[]);
    axis equal;
    axis([1 201 1 201]);

    nexttile;
    axis off; % 축 숨기기
    information = {...
        sprintf('%s', YPred_test(i)), ...
        sprintf('%s', isRain_categori(i)), ...
        };
    infoTextHandle = text(0.5, 0.5, information, ...
        'Units', 'normalized', ...
        'FontName', 'Courier New', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 16 ...
        );

    name = ['cnn_result_', dateStr, '.png'];
    saveas(gcf, name);
end