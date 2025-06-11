clear; close all; clc;

% for main = 10 : 12
% clearvars -except main


%% Search
%file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
%file_path = ['E:/png2019/',num2str(main), '/'];
%file_path = 'E:/png2019/10/';
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

energy_theta1 = deg2rad(100);
energy_theta2 = deg2rad(185);
energy_center = 1750 + 500/2;

%% 결과 저장용

tic

for i = 1 : length(file_list)
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

    %% 이미지 보기
    fig = figure(1);
    set(fig, "Position", [0, 0, 1500, 500]);
    tiledlayout(1, 3);

    nexttile;
    surf(img_energy(:, :, 1), 'EdgeAlpha', 0);
    xlim([1, 167]); ylim([1, 67]);
    view(0, 90);
    axis equal;
    axis off;
    xlabel('');
    ylabel('');
    title('LAND');

    nexttile;
    surf(img_surf(:, :, 1), 'EdgeAlpha', 0);
    xlim([1, 201]); ylim([1, 201]);
    view(0, 90);
    axis equal;
    axis off;
    xlabel('');
    ylabel('');
    title('SURF');

    nexttile;
    surf(img_wave(:, :, 1), 'EdgeAlpha', 0);
    xlim([1, 201]); ylim([1, 201]);
    view(0, 90);
    axis equal;
    axis off;
    xlabel('');
    ylabel('');
    title('WAVE');

    %% 이미지 저장
    name = ['./img3/img3_', dateStr, '.png'];
    saveas(fig, name);
end

toc