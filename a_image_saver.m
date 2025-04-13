clear; clc;

for main = 10 : 12
        clearvars -except main

%% Search
file_path = ['E:/png2019/',num2str(main), '/'];
%file_path = '/Users/limhyeonjong/Documents/Personal/GraduateProject/png2019/10/';
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
h_surf = 15;
h_wave = 30;

%% 결과 저장용
nFile = length(file_list);

Date  = NaT(nFile,1);

% main 루프 내, parfor 앞에서 출력 결과를 저장할 MAT 파일 생성
outFile = ['result_19', num2str(main), '.mat'];
mOut = matfile(outFile, 'Writable', true);

% 저장할 배열 사전 할당 (데이터 타입은 single)
% 결과 배열의 4번째 차원은 파일 인덱스(i)를 의미
mOut.img_surf = zeros(Nx, Ny, nFile, 'single');
mOut.img_wave = zeros(Nx, Ny, nFile, 'single');
mOut.img_energy = zeros(67, 167, nFile, 'single');

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

    %% 이미지 변환
    img_surf = rot90(flip(flip(img_surf), 2));
    img_wave = rot90(flip(flip(img_wave), 2));
    img_energy = rot90(flip(flip(img_energy), 2));

    %% 여기에서 이미지 저장
    Date(i)  = dateVal;

    mOut.img_surf(:, :, i) = img_surf(:, :, 1);
    mOut.img_wave(:, :, i) = img_wave(:, :, 1);
    mOut.img_energy(:, :, i) = img_energy(:, :, 1);


    %% 확인
    disp([num2str(i) ' / ' num2str(nFile) ' (폴더: ' num2str(main) ')']);

end

end
toc

system('shutdown -s -t 300');  % 60초 후 자동 종료