%% 라벨링 완료된 레이더 이미지를 이용하여 간단한 2D CNN 분류 모델을 학습하는 코드
clear; close all; clc;

%% 1) 이미지 추출
file_path = 'E:/png2019/10/';
file_list = dir([file_path, '*.png']);
nFile = 499;

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

img_surf = zeros(201, 201, 499, 'single');

parfor i = 1 : 499
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

    img_surf(:, :, i) = rot90(flip(single(png_surf(:, :, 1))));
    
    disp(i);
end

%% 2) 라벨링 데이터
load('labeledPngData.mat');

%% 3) CNN 학습
[H,W,N] = size(img_surf);
X = reshape(img_surf, [H, W, 1, N]);

layers = [
    imageInputLayer([201 201 1],"Name","input","Normalization","none")
    convolution2dLayer(3,8,"Padding","same","Name","conv1")
    batchNormalizationLayer("Name","bn1")
    reluLayer("Name","relu1")
    maxPooling2dLayer(2,"Stride",2,"Name","pool1")
    
    convolution2dLayer(3,16,"Padding","same","Name","conv2")
    batchNormalizationLayer("Name","bn2")
    reluLayer("Name","relu2")
    maxPooling2dLayer(2,"Stride",2,"Name","pool2")
    
    fullyConnectedLayer(2,"Name","fc")
    softmaxLayer("Name","softmax")
    classificationLayer("Name","classOut")
];

options = trainingOptions('adam', ...
    'MaxEpochs', 10, ...
    'MiniBatchSize', 16, ...
    'Plots','training-progress',...
    'Verbose',true);

net = trainNetwork(X, labelVec, layers, options);

%% 4) 결과
test_file_path = 'E:/png2019/10/';
test_list = dir([test_file_path '*.png']);
nTest = 1000;

img_surf_test = zeros(201,201,nTest,'single');

parfor j = 1:nTest
    png_path = fullfile(test_list(j).folder, test_list(j).name);
    dateStr  = test_list(j).name(5:end-4);
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

    img_surf_test(:, :, j) = rot90(flip(single(png_surf(:, :, 1))));
    
    disp(j);
end

X_test = reshape(img_surf_test, [201,201,1,nTest]);

% 라벨이 있다면(예: testLabelVec), 정확도 측정
YPred_test = classify(net, X_test);

% 만약 testLabelVec(길이 nTest)이 있다면
% accuracy = mean(YPred_test == testLabelVec);
% confusionchart(testLabelVec, YPred_test);

% 없으면 그냥 예측 결과만 확인
disp(YPred_test);