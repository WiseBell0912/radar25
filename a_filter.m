clear; close all; clc;

load('result_1910.mat');

%% 예제: 201x201 레이더 이미지에 필터 적용 및 비교

% 이미지가 이미 workspace에 201x201 행렬 (예: radarImage)로 로드되어 있다고 가정
radarImage = img_surf(:, :, 1779); % 필요시 이미지 로딩

% 윈도우 크기 설정 (예: 5x5)
window_size = 5;

%% Lee Filter 함수 정의
function filtered = leeFilter(img, window_size)
    % double 형으로 변환
    img = double(img);
    % 평균 필터 (정규화된 커널)
    h = ones(window_size) / (window_size^2);
    localMean = imfilter(img, h, 'symmetric');
    localVar = imfilter(img.^2, h, 'symmetric') - localMean.^2;
    % 노이즈 분산: 전체 이미지의 평균 분산으로 추정 (필요에 따라 다른 방법 적용 가능)
    noiseVar = mean(localVar(:));
    % 가중치 계산 (노이즈 분산과 지역분산을 이용)
    w = localVar ./ (localVar + noiseVar);
    % 보정된 결과
    filtered = localMean + w .* (img - localMean);
end

%% Frost Filter 함수 정의
function filtered = frostFilter(img, window_size, damping_factor)
    % img를 double 형으로 변환
    img = double(img);
    % 윈도우 중심 좌표
    halfWin = floor(window_size/2);
    [X, Y] = meshgrid(-halfWin:halfWin, -halfWin:halfWin);
    R = sqrt(X.^2 + Y.^2);
    % 지수 가중치: damping_factor는 조절 파라미터, 값이 클수록 중심에 가까운 픽셀에 높은 가중치 부여
    W = exp(-damping_factor * R);
    % 가중치 합계로 정규화
    W = W / sum(W(:));
    % 컨볼루션 적용
    filtered = conv2(img, W, 'same');
end

%% Kuan Filter 함수 정의
function filtered = kuanFilter(img, window_size)
    % double 형으로 변환
    img = double(img);
    h = ones(window_size) / (window_size^2);
    localMean = imfilter(img, h, 'symmetric');
    localVar = imfilter(img.^2, h, 'symmetric') - localMean.^2;
    % 노이즈 분산 추정 (전체 평균 분산 사용)
    noiseVar = mean(localVar(:));
    % 가중치 계산 (분모가 0인 경우를 방지)
    w = (localVar - noiseVar) ./ (localVar + eps);
    w(w < 0) = 0;  % 음수 값은 0으로 클리핑
    % 보정된 결과
    filtered = localMean + w .* (img - localMean);
end

%% 필터 적용 및 결과 비교
% Lee 필터 적용
lee_result = leeFilter(radarImage, window_size);

% Frost 필터 적용 (damping_factor는 예제에서는 0.5 정도로 설정, 필요에 따라 조정)
damping_factor = 0.5;
frost_result = frostFilter(radarImage, window_size, damping_factor);

% Kuan 필터 적용
kuan_result = kuanFilter(radarImage, window_size);

%% 결과 시각화
figure;
subplot(2,2,1);
surf(radarImage);
title('원본 레이더 이미지');

subplot(2,2,2);
surf(lee_result);
title('Lee Filter 결과');

subplot(2,2,3);
surf(frost_result);
title('Frost Filter 결과');

subplot(2,2,4);
surf(kuan_result);
title('Kuan Filter 결과');
