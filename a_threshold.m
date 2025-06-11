clear; close all; clc;

%% === 0. 데이터 불러오기
load("snr_y2020_0424.mat");
load("./ASOS.mat");
asos_Rain = asos_Precipitation; % 강우량 데이터
load("./AWS.mat");
asos_Rain = aws_RainVal; % 강우량 데이터

%% === 1. asos_Rain을 r_Date에 맞춰 보간 (Interpolation)
asos_Rain_interp = interp1(datenum(aws_Date), asos_Rain, datenum(r_Date), 'nearest', 'extrap');

%% === 2. 강우 여부 이진화
asos_RainTF = asos_Rain_interp > 0; % 강우량이 0mm 초과이면 비 왔다고 판단

%% === 3. Threshold 후보 생성
threshold_list = linspace(min(r_LandE), max(r_LandE), 100); % 100개 후보 생성

%% === 4. 성능 지표 저장용 배열 생성
nThreshold = length(threshold_list);
accuracy_list = zeros(nThreshold,1);
sensitivity_list = zeros(nThreshold,1);
specificity_list = zeros(nThreshold,1);

%% === 5. 각 Threshold에 대해 평가
for i = 1:nThreshold
    threshold = threshold_list(i);
    
    % 예측
    r_RainTF = r_LandE > threshold;
    
    % 성능 평가
    TP = sum((r_RainTF == 1) & (asos_RainTF == 1));
    TN = sum((r_RainTF == 0) & (asos_RainTF == 0));
    FP = sum((r_RainTF == 1) & (asos_RainTF == 0));
    FN = sum((r_RainTF == 0) & (asos_RainTF == 1));

    accuracy = (TP + TN) / (TP + TN + FP + FN);
    sensitivity = TP / (TP + FN);
    specificity = TN / (TN + FP);
    
    % 저장
    accuracy_list(i) = accuracy;
    sensitivity_list(i) = sensitivity;
    specificity_list(i) = specificity;
end

%% === 6. 최적 Threshold 찾기
% (1) Accuracy 기준
[max_accuracy, best_idx_acc] = max(accuracy_list);
best_threshold_acc = threshold_list(best_idx_acc);

% (2) Sensitivity + Specificity 기준
score_list = sensitivity_list + specificity_list;
[max_score, best_idx_score] = max(score_list);
best_threshold_score = threshold_list(best_idx_score);

%% === 7. 결과 출력
fprintf('[Accuracy 최대 기준]\n');
fprintf('최적 Threshold: %.4f\n', best_threshold_acc);
fprintf('Accuracy: %.2f%%\n', max_accuracy*100);
fprintf('Sensitivity: %.2f%%\n', sensitivity_list(best_idx_acc)*100);
fprintf('Specificity: %.2f%%\n', specificity_list(best_idx_acc)*100);

fprintf('\n[Sensitivity + Specificity 최대 기준]\n');
fprintf('최적 Threshold: %.4f\n', best_threshold_score);
fprintf('Accuracy: %.2f%%\n', accuracy_list(best_idx_score)*100);
fprintf('Sensitivity: %.2f%%\n', sensitivity_list(best_idx_score)*100);
fprintf('Specificity: %.2f%%\n', specificity_list(best_idx_score)*100);

%% === 8. ROC Curve 및 AUC 계산
FPR = 1 - specificity_list; % False Positive Rate
TPR = sensitivity_list;      % True Positive Rate

figure;
plot(FPR, TPR, '-o'); hold on;
plot([0 1], [0 1], 'k--'); % 랜덤 예측 기준선
xlabel('False Positive Rate (1 - Specificity)');
ylabel('True Positive Rate (Sensitivity)');
title('ROC Curve');
grid on;
legend('ROC Curve', 'Random Guess Line');
axis([0 1 0 1]);

% AUC 계산 (trapz로 근사)
AUC = trapz(FPR, TPR);
fprintf('\nROC Curve AUC: %.4f\n', AUC);

%% === 9. Threshold vs Performance 그래프
figure;
plot(threshold_list, accuracy_list, '-o', 'DisplayName', 'Accuracy'); hold on;
plot(threshold_list, sensitivity_list, '-x', 'DisplayName', 'Sensitivity');
plot(threshold_list, specificity_list, '-s', 'DisplayName', 'Specificity');
xlabel('Threshold (r\_LandE)');
ylabel('Performance');
title('Threshold vs Performance');
grid on;
legend show;