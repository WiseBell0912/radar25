clear; close all; clc;

load('kkkk.mat');
s_r_Date = r_Date;
s_r_LandE = r_LandE;
s_r_surf_Noise = r_surf_Noise1;
s_r_surf_Signal = r_surf_Signal1;
s_r_wave_Noise = r_wave_Noise1;
s_r_wave_Signal = r_wave_Signal1;

r_Date = s_r_Date;
r_LandE = s_r_LandE;
r_surf_Noise = s_r_surf_Noise;
r_surf_Signal = s_r_surf_Signal;
r_wave_Noise = s_r_wave_Noise;
r_wave_Signal = s_r_wave_Signal;
r_surf_SNR = r_surf_Signal ./ r_surf_Noise;
r_wave_SNR = r_wave_Signal ./ r_wave_Noise;

%% Bouy
load("./BOUY.mat", "b_Date", "b_SignificantWaveHeight");
b_Hs = b_SignificantWaveHeight;
valid_idx = isfinite(b_Hs) & ~isnat(b_Date);
b_Hs = interp1(b_Date(valid_idx), b_Hs(valid_idx), b_Date, 'linear');

clear b_SignificantWaveHeight

%% AWS
load("./AWS.mat");

%% ASOS
load("./ASOS.mat");
asos_Rain = asos_Precipitation;

%% Masking
common_dates = intersect(intersect(r_Date, asos_Date), b_Date);

[~, idx_r]   = ismember(common_dates, r_Date);
[~, idx_asos] = ismember(common_dates, asos_Date);
[~, idx_b]   = ismember(common_dates, b_Date);

m_Date          = common_dates;

m_aws_RainVal   = asos_Rain(idx_asos);

m_b_Hs          = b_Hs(idx_b);

m_r_LandE       = r_LandE(idx_r);

m_r_surf_Signal = r_surf_Signal(idx_r);
m_r_wave_Signal = r_wave_Signal(idx_r);

m_r_surf_Noise = r_surf_Noise(idx_r);
m_r_wave_Noise = r_wave_Noise(idx_r);

m_r_surf_SNR = r_surf_SNR(idx_r);
m_r_wave_SNR = r_wave_SNR(idx_r);

%% SNR-HS
rain_filter = (m_r_LandE < mean(r_LandE) + 0.5 * std(r_LandE));
rain_filter = (m_r_LandE > 0);

% 보간된 m_b_Hs 사용
y = m_b_Hs(rain_filter);

% ===== Surf SNR =====
x_surf = sqrt(m_r_surf_SNR(rain_filter));

% 회귀
surf_fit = fit(x_surf, y, 'poly1');
A_surf = surf_fit.p2;
B_surf = surf_fit.p1;

% 예측값
m_r_surf_Hs = A_surf + B_surf * x_surf;

% ===== Wave SNR =====
x_wave = sqrt(m_r_wave_SNR(rain_filter));

% 회귀
wave_fit = fit(x_wave, y, 'poly1');
A_wave = wave_fit.p2;
B_wave = wave_fit.p1;

% 예측값
m_r_wave_Hs = A_wave + B_wave * x_wave;

% ===== 결과 출력 =====
fprintf('Surf 근사식: Hs = %.4f + %.4f * sqrt(SNR)\n', A_surf, B_surf);
fprintf('Wave 근사식: Hs = %.4f + %.4f * sqrt(SNR)\n', A_wave, B_wave);

% ===== 시각화 =====
figure;
tiledlayout(1,2)

nexttile
scatter(x_surf, y, 10, 'k', 'filled'); hold on;
plot(x_surf, m_r_surf_Hs, 'r', 'LineWidth', 2);
title('Surf SNR vs Hs');
xlabel('sqrt(SNR)'); ylabel('Hs'); legend('Data','Prediction');

nexttile
scatter(x_wave, y, 10, 'k', 'filled'); hold on;
plot(x_wave, m_r_wave_Hs, 'r', 'LineWidth', 2);
title('Wave SNR vs Hs');
xlabel('sqrt(SNR)'); ylabel('Hs'); legend('Data','Prediction');

%% R2
% 결정계수 (R^2) 계산 함수
calc_r2 = @(y_true, y_pred) 1 - sum((y_true - y_pred).^2) / sum((y_true - mean(y_true)).^2);

% R^2 계산
r2_surf = calc_r2(m_b_Hs(rain_filter), m_r_surf_Hs);
r2_wave = calc_r2(m_b_Hs(rain_filter), m_r_wave_Hs);

% 시각화
figure;
tiledlayout(1,2)
axis equal

% Surf
nexttile
scatter(m_r_surf_Hs, m_b_Hs(rain_filter), 10, 'k', 'filled'); hold on;
plot([0, max([m_r_surf_Hs ; m_r_wave_Hs ; m_b_Hs])], [0, max([m_r_surf_Hs ; m_r_wave_Hs ; m_b_Hs])], 'r--');  % y=x 선
xlabel('Predicted Hs (Surf SNR)');
ylabel('Observed Hs');
title(sprintf('Surf Prediction\nR^2 = %.3f', r2_surf));
hold off;
xlim([0, max([m_r_surf_Hs ; m_r_wave_Hs ; m_b_Hs])]);
ylim([0, max([m_r_surf_Hs ; m_r_wave_Hs ; m_b_Hs])]);
grid on

% Wave
nexttile
scatter(m_r_wave_Hs, m_b_Hs(rain_filter), 10, 'k', 'filled'); hold on;
plot([0, max([m_r_surf_Hs ; m_r_wave_Hs ; m_b_Hs])], [0, max([m_r_surf_Hs ; m_r_wave_Hs ; m_b_Hs])], 'r--');  % y=x 선
xlabel('Predicted Hs (Wave SNR)');
ylabel('Observed Hs');
title(sprintf('Wave Prediction\nR^2 = %.3f', r2_wave));
hold off;
xlim([0, max([m_r_surf_Hs ; m_r_wave_Hs ; m_b_Hs])]);
ylim([0, max([m_r_surf_Hs ; m_r_wave_Hs ; m_b_Hs])]);
grid on

%% Time series Hs
% figure;
% hold on;
% plot(m_Date, m_b_Hs);
% plot(m_Date, m_r_surf_Hs);
% plot(m_Date, m_r_wave_Hs);
% hold off;
% legend('Bouy', 'Surf', 'Wave');
% title('Hs');
% xlabel('Time'); ylabel('Hs');

% % % %% Filter
% % % surf_flat = m_r_surf_K_max ./ m_r_surf_K_mean;
% % % wave_flat = m_r_wave_K_max ./ m_r_wave_K_mean;
% % % 
% % % criteria = mean(r_LandE) + 0.5 * std(r_LandE);
% % % filter = (m_r_LandE > criteria);
% % % 
% % % fig = figure;
% % % tiledlayout(2, 1)
% % % 
% % % nexttile;
% % % hold on;
% % % %plot(m_Date, m_r_surf_K_max);
% % % plot(m_Date(rain_filter), surf_flat(rain_filter));
% % % hold off;
% % % xlim([min(m_Date) max(m_Date)]);
% % % legend('surf', 'wave');
% % % 
% % % ax_img = nexttile;
% % % img_path   = '/Users/limhyeonjong/Documents/Personal/GraduateProject/img3/';
% % % img_prefix = 'img3_';
% % % img_ext    = '.png';
% % % 
% % % % ==== 클릭 이벤트 등록 (원래대로) ====
% % % fig.WindowButtonDownFcn = @click_callback;
% % % 
% % % % ==== 자동으로 데이터팁 모드 활성화 ====
% % % dcm = datacursormode(fig);
% % % dcm.Enable = 'on';
% % % dcm.SnapToDataVertex = 'on';   % 꼭 데이터점에만 달리도록
% % % dcm.DisplayStyle = 'window';   % 팝업 창으로 표시
% % % dcm.UpdateFcn = {@datatipUpdate, r_Date, img_path, img_prefix, img_ext, ax_img};
% % % 
% % % % ===== 4) UpdateFcn 콜백 함수 정의 =====
% % % function txt = datatipUpdate(~, event_obj, r_Date, img_path, img_pr, img_ext, ax_img)
% % % % 1) DataIndex로 원본 인덱스 추출
% % % idx = event_obj.DataIndex;
% % % sel_date = r_Date(idx);
% % % 
% % % % 2) y값(비율)도 바로 가져오기
% % % ratio = event_obj.Position(2);
% % % 
% % % % 3) 이미지 파일명 만들기
% % % ts = datestr(sel_date, 'yyyymmdd_HHMM');
% % % fname = fullfile(img_path, [img_pr, ts, img_ext]);
% % % 
% % % % 4) 이미지 로드/업데이트
% % % 
% % % img = imread(fname);
% % % imshow(img, 'Parent', ax_img);
% % % title(ax_img, ['Image: ', datestr(sel_date,'yyyy-mm-dd HH:MM')]);
% % % 
% % % 
% % % % 5) DataTip에 표시할 텍스트
% % % txt = {
% % %     ['Time: ', datestr(sel_date, 'yyyy-mm-dd HH:MM')]
% % %     ['Ratio: ', num2str(ratio)]
% % %     };
% % % end