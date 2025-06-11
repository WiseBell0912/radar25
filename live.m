clear; close all; clc;

load("snr_y2020_0423.mat");
load("./ASOS.mat");
asos_Rain = asos_Precipitation;
load("./AWS.mat");
load("./BOUY.mat", "b_Date", "b_SignificantWaveHeight");
b_Hs = b_SignificantWaveHeight;
valid_idx = isfinite(b_Hs) & ~isnat(b_Date);
b_Hs = interp1(b_Date(valid_idx), b_Hs(valid_idx), b_Date, 'linear');

clear b_SignificantWaveHeight

common_dates = intersect(intersect(r_Date, asos_Date), b_Date);

[~, idx_r]   = ismember(common_dates, r_Date);
[~, idx_asos] = ismember(common_dates, asos_Date);
[~, idx_b]   = ismember(common_dates, b_Date);

m_Date          = common_dates;

m_aws_RainVal   = asos_Rain(idx_asos);

m_b_Hs          = b_Hs(idx_b);

m_r_LandE       = r_LandE(idx_r);
m_r_surf_K_max  = r_surf_K_max(idx_r);
m_r_surf_K_mean = r_surf_K_mean(idx_r);
m_r_surf_Pdir   = r_surf_Pdir(idx_r);
m_r_surf_SNR    = r_surf_SNR(idx_r);
m_r_surf_Tp     = r_surf_Tp(idx_r);
m_r_surf_Ux     = r_surf_Ux(idx_r);
m_r_surf_Uy     = r_surf_Uy(idx_r);
m_r_surf_W_max  = r_surf_W_max(idx_r);
m_r_surf_W_mean = r_surf_W_mean(idx_r);
m_r_SurfE       = r_SurfE(idx_r);
m_r_wave_K_max  = r_wave_K_max(idx_r);
m_r_wave_K_mean = r_wave_K_mean(idx_r);
m_r_wave_Pdir   = r_wave_Pdir(idx_r);
m_r_wave_SNR    = r_wave_SNR(idx_r);
m_r_wave_Tp     = r_wave_Tp(idx_r);
m_r_wave_Ux     = r_wave_Ux(idx_r);
m_r_wave_Uy     = r_wave_Uy(idx_r);
m_r_wave_W_max  = r_wave_W_max(idx_r);
m_r_wave_W_mean = r_wave_W_mean(idx_r);
m_r_WaveE       = r_WaveE(idx_r);

%% SNR-HS
filter = (m_r_LandE < 30) & (m_r_wave_K_max > 0.1);
% 보간된 m_b_Hs 사용
y = m_b_Hs(filter);

% ===== Surf SNR =====
x_surf = sqrt(m_r_surf_SNR);
x_surf = x_surf(filter);

% 회귀
surf_fit = fit(x_surf, y, 'poly1');
A_surf = surf_fit.p2;
B_surf = surf_fit.p1;

% 예측값
m_r_surf_Hs = A_surf + B_surf * x_surf;

% ===== Wave SNR =====
x_wave = sqrt(m_r_wave_SNR);
x_wave = x_wave(filter);

% 회귀
wave_fit = fit(x_wave, y, 'poly1');
A_wave = wave_fit.p2;
B_wave = wave_fit.p1;

% 예측값
m_r_wave_Hs = A_wave + B_wave * x_wave;

fig = figure;
tiledlayout(2, 1)

nexttile;
hold on;
plot(m_Date, m_r_WaveE);
plot(m_Date, m_b_Hs);
hold off;
xlim([min(r_Date) max(r_Date)]);

ax_img = nexttile;
img_path   = '/Users/limhyeonjong/Documents/Personal/GraduateProject/img3/';
img_prefix = 'img3_';
img_ext    = '.png';

% ==== 클릭 이벤트 등록 (원래대로) ====
fig.WindowButtonDownFcn = @click_callback;

% ==== 자동으로 데이터팁 모드 활성화 ====
dcm = datacursormode(fig);
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'on';   % 꼭 데이터점에만 달리도록
dcm.DisplayStyle = 'window';   % 팝업 창으로 표시
dcm.UpdateFcn = {@datatipUpdate, r_Date, img_path, img_prefix, img_ext, ax_img};

% ===== 4) UpdateFcn 콜백 함수 정의 =====
function txt = datatipUpdate(~, event_obj, r_Date, img_path, img_pr, img_ext, ax_img)
% 1) DataIndex로 원본 인덱스 추출
idx = event_obj.DataIndex;
sel_date = r_Date(idx);

% 2) y값(비율)도 바로 가져오기
ratio = event_obj.Position(2);

% 3) 이미지 파일명 만들기
ts = datestr(sel_date, 'yyyymmdd_HHMM');
fname = fullfile(img_path, [img_pr, ts, img_ext]);

% 4) 이미지 로드/업데이트

img = imread(fname);
imshow(img, 'Parent', ax_img);
title(ax_img, ['Image: ', datestr(sel_date,'yyyy-mm-dd HH:MM')]);


% 5) DataTip에 표시할 텍스트
txt = {
    ['Time: ', datestr(sel_date, 'yyyy-mm-dd HH:MM')]
    ['Ratio: ', num2str(ratio)]
    };
end