clear; close all; clc;

%% ASOS
load("ASOS.mat");
asos_Date = asos_Date(1:end-180);
asos_Precipitation = asos_Precipitation(180:end);

%% Bouy
load("BOUY.mat");
b_Date = b_Date;
b_Hs = b_SignificantWaveHeight;
b_Wind = b_WindSpeed;
b_Pdir = b_WaveDirection;
b_Tp = b_MaximumWavePeriod;

b_Date = b_Date(~isnan(b_Hs));
b_Hs = b_Hs(~isnan(b_Hs));
b_Wind = b_Wind(~isnan(b_Hs));
b_Pdir = b_Pdir(~isnan(b_Hs));
b_Tp = b_Tp(~isnan(b_Hs));

b_Hs(b_Hs == 0) = NaN;
b_Hs = fillmissing(b_Hs, 'linear');

clear b_AirTemperature b_AtmosphericPressure b_CurrentDirection16Points b_CurrentDirectiondeg b_CurrentSpeed b_MaximumWaveHeight b_MaximumWavePeriod b_SignificantWaveHeight b_SignificantWavePeriod b_WaterTemperature b_WaveDirection b_WindDirection16Points b_WindDirectiondeg b_WindSpeed

%% Radar wave
load("./0327/snr_y1910_wave_0327.mat");
r_Date = Date;
r_wave_Pdir = Pdir;
r_wave_SNR = SNR;
r_wave_Tp = Tp;
r_wave_Ux = Ux;
r_wave_Uy = Uy;

load("./0327/snr_y1911_wave_0327.mat");
r_Date = [r_Date ; Date];
r_wave_Pdir = [r_wave_Pdir ; Pdir];
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_Tp = [r_wave_Tp ; Tp];
r_wave_Ux = [r_wave_Ux ; Ux];
r_wave_Uy = [r_wave_Uy ; Uy];

load("./0327/snr_y1912_wave_0327.mat");
r_Date = [r_Date ; Date];
r_wave_Pdir = [r_wave_Pdir ; Pdir];
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_Tp = [r_wave_Tp ; Tp];
r_wave_Ux = [r_wave_Ux ; Ux];
r_wave_Uy = [r_wave_Uy ; Uy];

r_wave_SNR = sqrt(r_wave_SNR);

clear Date Pdir SNR Tp Ux Uy

%% Radar surf
load("./0327/snr_y1910_surf_0327.mat");
r_Date = Date;
r_surf_Pdir = Pdir;
r_surf_SNR = SNR;
r_surf_Tp = Tp;
r_surf_Ux = Ux;
r_surf_Uy = Uy;

load("./0327/snr_y1911_surf_0327.mat");
r_Date = [r_Date ; Date];
r_surf_Pdir = [r_surf_Pdir ; Pdir];
r_surf_SNR = [r_surf_SNR ; SNR];
r_surf_Tp = [r_surf_Tp ; Tp];
r_surf_Ux = [r_surf_Ux ; Ux];
r_surf_Uy = [r_surf_Uy ; Uy];

load("./0327/snr_y1912_surf_0327.mat");
r_Date = [r_Date ; Date];
r_surf_Pdir = [r_surf_Pdir ; Pdir];
r_surf_SNR = [r_surf_SNR ; SNR];
r_surf_Tp = [r_surf_Tp ; Tp];
r_surf_Ux = [r_surf_Ux ; Ux];
r_surf_Uy = [r_surf_Uy ; Uy];

r_surf_SNR = sqrt(r_surf_SNR);

clear Date Pdir SNR Tp Ux Uy

%% Radar energy
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

%% Masking
mask = ismember(asos_Date, b_Date) & ismember(asos_Date, r_Date);
aasos_Date          = asos_Date(mask);
aasos_Precipitation = asos_Precipitation(mask);

mask = ismember(b_Date, asos_Date) & ismember(b_Date, r_Date);
bb_Date = b_Date(mask);
bb_Hs   = b_Hs(mask);

mask = ismember(r_Date, asos_Date) & ismember(r_Date, b_Date);
rr_Date     = r_Date(mask);
rr_LandE    = r_LandE(mask);
rr_SurfE    = r_SurfE(mask);
rr_WaveE    = r_WaveE(mask);
rr_surf_SNR = r_surf_SNR(mask);
rr_wave_SNR = r_wave_SNR(mask);

clear mask

%% Processing 1
rr_LandD = rr_LandE / (167*67*128);
rr_SurfD = rr_SurfE / (201*201*128);
rr_WaveD = rr_WaveE / (201*201*128);

%% Processing 2
modelfun = @(x, SNR) x(1) + x(2) * SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
x = lsqcurvefit(modelfun, initial_guess, rr_wave_SNR, bb_Hs, [], [], options);

rr_wave_Hs = x(1) + x(2) .* rr_wave_SNR;

modelfun = @(y, SNR) y(1) + y(2) * SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
y = lsqcurvefit(modelfun, initial_guess, rr_surf_SNR, bb_Hs, [], [], options);

rr_surf_Hs = y(1) + y(2) .* rr_surf_SNR;

clear modelfun initial_guess options x y

%% Figure 1 - PDF
pdf_data = rr_WaveD;

figure(1);
x_vals = linspace(min(pdf_data), max(pdf_data), 10000);
mu = mean(pdf_data);
sigma = std(pdf_data);
pdf_vals = normpdf(x_vals, mu, sigma);
plot(x_vals, pdf_vals, 'LineWidth', 2);
title('정규 분포 (PDF)');
xlabel('x'); ylabel('확률 밀도');
hold on;

% 평균선
xline(mu, '--k', 'Mean (μ)', 'LabelHorizontalAlignment','left');

% 표준편차 구간선
xline(mu - sigma, ':r', '-1σ');
xline(mu + sigma, ':r', '+1σ');
xline(mu - 2*sigma, ':g', '-2σ');
xline(mu + 2*sigma, ':g', '+2σ');
xline(mu - 3*sigma, ':b', '-3σ');
xline(mu + 3*sigma, ':b', '+3σ');

legend('PDF', 'Mean', '±1σ', '±2σ', '±3σ');
hold off;

%% Figure 2 - Energy
figure(2);
tiledlayout(2,1);
nexttile([1 1]);
hold on;

% 현재 y축의 한계를 얻기 위해 임의로 플롯 (아래에서 patch 적용 전에 ylimit이 설정되어야 함)
h1 = plot(rr_Date(rr_LandD < mean(rr_LandD) + 2 * std(rr_LandD)), pdf_data(rr_LandD < mean(rr_LandD) + 2 * std(rr_LandD)));
%h2 = plot(aasos_Date, ones(size(pdf_data))*( mean(pdf_data) + 0*std(pdf_data) ), 'r--');
%h3 = plot(aasos_Date, ones(size(pdf_data))*( mean(pdf_data) + 1*std(pdf_data) ), 'r--');
h4 = plot(aasos_Date, ones(size(pdf_data))*( mean(pdf_data) + 0*std(pdf_data) ), 'r--');
xlim([datetime(2019, 10, 1), datetime(2019, 10, 31)]);
title('Land Density');

% y축 한계를 가져옵니다.
yl = ylim;

% 0이 아닌 인덱스를 찾습니다.
nonzeroIdx = find(aasos_Precipitation ~= 0);

% 연속된 인덱스 구간(group) 찾기
groups = {};
if ~isempty(nonzeroIdx)
    group = nonzeroIdx(1);
    for k = 2:length(nonzeroIdx)
        % 만약 인덱스가 연속이면 group에 추가, 아니면 그룹을 저장하고 새 그룹 시작
        if nonzeroIdx(k) == nonzeroIdx(k-1) + 1
            group = [group; nonzeroIdx(k)];
        else
            groups{end+1} = group;
            group = nonzeroIdx(k);
        end
    end
    groups{end+1} = group;  % 마지막 그룹 저장
end

% 각 그룹별로 배경 채우기 (형광 노란색, 투명도 0.3)
% 시간 간격을 고려하여 각 구간의 시작과 끝을 구합니다.
% 여기서는 인접한 날짜간 간격이 일정하다고 가정하고,
% 양쪽에 절반 간격씩 확장하여 표시합니다.
if length(aasos_Date) > 1
    dt = median(diff(aasos_Date));  % 시간 간격 (날짜 형식이면 단위: days)
else
    dt = 0;
end

for i = 1:length(groups)
    idxGroup = groups{i};
    x_start = aasos_Date(idxGroup(1)) - dt/2;
    x_end   = aasos_Date(idxGroup(end))  + dt/2;
    
    % patch를 이용해 사각형 영역 채우기
    patch([x_start, x_end, x_end, x_start], [yl(1), yl(1), yl(2), yl(2)], ...
        [0.9 0.9 0.1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% patch가 배경에 깔리도록 플롯 순서를 조정
uistack(h1, 'top');
%uistack(h2, 'top');
%uistack(h3, 'top');
uistack(h4, 'top');

hold off;

% 데이터 커서 모드 활성화
dcm = datacursormode(gcf);
datacursormode on;

set(dcm, 'UpdateFcn', @(obj, event_obj) cursorCallback(obj, event_obj));

% 클릭 시 좌표를 출력하는 함수
function txt = cursorCallback(~, event_obj)
pos = event_obj.Position; % [x, y] 값
clickedDate = datetime(pos(1), 'ConvertFrom', 'datenum') + calyears(2019) + calmonths(9) + caldays(1); % x값을 datetime으로 변환

% 'yyyy-mm-dd HH:MM' 포맷으로 변환하여 출력
txt = {sprintf('.')}; % 추가 인수 사용

% 두 번째 타일에 이미지 업데이트
nexttile(2);
imshow(['/Users/limhyeonjong/Documents/Personal/GraduateProject/Image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
%imshow(['C:/Users/Hyeonjong Im/Documents/새 폴더/image/Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
end