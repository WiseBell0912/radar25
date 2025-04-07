clear; close all; clc;
%% radar_wave_quality.m
% 레이더 영상 기반 유의파고 추정의 품질 관리 스크립트
% 구현 내용
%   1) LandE 및 ASOS 강우 데이터를 이용한 강우 구간 판정
%   2) Surf / Wave 구역에서 레이더‑추정 파고의 신뢰도 평가
%   3) 위 두 결과를 종합한 품질 플래그 산출
%
% 필요 입력 (모두 동일 길이의 열 벡터, 시간 간격 동일)
%   LandE        - 육지 구역 픽셀 에너지 합계
%   SurfE        - Surf 구역 픽셀 에너지 합계
%   WaveE        - Wave 구역 픽셀 에너지 합계
%   b_Hs         - 부표 유의파고 [m]
%   r_surf_Hs    - Surf 구역 레이더 추정 유의파고 [m]
%   r_wave_Hs    - Wave 구역 레이더 추정 유의파고 [m]
%   asos_Precip  - ASOS 강우량 [mm/h] (0 = 무강우)

% Bouy
load("BOUY.mat");
b_Hs = b_SignificantWaveHeight;

b_Date = b_Date(~isnan(b_Hs));
b_Hs = b_Hs(~isnan(b_Hs));

b_Hs(b_Hs == 0) = NaN;
b_Hs = fillmissing(b_Hs, 'linear');

% ASOS
load("ASOS.mat");

% Rader E
load("snr_y1910_0406.mat");
r_Date = Date;
r_LandE = LandE; 
r_SurfE = SurfE;
r_WaveE = WaveE;
load("snr_y1911_0406.mat");
r_Date = [r_Date; Date];
r_LandE = [r_LandE; LandE];
r_SurfE = [r_SurfE; SurfE];
r_WaveE = [r_WaveE; WaveE];
load("snr_y1912_0406.mat");
r_Date = [r_Date; Date];
r_LandE = [r_LandE; LandE];
r_SurfE = [r_SurfE; SurfE];
r_WaveE = [r_WaveE; WaveE];

% Radar wave
load("./0327/snr_y1910_wave_0327.mat");
r_wave_SNR = SNR;
load("./0327/snr_y1911_wave_0327.mat");
r_wave_SNR = [r_wave_SNR ; SNR];
load("./0327/snr_y1912_wave_0327.mat");
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_SNR = sqrt(r_wave_SNR);

% Radar surf
load("./0327/snr_y1910_surf_0327.mat");
r_surf_SNR = SNR;
load("./0327/snr_y1911_surf_0327.mat");
r_surf_SNR = [r_surf_SNR ; SNR];
load("./0327/snr_y1912_surf_0327.mat");
r_surf_SNR = [r_surf_SNR ; SNR];
r_surf_SNR = sqrt(r_surf_SNR);

clear b_AirTemperature b_AtmosphericPressure b_CurrentDirection16Points b_CurrentDirectiondeg b_CurrentSpeed b_MaximumWaveHeight b_MaximumWavePeriod b_SignificantWaveHeight b_SignificantWavePeriod b_WaterTemperature b_WaveDirection b_WindDirection16Points b_WindDirectiondeg b_WindSpeed Date Pdir SNR surf_Pdir surf_SNR surf_Tp surf_Ux surf_Uy Tp Ux Uy wave_Pdir wave_SNR wave_Tp wave_Ux wave_Uy

% 처리
mask1 = ismember(b_Date, asos_Date);
mask2 = ismember(b_Date, r_Date);
mask = mask1 & mask2;

b_Date = b_Date(mask);
b_Hs = b_Hs(mask);

mask1 = ismember(asos_Date, b_Date);
mask2 = ismember(asos_Date, r_Date);
mask = mask1 & mask2;

asos_Date = asos_Date(mask);
asos_Precipitation = asos_Precipitation(mask);

mask1 = ismember(r_Date, b_Date);
mask2 = ismember(r_Date, asos_Date);
mask = mask1 & mask2;

r_Date = r_Date(mask);
r_LandE = r_LandE(mask);
r_surf_SNR = r_surf_SNR(mask);
r_SurfE = r_SurfE(mask);
r_wave_SNR = r_wave_SNR(mask);
r_WaveE = r_WaveE(mask);

clear mask mask1 mask2

% 개발
modelfun = @(x, SNR) x(1) + x(2) * SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
x = lsqcurvefit(modelfun, initial_guess, r_wave_SNR, b_Hs, [], [], options);

r_wave_Hs = x(1) + x(2) .* r_wave_SNR;

modelfun = @(y, SNR) y(1) + y(2) * SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
y = lsqcurvefit(modelfun, initial_guess, r_surf_SNR, b_Hs, [], [], options);

r_surf_Hs = y(1) + y(2) .* r_surf_SNR;

clear initial_guess modelfun options x y

% -------------------------------------------------------------------------

% ------------------------ 사용자 설정값 ----------------------------------
nPix_land   = 167*67*128;   % 육지 구역 픽셀 수 (필요 시 수정)
kRain       = 3;       % 고정 임계치: z‑score 배수
kMAD        = 3;       % 적응 임계치: median+MAD 배수
minHsTrain  = 0.2;     % 파고 학습에 포함할 최소 Hs (m)
minHsCalm   = 0.2;     % 잔잔한 바다로 간주할 최대 Hs (m)
corrWin     = 6;      % 이동 상관계수 창 길이 (샘플 수)
corrThr     = 0.5;     % 상관계수 임계치
plotFlag    = true;    % 진단 플롯 생성 여부
% -------------------------------------------------------------------------

% 0. 입력 벡터 길이 확인 ---------------------------------------------------
assert(isequal(size(r_LandE),size(r_SurfE),size(r_WaveE),size(b_Hs), ...
               size(r_surf_Hs),size(r_wave_Hs),size(asos_Precipitation)), ...
       '모든 입력 벡터의 크기가 동일해야 합니다.');

N = numel(r_LandE);
t = (1:N).';             % 시간 인덱스 (샘플 번호)

% 1. 정규화 ---------------------------------------------------------------
LandE_d = r_LandE ./ nPix_land;   % 육지 에너지 밀도
SurfR   = r_SurfE ./ r_LandE;       % Surf 에너지 / 육지 에너지 비
WaveR   = r_WaveE ./ r_LandE;       % Wave 에너지 / 육지 에너지 비

% 2. 강우 탐지 ------------------------------------------------------------
[isRain, thrRain, thrRainAdapt] = detectRain(LandE_d, asos_Precipitation, kRain, kMAD);

% 3. 파고 유효성 평가 ------------------------------------------------------
[isSurfValid, thrSurf, corrSurf] = detectWave(SurfR, b_Hs, isRain, ...
                                              minHsTrain, minHsCalm, ...
                                              corrWin, corrThr, r_surf_Hs);

[isWaveValid, thrWave, corrWave] = detectWave(WaveR, b_Hs, isRain, ...
                                              minHsTrain, minHsCalm, ...
                                              corrWin, corrThr, r_wave_Hs);

% 4. 통합 품질 플래그 ------------------------------------------------------
% 0 = 강우(사용불가), 1 = 의심, 2 = 양호
flag = zeros(N,1,'uint8');
flag(~isRain & isSurfValid & isWaveValid) = 2;   % 양호
flag(~isRain & ~(isSurfValid & isWaveValid)) = 1;% 의심

% 5. 진단 플롯 -------------------------------------------------------------
if plotFlag
    figure('Name','Radar‑wave QC','Color','w','Position',[100 100 900 700]);

    subplot(4,1,1)
    plot(t,LandE_d,'b'); hold on
    yline(thrRain,'r--','k = fixed');
    plot(t,thrRainAdapt,'m:');
    legend('LandE density','고정 임계치','적응 임계치','Location','best');
    title('강우 탐지');
    ylabel('LandE 밀도');
    grid on

    subplot(4,1,2)
    plot(t,SurfR,'b'); hold on
    yline(thrSurf,'g--','Surf thr');
    yyaxis right
    plot(t,corrSurf,'k'); yline(corrThr,'r--','corr thr');
    ylabel('상관계수');
    legend('SurfR','Surf 임계치','corr','corr 임계치','Location','best');
    title('Surf 구역 유효성');

    subplot(4,1,3)
    plot(t,WaveR,'b'); hold on
    yline(thrWave,'g--','Wave thr');
    yyaxis right
    plot(t,corrWave,'k'); yline(corrThr,'r--','corr thr');
    ylabel('상관계수');
    legend('WaveR','Wave 임계치','corr','corr 임계치','Location','best');
    title('Wave 구역 유효성');

    subplot(4,1,4)
    stairs(t,flag,'LineWidth',1.5); ylim([-0.5 2.5]);
    yticks([0 1 2]); yticklabels({'Rain','Doubt','Good'});
    xlabel('Time index');
    title('통합 품질 플래그');
    grid on
end

% -------------------------- 로컬 함수 ------------------------------------

function [isRain, thr_fixed, thr_adapt] = detectRain(LandE_d, precip, k, kMAD)
% detectRain  육지 구역 에너지를 이용해 강우 여부를 판정한다.
%
% 입력
%   LandE_d  - 육지 에너지 밀도 벡터
%   precip   - ASOS 강우량 (0 = 무강우)
%   k        - 고정 임계치용 z‑score 배수 (기본 3)
%   kMAD     - 적응 임계치용 MAD 배수 (기본 3)
%
% 출력
%   isRain     - 논리 벡터, true = 강우
%   thr_fixed  - 고정 임계치 값
%   thr_adapt  - 적응 임계치 벡터

if nargin < 3 || isempty(k),     k = 3;   end
if nargin < 4 || isempty(kMAD),  kMAD = 3; end

idxDry   = precip == 0;          % 무강우 구간 인덱스
mu       = mean(LandE_d(idxDry));
sigma    = std(LandE_d(idxDry));
thr_fixed = mu + k * sigma;      % 고정 임계치

% 적응 임계치: ±12시간(가정: 10분 간격) 이동창 median + k*MAD
win = 12*6;                      % 샘플 수
med  = movmedian(LandE_d, win, 'omitnan');
madv = movmad(LandE_d, win, 1);
thr_adapt = med + kMAD * madv;

% 두 임계치 중 하나라도 초과하면 강우로 판정
isRain = (LandE_d > thr_fixed) | (LandE_d > thr_adapt);

% 2샘플 이상 연속 검출 시에만 확정 (노이즈 완화)
isRain = conv(double(isRain), [1 1], 'same') >= 2;
isRain = logical(isRain);
end

% -------------------------------------------------------------------------
function [isValid, thr, corrVec] = detectWave(ERatio, b_Hs, isRain, ...
                                             minHsTrain, minHsCalm, ...
                                             corrWin, corrThr, r_est)
% detectWave  Surf/Wave 구역에서 레이더 파고 추정의 신뢰도를 평가한다.
%
% 입력
%   ERatio      - SurfR 또는 WaveR (E_zone / E_land)
%   b_Hs        - 부표 유의파고
%   isRain      - 강우 플래그
%   minHsTrain  - 학습에 포함할 최소 Hs (파 존재)
%   minHsCalm   - 잔잔한 바다로 간주할 최대 Hs
%   corrWin     - 이동 상관계수 창 길이
%   corrThr     - 상관계수 임계치
%   r_est       - 해당 구역 레이더 추정 Hs (없으면 [])
%
% 출력
%   isValid   - 논리 벡터, true = 유효
%   thr       - 에너지 임계치
%   corrVec   - 이동 상관계수 벡터 (r_est 없으면 NaN)

if nargin < 8, r_est = []; end
if nargin < 7 || isempty(corrThr), corrThr = 0.6; end
if nargin < 6 || isempty(corrWin), corrWin = 60;  end

N = numel(ERatio);
corrVec = nan(N,1);

% 에너지 기반 임계치 (방법 B)
idxCalm = ~isRain & (b_Hs < minHsCalm);
mu0  = mean(ERatio(idxCalm));
sig0 = std(ERatio(idxCalm));
thr  = mu0 + 2 * sig0;         % 2σ 규칙

isValidEnergy = (ERatio > thr) & ~isRain;

% 상관계수 기반 검증 (방법 A)
if ~isempty(r_est)
    corrVec = movingCorr(r_est, b_Hs, corrWin);
    isValidCorr = (corrVec > corrThr);
else
    isValidCorr = true(N,1);
end

% 두 기준 모두 만족해야 유효로 간주
isValid = isValidEnergy & isValidCorr;
end

% -------------------------------------------------------------------------
function c = movingCorr(x, y, win)
% movingCorr  단순 이동 Pearson 상관계수 계산.
%
% 입력
%   x, y - 열 벡터
%   win  - 창 길이 (샘플)
%
% 출력
%   c    - x/y와 동일 길이의 상관계수 벡터

N = numel(x);
c = nan(N,1);
half = floor(win/2);

for i = 1:N
    idx1 = max(1, i-half);
    idx2 = min(N, i+half);
    xx = x(idx1:idx2);
    yy = y(idx1:idx2);
    if numel(xx) >= 3 && ~all(isnan(xx)) && ~all(isnan(yy))
        r = corrcoef(xx,yy,'Rows','pairwise');
        c(i) = r(1,2);
    end
end
end

% -------------------------------------------------------------------------
% (선택) ROC 기반 강우 임계치 최적화 예시
% function thr_opt = optimizeRainThreshold(LandE_d, precip)
%     labels = precip > 0;
%     [~,~,~,thr_opt] = perfcurve(labels, LandE_d, true);
% end
