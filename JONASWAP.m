% 파라미터 설정
g = 9.81;           % 중력 가속도 (m/s^2)
alpha = 0.0081;     % 에너지 스케일링 상수
gamma = 3.3;        % 피크 강화 인자
f_p = 0.1;          % 피크 주파수 (Hz)

% 주파수 벡터 생성 (Hz)
f = linspace(0.01, 1, 1000);  % 예: 0.01Hz ~ 1Hz 범위

% JONSWAP 스펙트럼 계산
S = zeros(size(f));
for i = 1:length(f)
    if f(i) <= f_p
        sigma = 0.07;
    else
        sigma = 0.09;
    end
    S(i) = alpha * g^2 / (2*pi)^4 * f(i)^(-5) * exp(-5/4*(f_p/f(i))^4) * gamma^(exp(-((f(i)-f_p)^2)/(2*sigma^2*f_p^2)));
end

% S(f)를 플롯하여 확인
figure;
plot(f, S, 'LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('S(f) (m^2/Hz)');
title('JONSWAP Spectrum');
grid on;

% 이제 이 S(f)를 바탕으로 합성 해면을 생성할 수 있습니다.
% 예를 들어, 난수 위상(phi)을 부여하고 역 FFT를 수행하는 방법:
N = length(f);
phi = 2*pi*rand(1, N);  % 난수 위상
% 주파수 영역에서 복소수 계수를 구성
A = sqrt(2*S*f(2)-f(1)); % 스펙트럼의 크기를 적절히 조절 (이 부분은 보정 필요)
A = sqrt(S) .* exp(1i*phi);
% 역 FFT를 통해 시간 영역으로 변환 (간단 예시)
eta_t = ifft(A, 'symmetric');
figure;
plot(eta_t);
xlabel('Time Index');
ylabel('Surface Elevation');
title('Synthetic Sea Surface Time Series');
