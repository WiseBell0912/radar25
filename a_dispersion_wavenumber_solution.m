clear; close all; clc;

% 파라미터 설정
g = 9.81;                    % 중력가속도 [m/s^2]
d = 50;                      % 수심 [m]
f_th = 0.03;                 % 임계 주파수 [Hz]
omega_th = 2 * pi * f_th;    % 임계 각진동수 [rad/s]

% 초기 추정값 (얕은 수심 또는 깊은 수심 가정 기반으로)
k = omega_th^2 / g;          % 깊은 바다 가정의 초기값

% Newton-Raphson 반복 설정
tol = 1e-10;                 % 허용 오차
max_iter = 1000;
iter = 0;

while iter < max_iter
    % 함수값과 도함수 계산
    f = sqrt(g * k * tanh(k * d)) - omega_th;
    df = 0.5 * g * tanh(k * d) * (1 + k * d * (1 - tanh(k * d)^2)) / sqrt(g * k * tanh(k * d));
    
    % Newton-Raphson 업데이트
    k_new = k - f / df;
    
    % 수렴 검사
    if abs(k_new - k) < tol
        break;
    end
    
    k = k_new;
    iter = iter + 1;
end

% 결과 출력
fprintf('k_th = %.10f [rad/m] after %d iterations\n', k, iter);
