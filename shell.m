clear; close all; clc;

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

h = 50;

%% Frequency
kx = -pi/dx : 2*pi/Lx : pi/dx;
ky = -pi/dy : 2*pi/Ly : pi/dy;
w = -pi/dt : 2*pi/Lt : pi/dt;
w(65) = [];


[Kx, Ky, W] = meshgrid(ky, kx, w);
Kx = single(Kx);
Ky = single(Ky);
W = single(W);
K = sqrt(Kx.^2 + Ky.^2);

%% shell
dispersion = sqrt(g.*K.*tanh(K.*h));

surf(dispersion(:, :, 2), EdgeAlpha=0);