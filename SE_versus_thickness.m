% 参数设置
clear; clc;close all;

% 全局参数
% f = linspace(0, 18e9, 1000)'; % 频率区间，Hz
d = linspace(0, 15e-6, 1000)'; % 厚度区间，m
f = 10e9; % 频率，Hz
omega = 2 * pi * f; % 角频率
c0 = 3e8; % 真空光速 (m/s)
mu0 = 4*pi*1e-7; % 真空磁导率 (H/m)
epsilon0 = 8.854187817e-12; % 真空介电常数 (F/m)

% 材料参数
sigma = 5e5; % 电导率，单位S/m
% d = 10e-6;% 厚度，单位m

% 传输矩阵法
epsilon_air = 1; % 空气介电常数
epsilon_MXene = 1 + 1j *sigma./omega/epsilon0; % 复介电常数

k0 = omega/c0; % 真空波矢量

k_air = epsilon_air^0.5.*k0; % 空气波矢量，几乎同k_pz
k_MXene = (epsilon_MXene).^0.5.*k0; % MXene波矢量，几乎同k_qz

k_pz = sqrt(epsilon_air.*(k0.^2)-(pi./omega).^2);
k_qz = sqrt(epsilon_MXene.*(k0.^2)-(pi./omega).^2);

M11 = cos(k_qz.*d)+0.5j.*(k_qz./k_pz+k_pz./k_qz).*sin(k_qz.*d);
M12 = 0.5j.*(k_qz./k_pz-k_pz./k_qz).*sin(k_qz.*d);
M22 = cos(k_qz.*d)-0.5j.*(k_qz./k_pz+k_pz./k_qz).*sin(k_qz.*d);
M21 = -0.5j.*(k_qz./k_pz-k_pz./k_qz).*sin(k_qz.*d);   

S11 = -M21./M22;
S22 = M12./M22;
S12 = 1./M22;
S21 = M11-M12.*M21./M22;

T = abs(S12).^2;
R = abs(S11).^2;
A = 1-T-R;
SE_T = 10 * log10(1./T);

% 绘制结果
figure()
plot(d * 1e6,T,d * 1e6,R,d * 1e6,A);
legend("透射率","反射率","吸收率")

figure()
plot(d * 1e6, SE_T);
title("MXene薄膜在10 GHz下的电磁屏蔽效能")

% Export data to Excel
data = [d * 1e6, SE_T]; % Convert thickness to micrometers for output
filename = 'SE_MXene_10GHz.csv';
writematrix(data, filename);