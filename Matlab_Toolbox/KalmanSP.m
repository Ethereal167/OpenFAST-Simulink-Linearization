%% 卡尔曼滤波状态观测器
clc; clear;close all

% 系统模型
A = [1.2 0.5 0.3; 
     0.1 0.9 0.2; 
    -0.3 0.4 0.8];  % 状态转移矩阵
B = [0.2; 
     0.1; 
     0.5];    % 输入矩阵
H = [1, 0, 0;
     0, 0, 0;
     0, 0, 1];  % 量测状态矩阵

Q = 0.01 * eye(3);  % 过程噪声协方差矩阵
R = 0.01 * eye(3);           % 观测噪声协方差

% 初始化滤波器和观测器状态变量
x_real = zeros(3, 1);  % 初始真实状态
x_hat = 5*ones(3, 1);    % 初始的估计状态
P = eye(3);           % 估计协方差

% 生成观测数据
T = 15; % 时间步数
u = randn(T, 1); % 输入信号
y = zeros(T, 1); % 观测数据

for t = 2:T
    % 生成过程激励噪声和测量观测噪声
    w = sqrt(Q) * ones(3, 1) * randn;
    v = sqrt(R) * ones(3, 1) * randn;
    
    x_real(:, t) = A*x_real(:, t-1) + B*u(t-1) + w;
    z  = H * x_real(:, t) + v;

    % 时间更新（预测）
    x_hat_ = A * x_hat(:, t-1) + B*u(t-1);
    P_ = A * P * A' + Q;

    % 测量更新（校正）
    K = P_ * H' * (H * P_ * H' + R)^(-1);
    x_hat(:, t) = x_hat_ + K * (z - H * x_hat_);
    P = (eye(3) - K * H) * P_;

end

figure(1)
plot(x_hat(1, :))
hold on
plot(x_real(1, :))
legend('x-hat', 'x-real')

figure(2)
plot(x_hat(2, :))
hold on
plot(x_real(2, :))
legend('x-hat', 'x-real')

figure(3)
plot(x_hat(3, :))
hold on
plot(x_real(3, :))
legend('x-hat', 'x-real')