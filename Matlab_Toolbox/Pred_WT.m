% 风力机在工作点处的预测模型
%% 准备工作
clc
clear; close all;
load("StateSpaceList.mat")

Azimuth_list = StateSpaceList.Azimuth_list;
n_FixFrameInputs = StateSpaceList.n_FixFrameInputs;
n_RotTripletInputs = StateSpaceList.n_RotTripletInputs;
new_seq_inp = StateSpaceList.new_seq_inp;
n_FixFrameOutputs = StateSpaceList.n_FixFrameOutputs;
n_RotTripletOutputs = StateSpaceList.n_RotTripletOutputs;
new_seq_out = StateSpaceList.new_seq_out;

OP_num = 2;
A = StateSpaceList.A_total(:, :, OP_num);
B = StateSpaceList.B_total(:, :, OP_num);
Bd = StateSpaceList.Bd_total(:, :, OP_num);
C = StateSpaceList.C_total(:, :, OP_num);
D = StateSpaceList.D_total(:, :, OP_num);
Dd = StateSpaceList.Dd_total(:, :, OP_num);

% 定义步长相关变量
DT = 1e-3;  % 积分步长，一定要足够小，不然积分会发散
DT_gap = 200; % 实际步长为DT*DT_gap

WindSpeed = 14; % m/s
BldPitch = 8.47 * 1/(57.2957) *  ones(3, 1); %rad

xop = StateSpaceList.xop_total(:, :, OP_num);
uop = StateSpaceList.uop_total(:, :, OP_num);
wop = StateSpaceList.wop_total(:, :, OP_num);
yop_list = StateSpaceList.yop_List_total(:, :, OP_num);


yop = yop_list(:, 1);
u = uop;
w = wop;
w(1) = WindSpeed;
u(1) = BldPitch(1);
u(2) = BldPitch(2);
u(3) = BldPitch(3);
delta_u = u - uop;
delta_w = w - wop;

Azimuth = 0;
Tmbc  = mbc(Azimuth);
delta_w_dq = delta_w;
delta_xdq = xop .* 0;

Am = [DT*A+eye(size(A, 1)), zeros(size(A, 1), size(C, 1));
      C, eye(size(C, 1))];
Bm = [DT*B*Tmbc; D*Tmbc];
Bmd = [DT*Bd; Dd];
Cm = [zeros(size(C, 1), size(A, 1)), eye(size(C, 1))];

Np = 20;  % 预测时长为 Np*DT*DT_gap
Nc = 5;

% ---------------------------------- 计算矩阵P ----------------------------
for i = 1:Np
    if i == 1
        P = Cm*(Am^1);
    else
        P = [P; Cm*(Am^(1+DT_gap*(i-1)))];
    end
end

% ---------------------------------- 计算矩阵H ----------------------------
for i = 1:Np
    H_c = [];
    for j = 1:Nc
        if j==1
            H_c = Cm * Am^(DT_gap*(i-1)) * Bm;
        else
            if j > i
                H_c = [H_c, zeros(size(Cm, 1), size(Bm, 2))];
            else
                H_c = [H_c, Cm * Am^(DT_gap*(i-j)) * Bm];
            end
        end
    end
    if i == 1
        H = H_c;
    else
        H = [H; H_c];
    end
end

% ---------------------------------- 计算矩阵Hd ----------------------------
for i = 1:Np
    H_c = [];
    for j = 1:Np
        if j==1
            H_c = Cm * Am^(DT_gap*(i-1)) * Bmd;
        else
            if j > i
                H_c = [H_c, zeros(size(Cm, 1), size(Bmd, 2))];
            else
                H_c = [H_c, Cm * Am^(DT_gap*(i-j)) * Bmd];
            end
        end
    end
    if i == 1
        Hd = H_c;
    else
        Hd = [Hd; H_c];
    end
end
% ------------------------------- 初始状态 ---------------------------------------
xm = [delta_xdq; Transverse_y(Azimuth, n_FixFrameOutputs, n_RotTripletOutputs, new_seq_out, yop)];   % [delta_x; yop_dq]
u_initial = uop + delta_u;  % 初始控制动作

% ------------------------------- 风速序列 -------------------------------
delta_w_list = delta_w_dq;
for i = 2:Np
    delta_w_list = [delta_w_list; 0*delta_w_dq];  % 风速不变
end

% ------------------------------- 求解最优控制序列 ---------------------------------------
% 设置参考输出序列
Rs_1 = [1173.6; 12.1; 0; 0; 0]; Rs = zeros(5*Np, 1);
for i = 1:Np
    Rs(1+5*(i-1):5*i, 1) = Rs_1;
end
% 设置输出权重
Q_1 = eye(5); Q = eye(5*Np);
Q_1(1, 1) = 10000; Q_1(2, 2) = 0; Q_1(3, 3) = 10; Q_1(4, 4) = 10; Q_1(5, 5) = 10;
for i = 1:Np
    Q(1+5*(i-1):5*i, 1+5*(i-1):5*i) = Q_1;
end
% 设置输入权重
R_1 = eye(3); R = eye(3*Nc);
R_1(1, 1) = 100; R_1(2, 2) = 100; R_1(3, 3) = 100; 
for i = 1:Nc
    R(1+3*(i-1):3*i, 1+3*(i-1):3*i) = R_1;
end
% 控制输入的约束
delta_u_max = [2; 2; 2] * (pi / 180) * (DT * DT_gap); % 变桨速率约束的上限 rad
delta_u_min = delta_u_max; % 变桨速率约束的下限 rad
u_max = [90; 90; 90] * (pi / 180); % 变桨约束上限 rad
u_min = [0; 0; 0] * (pi / 180); % 变桨约束下限 rad
% 构造二次型的矩阵
% E
E = H' * Q * H + R;
% F
F = 2 * (P * xm + Hd * delta_w_list - Rs)' * Q * H;
% 约束相关矩阵
b = zeros(4*Nc*3, 1);
for i = 1:Nc
    b(1+3*(i-1):3*i, 1) = delta_u_max;
end
for i = 1:Nc
    b(1+3*Nc+3*(i-1):3*Nc+3*i, 1) = delta_u_min;
end
for i = 1:Nc
    b(1+6*Nc+3*(i-1):6*Nc+3*i, 1) = u_max-u_initial;
end
for i = 1:Nc
    b(1+9*Nc+3*(i-1):9*Nc+3*i, 1) = u_min+u_initial;
end
% Acons
Acons = zeros(12*Nc, 3*Nc);
Acons(1:3*Nc, :) = eye(3*Nc);
Acons(1+3*Nc:6*Nc, :) = -1*eye(3*Nc);
K = zeros(3*Nc, 3*Nc);
for i = 1:Nc
    for j = 1:Nc
        if i>=j
            K(1+3*(i-1):3*i, 1+3*(j-1):3*j) = eye(3);
        end
    end
end
Acons(1+6*Nc:9*Nc, :) = K;
Acons(1+9*Nc:12*Nc, :) = -1 * K;
% 得到二次型优化问题
% J = delta_u_list' * E * delta_u_list + F * delta_u_list
% Acons * delta_u_list <= b
options = optimoptions('quadprog','MaxIterations', 1000);
Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = [];
delta_u_list = quadprog(E, F', Acons, b, Aeq,beq,lb,ub,x0,options);

% delta_u_list = 0 .* delta_u_list;
% ------------------------------- 控制序列 ---------------------------------------
% delta_u_list = zeros(3*Nc, 1);
% delta_u_list(1:3, 1) = delta_u_dq;

% ------------------------------- 预测 ---------------------------------------
Y_dq_list = P * xm + H * delta_u_list + Hd * delta_w_list;

% ------------------------------- iMBC ---------------------------------------
Y_dq = reshape(Y_dq_list, [size(yop, 1), size(Y_dq_list, 1)/size(yop, 1)]);
delta_u_list = reshape(delta_u_list, [size(D, 2), size(delta_u_list, 1)/size(D, 2)]);

time_list = DT*DT_gap*(0:Np);

Y = zeros(size(D, 1), Np+1);
u_list = zeros(3, Np);

RotSpeed = Y_dq(2, :) * 2*pi/60;  % rad/s
Y(:, 1) = yop;
Y(:, 2) = iTransverse_y(Azimuth, n_FixFrameOutputs, n_RotTripletOutputs, new_seq_out, Y_dq(:, 1));

u_list(:, 1) = uop;
u_list(:, 2) = delta_u_list(:, 1) + uop;

for i = 3:Np+1
    if i > Nc
        u_list(:, i) = u_list(:, i-1);
    else
        u_list(:, i) = delta_u_list(:, i-1) + u_list(:, i-1);
    end
end
Azimuth_total = Azimuth;
for i = 3:Np+1
    Azimuth_total = Azimuth_total + 0.5*(RotSpeed(i-2)+RotSpeed(i-1))*(DT*DT_gap)*180/pi;
    Azimuth = rem(Azimuth_total, 360); % deg
    Y(:, i) = iTransverse_y(Azimuth, n_FixFrameOutputs, n_RotTripletOutputs, new_seq_out, Y_dq(:, i-1));
end

% ------------------------------- 画图 ---------------------------------------
figure(1)
plot(time_list, Y(1, :))
title('GenSpeed/rpm')

figure(2)
plot(time_list, Y(3, :))
hold on
plot(time_list, Y(4, :))
hold on
plot(time_list, Y(5, :))
title('RootMyc/kn-m')
legend('RootMyc1', 'RootMyc2', 'RootMyc3')



figure(3)
plot(time_list, rad2deg(u_list(1, :)))
hold on
plot(time_list, rad2deg(u_list(2, :)))
hold on
plot(time_list, rad2deg(u_list(3, :)))
title('BldPitch/deg')
legend('BldPitch1', 'BldPitch2', 'BldPitch3')



%% 用到的函数
function Tmbc  = mbc(Azimuth)
    % compute azimuth positions of blades:
    az = Azimuth*pi/180.0 + 2*pi/3* (0:(3-1)) ; % Eq. 1, azimuth in radians   3个叶片的方位角（弧度），叶片1的方位角为初始方位角
    
    % compute transformation matrices
    cos_col = cos(az(:));  % 列向量
    sin_col = sin(az(:));
    
    Tmbc = get_tt_inverse(sin_col, cos_col);     % inverse of tt (computed analytically in function below) 逆矩阵
end

function delat_u_dq  = Transverse(Azimuth, n_FixFrameInputs, n_RotTripletInputs, delta_u, new_seq_inp)
    % compute azimuth positions of blades:
    az = Azimuth*pi/180.0 + 2*pi/3* (0:(3-1)) ; % Eq. 1, azimuth in radians   3个叶片的方位角（弧度），叶片1的方位角为初始方位角
    
    % compute transformation matrices
    cos_col = cos(az(:));  % 列向量
    sin_col = sin(az(:));
    
    ttv = get_tt_inverse(sin_col, cos_col);     % inverse of tt (computed analytically in function below) 逆矩阵
    
    T1cv = zeros(n_FixFrameInputs+3*n_RotTripletInputs, n_FixFrameInputs+3*n_RotTripletInputs);
    
    T1cv(1:n_FixFrameInputs, 1:n_FixFrameInputs) = eye(n_FixFrameInputs);                % Eq. 21  T1c的逆矩阵
    
    for ii = 1:n_RotTripletInputs
        T1cv(n_FixFrameInputs+1+(ii-1)*3 : n_FixFrameInputs+ii*3, n_FixFrameInputs+1+(ii-1)*3 : n_FixFrameInputs+ii*3) = ttv;
    end
        
    
    %% 对u进行转换
    
    delat_u_dq = delta_u .* 0;
    delat_u_dq(new_seq_inp, :) = T1cv * delta_u(new_seq_inp);

end

function delat_u  = iTransverse(Azimuth, n_FixFrameInputs, n_RotTripletInputs, delta_u_dq, new_seq_inp)
    % compute azimuth positions of blades:
    az = Azimuth*pi/180.0 + 2*pi/3* (0:(3-1)) ; % Eq. 1, azimuth in radians   3个叶片的方位角（弧度），叶片1的方位角为初始方位角
    
    % compute transformation matrices
    cos_col = cos(az(:));  % 列向量
    sin_col = sin(az(:));
    
    tt  = [ones(3,1), cos_col, sin_col];        % Eq. 9, t_tilde
    
    T1c = zeros(n_FixFrameInputs+3*n_RotTripletInputs, n_FixFrameInputs+3*n_RotTripletInputs);
    
    T1c(1:n_FixFrameInputs, 1:n_FixFrameInputs) = eye(n_FixFrameInputs);                % Eq. 21  T1c的逆矩阵
    
    for ii = 1:n_RotTripletInputs
        T1c(n_FixFrameInputs+1+(ii-1)*3 : n_FixFrameInputs+ii*3, n_FixFrameInputs+1+(ii-1)*3 : n_FixFrameInputs+ii*3) = tt;
    end
        
    
    %% 对u进行转换
    
    delat_u = delta_u_dq .* 0;
    delat_u(new_seq_inp, :) = T1c * delta_u_dq(new_seq_inp, :);

end

function y  = iTransverse_y(Azimuth, n_FixFrameOutputs, n_RotTripletOutputs, new_seq_out, y_dq)
    % compute azimuth positions of blades:
    az = Azimuth*pi/180.0 + 2*pi/3* (0:(3-1)) ; % Eq. 1, azimuth in radians   3个叶片的方位角（弧度），叶片1的方位角为初始方位角
    
    % compute transformation matrices
    cos_col = cos(az(:));  % 列向量
    sin_col = sin(az(:));
    
    tt  = [ones(3,1), cos_col, sin_col];        % Eq. 9, t_tilde
    
    
    T1o = zeros(n_FixFrameOutputs+3*n_RotTripletOutputs, n_FixFrameOutputs+3*n_RotTripletOutputs);
    
    T1o(1:n_FixFrameOutputs, 1:n_FixFrameOutputs) = eye(n_FixFrameOutputs);              % Tlo (Eq. 23)  输出相关
    
    for ii = 1:n_RotTripletOutputs
        T1o(n_FixFrameOutputs+1+(ii-1)*3 : n_FixFrameOutputs+ii*3, n_FixFrameOutputs+1+(ii-1)*3 : n_FixFrameOutputs+ii*3) = tt;
    end
    
    y = y_dq .* 0;
    y(new_seq_out, :) = T1o * y_dq(new_seq_out);

end


function y_dq  = Transverse_y(Azimuth, n_FixFrameOutputs, n_RotTripletOutputs, new_seq_out, y)
    % compute azimuth positions of blades:
    az = Azimuth*pi/180.0 + 2*pi/3* (0:(3-1)) ; % Eq. 1, azimuth in radians   3个叶片的方位角（弧度），叶片1的方位角为初始方位角
    
    % compute transformation matrices
    cos_col = cos(az(:));  % 列向量
    sin_col = sin(az(:));

    ttv = get_tt_inverse(sin_col, cos_col);     % inverse of tt (computed analytically in function below) 逆矩阵
    
    T1ov = zeros(n_FixFrameOutputs+3*n_RotTripletOutputs, n_FixFrameOutputs+3*n_RotTripletOutputs);
    
    T1ov(1:n_FixFrameOutputs, 1:n_FixFrameOutputs) = eye(n_FixFrameOutputs);              % Tlo (Eq. 23)  输出相关
    
    for ii = 1:n_RotTripletOutputs
        T1ov(n_FixFrameOutputs+1+(ii-1)*3 : n_FixFrameOutputs+ii*3, n_FixFrameOutputs+1+(ii-1)*3 : n_FixFrameOutputs+ii*3) = ttv;
    end
    
    y_dq = y .* 0;
    y_dq(new_seq_out, :) = T1ov * y(new_seq_out);

end



function [ttv] = get_tt_inverse(sin_col, cos_col)
    c1 = cos_col(1);
    c2 = cos_col(2);
    c3 = cos_col(3);
    
    s1 = sin_col(1);
    s2 = sin_col(2);
    s3 = sin_col(3);

    
    ttv = [ c2*s3 - s2*c3,  c3*s1 - s3*c1, c1*s2 - s1*c2
               s2 - s3 ,       s3 - s1,       s1 - s2
               c3 - c2 ,       c1 - c3,       c2 - c1 ] / (1.5*sqrt(3));

    return
end