% 风力机在工作点处的线性化模型
%% 准备工作
clc
clear; close all
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

%% 迭代
DT = 1e-3;  % 积分步长，一定要足够小，不然积分会发散
Time_total = 10;  % s
WindSpeed = 14; % m/s
BldPitch = 8.47 * 1/(57.2957) *  ones(3, 1); %rad
xop = StateSpaceList.xop_total(:, :, OP_num);
uop = StateSpaceList.uop_total(:, :, OP_num);
wop = StateSpaceList.wop_total(:, :, OP_num);
yop_list = StateSpaceList.yop_List_total(:, :, OP_num);
y_out = zeros(6, Time_total/DT);

delta_xdq = xop .* 0;
delta_xdq(11) = 0;
Azimuth_total = 0;
time_now = 0;

for i = 1:(Time_total/DT)
    % 风速阶跃
    if time_now > 100
        WindSpeed = 15;
    end

    %% 生成控制动作u的偏差量
    u = uop;
    w = wop;
    w(1) = WindSpeed;
    u(1) = BldPitch(1);
    u(2) = BldPitch(2);
    u(3) = BldPitch(3);
    delta_u = u - uop;
    delta_w = w - wop;
    
    %% 生成当前的方位角
    RotSpeed = xop(11) + delta_xdq(11);  % rad/s
    Azimuth_total = Azimuth_total + RotSpeed * DT;
    Azimuth = rem(Azimuth_total, 2*pi) * (1/pi*180); % deg
    
    %% 根据方位角选择输出量的工作点
    yop  = SwitchOP(Azimuth, Azimuth_list, yop_list);
    
    %% 对控制动作u进行MBC变换
    delta_u_dq = Transverse(Azimuth, n_FixFrameInputs, n_RotTripletInputs, delta_u, new_seq_inp);
    delta_w_dq = delta_w;

    %% 利用状态空间进行迭代计算，默认工作点处delta_xdq的初始值为0
    delta_xdq_dot = A * delta_xdq + B * delta_u_dq + Bd * delta_w_dq;
    delta_y_dq = C * delta_xdq + D * delta_u_dq + Dd * delta_w_dq;
    delta_xdq = delta_xdq + DT * delta_xdq_dot;

    %% 对y_dq进行反MBC变换
    y  = iTransverse_y(Azimuth, n_FixFrameOutputs, n_RotTripletOutputs, new_seq_out, delta_y_dq) + yop;
    
    y_out(1, i) = time_now;
    y_out(2:6, i) = y;
    time_now = time_now + DT;
    
end
figure(1)
plot(y_out(1,:), y_out(4,:))
xlabel('Time/s')
ylabel('RootMyc1')

figure(2)
plot(y_out(1,:), y_out(2,:))
xlabel('Time/s')
ylabel('GenSpeed/rpm')







%% 用到的函数
function yop  = SwitchOP(Azimuth, Azimuth_list, yop_36)
    left = 1;
    right = length(Azimuth_list);
    
    while(left<right)
        mid = floor((left + right) / 2);
        if Azimuth_list(mid) == Azimuth
            index = mid;
            break
        elseif Azimuth_list(mid) < Azimuth
            left = mid + 1;
        else
            right = mid - 1;
        end
    end
    % If target is not found, return the index of the interval in which it falls
    if Azimuth < Azimuth_list(1)
        index = 1;
        % yop = yop_36(:, end) + (yop_36(:, 1)-yop_36(:, end))*(Azimuth+360-Azimuth_list(end))/(360/length(Azimuth_list));
    elseif Azimuth >= Azimuth_list(end)
        index = 36;
        % yop = yop_36(:, end) + (yop_36(:, 1)-yop_36(:, end))*(Azimuth-Azimuth_list(end))/(360/length(Azimuth_list));
    else
        index = left;
        % if index == 36
        %     index = index - 1;
        % end
        % yop = yop_36(:, index) + (yop_36(:, index+1)-yop_36(:, index))*(Azimuth-Azimuth_list(index))/(360/length(Azimuth_list));
    end
    yop = yop_36(:, index);

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