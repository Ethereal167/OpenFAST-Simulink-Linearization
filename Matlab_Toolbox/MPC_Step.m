function u_next = MPC_Step(A, B, C, D, Bd, Dd, xdq_init, u_init, Np, Nc, DT, DT_gap, delta_w_list, Rs, Q, R, Azimuth, RotSpeed)
    
    %% 差分方程的矩阵
    Am = [DT*A+eye(size(A, 1)), zeros(size(A, 1), size(C, 1));
          C, eye(size(C, 1))];
    Bm = [DT*B; D];
    Bmd = [DT*Bd; Dd];
    Cm = [zeros(size(C, 1), size(A, 1)), eye(size(C, 1))];
    
    %% 预测方程的矩阵
    % 计算矩阵Timbc
    % M = eye(size(D, 1));
    % M(3:5, 3:5) = imbc_in;
    % Timbc = zeros(size(D, 1)*Np, size(D, 1)*Np);
    % for i = 1:Np
    %     Timbc(1+size(D, 1)*(i-1):size(D, 1)*i, 1+size(D, 1)*(i-1):size(D, 1)*i) = M;
    % end
    % ---------------------------------- 计算矩阵P ----------------------------
    for i = 1:Np
        if i == 1
            P = Cm*(Am^1);
        else
            P = [P; Cm*(Am^(1+DT_gap*(i-1)))];
        end
    end
    % P = Timbc * P;
    % ---------------------------------- 计算矩阵H ----------------------------
    Azimuth_list = Azimuth * ones(1, Nc);
    for i = 2:Nc
        Azimuth_list(i) = rem(Azimuth_list(i-1) + rad2deg(RotSpeed * DT * DT_gap), 360);
    end
    
    for i = 1:Np
        H_c = [];
        for j = 1:Nc
            if j==1
                H_c = Cm * Am^(DT_gap*(i-1)) * Bm * mbc(Azimuth_list(j));
            else
                if j > i
                    H_c = [H_c, zeros(size(Cm, 1), size(Bm, 2))];
                else
                    H_c = [H_c, Cm * Am^(DT_gap*(i-j)) * Bm * mbc(Azimuth_list(j))];
                end
            end
        end
        if i == 1
            H = H_c;
        else
            H = [H; H_c];
        end
    end
    % H = Timbc*H;
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
    % Hd = Timbc*Hd;
    %% ------------------------------- 求解最优控制序列 ---------------------------------------
    % 控制输入的约束
    delta_u_max = [2; 2; 2] * (pi / 180) * (DT * DT_gap); % 变桨速率约束的上限 rad
    delta_u_min = delta_u_max; % 变桨速率约束的下限 rad
    u_max = [90; 90; 90] * (pi / 180); % 变桨约束上限 rad
    u_min = [0; 0; 0] * (pi / 180); % 变桨约束下限 rad
    % 构造二次型的矩阵
    % E
    E = H' * Q * H + R;
    % F
    F = 2 * (P * xdq_init + Hd * delta_w_list - Rs)' * Q * H;
    % 约束相关矩阵
    b = zeros(4*Nc*3, 1);
    for i = 1:Nc
        b(1+3*(i-1):3*i, 1) = delta_u_max;
    end
    for i = 1:Nc
        b(1+3*Nc+3*(i-1):3*Nc+3*i, 1) = delta_u_min;
    end
    for i = 1:Nc
        b(1+6*Nc+3*(i-1):6*Nc+3*i, 1) = u_max-u_init;
    end
    for i = 1:Nc
        b(1+9*Nc+3*(i-1):9*Nc+3*i, 1) = u_min+u_init;
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
    options = optimoptions('quadprog','MaxIterations', 1000, 'Display','off', 'Algorithm', 'active-set');
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    x0 = zeros(size(u_init, 1)*Nc, 1);
    [delta_u_list, ~, exitflag] = quadprog(E, F', Acons, b, Aeq,beq,lb,ub,x0,options);
    if exitflag==-2
        % 即无解
        u_next = u_init;
    else
        u_next = delta_u_list(1:3, 1) + u_init;
    end

end



function Tmbc  = mbc(Azimuth)
    % compute azimuth positions of blades:
    az = Azimuth*pi/180.0 + 2*pi/3* (0:(3-1)) ; % Eq. 1, azimuth in radians   3个叶片的方位角（弧度），叶片1的方位角为初始方位角
    
    % compute transformation matrices
    cos_col = cos(az(:));  % 列向量
    sin_col = sin(az(:));
    
    Tmbc = get_tt_inverse(sin_col, cos_col);     % inverse of tt (computed analytically in function below) 逆矩阵
end

function Timbc = imbc(Azimuth)
    az = Azimuth*pi/180.0 + 2*pi/3* (0:(3-1)) ; % Eq. 1, azimuth in radians   3个叶片的方位角（弧度），叶片1的方位角为初始方位角
    
    % compute transformation matrices
    cos_col = cos(az(:));  % 列向量
    sin_col = sin(az(:));
    
    Timbc  = [ones(3,1), cos_col, sin_col];        % Eq. 9, t_tilde
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