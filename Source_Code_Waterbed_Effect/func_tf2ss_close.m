function sys_cl = func_tf2ss_close(Gdev, Hnet)
% 假设 Gdev 和 Hnet 已定义为状态空间模型或传递函数矩阵
% 转换为状态空间模型
sys_Gdev = ss(Gdev);
sys_Hnet = ss(Hnet);

% 计算GdevHnet = Gdev * Hnet的串联系统
sys_GdevHnet = sys_Gdev * sys_Hnet;

% 提取状态空间矩阵
A = sys_GdevHnet.A;
B = sys_GdevHnet.B;
C = sys_GdevHnet.C;
D = sys_GdevHnet.D;

n = size(D, 1); % 系统维度
I = eye(n);

% 检查I + D是否可逆
if rank(I + D) < n
    error('I + D 矩阵不可逆，无法计算闭环系统。');
end

inv_ID = inv(I + D); % 计算逆矩阵

% 构造闭环系统的状态空间矩阵
A_cl = A - B * inv_ID * C;
B_cl = B * inv_ID;
C_cl = -inv_ID * C;
D_cl = inv_ID;

% 创建闭环状态空间模型
sys_cl = ss(A_cl, B_cl, C_cl, D_cl);
end