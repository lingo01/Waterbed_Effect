function sys_cl = func_tf2ss_close(Gdev, Hnet)
% Assume Gdev and Hnet are defined as state-space models or transfer-function matrices
% Convert to state-space form
sys_Gdev = ss(Gdev);
sys_Hnet = ss(Hnet);

% Series interconnection: GdevHnet = Gdev * Hnet
sys_GdevHnet = sys_Gdev * sys_Hnet;

% Extract state-space matrices
A = sys_GdevHnet.A;
B = sys_GdevHnet.B;
C = sys_GdevHnet.C;
D = sys_GdevHnet.D;

n = size(D, 1); % system dimension
I = eye(n);

% Check invertibility of (I + D)
if rank(I + D) < n
    error('Matrix (I + D) is singular; cannot form the closed-loop system.');
end

inv_ID = inv(I + D); % inverse of (I + D)

% Construct the closed-loop state-space matrices
A_cl = A - B * inv_ID * C;
B_cl = B * inv_ID;
C_cl = -inv_ID * C;
D_cl = inv_ID;

% Create closed-loop state-space model
sys_cl = ss(A_cl, B_cl, C_cl, D_cl);
end