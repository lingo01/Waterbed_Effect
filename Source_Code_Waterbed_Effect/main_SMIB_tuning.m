clear all; clc; close all;
%% system configuration
global M d tau Dq X R Z2 Vs As Tfl w;

M = 0.3;
d = 0.15;

tau = 0.05;
Dq = 0.8;

X = 0.30;
R = 0.20;
Z2 = X^2+R^2;

Vs = 1.03;
As = 21.3*pi/180;

Tfl = 0.01;

%% default
[eig_S, eig_matAc, eig_matAcn, w] = test();

figure(6524);
subplot(2,2,1); hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,1))), 'LineWidth', 3, 'Color', [0.6, 0.6, 0.6]);
subplot(2,2,2);  hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,2))), 'LineWidth', 3, 'Color', [0.6, 0.6, 0.6]);
subplot(2,2,3); hold on; grid off; box on;
plot(eig_matAc+1e-15j, 's', 'LineWidth', 1.2, 'Color', [0.6, 0.6, 0.6]);
subplot(2,2,4);  hold on; grid off; box on;
plot(eig_matAc+1e-15j, 's', 'LineWidth', 1.2, 'Color', [0.6, 0.6, 0.6]);

%% 2d
d = 2*d;
[eig_S, eig_matAc, eig_matAcn, w] = test();
d = d/2;

figure(6524);
subplot(2,2,1); hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,1))), '-', 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);
subplot(2,2,2);  hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,2))), '-', 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);
subplot(2,2,3); hold on; grid off; box on;
plot(eig_matAc+1e-15j, 'v', 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);

%% Dq/5
Dq = Dq/5;
[eig_S, eig_matAc, eig_matAcn, w] = test();
Dq = Dq*5;

figure(6524);
subplot(2,2,1); hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,1))), '--', 'LineWidth', 1.5, 'Color', [0.85, 0.3250, 0.0980]);
subplot(2,2,2);  hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,2))), '--', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
subplot(2,2,4); hold on; grid off; box on;
plot(eig_matAc+1e-15j, '^', 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);

%%
figure(6524);
subplot(2,2,1); 
xlim([0,1.5]); ylim([-0.5, 2.05])
legend({'$\rm default$', '$2\cdot d$', '$5\cdot k_q$'}, 'Interpreter', 'latex')
subplot(2,2,2); 
xlim([0,4]); ylim([-1, 0.15])
legend({'$\rm default$', '$2\cdot d$', '$5\cdot k_q$'}, 'Interpreter', 'latex')
subplot(2,2,3); 
xlim([-1,0]);
legend({'$\rm default$', '$2\cdot d$'}, 'Interpreter', 'latex')
subplot(2,2,4); 
xlim([-0.3,-0.2]);
legend({'$\rm default$', '$5\cdot k_q$'}, 'Interpreter', 'latex')

%%
function [eig_S, eig_matAc, eig_matAcn, w] = test()

global M d tau Dq X R Z2 Vs As Tfl;

% network transfer function
H_pa = Vs/Z2 * ( R*sin(As) + X*cos(As) );
H_pv = 1/Z2 * ( R*(2*Vs-cos(As)) + X*sin(As) );
H_qa = Vs/Z2 * ( X*sin(As) - R*cos(As) );
H_qv = 1/Z2 * ( X*(2*Vs-cos(As)) - R*sin(As) );

H_sys = [H_pa, H_pv; H_qa, H_qv];
Hn_sys = [H_pa, 0; 0, H_qv];

% state space model simulation
matA = [
0       1       0       0;
0       -d/M    0       0;
0       0       0       1;
0       0       -1/tau/Tfl      -(1/Tfl+1/tau);
];

matB = [
0       0;
1/M     0;
0       0;
0       Dq/tau/Tfl
];

matB2 = [
0       0;
-1/M    0;
0       0;
0       -Dq/tau/Tfl
];

matK2 = [
H_pa     0       H_pv       0;
H_qa     0       H_qv       0;
];

matK2n = [
H_pa     0       0          0;
0        0       H_qv       0;
];

% space space model: dx/dt = (A+B2K2)*x + B*u = Ac*x+B*u
matAc = matA + matB2 * matK2;
matAcn = matA + matB2 * matK2n;

eig_matAc = eig(matAc);
eig_matAcn = eig(matAcn);

% device transfer function
s = tf('s');

Kvp = 0.5;
Kvi = 5.0;
Cf = 1.59155e-6;
% Tf_sys = 1/(Tfl*s+1);
Tf_sys = (Kvp*s+Kvi) / (Cf*s^2+Kvp*s+Kvi);

G_pa = 1/(s*(M*s+d));
G_qa = 0;
G_pv = 0;
G_qv = Dq/(tau*s+1) * Tf_sys;

G_sys = [G_pa, G_qa; G_pv, G_qv];

w = 10.^[-8:0.001:8];

% sensitivity transfer function
S_sys = inv( eye(2) + G_sys * H_sys );

S_sys_arr = zeros(length(w), 2, 2);
for ii = 1:length(w)
    S_sys_arr(ii,:,:) = evalfr(S_sys, 1i*w(ii));
end

eig_S = zeros(length(w), 2);

for ii = 1:length(w)
    D_S = func_eig_modi('matS', squeeze(S_sys_arr(ii,:,:)));
    
    % eigenvalue
    eig_S(ii,:) = D_S';
end

% eigenvalue tracking
eig_S = func_track_eigenvalues(eig_S);

% directed area
area_eig_S_ang = func_area_bode(w, log10(abs(eig_S(:,1))));
fprintf('sensitivity integral: int log(eig_ang(S)) = %.6f\n', area_eig_S_ang);

max_eig_S = max(max(abs(eig_S)));
fprintf('max sensitivity S = %.4f\n', max_eig_S);

fprintf('\n')

end