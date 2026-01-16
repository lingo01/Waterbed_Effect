clear all; clc; close all;
%% system configuration
M = 0.3;
d = 0.15*1;

tau = 0.05;
Dq = 0.8;

X = 0.30;
R = 0.20;
Z2 = X^2+R^2;

Vs = 1.03;
As = 21.3*pi/180;

s = tf('s');

Tfl = 0.01;
Kvp = 0.5;
Kvi = 5.0;
Cf = 1.59155e-6;
% Tf_sys = 1/(Tfl*s+1);
Tf_sys = (Kvp*s+Kvi) / (Cf*s^2+Kvp*s+Kvi);

w = 10.^[-8:0.001:8];

kp = 1; ki = 50; Cf = 0.05;
G_vl = (kp + ki/s) * 1/(Cf * s);
G_vl = G_vl / (1+G_vl);

%% device transfer function

G_pa = 1/(s*(M*s+d));
G_qa = 0;
G_pv = 0;
G_qv = Dq/(tau*s+1) * Tf_sys * G_vl;

G_sys = [G_pa, G_qa; G_pv, G_qv];

%% network transfer function
H_pa = Vs/Z2 * ( R*sin(As) + X*cos(As) );
H_pv = 1/Z2 * ( R*(2*Vs-cos(As)) + X*sin(As) );
H_qa = Vs/Z2 * ( X*sin(As) - R*cos(As) );
H_qv = 1/Z2 * ( X*(2*Vs-cos(As)) - R*sin(As) );

H_sys = [H_pa, H_pv; H_qa, H_qv];
Hn_sys = [H_pa, 0; 0, H_qv];

%% transfer function model simulation
T_sys = inv(eye(2)+G_sys*H_sys) * G_sys;
Tn_sys = inv(eye(2)+G_sys*Hn_sys) * G_sys;

t_impulse = 0:0.001:20;
disturb_streng = 0.1;

disturb_direct = 80/180*pi;
disturb_vector = [cos(disturb_direct); sin(disturb_direct)];
A_impulse = impulse(disturb_vector(1)*T_sys(1,1)+disturb_vector(2)*T_sys(1,2), t_impulse);
V_impulse = impulse(disturb_vector(1)*T_sys(2,1)+disturb_vector(2)*T_sys(2,2), t_impulse);
An_impulse = impulse(disturb_vector(1)*T_sys(1,1), t_impulse);
Vn_impulse = impulse(disturb_vector(2)*T_sys(2,2), t_impulse);

%% state space model simulation
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

matAcn_a = matAcn(1:2, 1:2);
matAcn_v = matAcn(3:4, 3:4);

eig_matAc = eig(matAc);
eig_matAcn = eig(matAcn);

%% sensitivity transfer function
S_sys = inv( eye(2) + G_sys * H_sys );
Sn_sys = inv( eye(2) + G_sys * Hn_sys );

S_sys_arr = zeros(length(w), 2, 2);
Sn_sys_arr = zeros(length(w), 2, 2);
for ii = 1:length(w)
    S_sys_arr(ii,:,:) = evalfr(S_sys, 1i*w(ii));
    Sn_sys_arr(ii,:,:) = evalfr(Sn_sys, 1i*w(ii));
end

eig_S = zeros(length(w), 2);
eig_Sn = zeros(length(w), 2);

for ii = 1:length(w)
    D_S = func_eig_modi('matS', squeeze(S_sys_arr(ii,:,:)));
    D_Sn = func_eig_modi('matS', squeeze(Sn_sys_arr(ii,:,:)));
  
    % eigenvalue
    eig_S(ii,:) = D_S';
    eig_Sn(ii,:) = D_Sn';

    % determinant
    det_S(ii) = det(squeeze(S_sys_arr(ii,:,:)));
    det_Sn(ii) = det(squeeze(Sn_sys_arr(ii,:,:)));
end

% eigenvalue tracking
eig_S = func_track_eigenvalues(eig_S);
eig_Sn = func_track_eigenvalues(eig_Sn);

%% open-loop transfer function
L_sys = G_sys * H_sys;
Ln_sys =  G_sys * Hn_sys;

L_sys_arr = zeros(length(w), 2, 2);
Ln_sys_arr = zeros(length(w), 2, 2);

for ii = 1:length(w)
    L_sys_arr(ii,:,:) = evalfr(L_sys, 1i*w(ii));
    Ln_sys_arr(ii,:,:) = evalfr(Ln_sys, 1i*w(ii));
end

eig_L = zeros(length(w), 2);
eig_Ln = zeros(length(w), 2);

for ii = 1:length(w)
    eig_L(ii,:) = eig(squeeze(L_sys_arr(ii,:,:)))';
    eig_Ln(ii,:) = eig(squeeze(Ln_sys_arr(ii,:,:)))';
end

%% directed integral area
area_S = func_area_bode(w, log10(det_S));
area_Sn = func_area_bode(w, log10(det_Sn));

fprintf('sensitivity integral: int log(det(S)) = %.10f\n', area_S);
fprintf('sensitivity integral: int log(det(Sn)) = %.10f\n', area_Sn);

area_eig_S_1 = func_area_bode(w, log10(abs(eig_S(:,1))));
area_eig_S_2 = func_area_bode(w, log10(abs(eig_S(:,2))));

area_eig_Sn_1 = func_area_bode(w, log10(abs(eig_Sn(:,1))));
area_eig_Sn_2 = func_area_bode(w, log10(abs(eig_Sn(:,2))));

fprintf('sensitivity integral: int log(eig_1(S)) = %.6f\t log(eig_2(S)) = %.6f\n', area_eig_S_1, area_eig_S_2);
fprintf('sensitivity integral: int log(eig_1(Sn)) = %.6f\t log(eig_2(Sn)) = %.6f\n', area_eig_Sn_1, area_eig_Sn_2);

%% results plot
figure(1);
subplot(2,3,1); hold on;  grid off; box on;
plot(eig_matAc+1e-15j, '*', 'LineWidth', 1.2);
plot(eig(matAcn_a), 's', 'LineWidth', 1.2);
plot(eig(matAcn_v), 's', 'LineWidth', 1.2);
legend('H_{net}', 'H_{net}^{(n)}(A)', 'H_{net}^{(n)}(V)');

subplot(2,3,4);
hold on; box on;
plot(real(eig_L(:,1)), imag(eig_L(:,1)), 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);   
plot(real(eig_L(:,1)), -imag(eig_L(:,1)), 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);
plot(real(eig_L(:,2)), imag(eig_L(:,2)), 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);   
plot(real(eig_L(:,2)), -imag(eig_L(:,2)), 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);
plot(-1, 0, '*', 'LineWidth', 1.2)

plot(real(eig_Ln(:,1)), imag(eig_Ln(:,1)),'--', 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);    
plot(real(eig_Ln(:,1)), -imag(eig_Ln(:,1)),'--', 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);
plot(real(eig_Ln(:,2)), imag(eig_Ln(:,2)),'--', 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);    
plot(real(eig_Ln(:,2)), -imag(eig_Ln(:,2)),'--', 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);
plot(-1, 0, '*', 'LineWidth', 1.2)
ylim([-1.5, 1.5])

subplot(2,3,2); hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,1))), 'LineWidth', 1.2);
plot(log10(w), log(abs(eig_Sn(:,1))), 'LineWidth', 1.2);
legend('S', 'S_n');
xlabel('$\log(\omega)$', 'Interpreter','latex');  ylabel('$\log(|\sigma|)$', 'Interpreter','latex');
xlim([-1,4]);
ylim([-5, 3])

subplot(2,3,3);  hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,2))), 'LineWidth', 1.2);
plot(log10(w), log(abs(eig_Sn(:,2))), 'LineWidth', 1.2);
legend('S', 'S_n');
xlabel('$\log(\omega)$', 'Interpreter','latex');  ylabel('$\log(|\sigma|)$', 'Interpreter','latex');
xlim([-1,4])

subplot(2,3,5); 
hold on; grid off; box on;
plot(t_impulse, 180/pi * (As + disturb_streng * A_impulse), 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);
plot(t_impulse, 180/pi * (As + disturb_streng * An_impulse), 'LineWidth', 1.2, 'Color', [0.5, 0.5, 0.5]);
xlabel('t/s', 'FontName', 'Times New Roman')
ylabel('$\theta$/deg', 'FontName', 'Times New Roman', 'Interpreter','latex')

subplot(2,3,6); hold on; grid off; box on;
plot(t_impulse, Vs + disturb_streng * V_impulse, 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);
plot(t_impulse, Vs + disturb_streng * Vn_impulse, 'LineWidth', 1.2, 'Color', [0.5, 0.5, 0.5]);
xlabel('t/s', 'FontName', 'Times New Roman')
ylabel('V/pu', 'FontName', 'Times New Roman')
ylim([1.02, 1.05])