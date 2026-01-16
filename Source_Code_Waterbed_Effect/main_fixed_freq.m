clear all; clc; close all;
%%
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

w = 10.^[-8:0.001:8];

disturb_streng = 0.1;
disturb_direct = 80/180*pi;

%% default
[eig_S, eig_matAc, eig_matAcn, eig_Sn, t_impulse, A_impulse, V_impulse] = test(disturb_streng, disturb_direct);

figure(6532);
subplot(2,2,1); hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,1))), 'LineWidth', 3, 'Color', [0.6, 0.6, 0.6]);
subplot(2,2,2);  hold on; grid off; box on;
plot(log10(w), log(abs(eig_S(:,2))), 'LineWidth', 3, 'Color', [0.6, 0.6, 0.6]);
subplot(2,2,3);  hold on; grid off; box on;
plot(t_impulse, 180/pi * (As + disturb_streng * A_impulse), 'LineWidth', 1.2, 'Color', [0.6, 0.6, 0.6]);
subplot(2,2,4);  hold on; grid off; box on;
plot(t_impulse, Vs + disturb_streng * V_impulse, 'LineWidth', 1.2, 'Color', [0.6, 0.6, 0.6]);



%% scaled
d_seq = 10.^[0.5, 2, 4, 6];

linestyle_seq = {'-', '--', '.', ':.'};

for rr = 1:length(d_seq)
    d = d * d_seq(rr);
    [eig_S, eig_matAc, eig_matAcn, eig_Sn, t_impulse, A_impulse, V_impulse] = test(disturb_streng, disturb_direct);
    d = d / d_seq(rr);

    fprintf('max(abs(eig_S))=%.6f at %.2f * d\n', max(abs(eig_S(:,2))), d_seq(rr))

    if rr == 1
        w_plot = w;
    else
        w_plot = w(1:200:end);
        eig_S = eig_S(1:200:end, :);
    end

    figure(6532);
    subplot(2,2,1); hold on; grid off; box on;
    plot(log10(w_plot), log(abs(eig_S(:,1))), linestyle_seq{rr}, 'LineWidth', 1.5, 'Color', [0.85, 0.3250, 0.0980]);
    xlim([-4,3])

    subplot(2,2,2);  hold on; grid off; box on;
    plot(log10(w_plot), log(abs(eig_S(:,2))), linestyle_seq{rr}, 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
    xlim([-3,3])

    if rr == length(d_seq)
        subplot(2,2,3);  hold on; grid off; box on;
        plot(t_impulse, 180/pi * (As + disturb_streng * A_impulse), 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);
        xlim([0, 20]);

        subplot(2,2,4);  hold on; grid off; box on;
        plot(t_impulse, Vs + disturb_streng * V_impulse, 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);
        xlim([0, 20]); ylim([1.02, 1.04]);
    end

end

figure(6532);
subplot(2,2,1);
legend({'$\rm default$', '$10^{0.5}\cdot d$', '$10^{2}\cdot d$', '$10^{4}\cdot d$', '$10^{6}\cdot d$'}, 'Interpreter', 'latex', 'Location', 'southeast')
subplot(2,2,2);
legend({'$\rm default$', '$10^{0.5}\cdot d$', '$10^{2}\cdot d$', '$10^{4}\cdot d$', '$10^{6}\cdot d$'}, 'Interpreter', 'latex', 'Location', 'northwest')


%%
function [eig_S, eig_matAc, eig_matAcn, eig_Sn, t_impulse, A_impulse, V_impulse] = test(disturb_streng, disturb_direct)

global M d tau Dq X R Z2 Vs As Tfl w;

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
Tf_sys = 1/(Tfl*s+1);

G_pa = 1/(s*(M*s+d));
G_qa = 0;
G_pv = 0;
G_qv = Dq/(tau*s+1) * Tf_sys;

G_sys = [G_pa, G_qa; G_pv, G_qv];

% sensitivity transfer function
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

det_S = zeros(length(w), 1);
det_Sn = zeros(length(w), 1);

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

% time-domain simulation
t_impulse = 0:0.001:20;
disturb_vector = disturb_streng * [cos(disturb_direct); sin(disturb_direct)];

A_impulse = impulse(disturb_vector(1)*S_sys(1,1)+disturb_vector(2)*S_sys(1,2), t_impulse);
V_impulse = impulse(disturb_vector(1)*S_sys(2,1)+disturb_vector(2)*S_sys(2,2), t_impulse);


end