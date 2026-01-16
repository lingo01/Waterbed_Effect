clear all; clc; close all;
tic;
%% system configurations
vpa('1e400');

casename = 'case14';

w = 10.^[-9:0.002:10];

% network   
load(strcat('./data&figure/Hnet_info_', casename,'.mat'));
Hnet = Hnet_info.Hnet;
Hn_net = blkdiag(Hnet(1:size(Hnet,1)/2, 1:size(Hnet,1)/2), Hnet(size(Hnet,1)/2+1:end, size(Hnet,1)/2+1:end));

% device
dev_info = Hnet_info.dev_info;

s = tf('s');
Gdev = s*zeros(2*dev_info.dev.num, 2*dev_info.dev.num);

load(strcat('./data&figure/dev_info_', casename,'.mat'));

%% device transfer function model

Tfl = 0.01;
Kvp = 0.5;
Kvi = 5.0;
Cf = 1.59155e-6;
Tf_sys = (Kvp*s+Kvi) / (Cf*s^2+Kvp*s+Kvi);

for ii = 1:dev_info.dev.num
    if ~isempty(find(dev_info.dev.ID(ii) == dev_info.gen.ID))
        % M = abs(0.30 + 0.1 * randn()); d = abs(0.5 + 0.1 * randn());
        % Tfl = abs(0.02 + 0.002 * randn()); tau = abs(0.3 + 0.1 * randn()); Dq = abs(0.8 + 0.2 * randn());
        % kk = find(dev_info.dev.ID(ii) == dev_info.gen.ID);
        % dev_info.gen.M(kk) = M; dev_info.gen.d(kk) = d;
        % dev_info.gen.Dq(kk) = Dq; dev_info.gen.Tfl(kk) = Tfl; dev_info.gen.tau(kk) = tau;

        kk = find(dev_info.dev.ID(ii) == dev_info.gen.ID);
        M = dev_info.gen.M(kk); d = dev_info.gen.d(kk);
        Dq = dev_info.gen.Dq(kk); Tfl = dev_info.gen.Tfl(kk); tau = dev_info.gen.tau(kk);

        Gdev(ii, ii) = 1/( s*(M*s+d) );
        Gdev(ii+dev_info.dev.num, ii+dev_info.dev.num) = Tf_sys * Dq/(tau*s+1);

    elseif ~isempty(find(dev_info.dev.ID(ii) == dev_info.load.ID))
        % M = abs(0.30 + 0.1 * randn()); d = abs(0.5 + 0.1 * randn());
        % Tfl = abs(0.02 + 0.002 * randn()); tau = abs(0.3 + 0.1 * randn()); Dq = abs(0.8 + 0.2 * randn());
        % kk = find(dev_info.dev.ID(ii) == dev_info.load.ID);
        % dev_info.load.M(kk) = M; dev_info.load.d(kk) = d;
        % dev_info.load.Dq(kk) = Dq; dev_info.load.Tfl(kk) = Tfl; dev_info.load.tau(kk) = tau;

        kk = find(dev_info.dev.ID(ii) == dev_info.load.ID);
        M = dev_info.load.M(kk); d = dev_info.load.d(kk);
        Dq = dev_info.load.Dq(kk); Tfl = dev_info.load.Tfl(kk); tau = dev_info.load.tau(kk);

        Gdev(ii, ii) = 1/( s*(M*s+d) );
        Gdev(ii+dev_info.dev.num, ii+dev_info.dev.num) = Tf_sys * Dq/(tau*s+1);

    end
end
% save(strcat('./data&figure/dev_info_', casename,'_alter.mat'), 'dev_info');

%% state-space closed-loop system model
matA = zeros(4*dev_info.dev.num, 4*dev_info.dev.num);
matB2 = zeros(4*dev_info.dev.num, 2*dev_info.dev.num);

matK2 = zeros(2*dev_info.dev.num, 4*dev_info.dev.num);
matK2n = zeros(2*dev_info.dev.num, 4*dev_info.dev.num);

for ii = 1:dev_info.dev.num
    if ~isempty(find(dev_info.dev.ID(ii)==dev_info.gen.ID))
        kk = find(dev_info.dev.ID(ii)==dev_info.gen.ID);
        matA(0*dev_info.dev.num+ii, 1*dev_info.dev.num+ii) = 1;
        matA(1*dev_info.dev.num+ii, 1*dev_info.dev.num+ii) = -dev_info.gen.d(kk) / dev_info.gen.M(kk);
        matA(2*dev_info.dev.num+ii, 3*dev_info.dev.num+ii) = 1;
        matA(3*dev_info.dev.num+ii, 2*dev_info.dev.num+ii) = -1 / (dev_info.gen.tau(kk) *  dev_info.gen.Tfl(kk));
        matA(3*dev_info.dev.num+ii, 3*dev_info.dev.num+ii) = -(1/dev_info.gen.tau(kk) +  1/dev_info.gen.Tfl(kk));

        matB2(1*dev_info.dev.num+ii, 0*dev_info.dev.num+ii) = -1 / dev_info.gen.M(kk);
        matB2(3*dev_info.dev.num+ii, 1*dev_info.dev.num+ii) = -dev_info.gen.Dq(kk) / (dev_info.gen.tau(kk) *  dev_info.gen.Tfl(kk));
    elseif ~isempty(find(dev_info.dev.ID(ii)==dev_info.load.ID))
        kk = find(dev_info.dev.ID(ii)==dev_info.load.ID);
        matA(0*dev_info.dev.num+ii, 1*dev_info.dev.num+ii) = 1;
        matA(1*dev_info.dev.num+ii, 1*dev_info.dev.num+ii) = -dev_info.load.d(kk) / dev_info.load.M(kk);
        matA(2*dev_info.dev.num+ii, 3*dev_info.dev.num+ii) = 1;
        matA(3*dev_info.dev.num+ii, 2*dev_info.dev.num+ii) = -1 / (dev_info.load.tau(kk) *  dev_info.load.Tfl(kk));
        matA(3*dev_info.dev.num+ii, 3*dev_info.dev.num+ii) = -(1/dev_info.load.tau(kk) +  1/dev_info.load.Tfl(kk));

        matB2(1*dev_info.dev.num+ii, 0*dev_info.dev.num+ii) = -1 / dev_info.load.M(kk);
        matB2(3*dev_info.dev.num+ii, 1*dev_info.dev.num+ii) = -dev_info.load.Dq(kk) / (dev_info.load.tau(kk) *  dev_info.load.Tfl(kk));
    end
end

matK2 = [
    Hnet(1:dev_info.dev.num, 1:dev_info.dev.num), zeros(dev_info.dev.num, dev_info.dev.num), Hnet(1:dev_info.dev.num, dev_info.dev.num+1:end), zeros(dev_info.dev.num, dev_info.dev.num)
    Hnet(dev_info.dev.num+1:end, 1:dev_info.dev.num), zeros(dev_info.dev.num, dev_info.dev.num), Hnet(dev_info.dev.num+1:end, dev_info.dev.num+1:end), zeros(dev_info.dev.num, dev_info.dev.num)
    ];

matK2n = [
    Hnet(1:dev_info.dev.num, 1:dev_info.dev.num), zeros(dev_info.dev.num, dev_info.dev.num), zeros(dev_info.dev.num, dev_info.dev.num), zeros(dev_info.dev.num, dev_info.dev.num)
    zeros(dev_info.dev.num, dev_info.dev.num), zeros(dev_info.dev.num, dev_info.dev.num), Hnet(dev_info.dev.num+1:end, dev_info.dev.num+1:end), zeros(dev_info.dev.num, dev_info.dev.num)
    ];

matAc = matA + matB2 * matK2;
matAcn = matA + matB2 * matK2n;

eig_matAc = eig(matAc);
eig_matAcn = eig(matAcn);

matAcn_a = matAcn(1:size(matAcn,1)/2, 1:size(matAcn,1)/2);
matAcn_v = matAcn(size(matAcn,1)/2+1:end, size(matAcn,1)/2+1:end);

if max(real(eig_matAc)) >= 1e-8
    error('error: close-loop system is UNSTABLE');
end

%% compute matS^(-1) denoted as matR
matR_arr = zeros(length(w), size(Gdev, 1), size(Gdev, 2));
matRn_arr = zeros(length(w), size(Gdev, 1), size(Gdev, 2));

eig_R = zeros(length(w), size(Gdev,1));
eig_Rn = zeros(length(w), size(Gdev,1));

log_det_R = zeros(length(w), 1);
log_det_Rn = zeros(length(w), 1);

for ii = 1:length(w)
    matR_arr(ii, :, :) = eye(size(Gdev,1)) + evalfr(Gdev,1j*w(ii)) * Hnet;
    matRn_arr(ii, :, :) = eye(size(Gdev,1)) + evalfr(Gdev,1j*w(ii)) * Hn_net;
    fprintf('matS^(-1) computation: log10(w) = %.4f\n', log10(w(ii)));
end

%% compute matS
for ii = 1:length(w)
    D_R = func_eig_modi('matR', squeeze(matR_arr(ii,:,:)));
    D_Rn = func_eig_modi('matR', squeeze(matRn_arr(ii,:,:)));

    % eigenvalue
    eig_R(ii,:) = D_R';
    eig_Rn(ii,:) = D_Rn';

    % determinant
    log_det_R(ii) = sum(log10(abs(eig_R(ii,:))));
    log_det_Rn(ii) = sum(log10(abs(eig_Rn(ii,:))));

    fprintf('eigR computation: log10(w) = %.4f\n', log10(w(ii)));
end

min_log_eig_R = zeros(length(w), 1);
min_log_eig_Rn = zeros(length(w), 1);
for ii = 1:length(w)
    min_log_eig_R(ii) = log10(min(abs(eig_R(ii,:))));
    min_log_eig_Rn(ii) = log10(min(abs(eig_Rn(ii,:))));
end

%% eigenvalue classification and tracking
eig_R = func_track_eigenvalues(eig_R);
eig_Rn = func_track_eigenvalues(eig_Rn);

w_row = 1;
[eigA, eigAn, eigV, eigVn] = func_classify_singular(Gdev, Hnet, Hn_net, w(w_row), 0.01);
eigA = 1./eigA; eigAn = 1./eigAn; eigV = 1./eigV; eigVn = 1./eigVn;
[idx_A, idx_V, idx_An, idx_Vn] = func_classify_track(eig_R, eig_Rn, eigA, eigAn, eigV, eigVn, w_row);

%% directed area analysis

area_S_ang = 0;
area_S_vol = 0;
for ii = 1:size(eig_R, 2)
    if ~isempty(find(ii == idx_A))
        area_S_ang = area_S_ang - func_area_bode(w, log10(abs(eig_R(:,ii))));
    elseif ~isempty(find(ii == idx_V))
        area_S_vol = area_S_vol - func_area_bode(w, log10(abs(eig_R(:,ii))));
    end
end

fprintf('sensitivity integral of ANG-relative part: int log(S) = %.10f\n', area_S_ang);
fprintf('sensitivity integral of VOL-relative part: int log(S) = %.10f\n', area_S_vol);

%%

figure(1);
subplot(1,3,1); hold on; box on;
plot(eig_matAc+1e-15j, '*', 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]);
plot(eig_matAcn+1e-15j, 's', 'LineWidth', 1.5, 'Color', [0.85, 0.3250, 0.0980]);

subplot(1,3,2); hold on; box on;
plot(log10(w(1:2:end)), -log10(abs(eig_R(1:2:end,idx_A))), 'Color', [0.85, 0.3250, 0.0980]);
plot(log10(w(1:2:end)), -log10(abs(eig_R(1:2:end,idx_V))), 'Color', [0, 0.447, 0.741]);
plot(log10(w(1:2:end)), -log10(abs(eig_Rn(1:2:end,idx_An))), '--', 'LineWidth', 1.6, 'Color', [0.85, 0.3250, 0.0980]);
plot(log10(w(1:2:end)), -log10(abs(eig_Rn(1:2:end,idx_Vn))), '--', 'LineWidth', 1.6, 'Color', [0, 0.447, 0.741]);
xlim([-0.5, 2])

subplot(1,3,3); hold on; box on;
plot(log10(w), -log10(abs(eig_R(:,idx_A))), 'Color', [0.85, 0.3250, 0.0980]);
plot(log10(w), -log10(abs(eig_R(:,idx_V))), 'Color', [0, 0.447, 0.741]);
plot(log10(w), -log10(abs(eig_Rn(:,idx_An))), '--', 'LineWidth', 2, 'Color', [0.85, 0.3250, 0.0980]);
plot(log10(w), -log10(abs(eig_Rn(:,idx_Vn))), '--', 'LineWidth', 2, 'Color', [0, 0.447, 0.741]);
xlim([-0.5, 2])

%% generalized Nyquist plot

figure(2); hold on; box on;
plot(real(eig_R(:,idx_A))-1, imag(eig_R(:,idx_A)), '-', 'Color', [0.85, 0.3250, 0.0980], 'LineWidth', 1.2);
plot(real(eig_R(:,idx_V))-1, imag(eig_R(:,idx_V)), '-', 'Color', [0, 0.447, 0.741], 'LineWidth', 1.2);

plot(real(eig_R(:,idx_A))-1, -imag(eig_R(:,idx_A)), '-', 'Color', [0.85, 0.3250, 0.0980], 'LineWidth', 1.2);
plot(real(eig_R(:,idx_V))-1, -imag(eig_R(:,idx_V)), '-', 'Color', [0, 0.447, 0.741], 'LineWidth', 1.2);

ylim([-30,30]); xlim([-60, 60])

%%
toc;