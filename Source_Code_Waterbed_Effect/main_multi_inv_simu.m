clear all; clc; close all;
tic;
%% system configurations
vpa('1e400');     

casename = 'case14';

% network
load(strcat('./data&figure/Hnet_info_', casename,'.mat'));
Hnet = Hnet_info.Hnet;
Hn_net = blkdiag(Hnet(1:size(Hnet,1)/2, 1:size(Hnet,1)/2), Hnet(size(Hnet,1)/2+1:end, size(Hnet,1)/2+1:end));

% device
dev_info = Hnet_info.dev_info;

s = tf('s');
Gdev = s*zeros(2*dev_info.dev.num, 2*dev_info.dev.num);

load(strcat('./data&figure/dev_info_', casename,'.mat'));

for ii = 1:dev_info.dev.num
    if ~isempty(find(dev_info.dev.ID(ii) == dev_info.gen.ID))
        kk = find(dev_info.dev.ID(ii) == dev_info.gen.ID);
        M = dev_info.gen.M(kk); d = dev_info.gen.d(kk);
        Dq = dev_info.gen.Dq(kk); Tfl = dev_info.gen.Tfl(kk); tau = dev_info.gen.tau(kk);

        Gdev(ii, ii) = 1/( s*(M*s+d) );
        Gdev(ii+dev_info.dev.num, ii+dev_info.dev.num) = Dq /( (Tfl*s+1)*(tau*s+1) );

    elseif ~isempty(find(dev_info.dev.ID(ii) == dev_info.load.ID))
        kk = find(dev_info.dev.ID(ii) == dev_info.load.ID);
        M = dev_info.load.M(kk); d = dev_info.load.d(kk);
        Dq = dev_info.load.Dq(kk); Tfl = dev_info.load.Tfl(kk); tau = dev_info.load.tau(kk);

        Gdev(ii, ii) = 1/( s*(M*s+d) );
        Gdev(ii+dev_info.dev.num, ii+dev_info.dev.num) = Dq /( (Tfl*s+1)*(tau*s+1) );
 
    end
end

%% state-space closed-loop transfer function

% state-space function:
% dx/dt = Ax + B1*r + B2*n; y = Cx+w; n = Hnet*y
% equivalent form:
% dx/dt = (A+B2KC)*x + B1*r + B2*Kw
% y = C*x + w

matA = zeros(4*dev_info.dev.num, 4*dev_info.dev.num);
matB2 = zeros(4*dev_info.dev.num, 2*dev_info.dev.num);

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

matB1 = -matB2;

matC = zeros(size(Hnet,1), size(Hnet,1)*2);
matC(1:size(Hnet,1)/2, 1:size(Hnet,1)/2) = eye(size(Hnet,1)/2);
matC(size(Hnet,1)/2+1:end, size(Hnet,1)+1:size(Hnet,1)/2*3) = eye(size(Hnet,1)/2);

matAc = matA + matB2 * Hnet * matC;
matAcn = matA + matB2 * Hn_net * matC;

%% simluation: from Ps,Qs to V,theta
disturb_dir_seq = [10, 45, 80];
disturb_amp = 0.05;

for rr = 1:length(disturb_dir_seq)

dt_sim = 0.01;
t_sim = 0:dt_sim:10;

disturb_dir = disturb_dir_seq(rr);
disturb = disturb_amp * [cos(disturb_dir/180*pi)*ones(size(Hnet,1)/2,1); sin(disturb_dir/180*pi)*ones(size(Hnet,1)/2,1)];
disturb(1) = 0;

disturb_seq = ones(length(t_sim), 1) * disturb';

% time-domain simulation WITHOUT considering the coupling
sys_cl_n = ss(matAcn, matB2*Hn_net, matC, eye(size(Hnet))); % step: theta,V -> theta, V
[state_n, t_sim] = lsim(sys_cl_n, disturb_seq, t_sim);

ang_n = 180/pi*(state_n(:, 1:size(Hnet,1)/2) - state_n(:, 1));
vol_n = state_n(:, size(Hnet,1)/2+1:end);

figure(1+rr); 
subplot(2,2,1); hold on; box on;
plot(t_sim(1:10:end), ang_n(1:10:end,:), '--', 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);
title(['Angle  transient, phi = ', num2str(disturb_dir_seq(rr)), 'deg'])
subplot(2,2,2); hold on; box on;
plot(t_sim(1:10:end), vol_n(1:10:end,:), '--', 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);
title(['Voltage transient, phi = ', num2str(disturb_dir_seq(rr)), 'deg'])

% time-domain simulation WITH considering the coupling
sys_cl = ss(matAc, matB2*Hnet, matC, eye(size(Hnet))); % step: theta,V -> theta, V
[state, t_sim] = lsim(sys_cl, disturb_seq, t_sim);

ang = 180/pi*(state(:, 1:size(Hnet,1)/2) - state(:, 1));
vol = state(:, size(Hnet,1)/2+1:end);

figure(1+rr); 
subplot(2,2,3);
plot(t_sim(1:1:end), ang(1:1:end,:), 'LineWidth', 1.2, 'Color', [0.85, 0.3250, 0.0980]);
subplot(2,2,4);
plot(t_sim(1:1:end), vol(1:1:end,:), 'LineWidth', 1.2, 'Color', [0, 0.447, 0.741]);

subplot(2,2,1); ylim([min(min(min(ang_n)),min(min(ang))), max(max(max(ang_n)),max(max(ang)))]); xlim([0,5])
subplot(2,2,3); ylim([min(min(min(ang_n)),min(min(ang))), max(max(max(ang_n)),max(max(ang)))]); xlim([0,5])
subplot(2,2,2); ylim([0, 1.1*max(max(max(vol_n)),max(max(vol)))]); xlim([0,5])
subplot(2,2,4); ylim([0, 1.1*max(max(max(vol_n)),max(max(vol)))]); xlim([0,5])
set(gcf, 'Position', [100, 100, 500, 600]);
end