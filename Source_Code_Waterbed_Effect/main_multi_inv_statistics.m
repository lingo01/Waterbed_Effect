clear all; clc; close all; vpa('1e400');
tic;
%% system configuration
rng(2026);

bus_num = 14;

casename = 'case14';

w = 10.^[-9:0.002:10];

sample_num = 1000;

load_scaling_factor_seq = 1 + 0.2^2 * randn(bus_num, 2, sample_num);

area_S_ang_seq = zeros(sample_num, 1);
area_S_vol_seq = zeros(sample_num, 1);

%% device initialization
load(strcat('./data&figure/dev_info_', casename,'.mat'));

s = tf('s');
Gdev = s*zeros(2*dev_info.dev.num, 2*dev_info.dev.num);

Tfl = 0.01;
Kvp = 0.5;
Kvi = 5.0;
Cf = 1.59155e-6;
Tf_sys = (Kvp*s+Kvi) / (Cf*s^2+Kvp*s+Kvi);

for ii = 1:dev_info.dev.num
    if ~isempty(find(dev_info.dev.ID(ii) == dev_info.gen.ID))
        kk = find(dev_info.dev.ID(ii) == dev_info.gen.ID);
        M = dev_info.gen.M(kk); d = dev_info.gen.d(kk);
        Dq = dev_info.gen.Dq(kk); Tfl = dev_info.gen.Tfl(kk); tau = dev_info.gen.tau(kk);

        Gdev(ii, ii) = 1/( s*(M*s+d) );
        Gdev(ii+dev_info.dev.num, ii+dev_info.dev.num) = Tf_sys * Dq/(tau*s+1);

    elseif ~isempty(find(dev_info.dev.ID(ii) == dev_info.load.ID))
        kk = find(dev_info.dev.ID(ii) == dev_info.load.ID);
        M = dev_info.load.M(kk); d = dev_info.load.d(kk);
        Dq = dev_info.load.Dq(kk); Tfl = dev_info.load.Tfl(kk); tau = dev_info.load.tau(kk);

        Gdev(ii, ii) = 1/( s*(M*s+d) );
        Gdev(ii+dev_info.dev.num, ii+dev_info.dev.num) = Tf_sys * Dq/(tau*s+1);

    end
end

for rr = 1:sample_num
    load_scaling_factor = load_scaling_factor_seq(:, :, rr);
    %% network initialization
    [~, Hnet_info] = func_Hnet_generator(casename, 0, load_scaling_factor);
    Hnet = Hnet_info.Hnet;
    Hn_net = blkdiag(Hnet(1:size(Hnet,1)/2, 1:size(Hnet,1)/2), Hnet(size(Hnet,1)/2+1:end, size(Hnet,1)/2+1:end));
    
    %% stability check (time-domain eigenvalue analysis)
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
    
    matK2 = [
        Hnet(1:dev_info.dev.num, 1:dev_info.dev.num), zeros(dev_info.dev.num, dev_info.dev.num), Hnet(1:dev_info.dev.num, dev_info.dev.num+1:end), zeros(dev_info.dev.num, dev_info.dev.num)
        Hnet(dev_info.dev.num+1:end, 1:dev_info.dev.num), zeros(dev_info.dev.num, dev_info.dev.num), Hnet(dev_info.dev.num+1:end, dev_info.dev.num+1:end), zeros(dev_info.dev.num, dev_info.dev.num)
        ];
    
    matAc = matA + matB2 * matK2;
    
    eig_matAc = eig(matAc);
    
    if max(real(eig_matAc)) >= 1e-8
        % error('error: close-loop system is UNSTABLE');
        break;
    end
    
    %%
    matR_arr = zeros(length(w), size(Gdev, 1), size(Gdev, 2));
    matRn_arr = zeros(length(w), size(Gdev, 1), size(Gdev, 2));
    
    eig_R = zeros(length(w), size(Gdev,1));
    eig_Rn = zeros(length(w), size(Gdev,1));
    
    log_det_R = zeros(length(w), 1);
    log_det_Rn = zeros(length(w), 1);
    
    for ii = 1:length(w)
        matR_arr(ii, :, :) = eye(size(Gdev,1)) + evalfr(Gdev,1j*w(ii)) * Hnet;
        matRn_arr(ii, :, :) = eye(size(Gdev,1)) + evalfr(Gdev,1j*w(ii)) * Hn_net;
        % fprintf('matR computation: log10(w) = %.4f\n', log10(w(ii)));
    end
    
    %%
    for ii = 1:length(w)
        D_R = func_eig_modi('matR', squeeze(matR_arr(ii,:,:)));
        D_Rn = func_eig_modi('matR', squeeze(matRn_arr(ii,:,:)));
    
        % eigenvalue
        eig_R(ii,:) = D_R';
        eig_Rn(ii,:) = D_Rn';
    
        % determinant
        log_det_R(ii) = sum(log10(abs(eig_R(ii,:))));
        log_det_Rn(ii) = sum(log10(abs(eig_Rn(ii,:))));
    
        % fprintf('eigR computation: log10(w) = %.4f\n', log10(w(ii)));
    end
    
    %%
    eig_R = func_track_eigenvalues(eig_R);
    eig_Rn = func_track_eigenvalues(eig_Rn);
    
    %%
    w_row = 1;
    [eigA, eigAn, eigV, eigVn] = func_classify_singular(Gdev, Hnet, Hn_net, w(w_row), 0.01);
    eigA = 1./eigA; eigAn = 1./eigAn; eigV = 1./eigV; eigVn = 1./eigVn;
    [idx_A, idx_V, idx_An, idx_Vn] = func_classify_track(eig_R, eig_Rn, eigA, eigAn, eigV, eigVn, w_row);
    
    %%
    area_S_ang = 0;
    area_S_vol = 0;
    for ii = 1:size(eig_R, 2)
        if ~isempty(find(ii == idx_A))
            area_S_ang = area_S_ang - func_area_bode(w, log10(abs(eig_R(:,ii))));
        elseif ~isempty(find(ii == idx_V))
            area_S_vol = area_S_vol - func_area_bode(w, log10(abs(eig_R(:,ii))));
        end
    end
    
    fprintf('RD %d\t sensitivity integral of ANG/VOL-relative part: int log(S) = %.10f, %.10f\n', rr, area_S_ang, area_S_vol);

    area_S_ang_seq(rr) = area_S_ang;
    area_S_vol_seq(rr) = area_S_vol;
end

%%
figure(1); 
subplot(3,1,1);
histogram(area_S_ang_seq);
legend('$\mathcal{A}_{\mathcal{I}_A}$', '$\mathcal{A}_{\mathcal{I}_V}$', 'Interpreter', 'latex')
subplot(3,1,2);
histogram(area_S_vol_seq);
legend('$\mathcal{A}_{\mathcal{I}_V}$', '$\mathcal{A}_{\mathcal{I}_V}$', 'Interpreter', 'latex')
subplot(3,1,3);
histogram(1+area_S_vol_seq./area_S_ang_seq)

%%
toc;