function [Hnet, Hnet_info] = func_Hnet_generator(casename, saveflag, load_scaling_factor)

if nargin == 0
    casename = "case118";
    saveflag = 0;
    load_scaling_factor = 1;
elseif nargin == 1
    saveflag = 0;
    load_scaling_factor = 1;
elseif nargin == 2
    load_scaling_factor = 1;
end

mpc = loadcase(casename);

% Re-index buses in the case
bus_old = mpc.bus(:,1);
bus_new = (1:size(mpc.bus, 1))';

mpc.bus(:,1) = bus_new;
for ii = 1:size(mpc.branch, 1)
    mpc.branch(ii, 1) = bus_new( find(bus_old == mpc.branch(ii, 1)) );
    mpc.branch(ii, 2) = bus_new( find(bus_old == mpc.branch(ii, 2)) );
end
for ii = 1:size(mpc.gen, 1)
    mpc.gen(ii, 1) = bus_new( find(bus_old == mpc.gen(ii, 1)) );
end

% load scaling
mpc.bus(:,[3,4]) = load_scaling_factor .* mpc.bus(:,[3,4]);

% Solve AC power flow
mpc = runpf(mpc);
bus_Vs = mpc.bus(:, 8);
bus_As = mpc.bus(:, 9)*pi/180;

% Network admittance matrix (before reduction)
netY = full(makeYbus(mpc));
nNet = size(netY, 1);

% Net power estimate from branch flows
P_pure = zeros(nNet, 1);
Q_pure = zeros(nNet, 1);

for ii = 1:size(mpc.branch, 1)
    P_pure(mpc.branch(ii,1), 1) = P_pure(mpc.branch(ii,1), 1) + mpc.branch(ii, 14)/mpc.baseMVA;
    Q_pure(mpc.branch(ii,1), 1) = Q_pure(mpc.branch(ii,1), 1) + mpc.branch(ii, 15)/mpc.baseMVA;
    P_pure(mpc.branch(ii,2), 1) = P_pure(mpc.branch(ii,2), 1) + mpc.branch(ii, 16)/mpc.baseMVA;
    Q_pure(mpc.branch(ii,2), 1) = Q_pure(mpc.branch(ii,2), 1) + mpc.branch(ii, 17)/mpc.baseMVA;
end

%% Bus/device classification
dev.ID = [];
load.ID = [];
cnct.ID = [];

gen.ID = mpc.gen(:,1);
gen.ID = unique(gen.ID);

for ii = 1:nNet
    if (abs(P_pure(ii))>=1e-7) || (abs(Q_pure(ii))>=1e-7)
        dev.ID = [dev.ID, ii];
        if isempty(find(ii==gen.ID))
            load.ID = [load.ID, ii];
        end
    else
        cnct.ID = [cnct.ID, ii];
    end
end
dev.num = length(dev.ID);

%% Network reduction (Kron reduction)
netY_re = netY(dev.ID, dev.ID) - netY(dev.ID, cnct.ID) * inv(netY(cnct.ID, cnct.ID)) * netY(cnct.ID, dev.ID);
netG_re = real(netY_re);
netB_re = imag(netY_re);

%% Build network-side coupling Jacobian
syms As_sym [dev.num 1] real
syms Vs_sym [dev.num 1] real

% Define the power-flow mapping
power_flow = sym('power_flow', [2*dev.num, 1]);
for ii = 1:dev.num
    tmpP = 0;    tmpQ = 0;
    for jj = 1:dev.num
        tmpP = tmpP + Vs_sym(jj) * ( netG_re(ii,jj)*cos(As_sym(ii)-As_sym(jj)) + netB_re(ii,jj)*sin(As_sym(ii)-As_sym(jj)) );
        tmpQ = tmpQ + Vs_sym(jj) * ( netB_re(ii,jj)*cos(As_sym(ii)-As_sym(jj)) - netG_re(ii,jj)*sin(As_sym(ii)-As_sym(jj)) );
    end
    power_flow(ii) =  tmpP * Vs_sym(ii);
    power_flow(dev.num+ii) = -tmpQ * Vs_sym(ii);
end

% Compute the Jacobian matrix
Hnet_sym = jacobian(power_flow, [As_sym; Vs_sym]);

dev_As = bus_As(dev.ID);
dev_Vs = bus_Vs(dev.ID);

% Substitute the operating point and obtain the numeric matrix
Hnet = subs(Hnet_sym, [As_sym; Vs_sym], [dev_As; dev_Vs]);
Hnet = double(Hnet);

dev_info.dev = dev;
dev_info.gen = gen;
dev_info.load = load;
dev_info.As = dev_As;
dev_info.Vs = dev_Vs;
dev_info.mpc = mpc;

Hnet_info.Hnet = Hnet;
Hnet_info.Hnet_sym = Hnet_sym;

Hnet_info.dev_info = dev_info;

if saveflag == 1
    save(strcat('./data&figure/Hnet_info_', casename,'.mat'), 'Hnet_info');
end

% plot(eig(Hnet)+1e-10j, 'o', 'LineWidth', 2)
% 
% func_Hnet_generator("case5", 1);
% func_Hnet_generator("case9", 1);
% func_Hnet_generator("case14", 1);
% func_Hnet_generator("case30", 1);
% func_Hnet_generator("case39", 1);
% func_Hnet_generator("case57", 1);
% func_Hnet_generator("case118", 1);
% func_Hnet_generator("case145", 1);
% func_Hnet_generator("case300", 1);
end