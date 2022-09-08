%% Generate estimation based on experiment data
%% Setup
% load data
dir_name = '../data/covariance/coeff_4_cutted';
save_prefix = 'coeff_4_cutted';
data_path = dir(dir_name + "/*.mat");
n_data = length(data_path);
data = cell(n_data, 1);
state_exp = cell(n_data, 1);
frame_ts = cell(n_data);
for i = 1:n_data
    data{i} = load(fullfile(data_path(i).folder, data_path(i).name));
    frame_ts{i} = double(data{i}.frame_ts - data{i}.frame_ts(1));
    data{i}.torque = double(data{i}.torque);
    state_exp{i} = double([data{i}.feedback_pos, data{i}.feedback_vel]);
end
t_exp = frame_ts{1};
torque_exp = data{1}.torque;
% load model
n=14; %number of state
params = rob_model();
rob = modify_robot(importrobot("gen3.urdf"), params, n/2);
fs = params(:, end-2:end);
q=1e-6;       %std of process
Q=q^2*eye(n); % covariance of process
r_diag = [2.26300904e-04 1.37933811e-03 9.08683309e-05 4.36212329e-06 ...
          4.90829354e-05 5.36947627e-05 6.19376190e-05 ...
          5.86390899e-05 2.65790298e-03 6.49749867e-05 5.87998681e-05 ...
          2.87742624e-05 2.16910833e-05 2.80170833e-05];
% r_diag = r_diag;
R = diag(r_diag.^2);
h=@(x)(x);  % measurement equation

%% Simulation
start_time = 0;
time = 1;
step = 0.001;
N = ceil(time / step);
t = linspace(start_time, start_time + time, N);
torque = @(t) get_torque(t, t_exp, torque_exp);
amp = 2 * pi / 8;
% init state cov
P = cell(n_data); % cov
for i=1:n_data
    P{i} = zeros(n);
end
s = interp1(t_exp, state_exp{1}, start_time)';
[t_ode,sV] = ode45(@(t,y) deriv_state(y, rob, fs, torque(t)), [0,time], s);
N = length(t);
% state estimation
xV = cell(n_data);
for i=1:n_data
    xV{i} = zeros(n, N);
end
% init state
for i=1:n_data
    xV{i}(:, 1) = interp1(frame_ts{i}, state_exp{i}, start_time)' ...
                  + q * randn(n, 1);
end
%% EKF
zV = zeros(n, N);
sV = sV';
state_interp = zeros(n_data, n, N);
sim_interp = zeros(n, N);
for k=2:N
    f = @(x)[next_state(x, rob, fs, torque, k, t)];
    sim_interp(:, k) = interp1(t_ode, sV', t(k) - start_time);
    for i = 1:n_data
        state_interp(i, :, k) = interp1(frame_ts{i}, state_exp{i}, t(k));
        z = state_interp(1, :, k)';
        [xV{i}(:, k), P{i}] = ekf(f, xV{i}(:, k-1), P{i}, h, z, Q, R);
    end
    disp(sprintf('%dth iteration', k))
end

%% write back
for i = 1:n_data
    for j = 1:length(frame_ts{i})
        if frame_ts{i}(j) >= time - 0.01
            break
        end
        state = interp1(t, xV{i}', frame_ts{i}(j))
        data{i}.feedback_pos(j, :) = state(1:7);
        data{i}.feedback_vel(j, :) = state(8:14);
    end
    mkdir(save_prefix);
    d = data{i};
    save(sprintf('%s/%s_%d.mat', save_prefix, save_prefix, i), 'd');
end

%% Helper functions
function t = get_torque(time, t_exp, torque_exp)
    if time < t_exp(end)
        t = interp1(t_exp, torque_exp, time)';
    else
        t = zeros(7,1);
    end
end
