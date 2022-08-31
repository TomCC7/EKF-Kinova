%% Based on experiment data
%% Setup
% load data
data_path = dir('../data/covariance/coeff_1/*.mat');
n_data = length(data_path);
data = cell(n_data, 1);
state_exp = cell(n_data, 1);
for i = 1:n_data
    data{i} = load(fullfile(data_path(i).folder, data_path(i).name));
    data{i}.frame_ts = double(data{i}.frame_ts - data{i}.frame_ts(1));
    data{i}.torque = double(data{i}.torque);
    state_exp{i} = double([data{i}.feedback_pos, data{i}.feedback_vel]);
end
t_exp = data{1}.frame_ts;
torque_exp = data{1}.torque;
% load model
n=14; %number of state
params = rob_model();
rob = modify_robot(importrobot("gen3.urdf"), params, n/2);
fs = params(:, end-2:end);
q=1e-4;       %std of process
Q=q^2*eye(n); % covariance of process
r_diag = [2.26300904e-04 1.37933811e-03 9.08683309e-05 4.36212329e-06 ...
          4.90829354e-05 5.36947627e-05 6.19376190e-05 ...
          5.86390899e-05 2.65790298e-03 6.49749867e-05 5.87998681e-05 ...
          2.87742624e-05 2.16910833e-05 2.80170833e-05];
% r_diag = r_diag;
R = diag(r_diag.^2);
h=@(x)(x);  % measurement equation

%% Simulation
start_time = 0.5;
time = 1;
step = 0.001;
N = ceil(time / step);
t = linspace(start_time, start_time + time, N);
torque = @(t) get_torque(t, t_exp, torque_exp);
% init state
s = interp1(t_exp, state_exp{1}, start_time)';
x = s + q * randn(n, 1);
amp = 2 * pi / 8;
% init state cov
P = zeros(n); % cov
[t_ode,sV] = ode45(@(t,y) deriv_state(y, rob, fs, torque(t)), [0,time], s);
N = length(t);
% state estimation
xV = zeros(n, N);
% measurement
zV = zeros(n, N);
sV = sV';
state_interp = zeros(n_data, n, N);
sim_interp = zeros(n, N);
for k=2:N
    f = @(x)[next_state(x, rob, fs, torque, k, t)];
    for i = 1:n_data
        state_interp(i, :, k) = interp1(data{i}.frame_ts, state_exp{i}, t(k));
    end
    sim_interp(:, k) = interp1(t_ode, sV', t(k) - start_time);
    % z = h(sim_interp(:, k)) + r_diag * randn(n, 1); % measure
    z = state_interp(randi(n_data),:, k)';
    zV(:, k) = z;
    [x, P] = ekf(f, x, P, h, z, Q, R);
    xV(:, k) = x;
end

%% plot
figure(1);
a = subplot(7,2, 1);
for i=2:8
    a = subplot(7,2, 2*(i-1)-1);
    hold on;
    plot(t, zV(i-1,:), 'b');
    plot(t, xV(i-1,:), 'r', 'LineWidth',5);
    plot(t, sim_interp(i-1,:), 'k');
    for j=1:n_data
        plot(t, reshape(state_interp(j, i-1, :),1,[]), '--');
    end
    title("position "+ (i-1));
    hold off;
end

for i=9:15
    a = subplot(7,2, 2*(i-8));
    hold on;
    plot(t, zV(i-1,:), 'b');
    plot(t, xV(i-1,:), 'r', 'LineWidth',2);
    plot(t, sim_interp(i-1,:), 'k');
    for j=1:n_data
        plot(t, reshape(state_interp(j, i-1, :),1,[]), '--');
    end
    title("velocity "+ (i-8));
    hold off;
end
legend(a, ["raw measurement", "estimation", "simulation"]);
figure(2)
hold on;
plot(t, mean(xV - sim_interp));
plot(t, mean(zV - sim_interp));
legend(["error after filtering", "error before filtering"]);
title("error comparison")
hold off;

%% Helper functions
function t = get_torque(time, t_exp, torque_exp)
    if time < t_exp(end)
        t = interp1(t_exp, torque_exp, time)';
    else
        t = zeros(7,1);
    end
end
