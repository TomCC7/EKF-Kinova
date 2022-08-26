%% Based on experiment data
%% Setup
% load data
data = load('../data/covariance/12345/SinExperiment_1661463385.npz.mat');
t_exp = data.frame_ts - data.frame_ts(1);
torque_exp = data.torque;
state_exp = [data.feedback_pos, data.feedback_vel];
% load model
n=14; %number of state
params = rob_model();
rob = modify_robot(importrobot("gen3.urdf"), params, n/2);
fs = params(:, end-2:end);
q=5e-5;    %std of process
r=5e-5;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2*eye(n);        % covariance of measurement

h=@(x)(x);  % measurement equation

%% Simulation
time = 1;
step = 0.001;
N = ceil(time / step);
t = linspace(0, time, N);
torque = @(t) get_torque(t, t_exp, torque_exp);
% init state
s = zeros(n, 1);
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
sV_interp = zeros(n, N);
for k=2:N
    f = @(x)[next_state(x, rob, fs, torque, k, t)];
    sV_interp(:, k) = interp1(t_ode, sV', t(k));
    z = h(sV_interp(:, k)) + r * randn(n, 1); % measure
    zV(:, k) = z;
    [x, P] = ekf(f, x, P, h, z, Q, R);
    xV(:, k) = x;
end

%% plot
for i=2:15
figure(i);
hold on;
plot(t, zV(i-1,:));
plot(t, sV_interp(i-1,:));
plot(t, xV(i-1,:));
plot(t, interp1(t_exp, state_exp(:, i-1), t))
legend(["raw measurement", "ground truth", "estimation", "experiment measurement"]);
title("state "+ (i-1));
hold off;
end
figure(1);
hold on;
plot(t, mean(xV - sV_interp));
plot(t, mean(zV - sV_interp));
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
