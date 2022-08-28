%% Setup
% load model
n=14; %number of state
params = rob_model();
rob = modify_robot(importrobot("gen3.urdf"), params, n/2)
fs = params(:, end-2:end);
q=0.005;    %std of process
r=0.01;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2*eye(n);        % covariance of measurement

h=@(x)[x];  % measurement equation

%% Simulation
time = 4;
step = 0.01;
N = ceil(time / step);
t = linspace(0, time, N);
% init state
s = zeros(n, 1);
x = s + q * randn(n, 1);
amp = 2 * pi / 8;
torque = @(t) ones(7, 1) * sin(t) * amp;
% init state cov
P = zeros(n); % cov
[t_ode,sV] = ode45(@(t,y) deriv_state(y, rob, fs, torque(t)), [0,time], s);
N = length(t);
% state actual
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
figure(1);
hold on;
plot(t, mean(xV - sV_interp));
plot(t, mean(zV - sV_interp));
legend(["error after filtering", "error before filtering"]);
title("error comparison")
hold off;
for i=2:15
figure(i);
hold on;
plot(t, zV(i-1,:));
plot(t, sV_interp(i-1,:));
plot(t, xV(i-1,:));
legend(["raw measurement", "ground truth", "estimation"]);
title("state "+ (i-1));
hold off;
end
