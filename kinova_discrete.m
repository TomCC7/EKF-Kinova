%% Setup
% load model
rob = importrobot("gen3.urdf")
rob_param = rob_model();
rob.DataFormat = 'col';
rob.Gravity = [0;0;-9.81];
n=14;      %number of state
q=0.1;    %std of process
r=1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2*eye(n);        % covariance of measurement

h=@(x)[x];  % measurement equation

%% Simulation
step = 0.0001;
time = 0.5;
N = ceil(time / step);
% init state
s = zeros(n, 1);
x = s + q * randn(n, 1);
torque = @(t) ones(7, 1) * sin(t);
% init state cov
P = zeros(n); % cov
N = length(t);
% state estimation
xV = zeros(n, N);
% state actual
xV = zeros(n, N);
% measurement
zV = zeros(n, N);
t = linspace(0,time, N);
sV = zeros(n, N);
for k=2:N
    f = @(x)[next_state(x, rob, torque, k, t)];
    sV(:, k) = f(sV(:, k-1)) + q * randn(n, 1);
    % disp(sV(:,k))
    z = h(sV(:, k)) + r * randn(n, 1); % measure
    zV(:, k) = z;
    [x, P] = ekf(f, x, P, h, z, Q, R);
    xV(:, k) = x;
end

figure(1);
hold on;
plot(t, mean(xV - sV));
plot(t, mean(zV - sV));
legend(["error after filtering", "error before filtering"]);
title("error comparison")
hold off;
for i=2:15
figure(i);
hold on;
plot(t, zV(i-1,:));
plot(t, sV(i-1,:));
plot(t, xV(i-1,:));
legend(["raw measurement", "ground truth", "estimation"]);
title("state "+ (i-1));
hold off;
end
