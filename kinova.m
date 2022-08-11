%% Setup
% load model
rob = importrobot("gen3.urdf")
rob.DataFormat = 'col';
rob.Gravity = [0;0;-9.81];
n=14;      %number of state
q=0.1;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2*eye(n);        % covariance of measurement

h=@(x)[x];  % measurement equation

%% Simulation
time = 0.5;
% init state
s = zeros(n, 1);
x = s + q * randn(n, 1);
torque = zeros(7, 1);
% init state cov
P = eye(n); % cov
[t,sV] = ode45(@(t,y) deriv_state(y, rob, torque), [0,time], s);
f = @(x)[next_state(x, rob, torque, k, t)];
N = length(t);
% state estimation
xV = zeros(n, N);
% state actual
xV = zeros(n, N);
% measurement
zV = zeros(n, N);
sV = sV';
for k=2:N
    sV(:, k) = f(sV(:, k-1)) + q * randn(n, 1);
    disp(sV(:,k))
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
