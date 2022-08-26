function dydt = deriv_state(x, model, fs, torque)
num_joint = length(x) / 2;
q = x(1:num_joint);
qd = x(num_joint+1:end);
% disp(qd)
% disp(torque)
M = massMatrix(model, q); % (7, 7)
C = velocityProduct(model, q, qd); % (1, 7)
G = gravityTorque(model, q); % (1, 7)
% NOTE: warning for not differentiable dynamics because of dir
dir = (sigmf(qd,[10 0]) - 0.5) * 2; % 1 by 7
dydt = [
    qd;
    (M + fs(:, 3) .* eye(num_joint)) \ ...
          (- C - G ...
           + torque - fs(:, 1) .* dir - fs(:, 2) .* qd)
       ];
end
