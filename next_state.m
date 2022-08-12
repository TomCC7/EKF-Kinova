function y = next_state(x, model, torque, k, t)
vel = x(length(x)/2+1:end);
v = [
    vel;
    forwardDynamics(model, x(1:length(x)/2), vel, torque(t(k)))
    ];
y = v*(t(k) - t(k-1)) + x;
end
