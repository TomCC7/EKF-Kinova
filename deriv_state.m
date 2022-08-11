function dydt = deriv_state(x, model, torque)
vel = x(length(x)/2+1:end);
dydt = [
    vel;
    forwardDynamics(model, x(1:length(x)/2), vel, torque)
    ];
end
