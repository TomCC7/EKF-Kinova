function y = next_state(x, model, fs, torque, k, t)
v = deriv_state(x, model, fs, torque(t(k)));
y = v*(t(k) - t(k-1)) + x;
end
