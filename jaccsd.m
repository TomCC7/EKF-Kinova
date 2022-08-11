function [z,A]=jaccsd(fun,x)
% JACCSD Jacobian through complex step differentiation
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
z=fun(x);
n=numel(x);
m=numel(z);
A=zeros(m,n);
h=n*1e-8;
for k=1:n
    x1=x;
    % use real number first
    x1(k) = x1(k) + h;
    % x1(k)=x1(k)+h*1i;
    A(:,k)=(fun(x1) - fun(x))/h;
end
