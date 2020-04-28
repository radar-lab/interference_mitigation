% Function of the adaptive noise canceler
% Parameters:
% Iuput:    pri --- the primary channel input, should be M-by-1 vector
% Iuput:    ref --- the primary channel input, should be M-by-1 vector
% Iuput:    M   --- the length of the adaptive filter
% Iuput:    mu  --- the step size
% Output:   e   --- the estimation error
% Output:   wo  --- the adaptive filter tap-coefficient vector
function [e, wo] = af_with_ref_cx(pri, ref , M, mu)
% Find the data size
N = size(pri,1);
% Initialize output vector y to zero
y = zeros(N,1);
% Initialize error vector e to zero
e = zeros(N,1);
% Initialize weight vector to [1 0 0 ...]
wo = zeros(M,1);
wo(1) = 1;
% Initialize a vector for holding filter input
refvec = zeros(M,1);

for k = 1 : N
    % Form the reference input
    refvec = [ref(k) 
        refvec(1:M-1)];
    % Form the adaptive filter output
    y(k) = wo'*refvec;
    % Form the estimation error
    e(k) = pri(k) - y(k);
    % Update the tap-cofficient vector of the adaptive filter
    wo = wo + mu*refvec*conj(e(k));
end 