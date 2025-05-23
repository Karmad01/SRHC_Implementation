function R = reachability_matrix(A,B,k)
%REACHABILITY_MATRIX Summary of this function goes here
%   Detailed explanation goes here
R = zeros(height(B),k*width(B));
idx = 1;
for ii=1:k
    R(:,idx:idx+width(B)-1) = A^(k-ii)*B;
    idx = idx+width(B);
end
end

