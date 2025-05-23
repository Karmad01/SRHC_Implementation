function phi = phi_func(z,phi_max)
%PHI_FUNC Summary of this function goes here
%   Detailed explanation goes here
phi = max(min(z,phi_max),-phi_max);
end

