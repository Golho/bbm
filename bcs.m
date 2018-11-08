function [res] = bcs(xa,xb)
%bcs Return the boundary conditions for the flow for constitutive relation
% 2
res = [xa(2);
       xb(1)];
end

