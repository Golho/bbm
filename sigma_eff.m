function [sigma_e] = sigma_eff(sigma, mp, rotation)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
F = mp(1);
G = mp(2);
H = mp(3);
L = mp(4);
k = mp(5);
n = mp(6);

sigma_y0 = sqrt(3/(2*(F + G + H)));

Lmat = Lmatrix(rotation);
P = [F+G, -F,  -G,   0;
    -F,  F+H,  -H,   0;
    -G,   -H, G+H,   0;
     0,    0,   0, 2*L];
 
P_hat = Lmat'*P([1 2 4], [1 2 4])*Lmat;

sigma_e = sqrt(sigma_y0^2*sigma'*P_hat*sigma);
end

