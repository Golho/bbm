function [Dats] = alg_tan_stiff(sigma, dlambda, ep_eff, Dstar, mp, rotation )
% Algorithmic tangential stiffness matrix
%   Input:
%       sigma       - Stress at a point [s_11; s_22, s_12]
%       dlambda     - The change of the plastic multiplier from the LCPE
%       ep_eff       - The effective plastic strain
%       Dstar       - The linear elastic stiffness matrix [3x3] 
%       mp          - Material parameters [F, G, H, L, k, n]
%       rotation    - The rotation of the coordinate system
%   Output:
%       Dats        - Algoritmic tangential stiffness matrix [3x3]

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

sigma_eff = sqrt(sigma_y0^2*sigma'*P_hat*sigma);

df = sigma_y0^2/sigma_eff*P_hat*sigma;

D_a = inv(inv(Dstar) + dlambda*(sigma_y0^2/sigma_eff*P_hat-sigma_y0^4/sigma_eff^3*P_hat*(sigma*sigma')*P_hat'));

da = -(sigma_y0*n*k*ep_eff^(n-1)*(-1));

A_a = df'*D_a*df - (-1)*da;

Dats = D_a - 1/A_a*D_a*df*df'*D_a;

end

