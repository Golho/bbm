function [sigma, dlambda, ep_eff] = update_variables(sigma_old, ep_eff_old, delta_eps, Dstar, mp, rotation)
% UPDATE_VARIABLES Calculate the new stress for a elastoplastic strain change
% 
%     sigma_old     - The stress in the last calculated point in equilibrium (LCPE)
%       [s_11; s_22; s_12]
%     ep_eff_old    - The effective platic strain in the LCPE
%     delta_eps     - The change in total strain [d_eps_11; d_eps_22; d_eps_12]
%     Dstar         - The linear elastic stiffness matrix [3x3]
%     mp            - Material parameters [F, G, H, L, k, n]
% 
%     sigma         - The stress at the new point [s_11; s_22; s_12]
%     dlambda       - The change of the plastic multiplier from the LCPE
%     ep_eff        - The new effective plastic strain

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

% TRIAL STRESS
sigma_t = Dstar*delta_eps + sigma_old;

% INSERT TRIAL STRESS INTO YIELD FUNCTION
s_t_prim = Tmatrix*Lmat*sigma_t;
f_t = sigma_y0*(sqrt(s_t_prim'*P*s_t_prim) - (1+k*ep_eff_old^n));

function [f] = yield(dlambda)
%YIELD Yield function value dependent on dlambda (used to find dlambda when
% f = 0)
    sigma_y = sigma_y0*(1+k*(ep_eff_old + dlambda)^n);
    Amatrix = eye(3)+(sigma_y0^2/(sigma_y)*Dstar*Lmat'*Tmatrix'*P*Tmatrix*Lmat*dlambda);
    s_2 = Tmatrix*Lmat*(Amatrix\sigma_t);

    f = sigma_y0^2*s_2'*P*s_2 - sigma_y^2;
end

% USE VALUE OF YIELD FUNCTION TO DETERMINE THE RESPONSE
if f_t <= 0
    % ELASTIC RESPONSE
    sigma = sigma_t;
    dlambda = 0;
    ep_eff = ep_eff_old;
    return;
else
    % PLASTIC RESPONSE
    dlambda = fzero(@yield, [0 10]);
    ep_eff = ep_eff_old + dlambda;
    sigma_y = sigma_y0*(1+k*(ep_eff)^n);
    Amatrix = eye(3)+(sigma_y0^2/(sigma_y)*Dstar*Lmat'*Tmatrix'*P*Tmatrix*Lmat*dlambda);
    sigma = Amatrix\sigma_t;
    return;
end
end

