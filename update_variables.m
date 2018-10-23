function [sigma, dlambda, ep_eff] = update_variables(sigma_old, ep_eff_old, delta_eps, Dstar, mp, rotation)
% Calculate the new stress for a elastoplastic strain change
% 
% sigma_old - The stress in the last calculated point in equilibrium (LCPE)
% ep_eff_old - The effective platic strain in the LCPE
% delta_eps - The change in total strain
% Dstar - The linear elastic stiffness matrix
% mp - Material parameters [F, G, H, L, k, n]

F = mp(1);
G = mp(2);
H = mp(3);
L = mp(4);
k = mp(5);
n = mp(6);

sigma_y0 = sqrt(3/(2*(F + G + H)));

Lmat = Lmatrix(rotation);
P = [F+G, -F, -G, 0;
    -F, F+H, -H, 0;
    -G, -H, G+H, 0;
    0, 0, 0, 2*L];

sigma_t = Dstar*delta_eps + sigma_old;

% check if elastic of plastic response
s_t_prim = Tmatrix*Lmat*sigma_t;
f_t = sigma_y0*(sqrt(s_t_prim'*P*s_t_prim)-(1+k*ep_eff_old^n))
function [ f ] = lambda_zero( dlambda )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
hard = 1+k*(ep_eff_old+dlambda)^n;
Amatrix = eye(3)+(sigma_y0/(hard)*Dstar*Lmat'*P([1 2 4],:)*Tmatrix*Lmat*dlambda);
sigma_2 = Amatrix\sigma_t;

f = sigma_y0^2*(sigma_2'*Lmat'*Tmatrix'*P*Tmatrix*Lmat*sigma_2 - hard);

end

if f_t <= 0
    % ELASTIC RESPONSE
    sigma = sigma_t;
    dlambda = 0;
    ep_eff = ep_eff_old;
    return;
else
    % PLASTIC RESPONSE
    %dlambda = fzero(@lambda_zero, 0.0009);
    dlambda = 0.001;
    ep_eff = ep_eff_old + dlambda;
    hard = 1+k*(ep_eff)^n;
    Amatrix = eye(3)+(sigma_y0/(hard)*Dstar*Lmat'*Tmatrix'*P*Tmatrix*Lmat*dlambda);
    sigma = Amatrix\sigma_t;
end
end

