function [sigma, dlambda, ep_eff] = update_variables(sigma_old, ep_eff_old, delta_eps, Dstar, mp, rotation)
% Calculate the new stress for a elastoplastic strain change
% 
% sigma_old - The stress in the last calculated point in equilibrium (LCPE)
% ep_eff_old - The effective platic strain in the LCPE
% delta_eps - The change in total strain
% Dstar - The linear elastic stiffness matrix
% mp - Material parameters [F, G, H, L, k, n]
outputArg1 = inputArg1;
outputArg2 = inputArg2;

F = mp(1);
G = mp(2);
H = mp(3);
L = mp(4);
k = mp(5);
n = mp(6);

sigma_y0 = sqrt(3/(2*(F + G + H)));

P = [F+G, -F, -G, 0;
    -F, F+H, -H, 0;
    -G, -H, G+H, 0;
    0, 0, 0, 2*L];

end

