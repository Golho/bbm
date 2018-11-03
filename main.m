clear
close

MAX_LOAD = 8e8; % [Pa]
MAX_DISPLACEMENT = 0.2e-3; % [m]
nbrLoadSteps = 30; % number of load steps

load simple_mesh.mat
dtau_x = MAX_LOAD/nbrLoadSteps;
du = MAX_DISPLACEMENT/nbrLoadSteps;
theta = deg2rad(90); % [rad]
control = 0; % 0 - force control, 1 - displacement control

ptype = 1; % 1 - plane stress, 2 - plane strain
th = 1e-3; % [m]
ep = [ptype th];
NR_TOL = 10e-8; % Tolerance in the Newton Raphson loop

% Calculate the components in the P-matrix
sigma_y0 = [450 250 250]*1e6; % [Pa]
F = 1/2 * sum([1 1 -1]./sigma_y0.^2);
G = 1/2 * sum([1 -1 1]./sigma_y0.^2);
H = 1/2 * sum([-1 1 1]./sigma_y0.^2);
L = F + G + H;
E = 70e9; % [Pa]
nu = 0.32;
mp = [F G H L 5 0.45]; % F G H L k n

[bc, df, edof, dof, coord, enod] = TopConstMod(p, t, dtau_x, du, th, control);
plotDof = 7;
nbrDofs = numel(dof); % Number of degrees of freedom
nbrElems = size(edof, 1);
nbrNodes = size(dof, 1);
[Ex, Ey] = coordxtr(edof, coord, dof, 3);

a = zeros(nbrDofs, 1); % Displacement vector
K = zeros(nbrDofs); % Stiffness matrix
f = zeros(nbrDofs, 1); % External force vector
f_int = zeros(nbrDofs, 1); % Internal force vector
eps = zeros(nbrElems, 3); % Strain matrix
ep_eff = zeros(nbrElems, 1); % Effective plastic strain vector
stress = zeros(nbrElems, 3);

% INITIALIZE K
D_star = hooke(ptype, E, nu);
for e = 1:nbrElems
    ex = Ex(e, :);
    ey = Ey(e, :);
    index = edof(e, 2:end);
    Ke = plante(ex, ey, ep, D_star);
    K(index, index) = K(index, index) + Ke; 
end

if control == 1
    dbc = bc(:, 2)*[ones(1, nbrLoadSteps) -ones(1, nbrLoadSteps)];
    nbrLoadSteps = 2*nbrLoadSteps;
end
plot_dof = 3;
loads = zeros(1, nbrLoadSteps+1);
displacements = zeros(1, nbrLoadSteps+1);
displacements(1) = a(plot_dof);
loads(1) = f_int(plot_dof);
for loadStep = 1:(nbrLoadSteps*2)
    disp("Load step: " + loadStep);
    % APPLY LOAD
    if control == 0
        if loadStep > nbrLoadSteps
            f = f - df;
            K = zeros(nbrDofs);
            for e = 1:nbrElems
                ex = Ex(e, :);
                ey = Ey(e, :);
                index = edof(e, 2:end);
                Ke = plante(ex, ey, ep, D_star);
                K(index, index) = K(index, index) + Ke; 
            end
        else
            f = f + df;
        end
        bcIncr = bc;
    elseif control == 1
        bcIncr = [bc(:, 1) dbc(:, loadStep)];
    end
    
    % CALCULATE RESIDUAL
    resNorm = 1; % set residual norm to 1 to iterate at least once in the
    % NR loop
    res = f_int - f;
    c = 0;
    
    nr_stress = zeros(nbrElems, 3); % Temporary (not in equilibrium) stress
    nr_strain = zeros(nbrElems, 3); % Temporary (not in equilibrium) strain difference
    nr_ep_eff = zeros(nbrElems, 1);
    %% Newton-Raphson loop of iteration

    while resNorm > NR_TOL % NR LOOP
        c = c + 1;
        % Solve 
        da = solveq(K, -res, bcIncr);
        a = a + da;
        
        % UPDATE INTERNAL FORCE VECTOR AND GLOBAL STIFFNESS MATRIX
        elemDisplacements = extract(edof, da); % element displacement vector
        [~, tmp_strains] = plants(Ex, Ey, ep, eye(3), elemDisplacements);
        nr_strain = nr_strain + tmp_strains;
        % RESET THE VECTOR AND MATRIX
        K = zeros(nbrDofs);
        f_int = zeros(nbrDofs, 1);
        dlambda = zeros(nbrElems, 1);
        for e = 1:nbrElems
            ex = Ex(e, :);
            ey = Ey(e, :);
            index = edof(e, 2:end); % element degree of freedoms
            if t(4, e) == 1 % Sub domain with rotated orthotropic directions
                rotation = theta;
            else
                rotation = 0;
            end
            [tmp_stress, dlambda(e), nr_ep_eff(e)] = update_variables(stress(e, :)', ep_eff(e), nr_strain(e, :)', D_star, mp, rotation);
            nr_stress(e, :) = tmp_stress';
            D_te = alg_tan_stiff(nr_stress(e, :)', dlambda(e), ep_eff(e), D_star, mp, rotation);
            
            % ASSEMBLY STIFFNESS MATRIX
            Ke = plante(ex, ey, ep, D_te);
            K(index, index) = K(index, index) + Ke;
            
            % ASSEMBLY INTERNAL FORCE VECTOR
            f_inte = plantf(ex, ey, ep, nr_stress(e, :));
            f_int(index) = f_int(index) + f_inte';
        end
        % UPDATE RESIDUAL
        res = f_int - f;
        res(bc(:, 1)) = 0; % The Dirichlet BC nodes do not contribute to the residual
        resNorm = norm(res);
        disp("Residual norm: " + resNorm);
        disp("Iteration step: " + c);
        
        bcIncr(:, 2) = 0; % Set displacement increment to zero in NR loop
    end
    %% 
    % The Newton-Raphson loop ends and the calculated displacements and load 
    % is inserted into the memory vectors.
    % NOW IN EQULIBRIUM
    stress = nr_stress; % After NR loop, NR stresses are in equilibrium
    ep_eff = nr_ep_eff;
    % ACCEPT DISPLACEMENTS AND SAVE FOR PLOT DOF
    displacements(loadStep+1) = a(plotDof);
    loads(loadStep+1) = f_int(plotDof);
end
eldisp2(Ex, Ey, elemDisplacements, [3 2 0], 1);

figure
plot(loads, displacements);
if control == 0
    title("Force vs. displacements in dof " + plot_dof + " (load controlled)");
elseif control == 1
    title("Force vs. displacements in dof " + plot_dof + " (displacement controlled)");
end
xlabel("Force [N] in dof");
ylabel("Displacement [m] in dof");
grid on
grid minor

% stress = one row per element: [sigma_11, sigma_22, sigma_12]
elemEffStress = sqrt(sum(stress(:, 1:2).^2, 2) - stress(:, 1).*stress(:, 2) + ...
                     3*stress(:, 3).^2);

nodeEffStress = zeros(nbrElems, 3);
for node = 1:nbrNodes
    i = find(enod == node); % Find the elements where node is present
    [iRows, iCols] = ind2sub(size(enod), i); % Go from linear index to subscript index
    nodeEffStress(i) = sum(elemEffStress(iRows))/length(iRows); % extrapolate the stress from those elements
end

% CALC THE NEW COORDINATES
newEx = Ex + elemDisplacements(:, 1:2:end);
newEy = Ey + elemDisplacements(:, 2:2:end);

figure
axis equal
grid on
cb = colorbar;
ylabel(cb, "\sigma_e - [Pa]")
patch(newEx', newEy', nodeEffStress')
