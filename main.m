clear
close

MAX_LOAD = 8e8; % [Pa]
MAX_DISPLACEMENT = 0.2e-3; % [m]
nbrLoadSteps = 20; % number of load steps

ctrlCase = input("Skriv in lastfall (1-4):");
control = 0; % 0 - force control, 1 - displacement control

dtau_x = 0;
du = 0;

switch ctrlCase
    case 1
        control = 0;
        load simple_mesh.mat;
        theta = deg2rad(0); % [rad]
        MAX_LOAD = 4e8; % [Pa]
        mainDir = 2;
        plotDof = 7;
        dtau_x = MAX_LOAD/nbrLoadSteps;
        nbrLoadSteps = nbrLoadSteps*2;
    case 2
        control = 0;
        load simple_mesh.mat;
        theta = deg2rad(90);
        MAX_LOAD = 8e8;
        mainDir = 1;
        plotDof = 7;
        dtau_x = MAX_LOAD/nbrLoadSteps;
        nbrLoadSteps = nbrLoadSteps*2;
    case 3
        control = 1;
        load Geom_mesh_rough.mat;
        theta = deg2rad(90);
        plotDof = 3;
        du = MAX_DISPLACEMENT/nbrLoadSteps;
    case 4
        control = 1;
        load Geom_mesh_rough.mat;
        theta = deg2rad(60);
        plotDof = 3;
        du = MAX_DISPLACEMENT/nbrLoadSteps;
    case 5
        control = 1;
        load Geom_mesh_rough.mat;
        theta = deg2rad(0);
        plotDof = 3;
        du = MAX_DISPLACEMENT/nbrLoadSteps;
    otherwise
        error("Du har valt ett inkorrekt alternativ");
end



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

if ctrlCase == 5
    dbc = bc(:, 2)*[ones(1, nbrLoadSteps) -ones(1, nbrLoadSteps)];
    nbrLoadSteps = nbrLoadSteps*2;
end

nbrDofs = numel(dof); % Number of degrees of freedom
nbrElems = size(edof, 1);
nbrNodes = size(dof, 1);
[Ex, Ey] = coordxtr(edof, coord, dof, 3);

a = zeros(nbrDofs, 1); % Displacement vector
K = zeros(nbrDofs); % Stiffness matrix
f = zeros(nbrDofs, 1); % External force vector
f_int = zeros(nbrDofs, 1); % Internal force vector
eps = zeros(nbrElems, 3); % Strain matrix
ep_eff = zeros(nbrElems, 1);
ep_eff_memory = zeros(nbrElems, nbrLoadSteps);  % Effective plastic strain vector
stress = zeros(nbrElems, 3);

nr_stress = zeros(nbrElems, 3); % Temporary (not in equilibrium) stress
nr_strain = zeros(nbrElems, 3); % Temporary (not in equilibrium) strain difference
nr_ep_eff = zeros(nbrElems, 1);

% INITIALIZE K
D_star = hooke(ptype, E, nu);
for e = 1:nbrElems
    ex = Ex(e, :);
    ey = Ey(e, :);
    index = edof(e, 2:end);
    Ke = plante(ex, ey, ep, D_star);
    K(index, index) = K(index, index) + Ke; 
end

loads = zeros(1, nbrLoadSteps+1);
displacements = zeros(1, nbrLoadSteps+1);
displacements(1) = a(plotDof);
loads(1) = f_int(plotDof);
for loadStep = 1:nbrLoadSteps
    disp("Load step: " + loadStep);
    % APPLY LOAD
    if control == 0
        if loadStep > nbrLoadSteps/2
            f = f - df;
        else
            f = f + df;
        end
        bcIncr = bc;
    elseif control == 1
        if ctrlCase == 3 || ctrlCase == 4
            bcIncr = bc;
        elseif ctrlCase == 5
            bcIncr = [bc(:, 1) dbc(:, loadStep)];
        end
    end
    
    % CALCULATE RESIDUAL
    resNorm = 1; % set residual norm to 1 to iterate at least once in the
    % NR loop
    res = f_int - f;
    c = 0;

    nr_stress = nr_stress*0; 
    nr_strain = nr_strain*0; 
    nr_ep_eff = nr_ep_eff*0;
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
        K = K*0;
        f_int = f_int*0;
        
        for e = 1:nbrElems
            ex = Ex(e, :);
            ey = Ey(e, :);
            index = edof(e, 2:end); % element degree of freedoms
            if t(4, e) == 1 % Sub domain with rotated orthotropic directions
                rotation = theta;
            else
                rotation = 0;
            end
            [tmp_stress, dlambda, nr_ep_eff(e)] = update_variables(stress(e, :)', ep_eff(e), nr_strain(e, :)', D_star, mp, rotation);
            nr_stress(e, :) = tmp_stress';
            
            % Check if the response is elastic or plastic and set the
            % tangential stiffness matrix accordingly
            if dlambda ~= 0
                D_te = alg_tan_stiff(nr_stress(e, :)', dlambda, nr_ep_eff(e), D_star, mp, rotation);
            else
                D_te = D_star;
            end

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
    ep_eff_memory(:, loadStep) = nr_ep_eff;

    % ACCEPT DISPLACEMENTS AND SAVE FOR PLOT DOF
    displacements(loadStep+1) = a(plotDof);
    loads(loadStep+1) = f_int(plotDof);
end
elemDisplacements = extract(edof, a); % element displacement vector
eldisp2(Ex, Ey, elemDisplacements, [3 2 0], 1);

if ctrlCase == 1 || ctrlCase == 2
    figure;
    crossArea = 1*th;
    sigma = 2*loads / crossArea;
    epsilon = displacements;

    plot(epsilon, sigma, 'Color', 'k', 'LineWidth', 1);
    hold on
    plot(epsilon, sigma_y0(mainDir)*ones(size(epsilon)), 'Color', 'r');
    hold off
    if ctrlCase == 1
        legend(["Load path", "\sigma_{y0, CD}"]);
    else
        legend(["Load path", "\sigma_{y0, RD}"]);
    end
    title("Uniaxial stress vs. strain in dof " + plotDof + " (load controlled)");
    xlabel("Strain \epsilon");
    ylabel("Stress \sigma [Pa]");
    grid on
    grid minor
    hold off
end
% stress = one row per element: [sigma_11, sigma_22, sigma_12]
%% CALCULATING AND PLOTTING EFFECTIVE STRESS
elemEffStress = zeros(nbrElems, 1);
for e = 1:nbrElems
    if t(4, e) == 1 % Sub domain with rotated orthotropic directions
        rotation = theta;
    else
        rotation = 0;
    end
    elemEffStress(e) = sigma_eff(stress(e, :)', mp, rotation);
end

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

figure
axis equal
grid on
cb2 = colorbar;
patch(newEx', newEy', double(ep_eff_memory(:, end)~=0))
