function [bc,df,eDof,dof,coord,enod] = TopConstMod(p,t,dtau_x,du,th,control) 
% [bc,df,edof,dof,coord,enod,plot_dof] = TopConstMod(p,t,dtau_x,du,th,control) 
%-------------------------------------------------------------
% PURPOSE
% Extract bc,df,edof,dof,coord,enod,plot_dof from a mesh generated in
% pdetool. 
%
%
% INPUT:  p         Variable retrieved from pdetool contains 
%                   x- and y- coordinates of all nodes
%         t         Variable retrieved from pdetool contains information 
%                   of what nodes are connected to each element
%         dtau_x    Incremental traction stress applied at the left boundary
%         du        Incremental displacement on the left boundary
%         th        thickness
%         control   control=0 force controlled displacement
%                   control=1 displacement controlled displacement
%
% OUTPUT: bc        First column dof, second column known displacement
%         df        External applied force each time step (ndof X 1)
%         edof      Element topology matrix (how the dofs are conected to each element) (nelm X 7)
%         dof       Dof conexted to each node (nnod X 2)
%         coord     x- and y- coordinates of all nodes (nnod X 2)
%         enod      Element topology matrix (how the nodes are conected to each element) (nelm X 3)
%
%-------------------------------------------------------------
% LAST MODIFIED: Olov Günther-Hanssen 2018-11-03
%-------------------------------------------------------------

nen=3;                                  % number of nodes per element
enod = t(1:nen,:)';                     % 
nbrElem = size(enod,1);                 % number of elements  
coord = p';                             %
nbrDofNode=2;                           % number of dofs per node
nbrNodes = size(coord,1);               % number of nodes
dof = [(1:nbrNodes)',(nbrNodes+1:2*nbrNodes)'];       % give each dof a number
nbrDofs = max(max(dof));                % number of dofs  
eDof=zeros(nbrElem,nen*nbrDofNode+1);   % allocate space for edof

% Generate edof from enod and dof      
for ie = 1:nbrElem
   eDof(ie,:)=[ie dof(enod(ie,1),:), dof(enod(ie,2),:), dof(enod(ie,3),:)];
end

% Extract global element coordinates vectors and draw mesh
[ex,ey] = coordxtr(eDof,coord,dof,nen);
eldraw2(ex,ey,[1 2 0])

% Find the coordinates for the corners
minx = min(coord(:,1));                   % minimum x-coord. of the geometry
minXi = find(minx == coord(:, 1))
maxx = max(coord(:,1));                   % maximum x-coord. of the geometry
maxXi = find(maxx == coord(:, 1))
miny = min(coord(:,2));                   % minimum y-coord. of the geometry
minYi = find(miny == coord(:, 2))
maxy = max(coord(:,2));                   % maximum x-coord. of the geometry
maxYi = find(maxy == coord(:, 2));

% SET BOUNDARY CONDITIONS

% FIND DOFS CONNECTED TO COMMON BOUNDARY CONDITIONS
lowerLeftCornerNode = intersect(minXi, minYi);
lowerRightCornerNode = intersect(maxXi, minYi);
lowerLeftDof = dof(lowerLeftCornerNode, :);
lowerRightDofY = dof(lowerRightCornerNode, 2);

% Define df and bc
df = zeros(nbrDofs,1);                       % allocate space for df
bc = [lowerLeftDof', zeros(2,1);
      lowerRightDofY, 0];
if control == 1
    % FIND DOFS CONNECTED TO DISPLACEMENT CONTROLLED BOUNDARY CONDITIONS
    upperBoundaryDofs = dof(maxYi, 2);
    [~, I] = min(abs(coord(minYi, 1)-mean([maxx minx])));
    lowerMiddleBoundaryNode = minYi(I);
    lowerMiddleDofY = dof(lowerMiddleBoundaryNode, 2);
    
    bc = [bc;
          upperBoundaryDofs, du*ones(size(upperBoundaryDofs));
          lowerMiddleDofY, 0];
elseif control == 0
    % FIND DOFS CONNECTED TO FORCE CONTROLLED BOUNDARY CONDITIONS
    upperLeftCornerNode = intersect(minXi, maxYi);
    upperLeftCornerDofX = dof(upperLeftCornerNode, 1);
    bc = [bc;
          upperLeftCornerDofX, 0];

    [~, coordOrder] = sort(coord(maxYi, 1));
    upperNodesSorted = maxYi(coordOrder);
    upperDofsSorted = dof(upperNodesSorted, 2);
    L = diff(coord(upperNodesSorted, 1)); % Vector of all length between the nodes
    df(upperDofsSorted(1:end-1)) = dtau_x*L/2; % Divide half of the load to the "left" dofs
    df(upperDofsSorted(2:end)) = df(upperDofsSorted(2:end)) + dtau_x*L/2; % Divide of the load to the "right" dofs
end
 
 df = df*th;                        % multiply df with the thickness

% Below are optional plot routines where the added bc and df may be seen
% a=zeros(nbrDofs,1);
% a(bc(:,1))=1;
% ed=extract(eDof,a);
% figure(1)
% fill(ex',ey',ed(:,[1 3 5])');     %(markes all dofs with bc in x-dir)
% figure(2)
% fill(ex',ey',ed(:,[2 4 6])');     %(markes all dofs with bc in y-dir)
% 
% 
% a=zeros(nbrDofs,1);
% a(bc(:,1))=bc(:,2);
% ed=extract(eDof,a);
% figure(3)
% fill(ex',ey',ed(:,[1 3 5])');   %(markes all dofs with bc in x-dir)
% figure(4)
% fill(ex',ey',ed(:,[2 4 6])');   %(markes all dofs with bc in y-dir)
% 
% 
% 
% ed=extract(eDof,df);
% figure(5)
% fill(ex',ey',ed(:,[1 3 5])');   %(markes all dofs with fl in x-dir)
% figure(6)
% fill(ex',ey',ed(:,[2 4 6])');   %(markes all dofs with fl in y-dir)
