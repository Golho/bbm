% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
 % 1) Export the required variables from pdetool and create a MATLAB script
 %    to perform operations on these.
 % 2) Define the problem completely using a MATLAB script. See
 %    http://www.mathworks.com/help/pde/examples/index.html for examples
 %    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'PlotBoxAspectRatio',[613.11999999999989 575.40000000000009 16315.566110390817]);
set(ax,'XLimMode','auto');
set(ax,'YLim',[0 0.050000000000000003]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pdeellip(0.022499999999999999,0.025000000000000001,0.0025000000000000001,0.0025000000000000001,...
0,'C2');
pdeellip(0.02,0.014999999999999999,0.01,0.0050000000000000001,...
0,'P1');
pdepoly([ 0.01,...
 0.01,...
 0.02,...
 0.02,...
 0.029999999999999999,...
 0.029999999999999999,...
],...
[ 0.014999999999999999,...
 0.02,...
 0.029999999999999999,...
 0.040000000000000001,...
 0.029999999999999999,...
 0.014999999999999999,...
],...
 'P2');
pderect([0 0.040000000000000001 0.050000000000000003 0],'B1');
pdeellip(0.01,0.029999999999999999,0.0025000000000000001,0.0025000000000000001,...
0,'C1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','((P1-P2)+B1)-C1-C2')

% Boundary conditions:
pdetool('changemode',0)
pdetool('removeb',[18 ]);
pdetool('removeb',[9 ]);
pdesetbd(19,...
'dir',...
1,...
'1',...
'0')
pdesetbd(18,...
'dir',...
1,...
'1',...
'0')
pdesetbd(17,...
'dir',...
1,...
'1',...
'0')
pdesetbd(16,...
'dir',...
1,...
'1',...
'0')
pdesetbd(13,...
'dir',...
1,...
'1',...
'0')
pdesetbd(12,...
'dir',...
1,...
'1',...
'0')
pdesetbd(11,...
'dir',...
1,...
'1',...
'0')
pdesetbd(10,...
'dir',...
1,...
'1',...
'0')
pdesetbd(6,...
'dir',...
1,...
'1',...
'0')
pdesetbd(5,...
'dir',...
1,...
'1',...
'0')
pdesetbd(4,...
'dir',...
1,...
'1',...
'0')
pdesetbd(3,...
'dir',...
1,...
'1',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')
pdetool('refine')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','2568','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');