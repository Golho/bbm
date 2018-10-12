function inc = apply_load2(p,nstep,edof,ex,ey, inner_elements_right, inner_elements_left,ndof)
%
%Application of load for geometry 2 in Assignment 1.
%
%Input: p - total force applied to geometry
%       nstep - number of time steps
%       edof - element topology matrix
%       ex   - Elements x coordinates
%       ey   - Elements y coordiantes
%       inner_elements_right - Elements on the right side of the geometry
%       inner_elements_left  - Elements on the left side of the geometry
%       ndof - number of degrees of freedom
%
%Output: inc - incremental force vector 
%

inc = zeros(ndof,1);
for i=1:length(inner_elements_right)
    el = inner_elements_right(i);
    x1 = ex(el,1);
    x2 = ex(el,4);
    y1 = ey(el,1);
    y2 = ey(el,4);
    
    dx = abs(x1-x2);
    dy = abs(y1-y2);
    
    dofx1 = edof(el,2);
    dofx2 = edof(el,8);
    dofy1 = edof(el,3);
    dofy2 = edof(el,9);
    inc([dofx1;dofx2]) = inc([dofx1;dofx2]) + dy*p/nstep/2*[1;1];
    inc([dofy1;dofy2]) = inc([dofy1;dofy2]) + dx*p/nstep/2*[1;1];
end

for i=1:length(inner_elements_left)
    el = inner_elements_left(i);
    x1 = ex(el,1);
    x2 = ex(el,4);
    y1 = ey(el,1);
    y2 = ey(el,4);
    
    dx = abs(x1-x2);
    dy = abs(y1-y2);
    
    dofx1 = edof(el,2);
    dofx2 = edof(el,8);
    dofy1 = edof(el,3);
    dofy2 = edof(el,9);
    inc([dofx1;dofx2]) = inc([dofx1;dofx2]) - dy*p/nstep/2*[1;1];
    inc([dofy1;dofy2]) = inc([dofy1;dofy2]) + dx*p/nstep/2*[1;1];
end

end