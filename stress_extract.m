function stress=stress_extract(es_all,edof)
%-------------------------------------------------------------
% PURPOSE
%  Extract nodal effective stress (von Mises) from the element
%  stress matrix. Four node isoparametric element with 4 gauss points
%
% INPUT:  es = [ sigx sigy [sigz] tauxy    element stress matrix
%                  ......              ]   one row for each element Gauss
%                                          point
%         edof:  topology matrix
%
% OUTPUT: stress:   nodal effective stress per element
%
%-------------------------------------------------------------
%

% Element effective stress
seff_el = zeros(size(edof,1),4);
for el=1:size(edof,1)
    es = es_all(:,:,el);
    for gp=1:4
        tmp = sqrt((3/2)*(es(gp,1)^2+es(gp,2)^2+es(gp,3)^2+2*es(gp,4)^2 ...
              - (1/3)*(sum(es(gp,1:3)))^2));
        seff_el(el,gp) = tmp;
    end
    el = el + 1;
end
seff_av = zeros(size(seff_el));

for dof=1:max(max(edof(:,2:9)))
	[elm,pos] = find(edof(:,2:9)==dof);
	nconnel   = length(elm); % Number of elements connected to this dof
    seff_av   = 0;
	for el=1:nconnel
	   seff_av = seff_av + seff_el(elm(el),ceil(pos(el)/2));
	end
	seff_av = seff_av/nconnel;
	for el=1:nconnel
	 	stress(elm(el),ceil(pos(el)/2)) = seff_av;
	end
end