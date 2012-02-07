function dis = inner(p,g,w,price,alp,v,egr,s)
%this function maximizes utility given a coefficient guess g, and returns the distance between
%maximized utility coefficients and the guess.  The goal (equilibrium) is
%to reduce this distance to zero.

options=optimset('Display','iter','jacobian','off','TolFun',1e-1,'TolX',1e-1,'TolCon',1e-6,'DerivativeCheck','on',...
    'GradObj','off');%,'Algorithm','active-set');

%sobol points
s = egr(1) + (egr(end)-egr(1))*s;
s = sortrows(s);

scl = egr(end)^2;

init = g;

cp = fminunc(@(x) util_FOC(p,g,w,price,scl,x,alp,s,v,egr),init,options);

dis = norm(cp-g);
display(dis)

end

