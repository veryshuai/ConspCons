function dis = inner(p,g,w,price,alp,v,egr,s)
%this function maximizes utility given a coefficient guess g, and returns the distance between
%maximized utility coefficients and the guess.  The goal (equilibrium) is
%to reduce this distance to zero.

options=optimset('Display','iter','jacobian','on','TolFun',5e-2,'TolX',5e-2,'TolCon',1e-6,'DerivativeCheck','on',...
    'GradObj','on','MaxIter',20,'LargeScale','off');%,'Algorithm','active-set');

%sobol points
s = egr(1) + (egr(end)-egr(1))*s;
s = sortrows(s);

scl = egr(end)^2;

init = g;

cp = fminunc(@(x) util_FOC(p,x,w,price,scl,g,alp,s,v,egr),init,options);

dis = norm(cp-g)/norm(g);
if isnan(dis) == 1 || abs(dis) == inf
    dis = 1e4;
end
%display(dis)

end

