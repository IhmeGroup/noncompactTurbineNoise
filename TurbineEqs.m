function[res] = TurbineEqs(p)
%	A really ugly way of passing these parameters through
	global ess;
	global rho0;
	global rho1;
	global T0;
	global T1;
	global R;
	global M1;
	global U1;
	global gamma;
	global gm1;
	global gp1;
	global gm1o2;
	global gp1o2;
	global returnA;
	global returnTheta;
	A1 = returnA(0);

%	Unpack the solution vector
	rho = p(1);
	U 	= p(2);
	V 	= p(3);
	M	= p(4);
	T	= p(5);

	A = returnA(ess);
	theta = returnTheta(ess);

	res(1) = rho1*U1*A1 - rho*U*A;
	res(2) = U*atan(theta) - V;
	res(3) = M*M - (U*U + V*V)/(gamma*R*T);
	res(4) = rho - rho0*(1+gm1o2*M*M)^(-1/gm1);
	res(5) = T - T0*(1+gm1o2*M*M)^(-1);
end
