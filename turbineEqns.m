function[res] = TurbineEqs(p, stagnation_state, inlet_state, A, theta_bar)
%	Unpack the stagnation state
	gamma 		= stagnation_state(1);	
	R	 		= stagnation_state(2);	
	rho0 		= stagnation_state(3);	
	T0	 		= stagnation_state(4);	
	M1 			= stagnation_state(5);	

%	Pre-compute some quantities for brevity
	gm1 = gamma -1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;
	
%	Unpack the inlet state, b/c that is a thing that I need
	rho1		= inlet_state(1);
	T1			= inlet_state(2);
	a1			= inlet_state(3);
	U1			= inlet_state(4);
	V1			= inlet_state(5);
	A1			= inlet_state(6);

%	Unpack the solution vector
	rho = p(1);
	U 	= p(2);
	V 	= p(3);
	M	= p(4);
	T	= p(5);

	W1 = sqrt(U1*U1 + V1*V1);
	W = sqrt(U*U + V*V);
	res(1) = rho1*W1*A1 - rho*W*A;
	res(2) = U*atan(theta_bar) - V;
	res(3) = M*M - (U*U + V*V)/(gamma*R*T);
	res(4) = rho - rho0*(1+gm1o2*M*M)^(-1/gm1);
	res(5) = T - T0*(1+gm1o2*M*M)^(-1);
end
