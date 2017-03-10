function[baseFlow] = buildBaseFlow(boundaries, x_c, y_c, ytop, ybot, stagnation_state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	1. Unpack the boundaries vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	x_tl 		= boundaries(1);
	y_tl 		= boundaries(2);
	x_tr 		= boundaries(3);
	y_tr 		= boundaries(4);
	x_bl 		= boundaries(5);
	y_bl 		= boundaries(6);
	x_br 		= boundaries(7);
	y_br 		= boundaries(8);
	exit_slope 	= boundaries(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	2. Unpack the stagnation state vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	gamma 		= stagnation_state(1);	
	R	 		= stagnation_state(2);	
	rho0 		= stagnation_state(3);	
	T0	 		= stagnation_state(4);	
	M1 			= stagnation_state(5);	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	3. Precompute some stuff for speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	gm1 = gamma -1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;
	N = length(x_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	4. Compute the local normal vector of the cord and the arc length stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	First comput the local number
	dydx = zeros(N,1);
	dydx(1) = (y_c(2) - y_c(1))/(x_c(2) - x_c(1));
	for i = 2:N-1
		dydx(i) = (y_c(i+1) - y_c(i-1))	/(x_c(i+1) - x_c(i-1));
	end
	dydx(N) = (y_c(N) - y_c(N-1))/(x_c(N) - x_c(N-1));
%	Then compute the total cord length which will be used for normalization
	disp('cord length')
	L_c = trapz(x_c, sqrt(1 + dydx.^2))

%	Next, compute the instantaneous cord length to be used as a coordinate for the splines
	s = zeros(N,1);
	for i = 2:N
		s(i) = trapz(x_c(1:i), sqrt(1 + dydx(1:i).^2));
	end
%	Normalize this coordinate by the total cord length to paramaterize the arc
	s = s./L_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	5. Compute the X-sectional area of the geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	To compute the cord location, we need to determine where the normal to a point in the cord intersects the top and bottom boundaries
%	First, allocate memory to store these intersection locations
	xsect_top 	= zeros(N,1);
	ysect_top 	= zeros(N,1);
	xsect_bot 	= zeros(N,1);
	ysect_bot 	= zeros(N,1);
	A			= zeros(N,1);
	theta_bar	= zeros(N,1);

%	Loop over each point on the cord and find where the normal plane intersects the top and bottom geometry to determine the local X-sectional area
	for i = 1:N
%		If the slope is ill-conditioned (i.e. horizontal or vertical), we can do this analytically
		if ((abs(dydx(i)) <= 1E-3) || (abs(dydx(i)) > 1E6) || (isinf(dydx(i))))
			xsect_top(i) = x_c(i);
			xsect_bot(i) = x_c(i);
			[ysect_top(i)] = cleverYtop(xsect_top(i), ytop, x_tl, y_tl, x_tr, y_tr, dydx(i));
			[ysect_bot(i)] = cleverYbot(xsect_top(i), ybot, x_bl, y_bl, x_br, y_br, dydx(i));
%		If the slope and its inverse are finite, then we use a dumb solver to identify the intersection point. Sorry this is so slow, but better techniques seemed to struggle here. R.I.P. 1/11/17 
		else
			xsect_top(i) = x_c(i);
			xsect_bot(i) = x_c(i);
%			xsect_top(i) = ftopf();
			[ysect_top(i)] = cleverYtop(xsect_top(i), ytop, x_tl, y_tl, x_tr, y_tr, dydx(i));
%			xsect_bot(i) = fbotf();
			[ysect_bot(i)] = cleverYbot(xsect_top(i), ybot, x_bl, y_bl, x_br, y_br, dydx(i));
		end
		A(i) = sqrt((ysect_top(i) - ysect_bot(i)).^2 + (xsect_top(i) - xsect_bot(i)).^2);
		theta_bar(i) = atan2(dydx(i), 1);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	6. compute the inlet state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	rho1 = rho0*(1 + gm1o2*M1^2)^(-1/gm1);
	T1 = T0*(1 + gm1o2*M1^2)^(-1);
	a1 = sqrt(gamma*R*T1);
	U1 = M1*cos(atan2(dydx(1),1))*a1;
	V1 = M1*sin(atan2(dydx(1),1))*a1;
	inlet_state = [rho1, T1, a1, U1, V1, A(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	7. Compute the base flow and then spline it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Allocate memory for the solution for efficiency
	rho_bar = zeros(N,1);
	u_bar 	= zeros(N,1);
	v_bar 	= zeros(N,1);
	w_bar 	= zeros(N,1);
	M_bar 	= zeros(N,1);
	T_bar 	= zeros(N,1);
	p_bar 	= zeros(N,1);
	sos_bar = zeros(N,1);

%	Initialize the solver with the inflow state (This is not a very good guess)
	q = [rho1, U1, V1, M1, T1];%rho1, U1, V1, M1, T1


%	Loop over the n points to determine the base flow state 
	options = optimset('Display', 'off');%This turns off the convergence statement
	for i = 1:N
%		Doing the function-within-a-function thing to pass values to turbineEqs
		turbineEqnsCaller = @(q) turbineEqns(q, stagnation_state, inlet_state, A(i), theta_bar(i));
		ess = s(i);
		q = fsolve(turbineEqnsCaller, q, options);
		rho_bar(i) 		= q(1);	
		u_bar(i)  		= q(2);
		v_bar(i)		= q(3);
		w_bar(i) 		= sqrt(q(2).^2 + q(3).^2);
		T_bar(i) 		= q(5);
		p_bar(i) 		= rho_bar(i)*R*T_bar(i);
		sos_bar(i)		= sqrt(gamma*p_bar(i)/rho_bar(i));
		M_bar(i)		= w_bar(i)./sos_bar(i);
		psi_bar(i)		= 0;
	end

	baseFlow = BaseFlow(s, x_c, y_c, gamma, A, theta_bar, p_bar, rho_bar, T_bar, u_bar, v_bar, w_bar, M_bar, sos_bar, psi_bar);	
end
