function[extrema, x_c, y_c, ytop, ybot,stagnation_state] = nozzle()
	N = 101;
	rho0 = 1.2754;%Stagnation Condition
	T0 = 293;%Stagnation Condition
	M1 = 0.29;%Inlet Mach number
	M2 = 0.9;%Outlet Mach number
	R = 287.058;%Universal gas constant
	gamma = 1.4;%c_p/c_v
	delta = 1;
	L = 5;

%	Some parameters for convenience
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;

	x_l = sqrt(gp1o2*M1*M1/(1+gm1o2*M1*M1));
	x_r = sqrt(gp1o2*M2*M2/(1+gm1o2*M2*M2));
	x_bot = linspace(x_l, x_r, N)';
	x_top = x_bot;
	y_bot = zeros(N,1);
	y_top = zeros(N,1);
	for i = 1:N
		x2 = x_bot(i).^2;
		em = sqrt(2/gp1*x2/(1-gm1/gp1*x2));
		A = gp1o2.^(-gp1o2/gm1)*(1+gm1o2*em*em)^(gp1o2/gm1)/em;
		L = A;
		y_bot(i) = -L/2;
		y_top(i) = L/2;
	end
	x_top = x_top - x_l;
	x_top = x_top./max(x_top);
	x_bot = x_top;
	x_tl = x_top(1); y_tl = y_top(1); x_tr = x_top(end); y_tr = y_top(end);
	x_bl = x_bot(1); y_bl = y_bot(1); x_br = x_bot(end); y_br = y_bot(end);
	x_c = x_top;
	y_c = 0.5*(y_top + y_bot);
	exit_slope = 0;
	rslope = (y_top(end) - y_top(end-1))/(x_top(end) - x_top(end-1));
	lslope = (y_top(2) - y_top(1))/(x_top(2) - x_top(1));
	ytop = spline(x_top, [lslope; y_top; rslope]);
	ybot = spline(x_bot, [-lslope; y_bot; -rslope]);

	extrema = [x_tl; y_tl; x_tr; y_tr; x_bl; y_bl; x_br; y_br; exit_slope];
	stagnation_state = [gamma; R; rho0; T0; M1];
end
