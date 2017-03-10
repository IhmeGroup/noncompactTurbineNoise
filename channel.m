function[extrema, x_c, y_c, ytop, ybot,stagnation_state] = channel()
	N = 10;
	rho0 = 1.2754;%Stagnation Condition
	T0 = 293;%Stagnation Condition
	M1 = 0.3;%Inlet Mach number
	R = 287.058;%Universal gas constant
	gamma = 1.4;%c_p/c_v
	delta = 1;
	L = 5;

	x_bot = linspace(0, L, N)';
	x_top = x_bot;
	y_bot = zeros(N,1);
	y_top = y_bot + delta;
	ytop = spline(x_top, y_top);
	ybot = spline(x_bot, y_bot);
	x_tl = 0; y_tl = delta; x_tr = L; y_tr = delta;
	x_bl = 0; y_bl = 0; 	x_br = L; y_br = 0;
	x_c = x_top;
	y_c = 0.5*(y_top + y_bot);
	exit_slope = 0;

	extrema = [x_tl; y_tl; x_tr; y_tr; x_bl; y_bl; x_br; y_br; exit_slope];
	stagnation_state = [gamma; R; rho0; T0; M1];
end
