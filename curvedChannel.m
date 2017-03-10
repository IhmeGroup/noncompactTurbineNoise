function[extrema, x_c, y_c, ytop, ybot,stagnation_state] = curvedChannel()
	N = 100;
	phi_a = 1.01*pi;
	phi_b = 2*pi;
	R_outer = 2;%Outer radius (shoudln't matter)
	R_inner = 1;%Inner radius (shouldn't matter)
	rho0 = 1.2754;%Stagnation Condition
	T0 = 293;%Stagnation Condition
	M1 = 0.3;%Inlet Mach number
	R = 287.058;%Universal gas constant
	gamma = 1.4;%c_p/c_v
	
	phi = linspace(phi_a, phi_b,N);
	x_top = R_inner*cos(phi);
	y_top = R_inner*sin(phi);
	x_bot = R_outer*cos(phi);
	y_bot = R_outer*sin(phi);
	x_c	 = (R_inner + R_outer)/2*cos(phi);	
	y_c	 = (R_inner + R_outer)/2*sin(phi);
	ytop = spline(x_top, y_top);
	ybot = spline(x_bot, y_bot);
	x_tl = -R_outer - 1; y_tl = R_outer; x_tr = R_outer + 1; y_tr = R_outer;
	x_bl = -R_outer - 1; y_bl = R_outer; x_br = R_outer + 1; y_br = R_outer;
	exit_slope = 3;

	extrema = [x_tl; y_tl; x_tr; y_tr; x_bl; y_bl; x_br; y_br; exit_slope];
	stagnation_state = [gamma; R; rho0; T0; M1];
end
