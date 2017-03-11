function[baseFlow] = leyko(M_1, M_2, theta_1, theta_2)
	close all;

	N = 101;
	rho0 = 1.2754;%Stagnation Condition
	T0 = 293;%Stagnation Condition
	R = 287.058;%Universal gas constant
	p0 = rho0*R*T0;
	gamma = 1.4;%c_p/c_v
	delta = 1;
	L = 1;

%	Some parameters for convenience
	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;

	R_outer = 2;%Outer radius (shoudln't matter)
	R_inner = 1;%Inner radius (shouldn't matter)

	A_l = gp1o2.^(-gp1o2/gm1).*(1+gm1o2*M_1*M_1)^(gp1o2/gm1)/M_1;
	A_r = gp1o2.^(-gp1o2/gm1).*(1+gm1o2*M_2*M_2)^(gp1o2/gm1)/M_2;

	theta 	= linspace(theta_1, theta_2, N);
	M		= linspace(M_1, M_2, N);
	s		= linspace(0, 1, N);

	A = zeros(N,1);
	x_bot = zeros(N,1);
	x_c = zeros(N,1);
	x_top = zeros(N,1);
	y_bot = zeros(N,1);
	y_c = zeros(N,1);
	y_top = zeros(N,1);
	rho_bar = zeros(N,1);
	p_bar = zeros(N,1);
	T_bar = zeros(N,1);
	w_bar = zeros(N,1);
	u_bar = zeros(N,1);
	v_bar = zeros(N,1);
	M_bar = zeros(N,1);
	sos_bar = zeros(N,1);
	psi_bar = zeros(N,1);


	for i = 1:N
		A(i) 		= gp1o2.^(-gp1o2/gm1).*(1+gm1o2*M(i)*M(i))^(gp1o2/gm1)/M(i);
		r_inner 	= 1;
		r_center 	= r_inner + A(i)/2;
		r_outer		= r_inner + A(i);
		x_bot(i) 	= r_inner*cosd(theta(i));
		x_c(i)		= r_center*cosd(theta(i));
		x_top(i) 	= r_outer*cosd(theta(i));
		y_bot(i)	= r_inner*sind(theta(i));
		y_c(i)		= r_center*sind(theta(i));
		y_top(i)	= r_outer*sind(theta(i));
		rho_bar(i)	= rho0*(1+gm1o2*M(i)^2).^(-1/gm1);
		p_bar(i)	= p0*(1+gm1o2*M(i).^2).^(-gamma/gm1);
		T_bar(i)	= T0*(1+gm1o2*M(i).^2).^-1;
		w_bar(i)	= M(i)*sqrt(gamma*R*T_bar(i));
		u_bar(i)	= w_bar(i)*cosd(theta(i));
		v_bar(i)	= w_bar(i)*sind(theta(i));
		M_bar(i)	= M(i);
		sos_bar(i)	= sqrt(gamma*R*T_bar(i));
		psi_bar(i)	= 0;
	end


	x_tl = (1 + A_l/2)*cosd(theta_1);
	x_tr = (1 + A_r/2)*cosd(theta_2);
	y_tl = (1 + A_l/2)*sind(theta_1);
	y_tr = (1 + A_r/2)*sind(theta_2);
	x_bl = 1;
	x_br = 1;
	y_bl = 1;
	y_br = 1;

	theta = theta*2*pi/360;

	baseFlow = BaseFlow(s, x_c, y_c, gamma, A, theta, p_bar, rho_bar, T_bar, u_bar, v_bar, w_bar, M_bar, sos_bar, psi_bar);	
end

