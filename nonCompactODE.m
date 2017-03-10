function[dphi_dt] = nonCompactODE(t, phi, splines, Omega)
	rho 		= ppval(splines.rho_bar.raw, t);
	drho_dn		= ppval(splines.rho_bar.ddn, t);
	p			= ppval(splines.p_bar.raw, t);
	dp_dn		= ppval(splines.p_bar.ddn, t);
	w 			= ppval(splines.w_bar.raw, t);
	dw_dn		= ppval(splines.w_bar.ddn, t);
	Psi			= ppval(splines.psi_bar.raw, t);
	theta		= ppval(splines.theta_bar.raw, t);
	dtheta_dn	= ppval(splines.theta_bar.ddn, t);
	gamma		= splines.gamma;
	M			= ppval(splines.M_bar.raw, t);


	gm1 = gamma - 1;
	gp1 = gamma + 1;
	gm1og = gm1/gamma;
	gp1og = gp1/gamma;
	oog = 1/gamma;

	A = zeros(5,5);

	B = zeros(5,5);

	C = zeros(5,5);

	R = zeros(5,5);

%	1 - Continuity
%	2 - w' - transport
%	3 - \theta' - transport
%	4 - entropy transport
%	5 - mixture fraction transport

%					1						2			3			   4			5
%	q = [\frac{p'}{\gamma \bar{p}} \frac{w'}{\bar{w}} \theta' \frac{s'}{\bar{c}_p} Y_i']

%			1			2		3			4		5
	R = [	w			w		0			0			0;		%1
			w/M^2		w		0			0			0;		%2
			0			0		w			0			0;		%3
			0			0		0			w			0;		%4
			0			0		0			0			w];		%5

%									1		2	3		4	5
	A = 2*pi*sqrt(-1)*Omega*	[	-1 			0 		0 		0		0	;	%1
		 							0 			-1		0		0		0	;	%2
		 							0 			0		-1		0		0	;	%3
		 							0 			0		0		-1		0	;	%4
		 							0 			0		0		0		-1	];	%5

%		 			1			2				3					4				5
	B = [	0				0					w*dtheta_dn			0				0;					%1	
			gm1*dw_dn		-2*dw_dn			w*dtheta_dn			dw_dn			Psi*dw_dn;			%2
			0				-2*w*dtheta_dn 		-dw_dn				0				0;					%3
			0				0					0					0				0;					%4
			0				0					0					0				0];					%5


	dphi_dt = (R\(A + B))*phi;

end
