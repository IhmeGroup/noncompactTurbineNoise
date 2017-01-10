function[phidot] = nonCompactODE(s, phi)
	global gamma;
	global gm1; 
	global gp1;
	global gm1o2;
	global gp1o2;

	global omega;
	global areaSp;
	global thetSp;
	global xvelSp;
	global yvelSp;
	global machSp;
	global presSp;
	global densSp;
	global tempSp;

	global sos_Sp;
	global ddensdxSp;
	global ddensdySp;
	global dxveldxSp;
	global dxveldySp;
	global dyveldxSp;
	global dyveldySp;
	

	A_t = zeros(4);
	B_x = zeros(4);
	B_y = zeros(4);
	C 	= zeros(4);

	A_t = 2*pi*sqrt(-1)*omega*[	1 0 0 0; 
								0 1 0 0; 
								0 0 1 0; 
								0 0 0 1];

	a	= ppval(sos_Sp, s);
	rho = ppval(densSp, s);
	u 	= ppval(xvelSp, s);
	v 	= ppval(yvelSp, s);
	theta = ppval(thetSp, s);
	drhoudx = ppval(ddensdxSp, s)*u + ppval(dxveldxSp, s)*rho;
	drhovdy = ppval(ddensdySp, s)*v + ppval(dyveldySp, s)*rho;
	dudx = ppval(dxveldxSp, s);
	dudy = ppval(dxveldySp, s);
	dvdx = ppval(dyveldxSp, s);
	dvdy = ppval(dyveldySp, s);
	uoverv = u/(v + 1E-3);
	voveru = v/(u + 1E-3);


	B_x = [	u 		u 	0	0	;
			a^2/u	u	0	0	;
			0		0	u	0	;
			0		0	0	u];

	B_y = [ v		0 	v	0	;
			0		v	0	0	;
			a^2/v	0	v	0	;
			0		0	0	v	];

	C	= [				0					1/rho*drhoudx			1/rho*drhovdy				  0			;
			-gm1*(dudx + voveru*dudy)	2*dudx + voveru*dudy		voveru*dudy			-dudx - voveru*dudy	;
			-gm1*(uoverv*dvdx + dvdy)		uoverv*dvdx			uoverv*dvdx + 2*dvdy	-uoverv*dvdx - dvdy	;
						0						  0						  0						  0			];

	Astar = -(A_t + C);
	Bstar = B_x*cos(theta) + B_y*sin(theta);

	L = Bstar\Astar;
	phidot = L*phi;

end
