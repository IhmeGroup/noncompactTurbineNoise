function[areaSp, thetSp, xvelSp, yvelSp, machSp, presSp, densSp, tempSp, sos_Sp, ddensdxSp, ddensdySp, dxveldxSp, dxveldySp, dyveldxSp, dyveldySp] = buildMeanFlow(run, suppress)
	close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Set flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (suppress)
		plot_geom = false;
		plot_base = false;
		plot_checks = false;
	else
		plot_geom = false;
		plot_base = true;
		plot_checks = true;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Independent Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	A really ugly way of passing these parameters through to the mean flow solver, TurbineEqs
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
	global ddensdxSp;
	global ddensdySp;
	global dxveldxSp;
	global dxveldySp;
	global dyveldxSp;
	global dyveldySp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Independent Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (run == 1)%Linear nozzle
		N = 11;
		rho0 = 1.2754;%Stagnation Condition
		T0 = 293;%Stagnation Condition
		M1 = 0.3;%Inlet Mach number
		R = 287.058;%Universal gas constant
		gamma = 1.4;%c_p/c_v
	elseif (run == 2)%Curved channel
		N = 101;
		phi_a = 3*pi/2;%Entry angle
		phi_b = pi + 17*pi/18;%exit angle
		R_outer = 2;%Outer radius (shoudln't matter)
		R_inner = 1;%Inner radius (shouldn't matter)
		rho0 = 1.2754;%Stagnation Condition
		T0 = 293;%Stagnation Condition
		M1 = 0.3;%Inlet Mach number
		R = 287.058;%Universal gas constant
		gamma = 1.4;%c_p/c_v
	elseif (run == 3)%Single stator system
		N = 201;%# of points used to build interpolants
		delta = 40;%cord-cord separation
		rho0 = 1.2754;%Stagnation Condition
		T0 = 293;%Stagnation Condition
		M1 = 0.218;%Cascade entrance speed
		R = 287.058;%Universal gas constant
		gamma = 1.4;%c_p/c_v
	elseif (run == 4)%High-fidelity comparison
		N = 201;
		rho0 = 1.2754;%Stagnation Condition
		T0 = 293;%Stagnation Condition
		M1 = 0.3;%Inlet Mach number
		R = 287.058;%Universal gas constant
		gamma = 1.4;%c_p/c_v
	else
		error('Unrecognized run # detected');
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Dependent Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	gm1 = gamma -1;
	gp1 = gamma + 1;
	gm1o2 = gm1/2;
	gp1o2 = gp1/2;
	rho1 = rho0*(1 + gm1o2*M1^2)^(-1/gm1);
	T1 = T0*(1 + gm1o2*M1^2)^(-1);
	a1 = sqrt(gamma*R*T1);
	U1 = M1*a1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Load and interpret the blade geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (run == 1)%Linear nozzle

	elseif (run == 2)%Curved channel
%		Generate the channel geometry
		phi = linspace(phi_a, phi_b);
		xtop = R_inner*cos(phi);
		ytop = R_inner*sin(phi);
		xbot = R_outer*cos(phi);
		ybot = R_outer*sin(phi);
		[xmin, mindex] = min((xtop + xbot)/2);
		[xmax, maxdex] = max((xtop + xbot)/2);
	elseif (run == 3)%Real rotor
%		Load the stator geometry
		data = load('stator.txt');
		data = data(1:120,1:2); % Clip to insure that only one pass over the geometry is retained
%		Identify the extreme points to split the geometry into a "top" and "bottom"
		[xmax, maxdex] = max(data(:,1));
		[xmin, mindex] = min(data(:,1));	
		xtop = data(mindex:maxdex,1);
		ytop = data(mindex:maxdex,2);
%		B/c the splining flips out over the repeated point, remove the last point in the series
		xbot = [data(1:mindex,1); data(maxdex:end-4,1)];
		ybot = [data(1:mindex,2); data(maxdex:end-4,2)];
		ytopsp = spline(xtop, ytop);
		ybotsp = spline(xbot, ybot);

	elseif (run == 4)
	else
		error('Unrecognized run # detected');
	end

%	Build splines for the top and bottom section to allow an approximation of the chord
	top = spline(xtop, ytop);
	bot = spline(xbot, ybot);

	x = linspace(xmin, xmax, N);
	y1 = ppval(top, x);
	y2 = ppval(bot, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Establishing the arc and cord length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Parameterize the cord as c(x) = y
	c = (y1 + y2)/2;

%	Second, compute the local spatial gradient of the cord
	dydx = zeros(N,1);
	dydx(1) = (c(2) - c(1))/(x(2) - x(1));
	for i = 2:N-1
		dydx(i) = (c(i+1) - c(i-1))	/(x(i+1) - x(i-1));
	end
	dydx(N) = (c(N) - c(N-1))/(x(N) - x(N-1));

%	Then compute the total cord length which will be used for normalization
	disp('cord length')
	L_c = trapz(x, sqrt(1 + dydx.^2))

%	Next, compute the instantaneous cord length to be used as a coordinate for the splines
	s = zeros(N,1);
	for i = 2:N
		s(i) = trapz(x(1:i), sqrt(1 + dydx(1:i).^2));
	end

%	Lastly, normalize this coordinate by the total cord length to paramaterize the arc
	s = s./L_c;
	x_c = spline(s, x);
	y_c = spline(s, c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Parameterize geometry by the arc length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	First compute the thickness (Area) and direction (theta) values at each point 
%	t = abs(y2 - y1);
%	A = delta - t;
%	theta = atan2(dydx, 1);
%	IN THE FLOW PARALLEL DIRECTION!!!
	theta = atan2(dydx, 1);
	ex_top = zeros(length(s), 1);
	why_top = zeros(length(s), 1);
	ex_bot = zeros(length(s), 1);
	why_bot = zeros(length(s), 1);
%	for i = 1:length(s)
	for i = 101:101
		ex =  ppval(x_c, s(i));
		why = ppval(y_c, s(i));
		ytan = (why - (x - ex).*1./dydx(i));
		topsect = ytan - y2;
		botsect = ytan - y1;
		topsecter = spline(x, topsect);		
		botsecter = spline(x, botsect);
		toptop = ppval(topsecter, x);
		botbot = ppval(botsecter, x);
		toppy = @(x) ppval(topsecter, x);
		botty = @(x) ppval(botsecter, x);
		mean(x)
		[X1, Y1] = fzero(toppy, mean(x))
		[X2, Y2] = fzero(botty, mean(x))
%		ytan2 = (why - (x - ex).*1./dydx(i));
	end
		
	figure();
	plot(ex, why, 'rp');
	hold on;
	ex = ppval(x_c, s);
	why = ppval(y_c, s);
	plot(ex, why, 'k-');
	plot(x, ytan, 'b-p');
	plot(xtop, ytop, 'c--');
	plot(xbot, ybot, 'c--');
	plot(ex(i), why(i), 'rp');
	plot(x, c, 'm-');
	plot(X1, ppval(ybotsp, X1), 'rp');
	plot(X2, ppval(ytopsp, X2), 'gp');
	axis equal

	figure();
	plot(x, topsect,'k');
	hold on;
	plot(x, botsect,'b');
	plot(x, toptop, 'k--');
	plot(x, botbot, 'b--');
%	plot(ex_bot, why_bot, 'b-p');
		
%		X = [ex - 1; ex; ex + 1];
%		Y = [why - dydx; why; why + dydx];
%		normal = polyfit(X, Y, 1);
%		x
%		[ex_top(i), why_top(i)] = curveintersect(top, normal);
%		[ex_bot(i), why_bot(i)] = curveintersect(bot, normal);
%		[ex_top(i), why_top(i)] = intersections(xtop, ytop, X, Y, 1);
%		[ex_bot(i), why_bot(i)] = intersections(xbot, ybot, X, Y, 1);
%		A(i) = 
%	end


%	t = 

%	Next, compute splines parameterizing these quantities in terms of the arc length s
	areaSp 	= spline(s, A);
	thetSp = spline(s, theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Compute the base flow and then spline it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	1 Now that the geometry has been parameterized, replace s by a uniform mesh
	s = linspace(0,1,N);

%	2 While it's a bit clunky, define global functions to return A and theta from the splines
	global returnA;
	global returnTheta;
	returnA = @(s) ppval(areaSp,s);
	returnTheta = @(s) ppval(thetSp,s);
	
%	3 Allocate memory for the solution for efficiency
	rho = zeros(N,1);
	U = zeros(N,1);
	V = zeros(N,1);
	M = zeros(N,1);
	T = zeros(N,1);
	p = zeros(N,1);
	A = zeros(N,1);
	theta = zeros(N,1);
	sos = zeros(N,1);

%	4 Initialize the solver with the inflow state (This is not a very good guess)
	q = [rho1, U1, 0, M1, T1];%rho1, U1, V1, M1, T1

%	5 Loop over the n points to determine the base flow state 
	options = optimset('Display', 'off');%This turns off the convergence statement
	for i = 1:N
		ess = s(i);
		q = fsolve(@TurbineEqs, q, options);
		rho(i) 		= q(1);	
		U(i)  		= q(2);
		V(i)		= q(3);
		M(i)		= q(4);
		T(i) 		= q(5);
		p(i) 		= rho(i)*R*T(i);
		A(i) 		= returnA(ess);
		sos(i)		= sqrt(gamma*p(i)/rho(i));
		theta(i) 	= returnTheta(ess);
	end

%	6 Next, do the same thing for spatial gradient s
	drhodx = zeros(N,1);
	drhody = zeros(N,1);
	dudx = zeros(N,1);
	dudy = zeros(N,1);
	dvdx = zeros(N,1);
	dvdy = zeros(N,1);

	drhodx(1) 	= (rho(2) - rho(1))/(x(2) - x(1));
	drhody(1) 	= (rho(2) - rho(1))/(c(2) - c(1));
	dudx(1) 	= (U(2) - U(1))/(x(2) - x(1));
	dudy(1)		= (U(2) - U(1))/(c(2) - c(1));
	dvdx(1) 	= (V(2) - V(1))/(x(2) - x(1));
	dvdy(1)		= (V(2) - V(1))/(c(2) - c(1));

	for i = 2:N-1
		drhodx(i) 	= (rho(i+1) - rho(i-1))/(x(i+1) - x(i-1));
		drhody(i) 	= (rho(i+1) - rho(i-1))/(c(i+1) - c(i-1));
		dudx(i)		= (U(i+1) - U(i-1))/(x(i+1) - x(i-1));
		dudy(i)		= (U(i+1) - U(i-1))/(c(i+1) - c(i-1));
		dvdx(i)		= (V(i+1) - V(i-1))/(x(i+1) - x(i-1));
		dvdy(i)		= (V(i+1) - V(i-1))/(c(i+1) - c(i-1));
	end

	drhodx(N) 	= (rho(N) - rho(N-1))/(x(N) - x(N-1));
	drhody(N) 	= (rho(N) - rho(N-1))/(c(N) - c(N-1));
	dudx(N) 	= (U(N) - U(N-1))/(x(N) - x(N-1));
	dudy(N) 	= (U(N) - U(N-1))/(c(N) - c(N-1));
	dvdx(N) 	= (V(N) - V(N-1))/(x(N) - x(N-1));
	dvdy(N) 	= (V(N) - V(N-1))/(c(N) - c(N-1));
		

%	Lastly, spline everything so that the ODE solver is (relatively) fast
	densSp 	= spline(s, rho);
	xvelSp 	= spline(s, U);
	yvelSp	= spline(s, V);
	machSp	= spline(s, M);
	tempSp	= spline(s, T);
	presSp	= spline(s, p);
	areaSp	= spline(s, A);
	thetSp	= spline(s, theta);
	sos_Sp	= spline(s, sos);
	ddensdxSp	= spline(s, drhodx);
	ddensdySp	= spline(s, drhody);
	dxveldxSp	= spline(s, dudx);
	dxveldySp	= spline(s, dudy);
	dyveldxSp	= spline(s, dvdx);
	dyveldySp	= spline(s, dvdy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Plot the stator geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (plot_geom);
		h = figure();
%		First plot the raw data
		plot(data(:,1), data(:,2), 'k-o', 'LineWidth', 2);
		hold on;
%		Next, plot the inputs to the top and bottom splines to insure they're ok
		plot(xtop, ytop, 'r>', 'LineWidth', 2);
		plot(xbot, ybot, 'b>', 'LineWidth', 2);
%		Then plot the outputs from the splines to insure they're also reasonable
		plot(x, y1, 'r-', 'LineWidth', 2);
		plot(x, y2, 'b-', 'LineWidth', 2);
%		Plot the cord computed from the splines
		plot(x, c, 'g-', 'LineWidth', 2);
		plot(x, t, 'm-', 'LineWidth', 2);
		plot(x, -s, 'c-', 'LineWidth', 2);
	end%(plot_geom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Plot the base flow as a check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (plot_base)
		g = figure();
		set(g, 'Position', [0 0 1000 800]);

%		Plot hydrodynamic quantities
		subplot(3,3,1);
		plot(s, U, 'LineWidth', 2);
		xlabel('s');
		ylabel('$U(x)$', 'Interpreter', 'LaTeX');
		title('X Velocity');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		subplot(3,3,2);
		plot(s, V, 'LIneWidth', 2);
		xlabel('s');
		ylabel('$V(x)$', 'Interpreter', 'LaTeX');
		title('Y Velocity');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		subplot(3,3,3);
		plot(s, M, 'LineWidth', 2);
		xlabel('s');
		ylabel('$M(x)$', 'Interpreter', 'LaTeX');
		title('Mach #');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

%		Plot thermodynamic quantities
		subplot(3,3,4);
		plot(s, p, 'LineWidth', 2);
		xlabel('s');
		ylabel('$p(x)$', 'Interpreter', 'LaTeX');
		title('Pressure');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		subplot(3,3,5);
		plot(s, rho, 'LineWidth', 2);
		xlabel('s');
		ylabel('$\rho(x)$', 'Interpreter', 'LaTeX');
		title('Density');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		subplot(3,3,6);
		plot(s, T, 'LineWidth', 2);
		xlabel('s');
		ylabel('$T(x)$', 'Interpreter', 'LaTeX');
		title('Temperature');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

%		Plot geometric quantities
		subplot(3,3,7);
		plot(s, A, 'LineWidth', 2);
		xlabel('s');
		ylabel('$A(x)$', 'Interpreter', 'LaTeX');
		title('Area');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		subplot(3,3,8);
		plot(s, theta, 'LineWidth', 2);
		xlabel('x');
		ylabel('$\theta(x)$', 'Interpreter', 'LaTeX');
		title('Theta');
		set(gca, 'FontSize', 14, 'FontName', 'Times');
	end%(plot_base);

	if (plot_checks)
		figure();
		mdot = rho.*U.*A;
		plot(x, mdot, 'lineWidth', 2);
		xlabel('x');
		ylabel('mdot');
		title('Mass Flow Rate');
		set(gca, 'FontSize', 14, 'FontName', 'Times');

		figure();
		theta = atan2(V,U);
		plot(x, theta, 'lineWidth', 2);
		xlabel('x');
		ylabel('theta');
		title('Flow Direction');
		set(gca, 'FontSize', 14, 'FontName', 'Times');
	end%(plot_checks)
end
