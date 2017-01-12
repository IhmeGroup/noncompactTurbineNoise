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
		N = 10;
		rho0 = 1.2754;%Stagnation Condition
		T0 = 293;%Stagnation Condition
		M1 = 0.3;%Inlet Mach number
		R = 287.058;%Universal gas constant
		gamma = 1.4;%c_p/c_v
		delta = 1;
	elseif (run == 2)%Curved channel
		N = 100;
		phi_a = 1*pi;
		phi_b = 2*pi;
%		phi_a = 3*pi/2 + pi/18;%Entry angle
%		phi_b = pi + 17*pi/18;%exit angle
		R_outer = 2;%Outer radius (shoudln't matter)
		R_inner = 1;%Inner radius (shouldn't matter)
		rho0 = 1.2754;%Stagnation Condition
		T0 = 293;%Stagnation Condition
		M1 = 0.3;%Inlet Mach number
		R = 287.058;%Universal gas constant
		gamma = 1.4;%c_p/c_v
	elseif (run == 3)%Single stator system
		N = 500;%# of points used to build interpolants
		delta = 60;%cord-cord separation
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
		phi = linspace(phi_a, phi_b,N);
		x_top = R_inner*cos(phi);
		y_top = R_inner*sin(phi);
		x_bot = R_outer*cos(phi);
		y_bot = R_outer*sin(phi);
		x_c	 = (R_inner + R_outer)/2*cos(phi);
		y_c	 = (R_inner + R_outer)/2*sin(phi);
		ytop = spline(x_top, y_top);
		ybot = spline(x_bot, y_bot);
	elseif (run == 3)%Real rotor
%		Load the stator geometry
		data = load('stator.txt');
		data = data(1:120,1:2); % Clip to insure that only one pass over the geometry is retained
%		Identify the extreme points and split the geometry into a "top" and "bottom"
		[xmax, maxdex] = max(data(:,1));
		[xmin, mindex] = min(data(:,1));	
		x_top = data(mindex:maxdex,1);
		y_top = data(mindex:maxdex,2);
%		B/c the splining flips out over the repeated point, remove the last point in the series
		x_bot = [data(1:mindex,1); data(maxdex:end-4,1)];
		y_bot = [data(1:mindex,2); data(maxdex:end-4,2)];
%		Flip the geometry b/c the points are arranged clockwise
		y_bot = flipud(y_bot);
		x_bot = flipud(x_bot);
%		Offset the bottom layer to focus on the void instead of the blade
		y_bot = y_bot - delta;
%		Grab the corners of the top and bottom geometry
		x_tl = x_top(1); y_tl = y_top(1);
		x_tr = x_top(end); y_tr = y_top(end);
		x_bl = x_bot(1); y_bl = y_bot(1);
		x_br = x_bot(end); y_br = y_bot(end);
%		Compute the slope of the blade tail to determine the angle at which the flow leaves the cascade
		disp('slopes');
		[	x_top(end-3:end),	y_top(end-3:end), 	x_bot(end-3:end), 	y_bot(end-3:end)]
%		m_1 = (y_top(end) - y_top(end-1))/(x_top(end) - x_top(end-1))
		exit_slope = (y_bot(end) - y_bot(end-1))/(x_bot(end) - x_bot(end-1))
%		exit_slope = 0.5*(m_1 + m_2);
%		Spline the geometry so that it can be interpolated onto the same grid in terms of x or s or some coordinate
%			Specify the derivative of the spline at the endpoints so that it smoothly transitions into the thing we want it to
		xtop = spline(x_top, x_top);
		ytop = spline(x_top, [0; y_top; exit_slope]);
		xbot = spline(x_bot, x_bot);
		ybot = spline(x_bot, [0; y_bot; exit_slope]);
%		And evaluate the spline on a constant grid, first just over the cord
		x_c = linspace(xmin, xmax, N);
		x_top = ppval(xtop, x_c);
		y_top = ppval(ytop, x_c);
		x_bot = ppval(xbot, x_c);
		y_bot = ppval(ybot, x_c);
%		The cord is just the average of the top and bottom contours
		y_c = 0.5*(y_top + y_bot);

		

		figure();
		plot(x_c, y_top, 'b-', 'LineWidth', 2);
		hold on;
		plot(x_c, y_bot, 'm-', 'LineWidth', 2);
		plot(x_c(1), 		y_c(1), 	'bp');
		plot(x_c(end), 	y_c(end), 'bo');
		plot(x_c(1), 		y_c(1), 	'mp');
		plot(x_c(end), 	y_c(end), 'mo');
		plot(x_c, y_c, 'k-', 'LineWidth', 2)


%		x_left = linspace(x_top(1) - 30, x_top(1) - 0.001, 51)';
%		y_left = y_top(1)*ones(51,1);
%		x_right = linspace(x_top(end) + 0.001, x_top(end) + 30, 51)';
%		dydx_right = (y_bot(end) - y_bot(end-1))/(x_bot(end) - x_bot(end-1));
%		y_right = y_top(end) + (x_right - x_top(end))*dydx_right;
%		x_top = [x_left; x_top; x_right];
%		y_top = [y_left; y_top; y_right];
%		x_left = linspace(x_bot(1) - 30, x_bot(1) - 0.001, 51)';
%		y_left = y_bot(1)*ones(51,1);
%		x_right = linspace(x_bot(end) + 0.001, x_bot(end) + 30, 51)';
%		y_right = y_bot(end) + (x_right - x_bot(end))*dydx_right;
%		x_bot = [x_left; x_bot; x_right];
%		y_bot = [y_left; y_bot; y_right];
%		xmin = min(x_bot);
%		xmax = max(x_bot);
%		xtop = spline(x_top, x_top);
%		ytop = spline(x_top, y_top);
%		x_top = ppval(xtop, x_c);
%		y_top = ppval(ytop, x_c);
%		xbot = spline(x_bot, x_bot);
%		ybot = spline(x_bot, y_bot);
%		x_bot = ppval(xbot, x_c);
%		y_bot = ppval(ybot, x_c);
		y_c = (y_bot + y_top)/2;
	elseif (run == 4)
	else
		error('Unrecognized run # detected');
	end

	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Establishing the arc and cord length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Compute the local normal vector of the cord
%	N = length(x_c);
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
%	Compute the X-sectional area of the geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	To compute the cord location, we need to determine where the normal to a point in the cord intersects the top and bottom boundaries

%	First, allocate memory to store these intersection locations
	xsect_top = zeros(N,1);
	ysect_top = zeros(N,1);
	xsect_bot = zeros(N,1);
	ysect_bot = zeros(N,1);

%	Loop over each point on the cord
	global param
	b = [N/2:N, N/2-1:-1:1];
	for a = 1:N
		i = b(a)
%		Next, generate fns ftop and fbot which are equivalent to the 
%		difference between the geometry contour and the cord-normal curve passing through (x_loc, y_loc)
		if (abs(dydx(i)) <= 1E-6)
			xsect_top(i) = x_c(i);
			xsect_bot(i) = x_c(i);
			ysect_top(i) = ppval(ytop, xsect_top(i));
			ysect_bot(i) = ppval(ybot, xsect_bot(i));
		else
%			f_top = y_top - (y_c(i) + (x_top - x_c(i))*(1E-6 + -1/(dydx(i))));
%			f_bot = y_bot - (y_c(i) + (x_bot - x_c(i))*(1E-6 + -1/(dydx(i))));

%			Then run a spline through them
%			ftop = spline(x_top, f_top);
%			fbot = spline(x_bot, f_bot);

%			To help identify their zeros, define a function which just evaluates the splines
%			ftopf = @(x_top) ppval(ftop, x_top);
%			fbotf = @(x_bot) ppval(fbot, x_bot);

%			And lastly, use fzero to compute where these functions intersect the x axis
%			(i.e. where the normal line intersects the geometry)
%			if (i ~= N)
%				[xsect_top(i), ~, flag] = fzero(ftopf, xsect_top(i+1));
%			else
			m = -1/dydx(i);
			param = {x_tl, y_tl, x_tr, y_tr, m, x_c(i), y_c(i), ytop, dydx(i)};
			[xsect_top(i), ~, flag] = fzero(@ftopf, mean(x_c));
			ysect_top(i) = ppval(ytop, xsect_top(i));

			param = {x_bl, y_bl, x_br, y_br, m, x_c(i), y_c(i), ybot, dydx(i)};
			[xsect_bot(i), ~, flag] = fzero(@ftopf, mean(x_c));
			ysect_bot(i) = ppval(ybot, xsect_bot(i));
			
%			if (i ~= 1)
%				[xsect_bot(i), ~, flag] = fzero(fbotf, xsect_bot(i-1));
%			else
%				[xsect_bot(i), ~, flag] = fzero(fbotf, mean(x_c));
%			end
		end
%		if (flag  ~= 1) 
%			xsect_bot(i) = 0; 
%			ysect_bot(i) = 0;
%		else
			ysect_bot(i) = ppval(ybot, xsect_bot(i));
%		end
	end

%	Plot a bunch of thinks to insure the code is doing ok
	figure();
	plot(x_c, y_c, 'm-', 'LineWidth', 2);
	x_c(1:3)
	y_c(1:3)
	hold on;
	plot(x_top, y_top, 'k-o', 'LineWidth', 2);
	plot(x_bot, y_bot, 'k-o', 'LineWidth', 2);
	plot(x_top(1), y_top(1), 'bo')
	plot(x_bot(1), y_bot(1), 'ro')
	plot(x_top(N), y_top(N), 'b+')
	plot(x_bot(N), y_bot(N), 'r+')
	if (run == 3)
		plot(x_bot, y_bot + delta, 'k-+', 'LineWidth', 2);
		plot(x_top, y_top - delta, 'k-+', 'LineWidth', 2);
	end
	for i = 1:N
%		plot(x_loc(i), y_loc(i), 'rp');
		plot(xsect_top(i), ysect_top(i), 'bp');
		plot(xsect_bot(i), ysect_bot(i), 'mp');
		plot([xsect_top(i), xsect_bot(i)], [ysect_top(i), ysect_bot(i)], 'b-');
	end
	axis equal;
	
%	Next, compute splines parameterizing these quantities in terms of the arc length s
	areaSp 	= spline(s, A);
	thetSp = spline(s, theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Compute the base flow and then spline it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
	dxds_1 = (x_c(2) - x_c(1))/(s(2) - s(1))
	dyds_1 = (y_c(2) - y_c(1))/(s(2) - s(1))
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


function[err] = ftopf(x);%, ftop, xtl, ytl, xtr, ytr, m, x_c, y_c)
%need to use global structure to pass params in
	global param;
	x_tl 	= cell2mat(param(1));
	y_tl 	= cell2mat(param(2));
	x_tr 	= cell2mat(param(3));
	y_tr 	= cell2mat(param(4));
	m 		= cell2mat(param(5));
	x_c 	= cell2mat(param(6));
	y_c 	= cell2mat(param(7));
	ytop 	= cell2mat(param(8));
	dydx	= cell2mat(param(9));
	if (x <= x_tl)
		y_top = y_tl;
	elseif (x >= x_tr)
		y_top = y_tr + m*(x - x_tr);
	else
		y_top = ppval(ytop, x);
	end
	err = y_top - (y_c + (x - x_c))*(1E-6 + -1/dydx);
end




%		plot([x_bot(end), x_top(end)], [y_bot(end), y_top(end)], 'k-', 'LineWidth', 3)
%		plot([x_bot(1), x_top(1)], [y_bot(1), y_top(1)], 'k-', 'LineWidth', 3)

%		plot(x_top, y_top, 'm-', 'LineWidth', 3);
%		plot(x_bot, y_bot, 'b-', 'LineWidth', 3);

%		plot(x_c, y_c, 'k-o');

%	Lastly, spline the cord location by the normalized cord length s
%	dxds_1 = (x_c(2) - x_c(1))/(s(2) - s(1))
%	dyds_1 = (y_c(2) - y_c(1))/(s(2) - s(1))
%	dxds_2 = (x_c(end) - x_c(end-1))/(s(end) - s(end-1))
%	dyds_2 = (y_c(end) - y_c(end-1))/(s(end) - s(end-1))
%	dyds_1/dxds_1
%	dyds_2/dxds_2
%	dydx = spline(s, dydx);
%	xc = spline(s, x_c);
%	yc = spline(s, y_c);

%		First, evaluate the ?c splines to find the location of the ith point on the cord
%		x_loc(i) =  ppval(xc, s(i));
%		y_loc(i) = ppval(yc, s(i));

%	s = linspace(-0.05, 1.05, N);
%	xsect_top = mean(x_c)*ones(length(s), 1);
%	ysect_top = zeros(length(s), 1);
%	xsect_bot = mean(x_c)*ones(length(s), 1);
%	ysect_bot = zeros(length(s), 1);
%	x_c = ppval(xc, s);
%	y_c = ppval(yc, s);
%	dydx = ppval(dydx, s);
