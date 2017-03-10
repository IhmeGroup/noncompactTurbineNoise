function[extrema, x_c, y_c, ytop, ybot,stagnation_state] = rotor()
	N = 500;%# of points used to build interpolants
	delta = 60;%cord-cord separation
	rho0 = 1.2754;%Stagnation Condition
	T0 = 293;%Stagnation Condition
	M1 = 0.105;%Cascade entrance speed
	R = 287.058;%Universal gas constant
	gamma = 1.4;%c_p/c_v

%	Load the stator geometry
	data = load('stator.txt');
	data = data(1:120,1:2); % Clip to insure that only one pass over the geometry is retained
%	Identify the extreme points and split the geometry into a "top" and "bottom"
	[xmax, maxdex] = max(data(:,1));
	[xmin, mindex] = min(data(:,1));	
	x_top = data(mindex:maxdex,1);
	y_top = data(mindex:maxdex,2);
%	B/c the splining flips out over the repeated point, remove the last point in the series
	x_bot = [data(1:mindex,1); data(maxdex:end-4,1)];
	y_bot = [data(1:mindex,2); data(maxdex:end-4,2)];
%	Flip the geometry b/c the points are arranged clockwise
	y_bot = flipud(y_bot);
	x_bot = flipud(x_bot);
%	Clip some of the geometry so that the cord doesn't flip. 
	offset = 12;
	x_top = x_top(16:end);	
	y_top = y_top(16:end);	
	x_bot = x_bot(16:end);
	y_bot = y_bot(16:end);
%	Offset the bottom layer to focus on the void instead of the blade
	y_bot = y_bot - delta;
%	Grab the corners of the top and bottom geometry
	x_tl = x_top(1); y_tl = y_top(1);
	x_tr = x_top(end); y_tr = y_top(end);
	x_bl = x_bot(1); y_bl = y_bot(1);
	x_br = x_bot(end); y_br = y_bot(end);
%	Compute the slope of the blade tail to determine the angle at which the flow leaves the cascade
	exit_slope = (y_bot(end) - y_bot(end-1))/(x_bot(end) - x_bot(end-1))
%	Spline the geometry so that it can be interpolated onto the same grid in terms of x or s or some coordinate
%		Specify the derivative of the spline at the endpoints so that it smoothly transitions into the thing we want it to
	xtop = spline(x_top, x_top);
	ytop = spline(x_top, [0; y_top; exit_slope]);
	xbot = spline(x_bot, x_bot);
	ybot = spline(x_bot, [0; y_bot; exit_slope]);
%	And evaluate the spline on a constant grid, so that the top and bottom geometries (and thus the cord) are all on the same grid.
%		Add some "runway" on either side of the cord so that there is a region of homogenous flow
	x_c = linspace(xmin-10, xmax+10, N);
	x_top = x_c;
	x_bot = x_c;
	y_top = zeros(1,N);
	y_bot = zeros(1,N);
	dydx	= zeros(1,N);
	for i = 1:length(x_c)
%			xsect_top(i) = x_c(i);
%			xsect_bot(i) = x_c(i);
%			[ysect_top(i)] = cleverYtop(xsect_top(i), ytop, x_tl, y_tl, x_tr, y_tr, dydx(i));
%			[ysect_bot(i)] = cleverYbot(xsect_top(i), ybot, x_bl, y_bl, x_br, y_br, dydx(i));
%		param = {x_tl, y_tl, x_tr, y_tr, exit_slope, NaN, NaN, ytop, NaN};
%		y_top(i) = cleverYtop(x_c(i));
%		param = {x_bl, y_bl, x_br, y_br, exit_slope, NaN, NaN, ybot, NaN};
%		y_bot(i) = cleverYbot(x_c(i));
		xsect_top(i)	= x_c(i);
		xsect_bot(i)	= x_c(i);
		ysect_top(i)	= cleverYtop(xsect_top(i), ytop, x_tl, y_tl, x_tr, y_tr, dydx(i));
		ysect_bot(i)	= cleverYbot(xsect_top(i), ybot, x_bl, y_bl, x_br, y_br, dydx(i));
	end
%	The cord is just the average of the top and bottom contours
	y_c = 0.5*(y_top + y_bot);

	extrema = [x_tl; y_tl; x_tr; y_tr; x_bl; y_bl; x_br; y_br; exit_slope];
	stagnation_state = [gamma; R; rho0; T0; M1];
end
