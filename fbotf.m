function[x] = fbotf();%, ftop, xtl, ytl, xtr, ytr, m, x_c, y_c)
%need to use global structure to pass params in
	global param;
	x_bl 	= cell2mat(param(1));
	y_bl 	= cell2mat(param(2));
	x_br 	= cell2mat(param(3));
	y_br 	= cell2mat(param(4));
	m 		= cell2mat(param(5));
	x_c 	= cell2mat(param(6));
	y_c 	= cell2mat(param(7));
	ybot 	= cell2mat(param(8));
	dydx	= cell2mat(param(9));
	its = 0;
	tol = 1E-2;
	x = x_c;
	change = false;
	while (change == false)
		y_tan = y_c + (x-x_c)*(-1/dydx);
		if (x <= x_bl)
			y_bot = y_bl;
		elseif (x >= x_br)
			y_bot = y_br + m*(x - x_br);
		else
			y_bot = ppval(ybot, x);
		end
		if (y_bot < y_tan)
			x = x + sign(dydx)*tol;
		else
			change = true;
		end
		its = its + 1;
		if (its > 10000)
			disp('bottom solver timed out');
			break;
		end
	end
end
