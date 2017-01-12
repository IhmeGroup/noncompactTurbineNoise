function[x] = ftopf();%, ftop, xtl, ytl, xtr, ytr, m, x_c, y_c)
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

	its = 0;
	tol = 1E-2;
	x = x_c;
	change = false;
	while (change == false)
		y_tan = y_c + (x-x_c)*(-1/dydx);
		if (x <= x_tl)
			y_top = y_tl;
		elseif (x >= x_tr)
			y_top = y_tr + m*(x - x_tr);
		else
			y_top = ppval(ytop, x);
		end
		if (y_top > y_tan)
			x = x - sign(dydx)*tol;
		else
			change = true;
		end
		its = its + 1;
		if (its > 10000)
			disp('Top solver timed out');
			break;
		end
	end
end
