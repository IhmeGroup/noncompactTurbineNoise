function[y_top] = cleverYtop(x)
	global param;
	x_tl 	= cell2mat(param(1));
	y_tl 	= cell2mat(param(2));
	x_tr 	= cell2mat(param(3));
	y_tr 	= cell2mat(param(4));
	m 		= cell2mat(param(5));
%	x_c 	= cell2mat(param(6));
%	y_c 	= cell2mat(param(7));
	ytop 	= cell2mat(param(8));
%	dydx	= cell2mat(param(9));
	if (x <= x_tl)
		y_top = y_tl;
	elseif (x >= x_tr)
		y_top = y_tr + m*(x - x_tr);
	else
		y_top = ppval(ytop, x);
	end
end
