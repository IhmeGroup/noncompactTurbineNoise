function[y_bot] = cleverYbot(x)
	global param;
	x_bl 	= cell2mat(param(1));
	y_bl 	= cell2mat(param(2));
	x_br 	= cell2mat(param(3));
	y_br 	= cell2mat(param(4));
	m 		= cell2mat(param(5));
%	x_c 	= cell2mat(param(6));
%	y_c 	= cell2mat(param(7));
	ybot 	= cell2mat(param(8));
%	dydx	= cell2mat(param(9));
	if (x <= x_bl)
		y_bot = y_bl;
	elseif (x >= x_br)
		y_bot = y_br + m*(x - x_br);
	else
		y_bot = ppval(ybot, x);
	end
end
