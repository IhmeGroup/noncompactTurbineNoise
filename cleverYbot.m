function[y_bot] = cleverYbot(x, ybot, x_bl, y_bl, x_br, y_br, m)
	if (x <= x_bl)
		y_bot = y_bl;
	elseif (x >= x_br)
		y_bot = y_br + m*(x - x_br);
	else
		y_bot = ppval(ybot, x);
	end
end
