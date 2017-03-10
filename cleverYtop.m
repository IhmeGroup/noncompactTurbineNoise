function[y_top] = cleverYtop(x, ytop, x_tl, y_tl, x_tr, y_tr, m)
	if (x <= x_tl)
		y_top = y_tl;
	elseif (x >= x_tr)
		y_top = y_tr + m*(x - x_tr);
	else
		y_top = ppval(ytop, x);
	end
end
