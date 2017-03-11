function[res] = illPosedBCs(phi_a, phi_b, splines, type)
	res = zeros(5,1);
	M_1 = ppval(splines.M_bar.raw, 0);
	M_2 = ppval(splines.M_bar.raw, 1);
	res(1) = (phi_a(1) + phi_a(2)*M_1) - 0;% pi plus
	res(2) = (phi_b(1) - phi_b(2)*M_2) - 0;% pi minus
	res(3) = phi_a(3) - 0;% theta
	res(4) = phi_a(4) - 0;% entropy
	res(5) = phi_a(5) - 0;% Y

	if (type == 1)
		res(1) 	= (phi_a(1) + phi_a(2)*M_1) - 1;% pi plus
	elseif (type == 2)
		res(3) = phi_a(3) - 1;
	elseif (type == 3)
		res(4)	= phi_a(4) - 1;
	elseif (type == 4)
		res(5) 	= phi_a(5) - 1;
	elseif (type == 5)
		res(1) = phi_a(1) - 1;
	else
		error('Garbage forcing type picked');
	end
end
