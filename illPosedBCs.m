function[Res] = illPosedBCs(q_a, q_b)

	Res = zeros(4,1);

	Res(1) = q_a(1) - 1;
	Res(2) = q_a(2) - 1;
	Res(3) = q_a(3) - 1;
	Res(4) = q_a(4) - 1;

end
