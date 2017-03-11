classdef nonCompactSolution
	properties
		eta
		pi_plus
		pi_minus
		sigma
		theta
		Y
		p
		w
		s
		I_mass
		I_enthalpy
		I_entropy
		I_mixture
		I_theta
		M
		transfer
	end%properties
	methods
		function sol = nonCompactSolution(bvpsol, splines)
			gamma 	= splines.gamma;
			gm1 	= gamma - 1;
			gp1		= gamma + 1;
			gm1o2	= gm1/2;
			gp1o2	= gp1/2;

			eta 	= bvpsol.x(1:end);
			N 		= length(eta);
			p 		= bvpsol.y(1,1:end);
			w		= bvpsol.y(2,1:end);
			theta 	= bvpsol.y(3,1:end);
			s		= bvpsol.y(4,1:end);
			Y		= bvpsol.y(5,1:end);

%			build characteristics
			for i = 1:N
				M(i) 			= ppval(splines.M_bar.raw, eta(i));
				pi_plus(i) 		= (p(i) + w(i)*M(i)); 
				pi_minus(i) 	= (p(i) - w(i)*M(i));
				sigma(i)		= s(i);
			end%for i = 1:N

%			build invariants
			for i = 1:N
				alpha = 1./(1+gm1o2*M(i)*M(i));
				P = [	1 			1 						0 	-1 			0;
						alpha*gm1	M(i)*M(i)*alpha*gm1		0	alpha		0;
						0			0						0	1			0;
						0			0						0	0			1;
						0			0						1	0			0];
				x = [p(i); w(i); theta(i); s(i); Y(i)];
				I = P*x;
				I_mass(i) 		= I(1);
				I_enthalpy(i) 	= I(2);
				I_entropy(i) 	= I(3);
				I_mixture(i) 	= I(4);
				I_theta(i) 		= I(5);
			end%for i = 1:N

			sol.eta 		= eta;
			sol.p 			= p;
			sol.w			= w;
			sol.theta		= theta;
			sol.s			= s;
			sol.Y			= Y;
			sol.pi_plus		= pi_plus;
			sol.pi_minus	= pi_minus;
			sol.sigma		= sigma;
			sol.I_mass		= I_mass;
			sol.I_enthalpy 	= I_enthalpy;
			sol.I_entropy	= I_entropy;
			sol.I_mixture	= I_mixture;
			sol.I_theta		= I_theta;
			sol.M			= M;

			sol.transfer	= [pi_plus(1), pi_plus(end); pi_minus(1) pi_minus(end); sigma(1) sigma(end); theta(1) theta(end); Y(1) Y(end)];

		end%noncompactSolutionObject()

	end%methods
end%dclassdef
