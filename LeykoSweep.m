function[Theta_1 Theta_2 pi_plus, pi_minus, sigma, theta, Y] = LeykoSweep(N1, N2, type)
	M_1 = 0.12;
	M_2 = 0.12;

	theta_1 	= linspace(0, 90, N1);
	theta_2 	= linspace(1E-5, 360, N1);
	Theta_1 	= zeros(N1);
	Theta_2 	= zeros(N1);
	pi_plus 	= zeros(N1);
	pi_minus 	= zeros(N1);
	sigma		= zeros(N1);
	theta		= zeros(N1);
	Y			= zeros(N1);
	for i = 1:1
		for j = 1:N1
			leyko = SplineHolder(5, M_1, M_2, theta_1(i), theta_2(j));
			[transfer, sol, bvpsol] = nonCompactSolver(leyko, N2, 0, true, type)
			transfer
			pi_plus(i,j) 	= transfer(1,2);
			pi_minus(i,j)	= transfer(2,1);
			theta(i,j)		= transfer(4,2);
			sigma(i,j)		= transfer(3,2);
			Y(i,j)			= transfer(5,2);
			Theta_1(i,j) 	= theta_1(i);
			Theta_2(i,j) 	= theta_2(j);
		end
	end
	Theta_1
	Theta_2

	g = figure();
	for a = 1:5
		subplot(2,3,a)
		for i = 1:N1
			if (a == 1)
				plot(Theta_2(i,:), abs(pi_plus(i,:)), 'LineWidth', 2);
			elseif (a == 2)
				plot(Theta_2(i,:), abs(pi_minus(i,:)), 'LineWidth', 2);
			elseif (a == 3)
				plot(Theta_2(i,:), abs(sigma(i,:)), 'LineWidth', 2);
			elseif (a == 4)
				plot(Theta_2(i,:), abs(theta(i,:)), 'LineWidth', 2);
			elseif (a == 5)
				plot(Theta_2(i,:), abs(Y(i,:)), 'LineWidth', 2);
			end
			hold on;
		end
		xlabel('$\theta_2$', 'FontSize', 22, 'FontName', 'Times', 'Interpreter', 'LaTeX');
		if (a == 1)
			ylabel('$\pi^+$', 'FontSize', 22, 'FontName', 'Times', 'Interpreter', 'LaTeX');
%			legend('0', '10', '20', '30', '40', '50', '60', '70', '80', '90');
		elseif (a == 2)
			ylabel('$\pi^-$', 'FontSize', 22, 'FontName', 'Times', 'Interpreter', 'LaTeX');
		elseif (a == 3)
			ylabel('$\sigma$', 'FontSize', 22, 'FontName', 'Times', 'Interpreter', 'LaTeX');
		elseif (a == 4)
			ylabel('$\theta$', 'FontSize', 22, 'FontName', 'Times', 'Interpreter', 'LaTeX');
		elseif (a == 5)
			ylabel('$Y_I$', 'FontSize', 22, 'FontName', 'Times', 'Interpreter', 'LaTeX');
		end

end
