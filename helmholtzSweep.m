function[He, pi_plus, pi_minus, sigma, theta, Y] = helmholtzSweep(N, N2, splines)


	load('../composition_noise/entropicResponseFigure/entropicResponseData.mat');

	He = linspace(0,2,N);
	pi_plus 	= zeros(N,1);
	pi_minus 	= zeros(N,1);
	sigma 		= zeros(N,1);
	theta		= zeros(N,1);
	Y			= zeros(N,1);

	M_2 = 0.89;

	for i = 1:N
		[raw_transfer, sol] = nonCompactSolver(splines, N2, He(i), true);
		raw_transfer 	= abs(raw_transfer);
		pi_plus(i) 		= raw_transfer(1,2);
		pi_minus(i)		= raw_transfer(2,1);
		sigma(i)		= raw_transfer(3,2);
		theta(i)		= raw_transfer(4,2);
		Y(i)			= raw_transfer(5,2);
	end

	h = figure();
	set(h, 'Position', [0 0 1500 500]);
	subplot(2,3,1);
	plot(He,  abs(pi_plus), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
	hold on;
	plot(OMEGA, squeeze(abs(TRANS(1,:,1))), 'r-o', 'LineWidth', 2);
	xlabel('$he$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	ylabel('$\pi^+$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	title('Downstream', 'FontSize', 22, 'FontName', 'Times');
	set(gca, 'FontSize', 22, 'FontName', 'Times');

	subplot(2,3,2);
	plot(He, abs(pi_minus), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
	hold on;
	plot(OMEGA, squeeze(abs(TRANS(1,:,5))), 'r-o','LineWidth', 2);
	xlabel('$He$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	ylabel('$\pi^-$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	title('Upstream', 'FontSize', 22, 'FontName', 'Times');
	set(gca, 'FontSize', 22, 'FontName', 'Times');

	subplot(2,3,3);
	plot(He, abs(theta), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
	xlabel('$He$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	ylabel('$\theta$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	title('Theta', 'FontSize', 22, 'FontName', 'Times');
	set(gca, 'FontSize', 22, 'FontName', 'Times');

	subplot(2,3,4);
	plot(He, abs(sigma), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
	hold on;
	plot(OMEGA, squeeze(abs(TRANS(1,:,3))), 'r-o', 'LineWidth', 2);
	xlabel('$He$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	ylabel('$\frac{s}{\bar{c}_p}$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	title('Entropy', 'FontSize', 22, 'FontName', 'Times');
	set(gca, 'FontSize', 22, 'FontName', 'Times');

	subplot(2,3,5);
	plot(He,  abs(Y), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
	hold on;
	plot(OMEGA, squeeze(abs(TRANS(1,:,4))), 'r-o', 'LineWidth', 2);
	xlabel('$He$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	ylabel('$Y_i$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	title('Species', 'FontSize', 22, 'FontName', 'Times');
	set(gca, 'FontSize', 22, 'FontName', 'Times');

end
