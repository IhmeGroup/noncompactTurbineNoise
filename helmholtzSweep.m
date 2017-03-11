function[He, pi_plus, pi_minus, sigma, theta, Y] = helmholtzSweep(N, N2, splines, type)
	close all;


	if (type == 1)
		load('../composition_noise/acousticResponseFigure/acousticResponseData.mat');
	elseif (type == 3)
		load('../composition_noise/entropicResponseFigure/entropicResponseData.mat');
	elseif (type == 4)
		load('../composition_noise/compositionResponseFigure/compositionResponseData.mat');
	else
		error('Unrecognized Forcing type');
	end

	He = linspace(0,2,N);
	pi_plus 	= zeros(N,1);
	pi_minus 	= zeros(N,1);
	sigma 		= zeros(N,1);
	theta		= zeros(N,1);
	Y			= zeros(N,1);

	M_2 = 0.89;

	for i = 1:N
		[raw_transfer, sol] = nonCompactSolver(splines, N2, He(i), true, type);
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
	plot(OMEGA, squeeze(abs(TRANS(1,:,1))), 'ro', 'LineWidth', 2);
	hold on;
	plot(He,  abs(pi_plus), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
	xlabel('$He$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	ylabel('$\pi^+$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	title('Downstream', 'FontSize', 22, 'FontName', 'Times');
	set(gca, 'FontSize', 22, 'FontName', 'Times');

	subplot(2,3,2);
	plot(OMEGA, squeeze(abs(TRANS(1,:,5))), 'ro','LineWidth', 2);
	hold on;
	plot(He, abs(pi_minus), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
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
	plot(OMEGA, squeeze(abs(TRANS(1,:,3))), 'ro', 'LineWidth', 2);
	hold on;
	plot(He, abs(sigma), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
	xlabel('$He$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	ylabel('$\frac{s}{\bar{c}_p}$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	title('Entropy', 'FontSize', 22, 'FontName', 'Times');
	set(gca, 'FontSize', 22, 'FontName', 'Times');

	subplot(2,3,5);
	plot(OMEGA, squeeze(abs(TRANS(1,:,4))), 'ro', 'LineWidth', 2);
	hold on;
	plot(He,  abs(Y), 'LineWidth', 2);
	ylims = ylim();
	if (abs(ylims(2) - ylims(1)) < 1E-4)
		if (max(ylims) < 1E-4)
			ylim([0,1]);
		else
			ylim([0,2]);
		end
	end
	xlabel('$He$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	ylabel('$Y_i$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
	title('Species', 'FontSize', 22, 'FontName', 'Times');
	set(gca, 'FontSize', 22, 'FontName', 'Times');

end
