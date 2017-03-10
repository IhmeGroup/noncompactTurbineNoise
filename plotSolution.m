function[] = plotncsution(ncs)
	close all;
	plot_primitives 		= true;
	plot_invariants 		= true;
	plot_characteristics 	= true;

	if (plot_invariants)
		g = figure();
		set(g, 'Position', [0 0 1500 500]);
		subplot(2,3,1);
		plot(ncs.eta(:), abs(ncs.I_mass), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$I_m$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		set	(gca, 'FontSize', 22, 'FontName', 'Times');
		title('Mass', 'FontSize', 22, 'FontName', 'Times');

		subplot(2,3,2);
		plot(ncs.eta(:), abs(ncs.I_enthalpy), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$I_h$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
		title('Enthalpy', 'FontSize', 22, 'FontName', 'Times');

		subplot(2,3,3);
		plot(ncs.eta(:), abs(ncs.I_mixture), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$I_Y$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
		title('Mixture', 'FontSize', 22, 'FontName', 'times');
	
		subplot(2,3,4);
		plot(ncs.eta(:), abs(ncs.I_entropy), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$I_s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Entropy', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');

		subplot(2,3,5);
		plot(ncs.eta(:), abs(ncs.I_theta), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$I_\theta$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		set(gca, 'FontSize', 22, 'FontName', 'Times');

		subplot(2,3,6);
		plot(ncs.eta(:), ncs.M, 'LineWidth', 2);
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\bar{M}$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Mach', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	end%(plot_invariants)

	if (plot_characteristics)
		h = figure();
		set(h, 'Position', [0 0 1500 500]);
		subplot(2,3,1);
		plot(ncs.eta(:), abs(ncs.pi_plus), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\pi^+$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Downstream Acoustic', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	
		subplot(2,3,2);
		plot(ncs.eta(:), abs(ncs.pi_minus), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\pi^-$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Upstream Acoustic', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	
		subplot(2,3,3);
		plot(ncs.eta(:), abs(ncs.sigma), 'LineWidth', 2);
			ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\sigma$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Entropy', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	
		subplot(2,3,4);
		plot(ncs.eta(:), abs(ncs.theta), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\theta$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Theta', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');

		subplot(2,3,5);
		plot(ncs.eta(:), abs(ncs.Y), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$Y_i$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Species', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');

		subplot(2,3,6);
		plot(ncs.eta(:), ncs.M, 'LineWidth', 2);
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\bar{M}$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Mach', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	end%(plot_characteristics)
	if (plot_primitives)
		j = figure();
		set(j, 'Position', [0 0 1500 500]);
		subplot(2,3,1);
		plot(ncs.eta(:), abs(ncs.p), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\frac{p}{\gamma \bar{p}}$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Pressure', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	
		subplot(2,3,2);
		plot(ncs.eta(:), abs(ncs.w), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\frac{w}{\bar{w}}$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Velocity', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	
		subplot(2,3,3);
		plot(ncs.eta(:), abs(ncs.s), 'LineWidth', 2);
			ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\frac{s}{\bar{c}_p}$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Entropy', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	
		subplot(2,3,4);
		plot(ncs.eta(:), abs(ncs.theta), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\theta$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Theta', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');

		subplot(2,3,5);
		plot(ncs.eta(:), abs(ncs.Y), 'LineWidth', 2);
		ylims = ylim();
		if (abs(ylims(2) - ylims(1)) < 1E-4)
			if (max(ylims) < 1E-4)
				ylim([0,1]);
			else
				ylim([0,2]);
			end
		end
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$Y_i$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Species', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');

		subplot(2,3,6);
		plot(ncs.eta(:), ncs.M, 'LineWidth', 2);
		xlabel('$s$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		ylabel('$\bar{M}$', 'FontSize', 22, 'FontName', 'Times','Interpreter', 'LaTeX');
		title('Mach', 'FontSize', 22, 'FontName', 'Times');
		set(gca, 'FontSize', 22, 'FontName', 'Times');
	end%(plot_characteristics)
end
