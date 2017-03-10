function[] = splotMe(s, phi, xl, yl, tl, i)
	subplot(3,4,i)
	plot(s, phi, 'LineWidth', 2);
	xlabel(xl, 'Interpreter', 'LaTeX', 'FontSize', 14);
	ylabel(yl, 'Interpreter', 'LaTeX', 'FontSize', 14);
	title(tl, 'Interpreter', 'LaTeX', 'FontSize', 14);
end
