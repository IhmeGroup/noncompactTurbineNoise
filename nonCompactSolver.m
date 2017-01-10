function[] = nonCompactSolver()
	close all;
	N = 201;%number of spatial points

%	Constant parameters
	global gamma;
	global gm1; 
	global gp1;
	global gm1o2;
	global gp1o2;

%	Run-specific parameters
	global omega;

%	Thermodynamic and hydrodynamic splines
	global areaSp;
	global thetSp;
	global xvelSp;
	global yvelSp;
	global machSp;
	global presSp;
	global densSp;
	global tempSp;
	global sos_Sp;

%	Spatial gradient splines
	global ddensdxSp;
	global ddensdySp;
	global dxveldxSp;
	global dxveldySp;
	global dyveldxSp;
	global dyveldySp;
	
	suppress = true;
	[areaSp, thetSp, xvelSp, yvelSp, machSp, presSp, densSp, tempSp, sos_Sp, ddensdxSp, ddensdySp, dxveldxSp, dxveldySp, dyveldxSp, dyveldySp] = buildMeanFlow(suppress)

	omega = 0.1;
	phi0 = [1; 0; 0; 0];
	sol = bvpinit(linspace(0, 1, N), [0.5 0.5 0.5 0.5]);
	sol = bvp4c(@nonCompactODE, @illPosedBCs, sol);
	s = sol.x;
	phi = sol.y;
	size(s)
	size(phi)


%^	[s, phi] = ode45(@nonCompactODE, [0, 1], phi0);

	g = figure();
	set(g, 'Position', [0 0 800 600]);
	plot(s, abs(phi(1,:)), 'k-', 'LineWidth', 2);
	hold on;
	plot(s, abs(phi(2,:)), 'b-', 'LineWidth', 2);
	plot(s, abs(phi(3,:)), 'r-', 'LineWidth', 2);
	plot(s, abs(phi(4,:)), 'm-', 'LineWidth', 2);
	h = legend('$p \slash \gamma \bar{p}$', '$u \slash \bar{u}$', '$v \slash \bar{v}$', '$s \slash c_p$');
	set(h, 'Interpreter', 'LaTeX', 'FontSize', 14, 'FontName', 'Times');
	xlabel('$s$', 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	ylabel('$[-]$', 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'LaTeX');
	set(gca, 'FontSize', 14, 'FontName', 'Times');

end

