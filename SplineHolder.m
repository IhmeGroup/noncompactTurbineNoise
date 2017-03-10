classdef SplineHolder
	properties
		x_c
		y_c
		gamma
		A
		theta_bar
		p_bar
		rho_bar
		T_bar
		u_bar
		v_bar
		w_bar
		M_bar
		sos_bar
		psi_bar
	end%properties
	methods
		function splines = SplineHolder(run)

			if (run == 1)
				disp('Building base flow for straight channel');
				[boundaries, x_c, y_c, ytop, ybot, stagnation_state] = channel();

			elseif (run == 2)
				disp('Buildling base flow for linear nozzle');
				[boundaries, x_c, y_c, ytop, ybot, stagnation_state] = nozzle();

			elseif (run == 3)
				disp('Building base flow for curved channel');
				[boundaries, x_c, y_c, ytop, ybot, stagnation_state] = curvedChannel();

			elseif (run == 4)
				disp('Building base flow for single stator');
				[boundaries, x_c, y_c, ytop, ybot, stagnation_state] = rotor();

			else
				error('Run type not implemented yet');
			end

			[baseFlow] 	= buildBaseFlow(boundaries, x_c, y_c, ytop, ybot, stagnation_state);

			splines.x_c			= SplineField(baseFlow.s,baseFlow.x_c);
			splines.y_c			= SplineField(baseFlow.s,baseFlow.y_c);
			splines.gamma		= baseFlow.gamma;
			splines.A			= SplineField(baseFlow.s,baseFlow.A);
			splines.theta_bar 	= SplineField(baseFlow.s,baseFlow.theta_bar);
			splines.p_bar 		= SplineField(baseFlow.s,baseFlow.p_bar);
			splines.rho_bar 	= SplineField(baseFlow.s,baseFlow.rho_bar);
			splines.T_bar 		= SplineField(baseFlow.s,baseFlow.T_bar);
			splines.u_bar 		= SplineField(baseFlow.s,baseFlow.u_bar);
			splines.v_bar 		= SplineField(baseFlow.s,baseFlow.v_bar);
			splines.w_bar 		= SplineField(baseFlow.s,baseFlow.w_bar);
			splines.M_bar 		= SplineField(baseFlow.s,baseFlow.M_bar);
			splines.sos_bar 	= SplineField(baseFlow.s,baseFlow.sos_bar);
			splines.psi_bar 	= SplineField(baseFlow.s,baseFlow.psi_bar);
		end%SplineHolder()

		function checkRawSplines(splineHolder, N)
			s = linspace(0,1,N);

			theta	= zeros(N,1);
			p		= zeros(N,1);
			rho 	= zeros(N,1);
			T		= zeros(N,1);
			u		= zeros(N,1);
			v 		= zeros(N,1);
			w		= zeros(N,1);
			M		= zeros(N,1);
			sos		= zeros(N,1);

			x_c		= ppval(splineHolder.x_c.raw, s);
			y_c		= ppval(splineHolder.y_c.raw, s);
			A		= ppval(splineHolder.A.raw, s);
			theta 	= ppval(splineHolder.theta_bar.raw, s);
			p		= ppval(splineHolder.p_bar.raw, s);
			rho		= ppval(splineHolder.rho_bar.raw, s);
			T		= ppval(splineHolder.T_bar.raw, s);
			u		= ppval(splineHolder.u_bar.raw, s);
			v		= ppval(splineHolder.v_bar.raw, s);
			w		= ppval(splineHolder.w_bar.raw, s);
			M		= ppval(splineHolder.M_bar.raw, s);
			sos		= ppval(splineHolder.sos_bar.raw, s);
			psi		= ppval(splineHolder.psi_bar.raw, s);

			h = figure();
			set(h, 'Position', [0 0 1500 750]);
			splotMe(s, x_c, 's', '$x_c$', 'X-cord', 1);
			splotMe(s, y_c, 's', '$x_c$', 'Y-cord', 2);
			splotMe(s, A, 'A', '$A$', 'Area', 3);
			splotMe(s, theta, 's', '$\theta$', 'Theta', 4);
			splotMe(s, p, 's', '$p$', 'Pressure', 5);
			splotMe(s, rho, 's', '$\rho$', 'Density', 6);
			splotMe(s, T, 's', '$T$', 'Temp', 7);
			splotMe(s, psi, 's', '$\psi$', 'Chem Pot', 8);
			splotMe(s, u, 's', '$u$', 'X-vel', 9);
			splotMe(s, v, 's', '$v$', 'Y-Vel', 10);
			splotMe(s, w, 's', '$w$', 'Vel Mag', 11);
			splotMe(s, M, 's', '$M$', 'Mach', 12);

		end%checkSplines()

		function checkddxSplines(splineHolder, N)
			s = linspace(0,1,N);

			theta	= zeros(N,1);
			p		= zeros(N,1);
			rho 	= zeros(N,1);
			T		= zeros(N,1);
			u		= zeros(N,1);
			v 		= zeros(N,1);
			w		= zeros(N,1);
			M		= zeros(N,1);
			sos		= zeros(N,1);

			x_c		= ppval(splineHolder.x_c.ddx, s);
			y_c		= ppval(splineHolder.y_c.ddx, s);
			A		= ppval(splineHolder.A.ddx, s);
			theta 	= ppval(splineHolder.theta_bar.ddx, s);
			p		= ppval(splineHolder.p_bar.ddx, s);
			rho		= ppval(splineHolder.rho_bar.ddx, s);
			T		= ppval(splineHolder.T_bar.ddx, s);
			u		= ppval(splineHolder.u_bar.ddx, s);
			v		= ppval(splineHolder.v_bar.ddx, s);
			w		= ppval(splineHolder.w_bar.ddx, s);
			M		= ppval(splineHolder.M_bar.ddx, s);
			sos		= ppval(splineHolder.sos_bar.ddx, s);
			psi		= ppval(splineHolder.psi_bar.ddx, s);

			h = figure();
			splotMe(s, x_c, 	's', '$\frac{\partial x_c}{\partial x}$', 'X-cord', 1);
			splotMe(s, y_c, 	's', '$\frac{\partial x_c}{\partial x}$', 'Y-cord', 2);
			splotMe(s, A, 		's', '$\frac{\partial A}{\partial x}$', 'Area', 3);
			splotMe(s, theta, 	's', '$\frac{\partial \theta}{\partial x}$', 'Theta', 4);
			splotMe(s, p, 		's', '$\frac{\partial p}{\partial x}$', 'Pressure', 5);
			splotMe(s, rho, 	's', '$\frac{\partial \rho}{\partial x}$', 'Density', 6);
			splotMe(s, T, 		's', '$\frac{\partial T}{\partial x}$', 'Temp', 7);
			splotMe(s, psi, 	's', '$\frac{\partial \psi}{\partial x}$', 'Chem Pot', 8);
			splotMe(s, u, 		's', '$\frac{\partial u}{\partial x}$', 'X-vel', 9);
			splotMe(s, v, 		's', '$\frac{\partial v}{\partial x}$', 'Y-Vel', 10);
			splotMe(s, w, 		's', '$\frac{\partial w}{\partial x}$', 'Vel Mag', 11);
			splotMe(s, M, 		's', '$\frac{\partial M}{\partial x}$', 'Mach', 12);
			set(h, 'Position', [0 0 1500 750]);

		end

		function checkddySplines(splineHolder, N)
			s = linspace(0,1,N);

			theta	= zeros(N,1);
			p		= zeros(N,1);
			rho 	= zeros(N,1);
			T		= zeros(N,1);
			u		= zeros(N,1);
			v 		= zeros(N,1);
			w		= zeros(N,1);
			M		= zeros(N,1);
			sos		= zeros(N,1);

			x_c		= ppval(splineHolder.x_c.ddy, s);
			y_c		= ppval(splineHolder.y_c.ddy, s);
			A		= ppval(splineHolder.A.ddy, s);
			theta 	= ppval(splineHolder.theta_bar.ddy, s);
			p		= ppval(splineHolder.p_bar.ddy, s);
			rho		= ppval(splineHolder.rho_bar.ddy, s);
			T		= ppval(splineHolder.T_bar.ddy, s);
			u		= ppval(splineHolder.u_bar.ddy, s);
			v		= ppval(splineHolder.v_bar.ddy, s);
			w		= ppval(splineHolder.w_bar.ddy, s);
			M		= ppval(splineHolder.M_bar.ddy, s);
			sos		= ppval(splineHolder.sos_bar.ddy, s);
			psi		= ppval(splineHolder.psi_bar.ddy, s);

			h = figure();
			set(h, 'Position', [0 0 1500 750]);

			splotMe(s, x_c, 	's', '$\frac{\partial x_c}{\partial y}$', 'X-cord', 1);
			splotMe(s, y_c, 	's', '$\frac{\partial x_c}{\partial y}$', 'Y-cord', 2);
			splotMe(s, A, 		's', '$\frac{\partial A}{\partial y}$', 'Area', 3);
			splotMe(s, theta, 	's', '$\frac{\partial \theta}{\partial y}$', 'Theta', 4);
			splotMe(s, p, 		's', '$\frac{\partial p}{\partial y}$', 'Pressure', 5);
			splotMe(s, rho, 	's', '$\frac{\partial \rho}{\partial y}$', 'Density', 6);
			splotMe(s, T, 		's', '$\frac{\partial T}{\partial y}$', 'Temp', 7);
			splotMe(s, psi, 	's', '$\frac{\partial \psi}{\partial y}$', 'Chem Pot', 8);
			splotMe(s, u, 		's', '$\frac{\partial u}{\partial y}$', 'X-vel', 9);
			splotMe(s, v, 		's', '$\frac{\partial v}{\partial y}$', 'Y-Vel', 10);
			splotMe(s, w, 		's', '$\frac{\partial w}{\partial y}$', 'Vel Mag', 11);
			splotMe(s, M, 		's', '$\frac{\partial M}{\partial y}$', 'Mach', 12);
		end

		function checkDDtSplines(splineHolder, N)
			s = linspace(0,1,N);

			theta	= zeros(N,1);
			p		= zeros(N,1);
			rho 	= zeros(N,1);
			T		= zeros(N,1);
			u		= zeros(N,1);
			v 		= zeros(N,1);
			w		= zeros(N,1);
			M		= zeros(N,1);
			sos		= zeros(N,1);

			x_c		= ppval(splineHolder.x_c.DDt, s);
			y_c		= ppval(splineHolder.y_c.DDt, s);
			A		= ppval(splineHolder.A.DDt, s);
			theta 	= ppval(splineHolder.theta_bar.DDt, s);
			p		= ppval(splineHolder.p_bar.DDt, s);
			rho		= ppval(splineHolder.rho_bar.DDt, s);
			T		= ppval(splineHolder.T_bar.DDt, s);
			u		= ppval(splineHolder.u_bar.DDt, s);
			v		= ppval(splineHolder.v_bar.DDt, s);
			w		= ppval(splineHolder.w_bar.DDt, s);
			M		= ppval(splineHolder.M_bar.DDt, s);
			sos		= ppval(splineHolder.sos_bar.DDt, s);
			psi		= ppval(splineHolder.psi_bar.DDt, s);

			h = figure();
			set(h, 'Position', [0 0 1500 750]);

			splotMe(s, x_c, 	's', '$\frac{D x_c}{D t}$', 'X-cord', 1);
			splotMe(s, y_c, 	's', '$\frac{D y_c}{D t}$', 'Y-cord', 2);
			splotMe(s, A, 		's', '$\frac{DA}{Dt}$', 'Area', 3);
			splotMe(s, theta, 	's', '$\frac{D\theta}{Dt}$', 'Theta', 4);
			splotMe(s, p, 		's', '$\frac{Dp}{Dt}$', 'Pressure', 5);
			splotMe(s, rho, 	's', '$\frac{D\rho}{Dt}$', 'Density', 6);
			splotMe(s, T, 		's', '$\frac{DT}{Dt}$', 'Temp', 7);
			splotMe(s, psi, 	's', '$\frac{D\psi}{Dt}$', 'Chem Pot', 8);
			splotMe(s, u, 		's', '$\frac{Du}{Dt}$', 'X-vel', 9);
			splotMe(s, v, 		's', '$\frac{Dv}{Dt}$', 'Y-Vel', 10);
			splotMe(s, w, 		's', '$\frac{Dw}{Dt}$', 'Vel Mag', 11);
			splotMe(s, M, 		's', '$\frac{DM}{Dt}$', 'Mach', 12);
		end
	end%methods
end%class SplineHolder
