classdef BaseFlow
	properties
		s
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
		function BaseFlow =  BaseFlow(s, x_c, y_c, gamma, A, theta_bar, p_bar, rho_bar, T_bar, u_bar, v_bar, w_bar, M_bar, sos_bar, psi_bar);	
				rho_a 	= rho_bar(1);
				p_a		= p_bar(1);
				c_a		= sqrt(gamma*p_a/rho_a);

				BaseFlow.s 			= s;
				BaseFlow.x_c 		= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, x_c);
				BaseFlow.y_c 		= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, y_c);
				BaseFlow.gamma		= gamma;
				BaseFlow.A			= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, A);
				BaseFlow.theta_bar	= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, theta_bar);
				BaseFlow.p_bar		= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, p_bar./p_a);
				BaseFlow.rho_bar	= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, rho_bar./rho_a);
				BaseFlow.T_bar		= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, T_bar);
				BaseFlow.u_bar		= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, u_bar);
				BaseFlow.v_bar		= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, v_bar);
				BaseFlow.w_bar		= BaseFlowField(s, x_c, y_c, u_bar, w_bar, theta_bar, w_bar./c_a);
				BaseFlow.M_bar		= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, M_bar);
				BaseFlow.sos_bar	= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, sos_bar./c_a);
				BaseFlow.psi_bar	= BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta_bar, psi_bar);
		end%function
	end%methods
end%classdef
