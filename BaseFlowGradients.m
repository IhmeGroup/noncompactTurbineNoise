classdef BaseFlowGradients
	properties
		s
		x_c
		y_c
		A
		theta_bar
		p_bar
		rho_bar
		T_bar
		u_bar
		v_bar
		M_bar
		sos_bar
	end%properties
	methods
		function BaseFlowGradients =  BaseFlowGradients(BaseFlow)
			s 		= BaseFlow.s;
			x_c		= BaseFlow.x_c;
			y_c		= BaseFlow.y_c;
			A		= BaseFlow.A;

			BaseFlowGradients.p_bar = build
			BaseFlow.A			= A;
			BaseFlow.theta_bar	= theta_bar;
			BaseFlow.p_bar		= p_bar;
			BaseFlow.rho_bar	= rho_bar;
			BaseFlow.T_bar		= T_bar;
			BaseFlow.u_bar		= u_bar;
			BaseFlow.v_bar		= v_bar;
			BaseFlow.M_bar		= M_bar;
			BaseFlow.sos_bar	= sos_bar;
		end%function
	end%methods
end%classdef
