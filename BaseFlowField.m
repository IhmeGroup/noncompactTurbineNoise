classdef BaseFlowField
	properties
		raw
		ddx
		ddy
		ddn
		DDt
	end%properties
	methods
		function BaseFlowField =  BaseFlowField(s, x_c, y_c, u_bar, v_bar, theta, phi)
			N = length(phi);
			eps = 1E-10;

			ddx = zeros(N,1);
			ddy = zeros(N,1);

			dx = x_c(2) - x_c(1);
			dy = y_c(2) - y_c(1);
			ddn(1) = (phi(2) - phi(1))/sqrt(dx^2 + dy^2);
			ddx(1) = ddn(1)*cos(theta(1));
			ddy(1) = ddn(1)*sin(theta(1));

			for i = 2:N-1
				dx = x_c(i+1) - x_c(i-1);
				dy = y_c(i+1) - y_c(i-1);
				ddn(i) = (phi(i+1) - phi(i-1))/sqrt(dx^2 + dy^2);
				ddx(i) = ddn(i)*cos(theta(i));
				ddy(i) = ddn(i)*sin(theta(i));
			end

			dx = x_c(N) - x_c(N-1);
			dy = y_c(N) - y_c(N-1);
			ddn(N) = (phi(N) - phi(N-1))/sqrt(dx^2 + dy^2);
			ddx(N) = ddn(N)*cos(theta(N));
			ddy(N) = ddn(N)*sin(theta(N));

			BaseFlowField.raw = phi;
			BaseFlowField.ddx = ddx;
			BaseFlowField.ddy = ddy;
			BaseFlowField.ddn = ddn;
			BaseFlowField.DDt = u_bar.*ddx + v_bar.*ddy;
		end%function
	end%methods
end%classdef
