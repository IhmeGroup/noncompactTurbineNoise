classdef SplineField
	properties
		raw
		ddx
		ddy
		ddn
		DDt
	end%properties
	methods
		function SplineField =  SplineField(s, field)
			phi = field.raw;
			ddx = field.ddx;
			ddy = field.ddy;
			ddn = field.ddn;
			DDt = field.DDt;

			SplineField.raw = spline(s, phi);
			SplineField.ddx = spline(s, ddx);
			SplineField.ddy = spline(s, ddy);
			SplineField.ddn	= spline(s, ddn);
			SplineField.DDt = spline(s, DDt);
		end%function
	end%methods
end%classdef
