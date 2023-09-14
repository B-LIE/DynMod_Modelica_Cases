package ProjectHydroPower
	// Sept 2023
	model SimHydroPower
		//
		// Simulation experiment of hydro power system
		// Solution to project in course
		//		FM1015 Modeling of Dynamic Systems
		//		Fall semester 2016
		//
		//		Bernt Lie
		//		University College of Southeast Norway
		//		September 26, 2016
		//		modification of code by Liubomyr Vytvytskyi
		//
		ModHydroPower mhp;	// Instantiating model of hydro power system
		//
		Real _u_v;
		Real _Wd_e;
		Real _H_t;
		Real _p_i1, _p_m, _p_p2, _p_tr2, _p_d2;
		Real _w_g;
		Real _h_s, _Vd_i, _Vd_s, _Vd_p;
	equation
		_u_v = if time < 600 then 0.5 else 0.47;
		_Wd_e = if time < 1200 then 80 else 77;
		_H_t = if time < 1800 then 10 else 20;
		//
		mhp.u_v = _u_v;
		mhp.Wd_e = _Wd_e*1e6;
		mhp.H_t = _H_t;
		_p_i1 = mhp.p_i1/mhp.p_a;
		_p_m = mhp.p_m/mhp.p_a;
		_p_p2 = mhp.p_p2/mhp.p_a;
		_p_tr2 = mhp.p_tr2/mhp.p_a;
		_p_d2 = mhp.p_d2/mhp.p_a;
		_w_g = 4*mhp.w_a/2/mhp.PI;
		_h_s = mhp.h_s;
		_Vd_i = mhp.Vd_i;
		_Vd_s = mhp.Vd_s;
		_Vd_p = mhp.Vd_p;
	end SimHydroPower;
	//
	model ModHydroPower
		//
		// Model of hydro power system
		// Solution to project in course
		//		FM1015 Modeling of Dynamic Systems
		//		Fall semester 2016
		//
		//		Bernt Lie
		//		University College of Southeast Norway
		//		September 26, 2016
		//		modification of code by Liubomyr Vytvytskyi
		//
		// Constants
		constant Real g = 9.81 "Acceleration of gravity, m/s2";
		constant Real PI = 3.141592654 "pi, ";
		// General parameters
		parameter Real rho = 997 "Water density, kg/m3";
		parameter Real mu = 8.9e-4 "Viscosity of water, Pa.s";
		parameter Real epsilon = 1.5e-5 "Pipe roughness height, m";
		parameter Real H_r = 50 "Reservoir level above intake, m";
		parameter Real C_v = 6 "Guide vane capacity, m3/s";
		parameter Real p_a = 1.013e5 "Atmospheric pressure, Pa";
		// Rated values
		parameter Real Vd_0 = 20 "Nominal flow rate, m3/s";
		parameter Real w_a_0 = 50*2*PI/4 "Nominal angular velocity, rad/s";
		// Intake race
		parameter Real H_i = 25 "Vertical drop of intake race, m";
		parameter Real L_i = 6600 "Length of intake race, m";
		parameter Real D_i = 5.8 "Diameter of intake race, m";
		parameter Real A_i = PI*D_i^2/4 "Cross section area of intake race, m2";
		parameter Real V_i = A_i*L_i "Volume of intake race, m3";
		parameter Real m_i = rho*V_i "Mass of intake race, kg";
		parameter Real A_wi = PI*D_i*L_i "Wetting surface of intake race, m2";
		parameter Real p_i1 = p_a + rho*g*H_r "Intake race pressure, Pa";
		// Surge tank
		parameter Real H_s = 120 "Vertical drop of surge tank, m";
		parameter Real L_s = 140 "Length of surge tank, m";
		parameter Real D_s = 3.4 "Diameter of the surge tank, m";
		parameter Real A_s = PI*D_s^2/4 "Cross section area of surge tank, m2";
		parameter Real p_s2 = p_a "Outlet pressure of surge tank, Pa";
		// Penstock
		parameter Real H_p = 420 "Vertical drop of penstock, m";
		parameter Real L_p = 600 "Length of penstock, m";
		parameter Real D_p = 3.3 "Diameter of penstock, m";
		parameter Real A_p = PI*D_p^2/4 "Cross section area of penstock, m2";
		parameter Real V_p = A_p*L_p "Volume of penstock, m3";
		parameter Real m_p = rho*V_p "Mass of penstock, kg";
		parameter Real A_wp = PI*D_p*L_p "Wetting surface of penstock, m2";
		// Discharge
		parameter Real H_d = 5 "Vertical drop of discharge race, m";
		parameter Real L_d = 600 "Length of discharge race, m";
		parameter Real D_d = 5.8 "Diameter of discharge race, m";
		parameter Real A_d = PI*D_d^2/4 "Cross section area of discharge race, m2";
		parameter Real V_d = A_d*L_d "Volume of discharge race, m3";
		parameter Real m_d = rho*V_d "Mass of discharge race, kg";
		parameter Real A_wd = PI*D_d*L_d "Wetting surface of discharge race, m2";
		// Aggregate
		parameter Real J_a = 2e5 "Moment of inertia of aggregate, kg.m2";
		parameter Real k_ba = 1e3 "Friction factor in the aggregate bearing box, W-s3/rad3";
		parameter Real eta_e = 0.99 "Electricity generator efficiency";
		// Turbine
		parameter Real eta_h = 0.90 "Turbine hydraulic efficiency";
		// Variables
		Real m_s "Mass in surge tank, kg";
		Real V_s "Volume in surge tank, m3";
		Real ell_s "Length of water string in surge tank, m";
		Real h_s "Level in surge tank, m";
		Real A_ws "Wetting surface of surge tank string, m2";
		//
		Real K3p_i, K3p_s, K3p_p, K3p_d;// Volumetric kinetic energies
		Real md_i, md_s, md_p, md_d;	// mass flow rates
		Real Vd_i, Vd_s, Vd_p, Vd_d;	// Volumetric flow rates
		Real v_i, v_p, v_s, v_d;		// linear velocities
		Real N_Re_s, N_Re_i, N_Re_p, N_Re_d;	// Reynold numbers
		Real f_Di, f_Ds, f_Dp, f_Dd;	// Darcy friction coefficients
		Real F_fs, F_fi, F_fp, F_fd;	// Friction forces
		Real F_gs, F_gi, F_gp, F_gd;	// Gravity forces
		Real M_i, M_s, M_p, M_d;		// Linear momenta
		Real Md_s;						// Linear momentum flow into/out of surge tank
		Real p_m "Pressure in manifold, Pa";
		Real p_i2, p_s1, p_p1, p_p2, p_tr2, p_d1, p_d2;	// Pipe end pressures
		Real dp_tr "Pressure drop across turbine, Pa";
		Real K_a;						// Aggregate kinetic energy
		// Real Kd_p, Kd_tr2, Kd_d1;		// Kinetic flow rates
		Real Wd_fa, Wd_g, Wd_ts;		// Work powers
		Real w_a;						// Angular velocity of aggregate
		// Inputs
		input Real u_v "Input guide vane signal, -";
		input Real Wd_e "Electrical power usage, W";
		input Real H_t "Tail water level, m";
	initial equation
		Vd_i = Vd_0;
		Vd_p = Vd_0;
		w_a = w_a_0;
		h_s = H_r + H_i;
	equation
		// Surge tank liquid string
		m_s = rho*V_s;
		V_s = A_s*ell_s;
		h_s = ell_s*H_s/L_s;
		A_ws = PI*D_s*ell_s "Wetting surface of surge tank string, m2";
		// Discharge counter pressure
		p_d2 = p_a + rho * g * H_t;
		// Mass flow rates
		md_i = rho*Vd_i;
		md_s = rho*Vd_s;
		md_p = rho*Vd_p;
		md_d = rho*Vd_d;
		// Linear velocities
		Vd_i = A_i*v_i;
		Vd_s = A_i*v_s;
		Vd_p = A_p*v_p;
		Vd_d = A_d*v_d;
		// Manifold
		p_i2 = p_m;
		p_s1 = p_m;
		p_p1 = p_m;
		Vd_i = Vd_s + Vd_p;
		// Penstock-Discharge string
		Vd_d = Vd_p;
		// Momentums
		M_i = m_i*v_i;
		M_s = m_s*v_s;
		M_p = m_p*v_p;
		M_d = m_d*v_d;
		// Momentum flow rate into/out of surge tank
		Md_s = rho*Vd_s*abs(Vd_s)/A_s; 
		// Gravitational forces
		F_gi = m_i*g;
		F_gs = m_s*g;
		F_gp = m_p*g;
		F_gd = m_d*g;
		// Reynold's numbers
		N_Re_i = rho*abs(v_i)*D_i/mu;
		N_Re_s = rho*abs(v_s)*D_s/mu;
		N_Re_p = rho*abs(v_p)*D_p/mu;
		N_Re_d = rho*abs(v_d)*D_d/mu;
		// Darcy friction coefficients 
		f_Di = fDarcy(N_Re_i, D_i, epsilon);
		f_Ds = fDarcy(N_Re_s, D_s, epsilon);
		f_Dp = fDarcy(N_Re_p, D_p, epsilon);
		f_Dd = fDarcy(N_Re_d, D_d, epsilon);
		// Volumetric kinetic energies
		K3p_i = rho*v_i*abs(v_i)/2;
		K3p_s = rho*v_s*abs(v_s)/2;
		K3p_p = rho*v_p*abs(v_p)/2;
		K3p_d = rho*v_d*abs(v_d)/2;
		// Friction forces
		F_fi = K3p_i*A_wi*f_Di/4;
		F_fs = K3p_s*A_ws*f_Ds/4;
		F_fp = K3p_p*A_wp*f_Dp/4;
		F_fd = K3p_d*A_wd*f_Dd/4;
		// Turbine
		Vd_p = C_v*u_v*sqrt(dp_tr/p_a);
		dp_tr = p_p2 - p_tr2;
		Wd_ts = eta_h*dp_tr*Vd_p;
		// Aggregate
		K_a = J_a*w_a^2/2;
		Wd_e = eta_e*Wd_g;
		Wd_fa = k_ba*w_a^2/2;
		p_d1 = p_tr2;
		// Balance laws
		der(M_i) = (p_i1 - p_i2)*A_i + F_gi*H_i/L_i - F_fi;
		der(m_s) = md_s;
		der(M_s) = Md_s + (p_s1 - p_s2)*A_s - F_gs*H_s/L_s - F_fs;
		der(M_p) = (p_p1-p_p2)*A_p + F_gp*H_p/L_p - F_fp;
		der(M_d) = (p_d1-p_d2)*A_d + F_gd*H_d/L_d - F_fd;
		der(K_a) = Wd_ts - Wd_fa - Wd_g;
		//
	end ModHydroPower;
	//
	function fDarcy "Darcy friction factor"
		// Function for computing Darcy's friction factor
		// author:	Bernt Lie
		//			University College of Southeast Norway
		//			September 26, 2016
		//
		// Function input arguments
		input Real N_Re "Reynold number, -";
		input Real D "Pipe diameter, m";
		input Real epsilon "Pipe roughness height, m";
		// Function output (response) value
		output Real fD "Darcy friction factor, -";
		// Local (protected) quantities
	protected
		Real arg;
	// Algorithm for computing specific enthalpy
	algorithm
		arg := epsilon/3.7/D + 5.74/(N_Re + 1e-3)^0.9;
		if arg <= 0 then
			fD :=0;
		else
			fD := 1/(2*log10(arg))^2;
		end if;
	end fDarcy;
	//
end ProjectHydroPower;