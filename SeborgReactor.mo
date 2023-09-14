package SeborgReactor
	// Real time demo of Modelica
	// BL, November 4, 2016
	//
	//
	model SimReactor
		// Simulation set-up for SeborgReactor
		// BL, November 4, 2016
		//
		// Instantiate reactor model
		ModReactor reactor;
		ModReactor react1;
		ModReactor react2(V=101);
		//
		// Define variables
		Real _Vd1;	// input function
		Real _Vd2;
		Real _cAi;	// input function
		Real _Ti;	// input function
		Real _Tci;	// input function
		Real _cA;	// output
		Real _cA1;
		Real _cA2;
		//
	equation
		_Vd1 = if time < 5 then 100 else 90;
		_Vd2 = if time < 5 then 100 else 105;
		_cAi = 1;
		_Ti = 350;
		_Tci = 300;
		//
		reactor.Vd = _Vd1;
		reactor.cAi = _cAi;
		reactor.Ti = _Ti;
		reactor.Tci = _Tci;
		react1.Vd = _Vd2;
		react1.cAi = _cAi;
		react1.Ti = _Ti;
		react1.Tci = _Tci;
		react2.Vd = _Vd2;
		react2.cAi = _cAi;
		react2.Ti = _Ti;
		react2.Tci = _Tci;
		//
		_cA = reactor.cA;
		_cA1 = react1.cA;
		_cA2 = react2.cA;
	//
	end SimReactor;
	//
	model ModReactor
		// Model of SeborgReactor
		// BL, November 4, 2016
		//
		// Model parameters
		parameter Real V = 100 "Reactor volume, L";
		parameter Real k0 = exp(8750/350) "Preexponential factor, min-1";
		parameter Real EdR = 8750 "Activation 'temperature', K";
		parameter Real dHtr = -5e4 "Molar enthalpy of reaction, J/mol";
		parameter Real rho = 1000 "Density, g/L";
		parameter Real chp = 0.239 "Specific heat capacity, J/(g.K)";
		parameter Real UA = 5e4 "Heat transfer parameter, J/(min.K)";
		// Initial value parameters
		parameter Real cA0 = 0.5 "Initial concentration, mol/L";
		parameter Real nA0 = cA0*V "Initial amount of A, mol";
		parameter Real T0 = 350 "Initial reactor temperature, K";
		// Variables
		Real nA "Amount of A in reactor, mol";
		Real T "Reactor temperature, K";
		Real r "Reaction rate";
		Real Qd "Added heat flow";
		// Inputs
		input Real Vd "Influent volumetric flow rate, L/min";
		input Real cAi "Influent concentrationn, mol/L";
		input Real Ti "Influent reactor temperature, K";
		input Real Tci "Influent cooling temperature, K";
		// Outputs
		output Real cA "Reactor concentration of A, mol/L";
		//
	// Initialize model
	initial equation
		nA = nA0;	// Initializing nA
		T = T0;		// Initializing temperature
	// Equations
	equation
		der(nA) = Vd*(cAi - nA/V) - r*V;
		rho*V*chp*der(T) = rho*chp*Vd*(Ti-T) + (-dHtr)*r*V + Qd;
		r = k0*exp(-EdR/T)*nA/V;
		Qd = UA*(Tci-T);
		//
		cA = nA/V;		
	//
	end ModReactor;
//
end SeborgReactor;