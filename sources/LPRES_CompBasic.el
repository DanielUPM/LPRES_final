--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_CompBasic.el
-- DESCRIPTION: Defines basic components
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Abstract components
--------------------------------------------------------------------------------
-- Abstract component for definition of components with a fluid inlet and a
-- fluid outlet
--------------------------------------------------------------------------------
ABSTRACT COMPONENT FluidInFluidOut

  "Abstract component for definition of components with a fluid inlet and a fluid outlet"

   PORTS
      IN Fluid f_in						"Inlet fluid port"
      OUT Fluid f_out					"Outlet fluid port"

   CONTINUOUS	
		f_in.W = f_out.W
		f_in.fluid = f_out.fluid

END COMPONENT

--------------------------------------------------------------------------------
-- Abstract component for definition of gas turbomachinery components
--------------------------------------------------------------------------------
ABSTRACT COMPONENT GasTurbo IS_A FluidInFluidOut

   "Abstract component for definition of gas turbomachinery components"
	
	PORTS      
      IN Mechanical m				"Inlet mechanical port"

   DECLS
      REAL Power			UNITS u_W			"Mechanical power"	

   CONTINUOUS
	   Power = f_in.W * cp(f_in.fluid) * (f_in.Tt - f_out.Tt)
		
		m.Power = Power

END COMPONENT


-- Components
--------------------------------------------------------------------------------
-- Component that simulates a liquid discharge to the atmosphere
--------------------------------------------------------------------------------
COMPONENT Ambient (ENUM OnOffDesign Type = Design)

	"Liquid discharge to the atmosphere"
	
	PORTS
		IN Fluid l									"Inlet fluid port"	
		
	DATA
		REAL A = 0.01						UNITS u_m2		"Discharge area"
		
	DECLS		
		REAL p_amb 							UNITS u_Pa		"Ambient pressure"
		REAL A_d								UNITS u_m2		"Design discharge area"
		
	DISCRETE
		ASSERT (Liquid == State(l.fluid)) FATAL "ONLY LIQUIDS CAN GO THROUGH THE AMBIENT!"
						
	CONTINUOUS
		p_amb = ISA_Pressure(Altitude)
					
		IF (Type == Design) INSERT
			A_d = l.W / sqrt(2. * rho(l.fluid) * (l.pt - p_amb))	
		END IF
					
		IF (Type == Off_design) INSERT					
			l.pt = p_amb + (l.W / A)**2 / (2. * rho(l.fluid))	
						
			A_d = A
		END IF
			
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a turbine with a choked inlet
--------------------------------------------------------------------------------
COMPONENT Turbine_ch IS_A GasTurbo (ENUM Type_Turbines Type = Known_pi)

   "Turbine with a choked inlet"
		
	DATA
		REAL eta_d = 0.5			UNITS no_units		"Design efficiency"	
		REAL alpha_2 = 45.		UNITS u_deg			"Flow angle in section 2"	
		REAL A_in = 0.001			UNITS u_m2			"Input area"
		REAL rpm	= 30000.			UNITS "rpm"			"Design rotational speed [rpm]"
		REAL pi = 10.				UNITS no_units		"Design expansion ratio"
		REAL W = 5.					UNITS u_kg_s		"Design mass flow"	
		
	DECLS
		REAL A_in_d			UNITS u_m2			"Design input area"
		REAL eta				UNITS no_units		"Efficiency"
		REAL alpha			UNITS no_units		"Total temperature ratio"
		
	DISCRETE
		ASSERT (Gas == State(f_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE TURBINE!"
	
	CONTINUOUS		
		alpha = 1. - eta * (1. - (f_out.pt / f_in.pt)**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)))		
		alpha = f_out.Tt / f_in.Tt
		
		cos(alpha_2 * PI / 180.) * FGAMMA(f_in.fluid) = f_in.W * sqrt(f_in.Tt * R(f_in.fluid)) / (f_in.pt * A_in_d) 
		
		eta = eta_d		
		
		IF (Type == Off_design) INSERT
			A_in_d = A_in
		END IF
		
		IF (Type == Known_pi) INSERT
			m.N = rpm * 2. * PI / 60.
			pi = f_in.pt / f_out.pt		
		END IF
		
		IF (Type == Known_W) INSERT
			m.N = rpm * 2. * PI / 60.
			f_in.W = W
		END IF
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a turbine (the inlet may be choked or not)
--------------------------------------------------------------------------------
COMPONENT Turbine IS_A GasTurbo (ENUM Type_Turbines Type = Known_pi, ENUM AngCoef Type_AC = Coefficients)

   "Turbine (the inlet may be choked or not)"
		
	DATA
		REAL eta_d = 0.8					UNITS no_units		"Design efficiency"
		REAL phi_d = 0.05					UNITS no_units		"Design flow coefficient if Type_AC=Coefficients"
		REAL psi_d = 0.7					UNITS no_units		"Design loading coefficient if Type_AC=Coefficients"
		REAL alpha_2 = 45.				UNITS u_deg			"Flow angle in section 2"
		REAL alpha_4r = -30.				UNITS u_deg			"Relative flow angle in section 4 if Type_AC=Angles"
		REAL M = 1.							UNITS no_units		"Design Mach number in section 2 if Type_AC=Angles and Type=Design or initial Mach number in section 2 for iterative calculations"
		REAL A_in = 0.005					UNITS u_m2			"Input area"
		REAL r_m	= 0.01					UNITS u_m			"Average radius"
		REAL rpm = 30000.					UNITS "rpm"			"Design rotational speed [rpm]"
		REAL pi = 1.5						UNITS no_units		"Design expansion ratio"
		REAL W = 5.							UNITS u_kg_s		"Design mass flow"	
		REAL U_0 = 10000.					UNITS u_m_s			"Initial blade speed for iterative calculations"
		
	DECLS
		REAL A_in_d					UNITS u_m2			"Design input area"
		REAL r_m_d					UNITS u_m			"Design average radius"
		
		REAL eta						UNITS no_units		"Efficiency"
		REAL psi						UNITS no_units		"Loading coefficient"
		REAL phi						UNITS no_units		"Flow coefficient"
		
		ALG REAL U					UNITS u_m_s			"Blade speed"
		REAL alpha					UNITS no_units		"Total temperature ratio"	
		REAL tau						UNITS u_J_kg		"Specific work"
		
		REAL M_2						UNITS no_units		"Mach number in section 2"
		REAL V_z2					UNITS u_m_s			"Axial speed in section 2"
		REAL V_2						UNITS u_m_s			"Speed in section 2"
		
	INIT
		M_2 = M
		U = U_0
		
	DISCRETE
		ASSERT (Gas == State(f_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE TURBINE!"
	
	CONTINUOUS		
		alpha = 1. - eta * (1. - (f_out.pt / f_in.pt)**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)))		
		alpha = f_out.Tt / f_in.Tt		
		
		f_in.W * sqrt(f_in.Tt * R(f_in.fluid)) / (f_in.pt * A_in_d) = 	sqrt(gamma(f_in.fluid)) * min(M_2, 1.) * cos(alpha_2 * PI / 180.) * (1. + (gamma(f_in.fluid) - 1.) / 2. *  min(M_2, 1.)**2)**(- (gamma(f_in.fluid) + 1.) / 2. / (gamma(f_in.fluid) - 1.))
																						
		U = m.N * r_m_d
		psi = tau / U**2
		phi = V_z2 / U
		
		V_z2 = V_2 * cos(alpha_2 * PI / 180.)
		V_2 = M_2 * sqrt(gamma(f_in.fluid) * R(f_in.fluid) * f_in.Tt / (1. + (gamma(f_in.fluid) - 1.) / 2. * M_2**2))
		
		IF (Type_AC == Coefficients) INSERT
			psi = (1. + psi_d) / phi_d * phi - 1.
		ELSE
			psi = phi * (tan(alpha_2 * PI / 180.) - tan(alpha_4r * PI / 180.)) - 1.	
		END IF			
		eta = eta_d
		
		tau = Power / f_in.W
		
		IF (Type == Off_design) INSERT
			A_in_d = A_in
			r_m_d = r_m		
		END IF
		
		IF (Type == Known_pi) INSERT
			m.N = rpm * 2. * PI / 60.
			pi = f_in.pt / f_out.pt
		END IF
		
		IF (Type == Known_W) INSERT
			m.N = rpm * 2. * PI / 60.
			f_in.W = W
		END IF
		
		IF ((Type == Known_pi OR Type == Known_W) AND Type_AC == Angles) INSERT
			M_2 = min(M, 1.)
		ELSEIF ((Type == Known_pi OR Type == Known_W) AND Type_AC == Coefficients) INSERT
			phi = phi_d
		END IF		
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a compressor
--------------------------------------------------------------------------------
COMPONENT Compressor IS_A GasTurbo (ENUM OnOffDesign Type = Design, ENUM AngCoef Type_AC = Coefficients)

   "Compressor"
			
	DATA
		REAL eta_d = 0.8			UNITS no_units		"Design efficiency"
		REAL phi_d = 0.05			UNITS no_units		"Design flow coefficient"
		REAL psi_d = 0.7			UNITS no_units		"Design loading coefficient if Type_AC=Coefficients"
		REAL alpha_2r = -30.		UNITS u_deg			"Relative current angle in section 2 if Type_AC=Angles"
		REAL r_m = 0.3				UNITS u_m			"Average radius"
		REAL A_in = 0.35			UNITS u_m2			"Input area"
		REAL pi = 4.				UNITS no_units		"Design compression ratio"
		REAL U_0 = 500.			UNITS u_m_s			"Initial blade speed for iterative calculations"
		
	DECLS	
		REAL A_in_d					UNITS u_m2			"Design input area"
		REAL r_m_d					UNITS u_m			"Design average radius"
		
		REAL eta						UNITS no_units		"Efficiency"
		REAL psi						UNITS no_units		"Loading coefficient"
		REAL phi						UNITS no_units		"Flow coefficient"
		
		ALG REAL U					UNITS u_m_s			"Blade speed"
		REAL tau						UNITS u_J_kg		"Specific work"
		
		REAL rho_in					UNITS u_kg_m3		"Inlet density"
		ALG REAL M_in = 0.001	UNITS no_units		"Inlet Mach number"
		
	INIT
		U = U_0
		
	DISCRETE
		ASSERT (Gas == State(f_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE COMPRESSOR!"		
		ASSERT (M_in <= 1.) KILLPOINT "THE COMPRESSOR INLET CANNOT BE CHOKED!"
		
	CONTINUOUS
		(f_out.Tt / f_in.Tt)  = ((f_out.pt / f_in.pt)**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)) - 1.) / eta + 1.
		
		U = m.N * r_m_d		
		psi = tau / U**2
		phi = f_in.W /(A_in_d * U * rho_in)
		
		f_in.W * sqrt(R(f_in.fluid) * f_in.Tt) / A_in_d / f_in.pt = sqrt(gamma(f_in.fluid)) * M_in * (1. + (gamma(f_in.fluid) - 1.) / 2. * M_in**2)**(- (gamma(f_in.fluid) + 1.) / 2. / (gamma(f_in.fluid) - 1.))
		((f_in.pt / R(f_in.fluid) / f_in.Tt) / rho_in)**(gamma(f_in.fluid) - 1.) = 1. + (gamma(f_in.fluid) - 1.) / 2. * M_in**2		
		
		IF (Type_AC == Coefficients) INSERT
			psi = 1. - (1. - psi_d) / phi_d * phi
		ELSE
			psi = phi * tan(alpha_2r * PI / 180.) + 1.
		END IF			
		eta = eta_d
				
		tau = - Power / f_in.W
		
		IF (Type == Off_design) INSERT
			A_in_d = A_in
			r_m_d = r_m
		ELSE
			pi = f_out.pt / f_in.pt
			phi = phi_d
		END IF	
				
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a liquid turbine
--------------------------------------------------------------------------------
COMPONENT Turbine_liq IS_A FluidInFluidOut (ENUM Type_Turbine_liq Type = Known_dp)

	"Liquid turbine"
	
	PORTS      
      IN Mechanical m				"Inlet mechanical port"
		
	DATA
		REAL eta_d = 0.8			UNITS no_units		"Design efficiency"
		REAL phi_d = 0.05			UNITS no_units		"Design flow coefficient"
		REAL psi_d = 0.7			UNITS no_units		"Design loading coefficient"
		REAL A_in = 0.01			UNITS u_m2			"Input area"
		REAL r_m	= 0.05			UNITS u_m			"Average radius"
		REAL rpm = 30000.			UNITS "rpm"			"Design rotational speed [rpm]"
		REAL dp = 5000000.		UNITS u_Pa			"Design pressure decrease"
		REAL W = 5.					UNITS u_kg_s		"Design mass flow"	
		REAL U_0 = 500.			UNITS u_m_s			"Initial blade speed for iterative calculations"
		
	DECLS		
		REAL A_in_d					UNITS u_m2			"Design input area"
		REAL r_m_d					UNITS u_m			"Design average radius"
		
		REAL eta						UNITS no_units		"Efficiency"
		REAL psi						UNITS no_units		"Loading coefficient"
		REAL phi						UNITS no_units		"Flow coefficient"	
		
		ALG REAL U					UNITS u_m_s			"Blade speed"	
		REAL tau						UNITS u_J_kg		"Specific work"
      REAL Power					UNITS u_W			"Mechanical power"			
		
	INIT
		U = U_0
		
	DISCRETE
		ASSERT (Liquid == State(f_in.fluid)) FATAL "ONLY LIQUIDS CAN GO THROUGH THE LIQUID TURBINE!"
		
	CONTINUOUS
		tau = (f_in.pt - f_out.pt) * eta / rho(f_in.fluid)
		
		U = m.N * r_m_d		
		psi = tau / U**2
		phi = f_in.W /(A_in_d * U * rho(f_in.fluid))
			
		psi = (1. + psi_d) / phi_d * phi - 1.
		eta = eta_d
		
		tau = Power / f_in.W		
		
		IF (Type == Off_design) INSERT
			A_in_d = A_in
			r_m_d = r_m		
		END IF
		
		IF (Type == Known_dp) INSERT
			m.N = rpm * 2. * PI / 60.
			dp = f_in.pt - f_out.pt
			phi = phi_d
		END IF
		
		IF (Type == Known_W) INSERT
			m.N = rpm * 2. * PI / 60.
			f_in.W = W
			phi = phi_d
		END IF
		
	   (f_in.pt - f_out.pt) / rho(f_in.fluid) * (1. - eta) = cp(f_in.fluid) * (f_out.Tt - f_in.Tt)
		
		m.Power = Power
				
END COMPONENT	

--------------------------------------------------------------------------------
-- Component that represents a pump
--------------------------------------------------------------------------------
COMPONENT Pump IS_A FluidInFluidOut (ENUM OnOffDesign Type = Design)

	"Pump"
	
	PORTS      
      IN Mechanical m				"Inlet mechanical port"
		
	DATA
		REAL eta_d = 0.8			UNITS no_units		"Design efficiency"
		REAL phi_d = 0.05			UNITS no_units		"Design flow coefficient"
		REAL psi_d = 0.7			UNITS no_units		"Design loading coefficient"
		REAL A_in = 0.01			UNITS u_m2			"Input area"
		REAL r_m	= 0.05			UNITS u_m			"Average radius"
		REAL dp = 5000000.		UNITS u_Pa			"Design pressure increase"
		REAL U_0 = 500.			UNITS u_m_s			"Initial blade speed for iterative calculations"
		
	DECLS		
		REAL A_in_d					UNITS u_m2			"Design input area"
		REAL r_m_d					UNITS u_m			"Design average radius"
		
		REAL eta						UNITS no_units		"Efficiency"
		REAL psi						UNITS no_units		"Loading coefficient"
		REAL phi						UNITS no_units		"Flow coefficient"		
		
		ALG REAL U					UNITS u_m_s			"Blade speed"	
		REAL H						UNITS u_m			"Head"
		REAL tau						UNITS u_J_kg		"Specific work"
      REAL Power					UNITS u_W			"Mechanical power"			
		
	INIT
		U = U_0
		
	DISCRETE
		ASSERT (Liquid == State(f_in.fluid)) FATAL "ONLY LIQUIDS CAN GO THROUGH THE PUMP!"
		
	CONTINUOUS
		H = (f_out.pt - f_in.pt) / (g_0 * rho(f_in.fluid))
		tau = (f_out.pt - f_in.pt) / (eta * rho(f_in.fluid))
		
		U = m.N * r_m_d		
		psi = tau / U**2
		phi = f_in.W /(A_in_d * U * rho(f_in.fluid))
			
		psi = 1. - (1. - psi_d) / phi_d * phi
		eta = eta_d
		
		tau = - Power / f_in.W	
		
		IF (Type == Off_design) INSERT
			A_in_d = A_in
			r_m_d = r_m
		ELSE
			dp = f_out.pt - f_in.pt
			phi = phi_d
		END IF		
		
	   (f_out.pt - f_in.pt) / rho(f_in.fluid) * (1. / eta - 1.) = cp(f_in.fluid) * (f_out.Tt - f_in.Tt)
		
		m.Power = Power
				
END COMPONENT	

--------------------------------------------------------------------------------
-- Component that represents a liquid pipe with pressure drop
--------------------------------------------------------------------------------
COMPONENT Pipe IS_A FluidInFluidOut

   "Liquid pipe with pressure drop"
	
	DATA	
		REAL L = 1.				UNITS u_m 			"Length"	
		REAL D = 0.1			UNITS u_m 			"Diameter"	
		REAL K = 5				UNITS no_units		"Additional pressure losses"
		REAL rug	= 1.5e-6		UNITS u_m			"Absolute rugosity"
		
	DECLS
		REAL v					UNITS u_m_s				"Liquid speed"
		REAL f					UNITS no_units			"Darcy friction factor"		
		REAL Re					UNITS no_units			"Reynolds number"
		
	DISCRETE
		ASSERT (Liquid == State(f_in.fluid)) FATAL "ONLY LIQUIDS CAN GO THROUGH THE PIPE!"
		
	CONTINUOUS
		f_out.pt = f_in.pt - (f * L / D + K) * 0.5 * rho(f_in.fluid) * v**2
		f = hdc_fric(D, rug, Re)
		
		v = f_in.W / (PI * D**2 / 4 * rho(f_in.fluid))
		Re = rho(f_in.fluid) * v * D / visc(f_in.fluid)
		
		f_out.Tt = f_in.Tt
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a flow splitter with pressure drop
--------------------------------------------------------------------------------
COMPONENT SplitFrac

   "Flow splitter with pressure drop"

	PORTS
      IN Fluid f_in                 				"Inlet fluid port"
      OUT Fluid f_out               				"Outlet fluid port"
      OUT Fluid f_b           						"Branch fluid port"
		
	DATA	
		REAL TPL = 0.9			UNITS no_units			"Total pressure loss"

   CONTINUOUS	
		f_in.W = f_out.W + f_b.W
		
		f_in.Tt = f_out.Tt
		f_in.Tt = f_b.Tt
		
		f_in.pt = f_out.pt / TPL
		f_in.pt = f_b.pt / TPL 
		
		f_in.fluid = f_out.fluid
		f_in.fluid = f_b.fluid
	
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a junction with pressure drop
--------------------------------------------------------------------------------
COMPONENT Junction

   "Junction with pressure drop"

	PORTS
      IN Fluid f_in1                 				"Inlet fluid port"
		IN Fluid f_in2                				"Inlet fluid port"
      OUT Fluid f_out               				"Outlet fluid port"
		
	DATA	
		REAL TPL = 0.9			UNITS no_units			"Total pressure loss"

   CONTINUOUS	
		f_in1.W + f_in2.W = f_out.W
		
		(f_in1.W / f_out.W) * cp(f_in1.fluid) * (f_out.Tt - f_in1.Tt) + (f_in2.W / f_out.W) * cp(f_in2.fluid) * (f_out.Tt - f_in2.Tt) = 0
		
		f_in1.pt = f_out.pt / TPL
		f_in2.pt = f_out.pt / TPL
		
		f_out.fluid = f_in1.fluid
		f_out.fluid = f_in2.fluid
	
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a pressure regulator
--------------------------------------------------------------------------------
COMPONENT Regulator IS_A FluidInFluidOut (ENUM Type_Regulator Type = Design)

   "Pressure regulator"
		
	DATA
		REAL dp_min = 1500. 			UNITS u_Pa			"Minumum total pressure drop"
		REAL pt_out = 1200000.		UNITS u_Pa			"Outlet total pressure"
		
		REAL dp = 100000.				UNITS u_Pa			"Imposed total pressure drop"
		
		REAL dpr = 0.1					UNITS no_units		"Imposed total pressure drop ratio"
		
		--REAL Ag  = 0.001				UNITS u_m2			"Area in gas or liquid valve"	
		
	DECLS
		REAL dp_d						UNITS u_Pa			"Design total pressure drop"

		--REAL M_g							UNITS no_units		"Mach number in gas valve"
		--REAL rho_g						UNITS	u_kg_m3		"Gas density in gas valve"
		
   CONTINUOUS	
		f_out.Tt = f_in.Tt
		
		IF (Type == Known_pt_out) INSERT
			--M_g=0
			--rho_g=0
			f_out.pt = min(pt_out, f_in.pt - dp_min)
		END IF
		
		IF (Type == Known_dp) INSERT
			dp_d = dp
			
			--M_g=0
			--rho_g=0
			--f_out.pt = f_in.pt - dp
	
		END IF
		
		IF (Type == Known_dpr) INSERT
			--M_g=0
			--rho_g=0
			f_out.pt = f_in.pt * dpr
		END IF		
		
		--IF (Type == Gas_valve) INSERT
		--	f_in.W * sqrt(R(f_out.fluid) * f_in.Tt) / Ag / f_in.pt = sqrt(gamma(f_out.fluid)) * M_g * (1. + (gamma(f_out.fluid) - 1.) / 2. * M_g**2)**(- (gamma(f_out.fluid) + 1.) / 2. / (gamma(f_out.fluid) - 1.))
		--	((f_in.pt / R(f_in.fluid) / f_in.Tt) / rho_g)**(gamma(f_in.fluid) - 1.) = 1. + (gamma(f_in.fluid) - 1.) / 2. * M_g**2
		--	(f_in.pt / f_out.pt)**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)) = 1. + (gamma(f_in.fluid) - 1.) / 2. * M_g**2	
		--END IF
		f_out.pt = f_in.pt - dp_d
	
END COMPONENT

--------------------------------------------------------------------------------
-- Component that specifies the conditions of a fluid inlet in the system
--------------------------------------------------------------------------------
COMPONENT Inlet (ENUM Type_Inlet Type = All)

	"Conditions of a fluid inlet in the system"
	
	PORTS
		OUT Fluid f						"Outlet fluid port"
	
	DATA
		REAL Tt = 288.15						UNITS u_K				"Total temperature" 
		REAL pt = 101325.						UNITS u_Pa				"Total pressure"	
		REAL W = 1.								UNITS u_kg_s			"Mass flow"	
		ENUM LiquidsGases fluid = LOX									"Working fluid name"

	INIT PRIORITY 100
		Init_fluid(fluid, f.fluid)

	CONTINUOUS
		Init_fluid(fluid, f.fluid)
		f.pt = pt
		f.Tt = Tt
		
		IF (Type == All) INSERT
			f.W = W
		END IF

END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a liquid tank
--------------------------------------------------------------------------------
COMPONENT Tank (ENUM OnOffDesign Type = Design)

	"Liquid tank"
	
	PORTS
		IN Fluid g						"Inlet fluid port"
		OUT Fluid l						"Outlet fluid port"
	
	DATA		
		ENUM Liquids fluid_l	= LOX									"Working liquid name"
		REAL T_d = 288.15					UNITS u_K				"Tank temperature"		
		REAL A_g = 0.001					UNITS u_m2				"Gas input area"
		REAL p_d	= 1000000				UNITS u_Pa				"Design tank pressure"
		
	DECLS		
		REAL p_g							UNITS u_Pa			"Pressurisation gas pressure"
		REAL rho_g						UNITS u_kg_m3		"Pressurisation gas density"
		ALG REAL M_g = 0.1			UNITS no_units		"Pressurisation gas Mach number"
		REAL A_g_d						UNITS u_m2			"Design pressurisation gas input area"
		
	INIT PRIORITY 100
		Init_fluid(fluid_l, l.fluid)
		
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN PRESSURISE THE TANK!"
		ASSERT (M_g <= 1.) KILLPOINT "THE PRESSURISATION GAS INLET CANNOT BE CHOKED!"
			
	CONTINUOUS
		Init_fluid(fluid_l, l.fluid)
		l.Tt = T_d
		
		l.pt = p_g
		
		l.W / rho(l.fluid) = g.W / rho_g
		
		g.W * sqrt(R(g.fluid) * g.Tt) / A_g_d / g.pt = sqrt(gamma(g.fluid)) * M_g * (1. + (gamma(g.fluid) - 1.) / 2. * M_g**2)**(- (gamma(g.fluid) + 1.) / 2. / (gamma(g.fluid) - 1.))
		((g.pt / R(g.fluid) / g.Tt) / rho_g)**(gamma(g.fluid) - 1.) = 1. + (gamma(g.fluid) - 1.) / 2. * M_g**2
		(g.pt / p_g)**((gamma(g.fluid) - 1.) / gamma(g.fluid)) = 1. + (gamma(g.fluid) - 1.) / 2. * M_g**2	
		
		IF (Type == Off_design) INSERT
			A_g = A_g_d
		ELSE
			p_g = p_d
		END IF
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a liquid tank that is pressurised by the atmosphere
--------------------------------------------------------------------------------
COMPONENT TankOpen

	"Liquid tank that is pressurised by the atmosphere"
	
	PORTS
		OUT Fluid l				"Outlet fluid port"
	
	DATA
		ENUM Liquids fluid = LOX									"Working liquid name"
		REAL T_d = 288.15					UNITS u_K				"Tank temperature"
		
	DECLS		
		REAL p_d						UNITS u_Pa		"Tank pressure"
		
	INIT PRIORITY 100
		Init_fluid(fluid, l.fluid)
			
	CONTINUOUS
		Init_fluid(fluid, l.fluid)
		l.Tt = T_d
		
		l.pt = p_d
		
		p_d = ISA_Pressure(Altitude)

END COMPONENT


--------------------------------------------------------------------------------
-- Component that represents a mechanical shaft without acceleration
--------------------------------------------------------------------------------
COMPONENT Shaft

   "Mechanical shaft without acceleration"

   PORTS
      OUT Mechanical m_1				"Outlet mechanical port"
      OUT Mechanical m_2				"Outlet mechanical port"
		
	DATA
		REAL eta	= 1.	UNITS no_units		"Efficiency"
	
	DECLS
		REAL rpm			UNITS "rpm"			"Rotational speed [rpm]"

   CONTINUOUS
      m_1.N = m_2.N		
		m_1.N = rpm * 2. * PI / 60.
		
		m_1.Power = ZONE (m_2.Power > 0)	- m_2.Power * eta
						OTHERS					- m_2.Power / eta

END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a mechanical shaft with a gearbox and without 
-- acceleration
--------------------------------------------------------------------------------
COMPONENT Gearbox

   "Mechanical shaft with a gearbox and without acceleration"

   PORTS
      OUT Mechanical m_in				"Outlet mechanical port"
      OUT Mechanical m_out				"Outlet mechanical port"
		
	DATA
		REAL GR = 1.	UNITS no_units		"Gear ratio"	
		REAL eta	= 1.	UNITS no_units		"Efficiency"

   CONTINUOUS
      m_out.N = GR * m_in.N
		
      m_out.Power = 	ZONE (m_in.Power > 0) 	- m_in.Power * eta
							OTHERS					 	- m_in.Power / eta

END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a convergent nozzle with the ambient already 
-- connected
--------------------------------------------------------------------------------
COMPONENT Nozzle (ENUM OnOffDesign Type = Design)

	"Convergent nozzle with the ambient already connected"
	
	PORTS
		IN Fluid g						"Inlet fluid port"
		OUT Info (n = 1) i			"Outlet information port"
		
	DATA
		REAL A = 0.02			UNITS u_m2			"Discharge area"
		
	DECLS
		REAL p_amb 				UNITS u_Pa			"Ambient pressure"
		
		REAL A_d					UNITS u_m2			"Design discharge area"
		REAL A_sl				UNITS u_m2			"Sonic lock area"		
		
		REAL PR					UNITS no_units		"Pressure ratio"
		REAL PR_sl				UNITS no_units		"Sonic lock pressure ratio"		
		REAL p_out				UNITS u_Pa			"Outlet pressure"	
		REAL p_out_ch			UNITS u_Pa			"Choked outlet pressure"
			
		REAL M_out				UNITS no_units		"Outlet Mach number"	
		REAL T_out			   UNITS u_K			"Outlet temperature"
		REAL v_out				UNITS u_m_s			"Outlet speed"
		REAL Thrust				UNITS u_N			"Thrust"
		
	INIT
		g.pt = 10000000
		
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE NOZZLE!"	
						
	CONTINUOUS
		p_amb = ISA_Pressure(Altitude)
		
		PR = g.pt / p_amb
		PR_sl = ((gamma(g.fluid) + 1.) / 2.)**(gamma(g.fluid) / (gamma(g.fluid) - 1.))
		
		p_out_ch = g.pt / PR_sl
		
		g.W = A_sl * FGAMMA(g.fluid) * g.pt / sqrt(g.Tt * R(g.fluid))	
		
		A_sl = A_d * M_out / ((2. + (gamma(g.fluid) - 1.) * M_out**2) / (gamma(g.fluid) + 1.))**((gamma(g.fluid) + 1.) / (2. * (gamma(g.fluid) - 1.)))
		
		M_out =	ZONE (PR < PR_sl)		sqrt(2. * (PR**((gamma(g.fluid) - 1.) / gamma(g.fluid)) - 1.) / (gamma(g.fluid) - 1.))
					OTHERS					1.
					
		p_out = max(p_amb, p_out_ch)
		
		IF (Type == Off_design) INSERT		
			A_d = A						
		END IF
		
		T_out = g.Tt / (1. + (gamma(g.fluid) - 1.) / 2. * M_out**2)
		v_out = M_out * sqrt(gamma(g.fluid) * R(g.fluid) * T_out)
		
		Thrust = g.W * v_out + A_d * (p_out - p_amb)
		
		i.Data[1] = Thrust
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a choked convergent-divergent nozzle
--------------------------------------------------------------------------------
COMPONENT NozzleConDiv (ENUM OnOffDesign Type = Design)

   "Choked Convergent-divergent nozzle"
	
	PORTS
		IN Fluid g_in						"Inlet fluid port"
		OUT GasNozzle g_out				"Outlet gas port through a nozzle"
		
	DATA
		REAL AR = 10.					UNITS no_units		"Area ratio"
		REAL A_th = 0.05				UNITS u_m2			"Throat area"
		REAL W = 100.					UNITS u_kg_s		"Mass flow"
		
	DECLS
		REAL A_out						UNITS u_m2			"Output area"	
		REAL A_th_d						UNITS u_m2			"Design throat area"
				
		ALG REAL p_out_ch = 100.	UNITS u_Pa			"Choked outlet pressure"
		
	DISCRETE
		ASSERT (Gas == State(g_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE CHOKED CONVERGENT-DIVERGENT NOZZLE!"	
		
	CONTINUOUS
		g_in.fluid = g_out.fluid
		g_in.W = g_out.W
		
		g_out.A_out = A_out 
		AR = A_out / A_th_d
		
		FGAMMA(g_in.fluid) = g_in.W * sqrt(g_in.Tt * R(g_in.fluid)) / (g_in.pt * A_th_d) 
		
		IF (Type == Off_design) INSERT
			A_th = A_th_d
		ELSE
			g_in.W = W
		END IF
		
		AR = FGAMMA(g_out.fluid) / ((p_out_ch / g_in.pt)**(1. / gamma(g_out.fluid)) * sqrt(2 * gamma(g_out.fluid) * (1 - (p_out_ch / g_in.pt)**((gamma(g_out.fluid) - 1.) / gamma(g_out.fluid))) / (gamma(g_out.fluid) - 1)))
		
		g_out.pt = g_in.pt 
		g_out.Tt = g_in.Tt
	
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a nozzle extension
--------------------------------------------------------------------------------
COMPONENT NozzleExt
	
	"Nozzle extension"
	
	PORTS
		IN GasNozzle g_in				"Inlet gas port through a nozzle"
		OUT GasNozzle g_out			"Outlet gas port through a nozzle"
			
	DATA
		REAL AR = 2						UNITS no_units		"Extension area ratio"
		
	DECLS
		REAL A_in						UNITS u_m2			"Extension input area"
		REAL A_out 						UNITS u_m2			"Extension output area"
		
	DISCRETE
		ASSERT (Gas == State(g_in.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE NOZZLE EXTENSION!"
		
	CONTINUOUS
		g_in.fluid = g_out.fluid
		g_in.W = g_out.W

		g_in.A_out = A_in
		g_out.A_out = A_out
		AR = A_out / A_in
		
		g_out.pt = g_in.pt 
		g_out.Tt = g_in.Tt 
		
END COMPONENT


-- Measuring components
--------------------------------------------------------------------------------
-- Component that represents a monitor to measure the thrust
--------------------------------------------------------------------------------
COMPONENT ThrustMonitor

	"Monitor to measure the thrust"
	
	PORTS
		IN GasNozzle g						"Inlet gas port through a nozzle"
		OUT Info (n = 1) i				"Outlet information port"
		
	DECLS
		REAL Thrust						UNITS u_N			"Thrust"
		
		REAL p_amb 						UNITS u_Pa			"Ambient pressure"
		
		REAL T_out						UNITS u_K			"Outlet temperature"			
		REAL p_out						UNITS u_Pa			"Outlet pressure"	
		REAL v_out						UNITS u_m_s			"Outlet speed"
		REAL A_out						UNITS u_m2			"Output area"
		ALG REAL M_out = 100.		UNITS no_units		"Gas Mach number"
		
	DISCRETE
		ASSERT (Gas == State(g.fluid)) FATAL "ONLY GASES CAN GO THROUGH THE THRUST MONITOR!"
		ASSERT (M_out >= 1.) KILLPOINT "THE CONVERGENT-DIVERGENT NOZZLE MUST BE CHOKED!"
		ASSERT (p_out / p_amb > (1.88 * M_out - 1)**(-0.64)) WARNING "ACCORDING TO SCHMUCKER CRITERION, THE NOZZLE HAS A REGION OF DETACHMENT!"
						
	CONTINUOUS
		g.A_out = A_out
		p_amb = ISA_Pressure(Altitude)
		
		g.W * sqrt(R(g.fluid) * g.Tt) / A_out / g.pt = sqrt(gamma(g.fluid)) * M_out * (1. + (gamma(g.fluid) - 1.) / 2. * M_out**2)**(- (gamma(g.fluid) + 1.) / 2. / (gamma(g.fluid) - 1.))
		g.Tt / T_out = 1. + (gamma(g.fluid) - 1.) / 2. * M_out**2
		(g.pt / p_out)**((gamma(g.fluid) - 1.) / gamma(g.fluid)) = 1. + (gamma(g.fluid) - 1.) / 2. * M_out**2
		
		v_out = M_out * sqrt(gamma(g.fluid) * R(g.fluid) * T_out)	
		
		Thrust = g.W * v_out + A_out * (p_out - p_amb)
		
		i.Data[1] = Thrust

END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a mass flow meter
--------------------------------------------------------------------------------
COMPONENT FlowMeter IS_A FluidInFluidOut

	"Mass flow meter"
	
	PORTS
		OUT Info (n = 1) i							"Outlet information port"
						
	CONTINUOUS
		f_in.pt = f_out.pt
		f_in.Tt = f_out.Tt
		
		i.Data[1] = abs(f_in.W)

END COMPONENT

--------------------------------------------------------------------------------
-- Component that performs the calculations with the measurements done by other
-- components
--------------------------------------------------------------------------------
COMPONENT ControlPanel

	"Calculations with the measurements done by other components"
	
	PORTS
		IN Info (n = 1) i_W				"Inlet information port"
		IN Info (n = 1) i_Thrust		"Inlet information port"
		IN Info (n = 1) i_Comb			"Inlet information port"
		
	DECLS
		REAL Thrust				UNITS u_N			"Thrust"
		REAL W_tot 				UNITS u_kg_s		"Total mass flow"
		REAL c_star				UNITS u_m_s			"Characteristic exhaust velocity"
		
		REAL Isp					UNITS u_m_s			"Specific impulse [m/s]"
		REAL Isp_0				UNITS u_s			"Specific impulse [s]"
		REAL C_E					UNITS no_units		"Thrust coefficient"
						
	CONTINUOUS		
		i_Thrust.Data[1] = Thrust
		i_W.Data[1] = W_tot
		i_Comb.Data[1] = c_star
		
		Isp = Thrust / W_tot
		Isp = C_E * c_star
		Isp = g_0 * Isp_0

END COMPONENT