--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_CompHeat.el
-- DESCRIPTION: Defines components in which there are heat exchanges
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 19/02/2015
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Components
--------------------------------------------------------------------------------
-- Component that represents a cooling jacket
--------------------------------------------------------------------------------
COMPONENT CoolingJacket (ENUM Type_Cooling Type = Darcy)

  "Cooling jacket"

   PORTS
      IN Fluid l							"Inlet fluid port"
      OUT Fluid g							"Outlet fluid port"
		IN Heat h							"Inlet heat exchange port"
		
	DATA
		REAL L = 1.				UNITS u_m 			"Length of the channels" 
		REAL a = 0.002			UNITS u_m 			"Average width of the channels" 
		REAL b = 0.004			UNITS u_m 			"Average height of the channels" 
		INTEGER N = 100		UNITS no_units		"Number of channels"
		REAL rug = 5.e-5		UNITS u_m			"Absolute rugosity"
		REAL k_w = 370			UNITS	u_W_mK		"Thermal conductivity of the nozzle wall material"
		REAL t_w = 0.003		UNITS u_m 			"Nozzle wall thickness" 
		REAL dp = 1.e5			UNITS u_Pa			"Imposed total pressure drop if Type=Known_dp"
		
		REAL T_ch=3180.93		UNITS u_K			"Characteristic Temperature"
		REAL rho_ch=0.567		UNITS u_kg_m3		"Characteristic Density"
		REAL mfr_ch=19.05		UNITS u_kg_s		"Characteristic Mass flow"
		REAL temp_ch=0.032	UNITS u_s			"Characteristic filling Time"
		
	DECLS
		REAL Q							UNITS u_W 				"Heat flux"
		REAL T_w_hot					UNITS u_K				"Wall temperature of the combustion gases side"				
		REAL T_w_cold					UNITS u_K				"Wall temperature of the cooling liquid side"	
		REAL A_wet_cooling			UNITS u_m2				"Cooling wet area"
		REAL A_wet_nozzle				UNITS u_m2				"Nozzle wet area"
		REAL h_l							UNITS u_W_m2K			"Cooling liquid heat transfer coefficient"
		REAL Pr							UNITS no_units			"Cooling liquid Prandtl number"
		REAL Nu							UNITS no_units			"Cooling liquid Nusselt number"
		REAL Re							UNITS no_units			"Cooling liquid Reynolds number"
		REAL v							UNITS u_m_s				"Cooling liquid speed"
		REAL f							UNITS no_units			"Darcy friction factor"
		REAL D_eq						UNITS u_m				"Circular equivalent diameter of a rectangular duct for equal friction and flow capacity"
		REAL D_hy						UNITS u_m				"Hydraulic diameter"
		
		DISCR REAL k_1					UNITS no_units			"Inverso del gasto característico"
		DISCR REAL k_2					UNITS no_units			"Inverso de la entalpía característica"
		REAL rho_trans					UNITS u_kg_m3			"Gas density calculated in transient model"
		
	INIT PRIORITY 90
		Init_fluid(Vaporisation(l.fluid), g.fluid)
	
		T_w_cold = 500.
		
		k_1 = 1./mfr_ch
		k_2 = 1./(mfr_ch*T_ch)
		
	DISCRETE
		ASSERT (Liquid == State(l.fluid)) FATAL "ONLY LIQUIDS CAN ENTER TO THE COOLING JACKET!"

	CONTINUOUS
		Init_fluid(Vaporisation(l.fluid), g.fluid)
		--g.W = l.W	
		

		(temp_ch/(rho_ch*k_1)) * rho_trans' = l.W - g.W 
				
		h.Q = Q
		h.T = T_w_hot
		h.A = A_wet_nozzle
		
		A_wet_cooling = N * 2 * (a + b) * L
		
		Q = h_l * (T_w_cold - l.Tt) * A_wet_cooling
		Q = k_w / t_w * (T_w_hot - T_w_cold) * A_wet_nozzle
		--Q = l.W * cp(g.fluid) * (g.Tt - l.Tt)
		
		(temp_ch/(rho_ch*k_2*T_ch)) * (rho_trans*cv(g.fluid)*g.Tt'+g.Tt*cv(g.fluid)*rho_trans') = Q + cp(g.fluid) * (l.W * l.Tt - g.W * g.Tt)
		g.pt=rho_trans*R(g.fluid)*g.Tt
				
		h_l = Nu * cond(l.fluid) / D_hy
		Re = rho(l.fluid) * v * D_hy / visc(l.fluid)
		Pr = visc(l.fluid) * cp(l.fluid) / cond(l.fluid)
		
		-- Sieder-Tate correlation for turbulent flow
		Nu = 0.027 * Re**0.8 * Pr**0.33	
		
		v = l.W / (a * b * rho(l.fluid)) / N
		
		IF (Type == Known_dp) INSERT
			g.pt = l.pt - dp
		ELSE
			g.pt = l.pt - f * L / D_eq * 0.5 *  rho(l.fluid) * v**2
		END IF
		f = hdc_fric(D_eq, rug, Re)
		
		D_hy = 2. * a * b / (a + b)
		
		-- Huebscher (1948)
		D_eq = 1.3 * (a * b)**0.625 / (a + b)**0.25
		
END COMPONENT

