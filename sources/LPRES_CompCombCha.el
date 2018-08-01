--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_CompGasGen.el
-- DESCRIPTION: Defines combustion chamber type components with their inyectors
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Components
--------------------------------------------------------------------------------
-- Component that represents an injector
--------------------------------------------------------------------------------
COMPONENT Injector (ENUM OnOffDesign Type = Design)

   "Injector"
	
	PORTS
      IN Fluid f_in						"Inlet fluid port"
      OUT FluidInj f_out				"Outlet fluid injection port"
	
	DATA
		REAL C_D = 0.5			UNITS no_units    	"Discharge coefficient (used only in liquid state)"
		REAL A = 0.05		 	UNITS u_m2				"Output area"
		REAL W = 100.			UNITS u_kg_s			"Mass flow"
			
	DECLS
		REAL A_d					UNITS u_m2			"Design output area"
		REAL A_sl				UNITS u_m2			"Sonic lock area (calculated only for gases)"
		
		REAL PR = 10			UNITS no_units		"Pressure ratio"
		REAL PR_sl				UNITS no_units		"Sonic lock pressure ratio (calculated only for gases)"
			
		REAL M_out				UNITS no_units		"Outlet Mach number (calculated only for gases)"		
		
		REAL p_out_ch			UNITS u_Pa			"Choked outlet pressure (calculated only for gases)"		
		
		REAL v_ideal			UNITS u_m_s			"Ideal outlet speed"
		REAL Re					UNITS no_units		"Outlet Reynolds number (calculated only for liquids)"
		

	CONTINUOUS
	
		f_in.fluid = f_out.fluid
		f_in.W = f_out.W
				
		
		PR = f_in.pt / f_out.p_c
		PR_sl = 		IF (State(f_in.fluid) == Gas) 	((gamma(f_in.fluid) + 1.) / 2.)**(gamma(f_in.fluid) / (gamma(f_in.fluid) - 1.))
						ELSE										0
					
		p_out_ch = 	IF (State(f_in.fluid) == Gas) 	f_in.pt / PR_sl	
						ELSE	 									0
						
		A_sl = 		IF (State(f_in.fluid) == Gas) 	sqrt(f_in.Tt * R(f_in.fluid)) * f_in.W / (FGAMMA(f_in.fluid) * f_in.pt)
						ELSE	 									0	
				
		M_out =		ZONE (State(f_in.fluid) == Gas AND PR < PR_sl)		sqrt(2. * (PR**((gamma(f_in.fluid) - 1.) / gamma(f_in.fluid)) - 1.) / (gamma(f_in.fluid) - 1.))
						ZONE (State(f_in.fluid) == Gas AND PR >= PR_sl)		1.
						OTHERS 									0		
					
		f_out.p = 	IF (State(f_in.fluid) == Gas) 	max(f_out.p_c, p_out_ch)
						ELSE	 									f_out.p_c
						
		f_out.T = 	IF (State(f_in.fluid) == Gas) 	f_in.Tt / (1. + (gamma(f_in.fluid) - 1.) / 2. * M_out**2)
						ELSE	 									f_in.Tt
						
		IF (Type == Design) INSERT
		
		f_in.W = W
		
		A_d = 		IF (State(f_in.fluid) == Gas) 	A_sl * ((2. + (gamma(f_in.fluid) - 1.) * M_out**2) / (gamma(f_in.fluid) + 1.))**((gamma(f_in.fluid) + 1.) / (2. * (gamma(f_in.fluid) - 1.))) / M_out
						ELSE 										f_in.W / C_D / sqrt(2. * rho(f_in.fluid) * (f_in.pt - f_out.p))	
		
		ELSE
		
		A_d = A
		
		f_in.W = 	IF (State(f_in.fluid) == Gas) 	A * FGAMMA(f_in.fluid) * f_in.pt / sqrt(f_in.Tt * R(f_in.fluid)) / (((2. + (gamma(f_in.fluid) - 1.) * M_out**2) / (gamma(f_in.fluid) + 1.))**((gamma(f_in.fluid) + 1.) / (2. * (gamma(f_in.fluid) - 1.))) / M_out)
						ELSE	 									A_d * C_D * sqrt(2. * rho(f_in.fluid) * (f_in.pt - f_out.p))	
						
		END IF
		
		v_ideal =  	IF (State(f_in.fluid) == Gas) 	M_out * sqrt(gamma(f_in.fluid) * R(f_in.fluid) * f_out.T)
						ELSE	 									f_in.W / (rho(f_in.fluid) * A_d * C_D)
						
		Re =  		IF (State(f_in.fluid) == Gas) 	0
						ELSE										rho(f_in.fluid) * v_ideal * C_D * sqrt(4. * A_d / PI) / visc(f_in.fluid)
		
END COMPONENT

--------------------------------------------------------------------------------
-- Component that represents a combustion chamber
--------------------------------------------------------------------------------
COMPONENT CombCha (ENUM OnOffDesign Type = Design, ENUM YesNo Cooled = No)

	"Combustion chamber"

	PORTS
      IN FluidInj f_O                			"Inlet fluid injection port"
      IN FluidInj f_F              				"Inlet fluid injection port"
		OUT GasNozzle g								"Outlet gas port through a nozzle"
		OUT Info (n = 1) i							"Outlet information port"
		OUT Heat h										"Outlet heat exchange port"
		
	DATA
		REAL eta_d = 0.9						UNITS no_units		"Design combustion efficiency"
		REAL OF_st = 8.						UNITS no_units		"Stoichiometric mixture ratio"
		REAL Q_comb	= 2000000.				UNITS u_J_kg		"Heat of combustion per oxidant mass flow unit"
		REAL cp_P = 4182.						UNITS u_J_kgK		"Specific heat at constant pressure of the products using a stoichiometric mixture"
		REAL M_P	= 32.							UNITS u_g_mol		"Molar mass of the products using a stoichiometric mixture"			
		REAL AR = 10.							UNITS no_units		"Area ratio"	
		REAL A_th = 0.05						UNITS u_m2			"Throat area"
		REAL p_c	= 5000000.					UNITS u_Pa			"Design combustion pressure"
		REAL AR_r = 10. / 2.					UNITS no_units		"Area at the characteristic section of heat exchange divided by the throat area"			
		REAL A_wet = 1.						UNITS u_m2			"Nozzle wet area of the cooled zone"
		ENUM ConDiv Zone = Divergent								"Convergent if the characteristic section of heat exchange is placed in the convergent zone of the nozzle, Divergent if it is placed in the divergent zone"
		REAL p_c0 = 5000000.					UNITS u_Pa			"Initial combustion pressure for iterative calculations"
		REAL T_c0 = 4000.						UNITS u_K			"Initial combustion temperature for iterative calculations"
		REAL W_F0 = 1000.						UNITS u_kg_s		"Initial fuel mass flow for iterative calculations"
		
		REAL M_oxid = 31.9988				UNITS u_g_mol		"Masa molar oxidante escogido"				
		REAL M_fuel = 2.01594				UNITS u_g_mol		"Masa molar fuel escogido"
				
		REAL T_ch=3000.						UNITS u_K			"Characteristic Temperature"
		REAL rho_ch=1.							UNITS u_kg_m3		"Characteristic Density"
		REAL mfr_ch=10.							UNITS u_kg_s		"Characteristic Mass flow"
		REAL temp_ch=0.032					UNITS u_s			"Characteristic filling Time"
				
	DECLS
		REAL T_in								UNITS u_K			"Average inlet temperature"
		REAL T_c									UNITS u_K			"Combustion temperature"
	
		REAL eta									UNITS no_units		"Combustion efficiency"
		REAL OF									UNITS no_units		"Mixture ratio"
		REAL phi									UNITS no_units		"Equivalence ratio"
		REAL W_F_st								UNITS u_kg_s		"Fuel mass flow that provides a stoichiometric mixture when combined with the oxidant mass flow"
		ALG REAL W_F							UNITS u_kg_s		"Fuel mass flow"
		REAL W_O									UNITS u_kg_s		"Oxidant mass flow"
		REAL W_IF								UNITS u_kg_s		"Inert mass flow through the fuel port"
		REAL W_IO								UNITS u_kg_s		"Inert mass flow through the oxidant port"
		
		HIDDEN REAL fluid_O[ChemName]		UNITS no_units		"Oxidant fluid"		
		HIDDEN REAL fluid_F[ChemName]		UNITS no_units		"Fuel fluid"
		HIDDEN REAL fluid_P[ChemName]		UNITS no_units		"Combustion products fluid"
		REAL cp_R								UNITS u_J_kgK		"Reactant specific heat at constant pressure"
		
		BOOLEAN Combustion											"TRUE if there is combustion and FALSE if either oxidant or fuel is lacking"	
		REAL Q_comb_effective				UNITS u_J_kg		"Effective heat of combustion per oxidant mass flow unit"
				
		REAL A_out								UNITS u_m2			"Output area"
		REAL A_th_d								UNITS u_m2			"Design throat area"
				
		ALG REAL p_out_ch = 100.			UNITS u_Pa			"Choked outlet pressure"	
		
		REAL c_star								UNITS u_m_s			"Characteristic exhaust velocity"	
				
		REAL A_r									UNITS u_m2			"Area at the characteristic section of heat exchange"
		REAL M_r									UNITS no_units		"Mach number at the characteristic section of heat exchange"
		REAL Pr_r								UNITS no_units		"Prandtl number at the characteristic section of heat exchange"	
		REAL visc_r								UNITS u_Pas			"Dynamic viscosity at the characteristic section of heat exchange"
		REAL T_aw								UNITS u_K			"Adiabatic wall temperature of combustion gases"
		REAL h_g									UNITS u_W_m2K		"Combustion gases heat transfer coefficient"
		
		DISCR REAL k_1									UNITS no_units		"Inverso del gasto característico"
		DISCR REAL k_2									UNITS no_units		"Inverso de la entalpía característica"
		REAL rho_trans							UNITS u_kg_m3		"Density calculated in transient model"
		
	INIT PRIORITY 100
		ASSERT (AR_r >= 1) KILLPOINT "AR_r CAN NOT BE LOWER THAN 1!"
		
		Init_fluid(Comb_prod, g.fluid)		
		
		Combustion = TRUE
		
		phi = 1.
--		IF (Type == Off_design) THEN
--			g.pt = p_c0
--		END IF
--		g.Tt = T_c0
		W_F = W_F0
		
		IF (Zone == Convergent) THEN
			M_r = 0.001
		ELSE
			M_r = 100
		END IF
		
		k_1 = 1./mfr_ch
		k_2 = 1./(mfr_ch*T_ch)

	DISCRETE
		WHEN (SUM (i IN LiquidsGases ; fluid_O[i] * fluid_F[i]) != 0) THEN
			Combustion = FALSE
		END WHEN
		
		WHEN (SUM (i IN LiquidsGases ; fluid_O[i] * fluid_F[i]) == 0) THEN
			Combustion = TRUE
		END WHEN
		
	CONTINUOUS	
		g.A_out = A_out
		AR = A_out / A_th_d		
		
		eta = eta_d		
		
		g.pt = f_O.p_c
		g.pt = f_F.p_c
      		
		IF (Type == Design) INSERT
			g.pt = p_c
		END IF
		
		--f_O.W + f_F.W = g.W
		
		(temp_ch/(rho_ch*k_1)) * rho_trans' = f_O.W + f_F.W - g.W
		
		W_O = f_O.W * (1 - f_O.fluid[Comb_prod])
		W_F = f_F.W * (1 - f_F.fluid[Comb_prod])
		W_IO = f_O.W - W_O
		W_IF = f_F.W - W_F
		OF = W_O / W_F
		phi = OF_st / OF
		W_F_st = W_F / phi
		
		Q_comb_effective = 	ZONE (Combustion)		Q_comb
									OTHERS					0
				
		eta * Q_comb_effective = ((1. + OF) / min(OF,OF_st)) * (cp(fluid_P) * (T_c - T_ref) - cp_R * (T_in - T_ref))
	
		
		g.pt=rho_trans*R(fluid_P)*g.Tt
				
		(1 + phi / OF_st) * cp_R * T_in = cp(fluid_O) * f_O.T + phi / OF_st * cp(fluid_F) * f_F.T		
		
		--((W_O + W_F) / g.W) * cp(fluid_P) * (g.Tt - T_c) + (W_IO / g.W) * f_O.fluid[Comb_prod_cp] * (g.Tt - f_O.T) + (W_IF / g.W) * f_F.fluid[Comb_prod_cp] * (g.Tt - f_F.T) = 0
		
		((W_O + W_F)/ g.W) * cp(fluid_P) * (g.Tt - T_c) + (W_IO / g.W) * f_O.fluid[Comb_prod_cp] * (g.Tt - f_O.T) + (W_IF / g.W) * f_F.fluid[Comb_prod_cp] * (g.Tt - f_F.T) = -(temp_ch/(g.W*rho_ch*k_2*T_ch)) * (rho_trans*cv(fluid_P)*g.Tt'+g.Tt*cv(fluid_P)*rho_trans')
		
		EXPAND (i IN LiquidsGases) fluid_O[i] = f_O.fluid[i] / (1 - f_O.fluid[Comb_prod])
		fluid_O[Comb_prod] = 0
		fluid_O[Comb_prod_M] = 0
		fluid_O[Comb_prod_cp] = 0
		
		EXPAND (i IN LiquidsGases) fluid_F[i] = f_F.fluid[i] / (1 - f_F.fluid[Comb_prod])
		fluid_F[Comb_prod] = 0	
		fluid_F[Comb_prod_M] = 0
		fluid_F[Comb_prod_cp] = 0			
				
		cp_R = (W_O * cp(fluid_O) + W_F * cp(fluid_F)) / (W_O + W_F)
		
		EXPAND (i IN LiquidsGases) fluid_P[i] = 	IF (Combustion)	(fluid_O[i] * max(1 - phi,0) * W_O + fluid_F[i] * max(phi - 1,0) * W_F_st) / (W_O + W_F)
																ELSE					(fluid_O[i] * W_O + fluid_F[i] * W_F) / (W_O + W_F)
		fluid_P[Comb_prod] = 	IF (Combustion)	(fluid_O[Comb_prod] * max(1 - phi,0) * W_O + fluid_F[Comb_prod] * max(phi - 1,0) * W_F_st + (1 - max(1 - phi,0)) * (W_O + W_F_st)) / (W_O + W_F)
										ELSE					0
		fluid_P[Comb_prod_M] = 	IF (Combustion)	M_P
										ELSE					0
		fluid_P[Comb_prod_cp] = IF (Combustion)	cp_P
										ELSE					0
		
		EXPAND (i IN LiquidsGases) g.fluid[i] = fluid_P[i] * (W_O + W_F) / g.W
		g.fluid[Comb_prod] = (fluid_P[Comb_prod] * (W_O + W_F) + (W_IO + W_IF)) / g.W		
		g.fluid[Comb_prod_M] = fluid_P[Comb_prod_M] -- (M_P * fluid_P[Comb_prod] * (W_O + W_F) + f_O.fluid[Comb_prod_M] * f_O.fluid[Comb_prod] * f_O.W + f_F.fluid[Comb_prod_M] * f_F.fluid[Comb_prod] * f_F.W) / (fluid_P[Comb_prod] * (W_O + W_F) + f_O.fluid[Comb_prod] * f_O.W + f_F.fluid[Comb_prod] * f_F.W)
		g.fluid[Comb_prod_cp] = fluid_P[Comb_prod_cp] -- (cp_P * fluid_P[Comb_prod] * (W_O + W_F) + f_O.fluid[Comb_prod_cp] * f_O.fluid[Comb_prod] * f_O.W + f_F.fluid[Comb_prod_cp] * f_F.fluid[Comb_prod] * f_F.W) / (fluid_P[Comb_prod] * (W_O + W_F) + f_O.fluid[Comb_prod] * f_O.W + f_F.fluid[Comb_prod] * f_F.W)
		
		g.W = g.pt * A_th_d / (sqrt(R(g.fluid) * g.Tt)/FGAMMA(g.fluid))
		--g.W = g.pt * A_th_d / (c_star)
		
		IF (Type == Off_design) INSERT
			A_th = A_th_d
		END IF
		
		A_out / A_th_d = FGAMMA(g.fluid) / ((p_out_ch / g.pt)**(1. / gamma(g.fluid)) * sqrt(2 * gamma(g.fluid) * (1 - (p_out_ch / g.pt)**((gamma(g.fluid) - 1.) / gamma(g.fluid))) / (gamma(g.fluid) - 1)))
		
		c_star = sqrt(R(g.fluid) * g.Tt) / FGAMMA(g.fluid)
		i.Data[1] = c_star		
		
		
		-- Cooling of the nozzle
		
		T_aw = g.Tt * ((1. + Pr_r**0.33 * (gamma(g.fluid) - 1.) * M_r**2 / 2.) / (1. + (gamma(g.fluid) - 1.) * M_r**2 / 2.))
		Pr_r = 4. * gamma(g.fluid) / (9. * gamma(g.fluid) - 5.)
		visc_r = 1.184e-7 * M(g.fluid)**0.5 * T_aw**0.6		
		
		IF (Cooled == No) INSERT
			h.A = 0
			A_r = A_th_d
			M_r = 1.
			
			h_g = 0
			h.T = T_aw
			h.Q = 0
			
		ELSE		
			h.A = A_wet
			A_r = A_th_d * AR_r
			AR_r = 1 / M_r * FGAMMA(g.fluid) / sqrt(gamma(g.fluid)) * (1. + (gamma(g.fluid) - 1.) * M_r**2 / 2.)**((gamma(g.fluid) + 1.) / (2. * (gamma(g.fluid) - 1.)))
		
			-- Bartz correlation
			h_g = 0.026 / sqrt(4 * A_th_d / PI)**0.2 * (visc_r**0.2 * cp(g.fluid) / Pr_r**0.6) * (g.pt / c_star)**0.8 * (A_th_d / A_r)**0.9
			h.Q = h_g * (T_aw - h.T) * h.A
			
		END IF		
		
END COMPONENT