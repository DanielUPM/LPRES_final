--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_CompGasGen.el
-- DESCRIPTION: Defines gas generator type components
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Components
--------------------------------------------------------------------------------
-- Component that represents a gas generator
--------------------------------------------------------------------------------
COMPONENT GasGen (ENUM OnOffDesign Type = Off_design)

	"Gas generator"

	PORTS
      IN FluidInj f_O                			"Inlet fluid injection port"
      IN FluidInj f_F              				"Inlet fluid injection port"
      OUT Fluid g          						"Outlet fluid port"
		OUT Info (n = 1) i							"Outlet information port"
		
	DATA
		REAL eta_d = 0.9						UNITS no_units		"Design combustion efficiency"
		REAL OF_st = 8.						UNITS no_units		"Stoichiometric mixture ratio"
		REAL Q_comb	= 2000000.				UNITS u_J_kg		"Heat of combustion per oxidant mass flow unit"
		REAL cp_P = 4182.						UNITS u_J_kgK		"Specific heat at constant pressure of the products using a stoichiometric mixture"
		REAL M_P	= 32.							UNITS u_g_mol		"Molar mass of the products using a stoichiometric mixture"			
		REAL p_c	= 5000000.					UNITS u_Pa			"Design combustion pressure"
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
					
		REAL c_star								UNITS u_m_s			"Characteristic exhaust velocity"	
				
		DISCR REAL k_1									UNITS no_units		"Inverso del gasto característico"
		DISCR REAL k_2									UNITS no_units		"Inverso de la entalpía característica"
		REAL rho_trans							UNITS u_kg_m3		"Density calculated in transient model"
		
	INIT PRIORITY 100
			
		Init_fluid(Comb_prod, g.fluid)		
		
		Combustion = TRUE
		
		phi = 1.
--		IF (Type == Off_design) THEN
--			g.pt = p_c0
--		END IF
--		g.Tt = T_c0
		W_F = W_F0
			
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
		
		--(temp_ch/(rho_ch*k_2*T_ch)) * rho_trans*cv_R*T_trans' = ((1. + OF) / min(OF,OF_st)) * (cp(fluid_P) * (T_c - T_ref) - cp_R * (T_in - T_ref)) - eta * Q_comb_effective
		
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
		
		c_star = sqrt(R(g.fluid) * g.Tt) / FGAMMA(g.fluid)
		i.Data[1] = c_star

END COMPONENT