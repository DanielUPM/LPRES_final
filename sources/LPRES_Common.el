--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_Common.el
-- DESCRIPTION: Defines basic components
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Global constants
CONST REAL p_std = 101325.						UNITS u_Pa			"Standard pressure"
CONST REAL T_std = 288.15 						UNITS u_K			"Standard temperature"

CONST REAL g_0 = 9.80665    					UNITS u_m_s2		"Gravity acceleration at the Earth surface"

CONST REAL R_u = 8314.							UNITS u_J_kmolK	"Universal Gas Constant"

CONST REAL T_ref = 298.15						UNITS u_K			"Reference temperature"


-- Global variables
BOUND REAL Altitude = 0.						UNITS u_m			"Geometric altitude"


-- Global enumerations
ENUM ChemName = {LOX, NTO, H2O2, HNO3, LF2, RP_1, LCH4, LH2, N2H4, UDMH, MMH, JP_10, Kerox, Oil, H2O, IPA, Air, Ar, CH4, CO, CO2, H2, He, N2, O2, MMH_vapour, Comb_prod, Comb_prod_M, Comb_prod_cp}		\
								"Names of available chemicals"

SET_OF(ChemName) Liquids = {LOX, NTO, H2O2, HNO3, LF2, RP_1, LCH4, LH2, N2H4, UDMH, MMH, JP_10, Kerox, Oil, H2O, IPA}   "Names of available liquids"
SET_OF(ChemName) Gases = {Air, Ar, CH4, CO, CO2, H2, He, N2, O2, MMH_vapour, Comb_prod}            "Names of available gases"  -- añadido "Comb_prod" por problemas de inicialización
SET_OF(ChemName) LiquidsGases = {LOX, NTO, H2O2, HNO3, LF2, RP_1, LCH4, LH2, N2H4, UDMH, MMH, JP_10, Kerox, Oil, H2O, IPA, Air, Ar, CH4, CO, CO2, H2, He, N2, O2, MMH_vapour}			\
								"Names of available chemicals except Comb_prod"
								
SET_OF(ChemName) LV = {LOX, LCH4, LH2, MMH}			"Liquids that can be vaporised"

ENUM ChemState = {Liquid, Gas}							"State of available chemicals"

ENUM ConDiv = {Convergent, Divergent}
ENUM YesNo = {Yes, No}
ENUM Type_Inlet = {All, Unknown_W}
ENUM Type_All = {Design, Known_pi, Known_W, Off_design, Known_pt_out, Known_dp, Known_dpr, Darcy}
SET_OF(Type_All) OnOffDesign = {Design, Off_design}
SET_OF(Type_All) Type_Turbines = {Known_pi, Known_W, Off_design}
SET_OF(Type_All) Type_Turbine_liq = {Known_dp, Known_W, Off_design}
SET_OF(Type_All) Type_Regulator = {Design, Known_pt_out, Known_dp, Known_dpr}
SET_OF(Type_All) Type_Cooling = {Darcy, Known_dp}
ENUM AngCoef = {Angles, Coefficients}

