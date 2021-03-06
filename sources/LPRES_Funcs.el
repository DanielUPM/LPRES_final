--------------------------------------------------------------------------------
-- EcosimPro Simulation Source Code
-- FILE NAME: LPRES_Funcs.el
-- DESCRIPTION: Defines functions for thermodynamic properties calculation
-- NOTES:
-- AUTHOR: Pablo Sierra Heras
-- CREATION DATE: 05/12/2014
--------------------------------------------------------------------------------

-- Libraries
USE MATH VERSION "3.1"


-- Functions
--------------------------------------------------------------------------------
-- Purpose: To initialise the value of the variable fluid that goes through the 
-- ports
--------------------------------------------------------------------------------
FUNCTION NO_TYPE Init_fluid
   (
	IN ENUM ChemName fluid_name  										"Working fluid name",
	OUT REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Initialises the value of the variable fluid that goes through the ports"

   BODY	
		FOR (i IN ChemName)
    		fluid[i] = 0
		END FOR
		
		fluid[fluid_name] = 1.
		
		RETURN

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To know the chemical name of the fluid by means of the variable 
-- fluid that goes through the ports
--------------------------------------------------------------------------------
FUNCTION ENUM ChemName Know_fluid
   (
	IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Returns the chemical name of the fluid by means of the variable fluid that goes through the ports"
	
	DECLS
		ENUM ChemName fluid_name 			"Working fluid name"

   BODY	
		FOR(i IN LiquidsGases)
    		IF (fluid[i] != 0) THEN
				fluid_name = i
			END IF
		END FOR
		
    	IF (fluid[Comb_prod] != 0) THEN
			fluid_name = Comb_prod
		END IF
		
		RETURN fluid_name

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To know the state of a fluid
--------------------------------------------------------------------------------
FUNCTION ENUM ChemState State
   (
	IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"	
   )
	
	"Returns the state of a fluid"
	
	DECLS
		ENUM ChemState fluid_state			"Working fluid state"
		ENUM ChemName fluid_name 			"Working fluid name"

   BODY
		fluid_name = Know_fluid(fluid)
		
		IF (setofPos(Liquids,fluid_name) != 0) THEN	
      	fluid_state = Liquid
		ELSE
			fluid_state = Gas
		END IF

		RETURN fluid_state

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To know the gas to which a liquid vaporises
--------------------------------------------------------------------------------
FUNCTION ENUM Gases Vaporisation
   (
	IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Returns the gas to which a liquid vaporises"
	
	DECLS
		ENUM Gases Vapour[LV] = {O2, CH4, H2, MMH_vapour}								"Gases name to which each liquid vaporises"
		ENUM Liquids liquid_name  															"Working liquid name"

   BODY		
		ASSERT (Liquid == State(fluid)) FATAL "ONLY LIQUIDS CAN BE USED IN THIS FUNCTION!"
		
		liquid_name = Know_fluid(fluid)		
		
		ASSERT (setofPos(LV,liquid_name) != 0) FATAL "ONLY LIQUIDS IN ENUM LV CAN BE VAPORISED!"
		
		RETURN Vapour[liquid_name]

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate molar mass
--------------------------------------------------------------------------------
FUNCTION REAL M
   (
	IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"	
   )
	
	"Calculates molar mass"
	
	DECLS
		REAL Chem_M[LiquidsGases] = {31.9988, 92.011, 34.0147, 63.01, 38, 172, 16.0426, 2.01594, 32.04516, 60.1, 46.07, 136.234, 170.34, 60, 18, 60.1, 28.958538, 39.948, 16.0426, 28.0104, 44.0098, 2.01594, 4.0026, 28.01348, 31.9988, 46.07}		\
															UNITS u_g_mol		"Molar mass of each chemical"

   BODY
		RETURN  1 / (SUM (i IN LiquidsGases ; fluid[i] / Chem_M[i]) + fluid[Comb_prod] / fluid[Comb_prod_M]) 

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate gas constant
--------------------------------------------------------------------------------
FUNCTION REAL R
   (
	IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Calculates gas constant"

   BODY
		ASSERT (Gas == State(fluid)) FATAL "ONLY GASES CAN BE USED IN THIS FUNCTION!"
		
      RETURN R_u / M(fluid)

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate specific heat at constant pressure
--------------------------------------------------------------------------------
FUNCTION REAL cp
   (
   IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Calculates specific heat at constant pressure"
	
	DECLS
		REAL Chem_cp[LiquidsGases] = {1680, 1600, 2629, 1720, 1640, 2093, 3480, 7320, 3080, 2720, 2840, 1675, 1800, 1800, 4182, 2400, 1004, 520, 3000, 1080, 1000, 15800, 5193, 1039, 1000, 1600}		\
															UNITS u_J_kgK		"Specific heat at constant pressure of each chemical"

   BODY
      RETURN SUM (i IN LiquidsGases ; Chem_cp[i] * fluid[i]) + fluid[Comb_prod_cp] * fluid[Comb_prod] 

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate specific heat at constant volume
--------------------------------------------------------------------------------
FUNCTION REAL cv
   (
   IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Calculates specific heat at constant volume"

   BODY
      RETURN cp(fluid) - R(fluid)

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate ratio of specific heats
--------------------------------------------------------------------------------
FUNCTION REAL gamma
   (
   IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Calculates ratio of specific heats"

   BODY
      RETURN cp(fluid) / cv(fluid)

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate the value of the Gamma function
--------------------------------------------------------------------------------
FUNCTION REAL FGAMMA
   (
   IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Calculates the value of the Gamma function"

   BODY
      RETURN sqrt(gamma(fluid)) * (2. / (gamma(fluid) + 1.))**((gamma(fluid) + 1.) / (2. * (gamma(fluid) - 1.)))

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate the density of a fluid in liquid state
--------------------------------------------------------------------------------
FUNCTION REAL rho
   (
   IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Calculates the density of a fluid in liquid state"
	
	DECLS
		REAL Chem_rho[Liquids] = {1200, 1477, 1450, 1560, 1500, 809, 422, 73, 1008, 783, 874, 920, 749, 950, 1000, 786}		UNITS u_kg_m3		"Liquid state density of each liquid"

   BODY
		ASSERT (Liquid == State(fluid)) FATAL "ONLY LIQUIDS CAN BE USED IN THIS FUNCTION!"
		
      RETURN SUM (i IN Liquids ; Chem_rho[i] * fluid[i])

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate thermal conductivity
--------------------------------------------------------------------------------
FUNCTION REAL cond
   (
   IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Calculates thermal conductivity"
	
	DECLS
		REAL Chem_cond[Liquids] = {0.149, 0.131, 0.353, 0.294, 0.1, 0.137, 0.15, 0.117, 0.488, 0.035, 0.246, 0.11, 0.17, 0.12, 0.607, 0.145}		UNITS	u_W_mK		"Thermal conductivity of each liquid"

   BODY      
		ASSERT (Liquid == State(fluid)) FATAL "ONLY LIQUIDS CAN BE USED IN THIS FUNCTION!"
		
      RETURN SUM (i IN Liquids ; Chem_cond[i] * fluid[i])

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate dynamic viscosity
--------------------------------------------------------------------------------
FUNCTION REAL visc
   (
   IN REAL fluid[ChemName] 				UNITS no_units			"Working fluid"
   )
	
	"Calculates dynamic viscosity"
	
	DECLS
		REAL Chem_visc[Liquids] = {0.00019, 0.000423, 0.00127, 0.001, 0.000079, 0.00021, 0.00002, 0.00002, 0.00097, 0.000754, 0.0004, 0.0028, 0.0045, 0.0045, 0.001, 0.002}			UNITS u_Pas			"Dynamic viscosity of each liquid"

   BODY      
		ASSERT (Liquid == State(fluid)) FATAL "ONLY LIQUIDS CAN BE USED IN THIS FUNCTION!"
		
      RETURN SUM (i IN Liquids ; Chem_visc[i] * fluid[i])

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate the friction factor
--------------------------------------------------------------------------------
FUNCTION REAL hdc_fric
	(
	IN REAL D					UNITS u_m			"Hydraulic diameter",
	IN REAL rug					UNITS u_m			"Absolute rugosity",
	IN REAL Re					UNITS no_units		"Reynolds number"
   )
	 
	"Calculates the friction factor"
	 
   DECLS
   	REAL rey					UNITS no_units		"Auxiliar variable"
		REAL fric				UNITS no_units		"Friction factor"
		REAL a					UNITS no_units		"Auxiliar variable"
		REAL b					UNITS no_units		"Auxiliar variable"
		
   BODY
		rey = max(abs(Re), 1e-5)

		a = (2.457 * log(1. / ((7. / Re)**.9 + .27 * rug / D)))**16
		b = (37530. / rey)**16

		fric = 8. * ((8. / rey)**12 + 1. / ((a + b)**1.5))**.0833333333
		  
		RETURN fric
		  
END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate geopotential altitude
--------------------------------------------------------------------------------
FUNCTION REAL Geopotential_Altitude
   (
   IN REAL z			UNITS u_m			"Geometric altitude"
   )
	
	"Calculates geopotential altitude"

   BODY
      RETURN z * 6371000. / (z + 6371000.)

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate International Standard Atmosphere (ISA) temperature
--------------------------------------------------------------------------------
FUNCTION REAL ISA_Temperature
   (
   IN REAL z			UNITS u_m			"Geometric altitude"
   )

   "Calculates International Standard Atmosphere (ISA) temperature"
	
   DECLS
		REAL T_amb = T_std																						UNITS u_K			"Ambient temperature"
		REAL a[8] = { -6.5 , 0. , 1. , 2.8 , 0. , -2.8 , -2. , 0. }									UNITS "K/km"		"Thermal gradient in each atmosphere layer"
		REAL T_0	= T_std																							UNITS u_K			"Base temperature of the atmosphere layer"
		REAL h_0	= 0.																								UNITS u_m			"Base geopotential altitude of the atmosphere layer"
		REAL h_max[7] = { 11000. , 20000. , 32000. , 47000. , 51000. , 71000. , 84852. }		UNITS u_m			"Maximum geopotential altitude in each atmosphere layer"
		INTEGER i																															"Variable used as an index"	
		INTEGER j = 1																														"Variable used as an index"			
		REAL h																										UNITS u_m			"Geopotential altitude"					

   BODY		
		h = Geopotential_Altitude(z)
						
		FOR (i IN 1,7) 
      	IF (h > h_max[i]) THEN				
				T_0 = T_0 + a[i] * (h_max[i] - h_0) / 1000.
				h_0 = h_max[i]
				j = i + 1
			END IF
		END FOR
		
		T_amb = T_0 + a[j] * (h - h_0) / 1000.
		
		RETURN T_amb

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate International Standard Atmosphere (ISA) pressure
--------------------------------------------------------------------------------
FUNCTION REAL ISA_Pressure
   (
   IN REAL z			UNITS u_m			"Geometric altitude"
   )

   "Calculates International Standard Atmosphere (ISA) pressure"
	
	DECLS
		REAL p_amb = p_std																						UNITS u_Pa			"Ambient pressure"
		REAL a[8] = { -6.5 , 0. , 1. , 2.8 , 0. , -2.8 , -2. , 0. }									UNITS "K/km"		"Thermal gradient in each atmosphere layer"
		BOOLEAN b[8] = { FALSE , TRUE , FALSE , FALSE , TRUE , FALSE , FALSE , TRUE }									"Variable that shows if the value of the temperature gradient is 0 (TRUE) or not (FALSE) in each atmosphere layer"
		REAL T_0	= T_std																							UNITS u_K			"Base temperature of the atmosphere layer"
		REAL p_0	= p_std																							UNITS u_Pa			"Base pressure of the atmosphere layer"
		REAL h_0	= 0.																								UNITS u_m			"Base geopotential altitude of the atmosphere layer"
		REAL h_max[7] = { 11000. , 20000. , 32000. , 47000. , 51000. , 71000. , 84852. }		UNITS u_m			"Maximum geopotential altitude in each atmosphere layer"
		INTEGER i																															"Variable used as an index"	
		INTEGER j = 1																														"Variable used as an index"	
		REAL h																										UNITS u_m			"Geopotential altitude"

   BODY		
		h = Geopotential_Altitude(z)
						
		FOR (i IN 1,7) 
      	IF (h > h_max[i]) THEN
				IF (b[i]) THEN
					p_0 = p_0 * exp(- g_0 * (h_max[i] - h_0) / (T_0 * 287.))
				ELSE
					p_0 = p_0 * ((T_0 + a[i] * (h_max[i] - h_0) / 1000.) / T_0)**(- g_0 / (a[i] * 287. / 1000.))
				END IF
				T_0 = T_0 + a[i] * (h_max[i] - h_0) / 1000.				
				h_0 = h_max[i]
				j = i + 1
			END IF
		END FOR
		
		IF (b[j]) THEN
			p_amb = p_0 * exp(- g_0 * (h - h_0) / (T_0 * 287.))
		ELSE
			p_amb = p_0 * ((T_0 + a[j] * (h - h_0) / 1000.) / T_0)**(- g_0 / (a[j] * 287. / 1000.))
		END IF
		
		RETURN p_amb

END FUNCTION

--------------------------------------------------------------------------------
-- Purpose: To calculate International Standard Atmosphere (ISA) density
--------------------------------------------------------------------------------
FUNCTION REAL ISA_Density
   (
   IN REAL z			UNITS u_m			"Geometric altitude"
   )

   "Calculates International Standard Atmosphere (ISA) density"
	
   DECLS
		REAL rho_amb = p_std / (T_std * 287.)																UNITS u_kg_m3		"Ambient density"
		REAL a[8] = { -6.5 , 0. , 1. , 2.8 , 0. , -2.8 , -2. , 0. }									UNITS "K/km"		"Thermal gradient in each atmosphere layer"
		BOOLEAN b[8] = { FALSE , TRUE , FALSE , FALSE , TRUE , FALSE , FALSE , TRUE }									"Variable that shows if the value of the temperature gradient is 0 (TRUE) or not (FALSE) in each atmosphere layer"
		REAL T_0	= T_std																							UNITS u_K			"Base temperature of the atmosphere layer"
		REAL rho_0 = p_std / (T_std * 287.)																	UNITS u_kg_m3		"Base pressure of the atmosphere layer"
		REAL h_0	= 0.																								UNITS u_m			"Base geopotential altitude of the atmosphere layer"
		REAL h_max[7] = { 11000. , 20000. , 32000. , 47000. , 51000. , 71000. , 84852. }		UNITS u_m			"Maximum geopotential altitude in each atmosphere layer"
		INTEGER i																															"Variable used as an index"	
		INTEGER j = 1																														"Variable used as an index"	
		REAL h																										UNITS u_m			"Geopotential altitude"
		
   BODY		
		h = Geopotential_Altitude(z)
						
		FOR (i IN 1,7) 
      	IF (h > h_max[i]) THEN
				IF (b[i]) THEN
					rho_0 = rho_0 * exp(- g_0 * (h_max[i] - h_0) / (T_0 * 287.))
				ELSE
					rho_0 = rho_0 * ((T_0 + a[i] * (h_max[i] - h_0) / 1000.) / T_0)**(- g_0 / (a[i] * 287. / 1000.) - 1.)
				END IF
				T_0 = T_0 + a[i] * (h_max[i] - h_0) / 1000.
				h_0 = h_max[i]
				j = i + 1
			END IF
		END FOR
		
		IF (b[j]) THEN
			rho_amb = rho_0 * exp(- g_0 * (h - h_0) / (T_0 * 287.))
		ELSE
			rho_amb = rho_0 * ((T_0 + a[j] * (h - h_0) / 1000.) / T_0)**(- g_0 / (a[j] * 287. / 1000.) - 1.)
		END IF
		
		RETURN rho_amb

END FUNCTION