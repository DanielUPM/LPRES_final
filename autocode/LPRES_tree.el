-- Generated automatically by - PROOSIS - 3.8.1 




-- EL code of the schematic.
-- The COMPONENT definition lines are blocked for edition.
-- You can edit the parameters clicking over them.

COMPONENT LPRES_tree 

PORTS
   IN LPRES.Fluid Fluid

   IN LPRES.FluidInj FluidInj

   IN LPRES.GasNozzle GasNozzle

   IN LPRES.Heat Heat

   IN LPRES.Info(
      n = XXX_EMPTY_VALUE_XXX) Info

   IN LPRES.Mechanical Mechanical



TOPOLOGY
   LPRES.Ambient(
      Type = Design) Ambient(
      A = 0.01)

   LPRES.CombCha(
      Type = Design,
      Cooled = No) CombCha(
      eta_d = 0.9,
      OF_st = 8,
      Q_comb = 2000000,
      cp_P = 4182,
      M_P = 32,
      AR = 10,
      A_th = 0.05,
      p_c = 5000000,
      AR_r = 10 / 2,
      A_wet = 1,
      Zone = Divergent,
      p_c0 = 5000000,
      T_c0 = 4000,
      W_F0 = 1000)

   LPRES.Compressor(
      Type = Design,
      Type_AC = Coefficients) Compressor(
      eta_d = 0.8,
      phi_d = 0.05,
      psi_d = 0.7,
      alpha_2r = -30,
      r_m = 0.3,
      A_in = 0.35,
      pi = 4,
      U_0 = 500)

   LPRES.ControlPanel ControlPanel

   LPRES.CoolingJacket(
      Type = Darcy) CoolingJacket(
      L = 1,
      a = 0.002,
      b = 0.004,
      N = 100,
      rug = 5e-05,
      k_w = 370,
      t_w = 0.003,
      dp = 100000)

   LPRES.FlowMeter FlowMeter

   LPRES.GasGen GasGen(
      TPL_d = 0.9,
      eta_d = 0.9,
      OF_st = 8,
      Q_comb = 2000000,
      cp_P = 4182,
      M_P = 32,
      OF = 4,
      W_F0 = 100)

   LPRES.Gearbox Gearbox(
      GR = 1,
      eta = 1)

   LPRES.Injector(
      Type = Design) Injector(
      C_D = 0.5,
      A = 0.05,
      W = 100)

   LPRES.Inlet(
      Type = All) Inlet(
      Tt = 288.15,
      pt = 101325,
      W = 1,
      fluid = LOX)

   LPRES.Junction Junction(
      TPL = 0.9)

   LPRES.Nozzle(
      Type = Design) Nozzle(
      A = 0.02)

   LPRES.NozzleConDiv(
      Type = Design) NozzleConDiv(
      AR = 10,
      A_th = 0.05,
      W = 100)

   LPRES.NozzleExt NozzleExt(
      AR = 2)

   LPRES.Pipe Pipe(
      L = 1,
      D = 0.1,
      K = 5,
      rug = 1.5e-06)

   LPRES.Pump(
      Type = Design) Pump(
      eta_d = 0.8,
      phi_d = 0.05,
      psi_d = 0.7,
      A_in = 0.01,
      r_m = 0.05,
      dp = 5000000,
      U_0 = 500)

   LPRES.Regulator(
      Type = Design) Regulator(
      dp_min = 1500,
      pt_out = 1200000,
      dp = 100000,
      dpr = 0.1)

   LPRES.Shaft Shaft(
      eta = 1)

   LPRES.SplitFrac SplitFrac(
      TPL = 0.9)

   LPRES.Tank(
      Type = Design) Tank(
      fluid_l = LOX,
      T_d = 288.15,
      A_g = 0.001,
      p_d = 1000000)

   LPRES.TankOpen TankOpen(
      fluid = LOX,
      T_d = 288.15)

   LPRES.ThrustMonitor ThrustMonitor

   LPRES.Turbine(
      Type = Known_pi,
      Type_AC = Coefficients) Turbine(
      eta_d = 0.8,
      phi_d = 0.05,
      psi_d = 0.7,
      alpha_2 = 45,
      alpha_4r = -30,
      M = 1,
      A_in = 0.005,
      r_m = 0.01,
      rpm = 30000,
      pi = 1.5,
      W = 5,
      U_0 = 10000)

   LPRES.Turbine_ch(
      Type = Known_pi) Turbine_ch(
      eta_d = 0.5,
      alpha_2 = 45,
      A_in = 0.001,
      rpm = 30000,
      pi = 10,
      W = 5)

   LPRES.Turbine_liq(
      Type = Known_dp) Turbine_liq(
      eta_d = 0.8,
      phi_d = 0.05,
      psi_d = 0.7,
      A_in = 0.01,
      r_m = 0.05,
      rpm = 30000,
      dp = 5000000,
      W = 5,
      U_0 = 500)


END COMPONENT