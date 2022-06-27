# --- set up EV generator data ---
ev_idx = 'PV_10'
ssa.PV.set(src='p0', idx=ev_idx, attr='v', value=sse.Ptc/ssa.config.mva)
ssa.PV.set(src='pmax', idx=ev_idx, attr='v', value=2*sse.Pu/ssa.config.mva)
ssa.PV.set(src='pmin', idx=ev_idx, attr='v', value=2*sse.Pl/ssa.config.mva)

# --- setup pandapower: ssp ---
ssp = to_pandapower(ssa)

# set EV generator as uncontrollable
ssp.gen.controllable.iloc[9] = False

# add gen cost, unit: $/MWh; G1-11. G10 EV, G11 Slack; 
# set EV as -10, for the cost of SFR mileage
c1 = [20, 20, 20, 20, 20, 20, 20, 20, 20, -10, 20]
c0 = [500, 380, 42, 380, 295, 400, 350, 330, 490, 0000, 550]
c2 = [0.014, 0.020, 0.194, 0.020, 0.0255, 0.0210, 0.230, 0.0222, 0.0150, 0.000, 0.0300]
gen_cost = np.array([[2., 0., 0., 3., 0., 0., 0.]] * ssp.gen.shape[0])
gen_cost[:, 4] = np.array(c2) / ssa.config.mva
gen_cost[:, 5] = np.array(c1) / ssa.config.mva
gen_cost[:, 6] = np.array(c0) / ssa.config.mva

add_gencost(ssp, gen_cost)

# --- setup RTED: ssd ---
ssd = rted3(name='RTED', OutputFlag=0)
ssd.from_andes(ssa)
ssd.build()

# set EV generator as uncontrollable
ssd.gen.ctrl.loc[ev_idx] = 0

# set EV geenrator as type2
prumax = sse.g_frc()[0]
prdmax = sse.g_frc()[1]
ssd.def_type2([ev_idx], [prumax], [prdmax])

# Case data comes from a MPCE
# https://ieeexplore.ieee.org/document/9018441
# set ramp_5
ramp_15 = [156, 120, 130, 110, 80, 105, 90, 90, 150, 999, 200]
# ramp_hour = [80, 80, 80, 50, 50, 50, 30, 30, 30, 999, 30]
ssd.gen['ramp_5'] = np.array(ramp_15) / 3 / ssd.mva

# set cost
ssd.cost['c1'] = c1
ssd.cost['c2'] = c2
ssd.cost['c0'] = c0

# adjust SFR cost of EV lower than SynGen
ssd.cost['cru'] = [0] * ssd.gen.shape[0]
ssd.cost['cru'].loc[ev_idx] = - 0.000001
ssd.cost['crd'] = ssd.cost.cru

# --- benchmark ssd with ssp using DCOPF ---
pp.rundcopp(ssp)
dc_comp = ssp.res_gen.copy()
gb_res = ssd.solve(disable_sfr=True, disable_ramp=True, info=False)
dc_comp['p_mw(GB)'] = ssd.mva * gb_res['pg'].values
print(f"pp cost={ssp.res_cost}, gb cost={ssd.res_cost}")
