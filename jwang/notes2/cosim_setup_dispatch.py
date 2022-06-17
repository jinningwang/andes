# --- set up EV generator data ---
ev_idx = 'PV_10'
ssa.PV.set(src='p0', idx=ev_idx, attr='v', value=sse.Ptc)
ssa.PV.set(src='pmax', idx=ev_idx, attr='v', value=sse.Pu)
ssa.PV.set(src='pmin', idx=ev_idx, attr='v', value=sse.Pl)

# --- setup pandapower: ssp ---
ssp = to_pandapower(ssa)

# set EV generator as uncontrollable
ssp.gen.controllable.iloc[9] = False

# add gen cost, G1-11. G10 EV, G11 Slack
linearcost = [20.14, 20.20, 20.19, 20.20, 20.25, 20.21, 20.23, 20.22, 20.15, 0, 20.30]

gen_cost = np.array([[2., 0., 0., 3., 0., 0., 0.]] * ssp.gen.shape[0])
gen_cost[:, 5] = linearcost  # c1

add_gencost(ssp, gen_cost)

# --- setup RTED: ssd ---
ssd = rted2(name='RTED', OutputFlag=0)
ssd.from_andes(ssa)
ssd.build()

# set EV generator as uncontrollable
ssd.gen.ctrl.loc[ev_idx] = 0

# set EV geenrator as type2
prumax = sse.g_frc()[0]
prdmax = sse.g_frc()[1]
ssd.def_type2([ev_idx], [prumax], [prdmax])

# set ramp5
ramp_15 = [156, 120, 130, 110, 80, 105, 90, 90, 150, 999, 200]
# ramp_hour = [80, 80, 80, 50, 50, 50, 30, 30, 30, 999, 30]
ssd.gen['ramp_5'] = np.array(ramp_15) / 3 / ssd.mva

# set cost
ssd.cost['c1'] = linearcost

# adjust SFR cost of EV lower than SynGen
ssd.cost['cru'] = [0] * ssd.gen.shape[0]
ssd.cost['cru'].loc[ev_idx] = - 0.001
ssd.cost['crd'] = ssd.cost.cru

# --- benchmark ssd with ssp using DCOPF ---
pp.rundcopp(ssp)
ppres = ssp.res_gen.copy()
gb_res = ssd.solve(disable_sfr=True, info=False)
ppres['p_mw(GB)'] = ssd.mva * gb_res['pg'].values
print(f"pp cost={ssp.res_cost}, gb cost={ssd.res_cost}")
