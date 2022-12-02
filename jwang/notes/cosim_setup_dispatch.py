# --- set up EV generator data ---
if case_name == 'npcc':
    ev_gen = 'PV_49'
    ssp_tol = 1e-1
    opf_alt = 0
elif case_name == 'ieee39':
    ev_gen = 'PV_10'
    ssp_tol = 1e-3
    opf_alt = 1

ev_idx = ssa.DG.find_idx(keys='gen', values=[ev_gen])
ssa.PV.set(src='p0', idx=ev_gen, attr='v', value=30 / ssa.config.mva)
ssa.PV.set(src='pmax', idx=ev_gen, attr='v', value=99999 / ssa.config.mva)
ssa.PV.set(src='pmin', idx=ev_gen, attr='v', value=-99999 / ssa.config.mva)


# --- setup pandapower: ssp ---
ssp = to_pandapower(ssa, tol=ssp_tol)
# NOTE: only effective for a single EV
ev_rid = ssp.gen[ssp.gen['name'] == ev_gen].index[0]
# ssp.gen.loc[ev_rid, 'controllable'] = False
ssp.gen.loc[ev_rid, 'max_p_mw'] = 9999
ssp.gen.loc[ev_rid, 'min_p_mw'] = -9999

if case_name == 'npcc':
    # --- cost ---
    # add gen cost, unit: $/MWh; EV1, EV2, EV; 27 Slack;
    # set EV as -10, for the cost of SFR mileage
    ncost = pd.read_csv(dir_path + '/case/npcc_cost.csv', index_col=0)
    c1 = list(ncost['Incremental cost'].values)
    c1.append(-10)
    c0 = list(ncost['Fixed cost'].values)
    c0.append(0)
    c2 = [0] * ssp.gen.shape[0]
    gen_cost = np.array([[2., 0., 0., 3., 0., 0., 0.]] * ssp.gen.shape[0])
    gen_cost[:, 5] = np.array(c1)
    gen_cost[:, 6] = np.array(c0)
    add_gencost(ssp, gen_cost)
elif case_name == 'ieee39':
    # --- cost ---
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
ssd.gen.loc[ev_gen, 'ctrl'] = 0
ssd.build()

# set EV generator as uncontrollable
prumax = 0
prdmax = 0
# set EV geenrator as type2
ssd.def_type2([ev_gen], [prumax], [prdmax])
# --- ramp ---
if case_name == 'npcc':
    ramp_15 = [99999] * ssp.gen.shape[0]
    # ramp_hour = [80, 80, 80, 50, 50, 50, 30, 30, 30, 999, 30]
    ssd.gen['ramp_5'] = np.array(ramp_15) / 3 / ssd.mva
elif case_name == 'ieee39':
    # --- ramp ---
    # Case data comes from a MPCE
    # https://ieeexplore.ieee.org/document/9018441
    # set ramp_5
    ramp_15 = [156, 120, 130, 110, 80, 105, 90, 90, 150, 99999, 200]
    # ramp_hour = [80, 80, 80, 50, 50, 50, 30, 30, 30, 999, 30]
    ssd.gen['ramp_5'] = 10 * np.array(ramp_15) / 3 / ssd.mva

# --- SFR cost ---
# adjust SFR cost of EV lower than SynGen
ssd.cost['cru'] = [0] * ssd.gen.shape[0]
ssd.cost['cru'].loc[ev_gen] = -1
ssd.cost['crd'] = ssd.cost.cru

# set cost
ssd.cost['c1'] = c1
ssd.cost['c2'] = c2
ssd.cost['c0'] = c0

# --- link_table ---
ssa_key = make_link_table(ssa)

# --- add controllable in the link table ---
ssa_bus = ssa.Bus.as_df().reset_index().rename(columns={'uid': 'bus_pp', 'idx': 'bus_idx'})
ssp_gen = ssp.gen.reset_index().rename(columns={'index': 'gen_pp', 'name': 'stg_idx', 'controllable': 'ctrl'})
ssa_key2 = pd.merge(left=ssa_key,
                    right=ssp_gen[['stg_idx', 'gen_pp', 'ctrl']],
                    on='stg_idx', how='left')
# --- gen ctrl ---
# update ctrl by existance of governor
ssa_key2['ctrl'] = ssa_key2['gov_idx'].astype(bool)
# ev_gen ctrl
ssa_key2.loc[ssa_key2[ssa_key2['stg_idx'] == ev_gen].index[0], 'ctrl'] = False
# NOTE: for NPCC case, set slack gen as uncontrollable
if case_name == 'npcc':
    ssa_key2['ctrl'] = False
    gen_op_idx = ['PV_1', 'PV_2', 'PV_3', 'PV_4', 'PV_5', 'PV_6', 'PV_7',
                  'PV_8', 'PV_9', 'PV_10', 'PV_11', 'PV_12', 'PV_13',
                  'PV_14', 'PV_16', 'PV_17', 'PV_18', 'PV_19',
                  'PV_20', 'PV_21', 'PV_22', 'Slack_27',
                  'PV_28', 'PV_29', 'PV_30', 'PV_31', 'PV_35', 'PV_36']
    ssa_key2.loc[ssa_key2[ssa_key2['stg_idx'].isin(gen_op_idx)].index, 'ctrl'] = True
# ssp
sspc_df = ssa_key2[['stg_name', 'ctrl']]
sspc_df.columns = ['name', 'controllable']
sspc_df = ssp.gen[['name']].merge(right=sspc_df, on='name', how='left')
ssp.gen['controllable'] = sspc_df['controllable']

gen_rid = ssp.gen[ssp.gen['controllable'] == False].index
ssp.gen.loc[gen_rid, 'max_p_mw'] = ssp.gen.loc[gen_rid, 'p_mw']
ssp.gen.loc[gen_rid, 'min_p_mw'] = ssp.gen.loc[gen_rid, 'p_mw']
ssp.gen['controllable'] = True
ssp.gen.loc[ev_rid, 'controllable'] = False

# ssd
ssdc_df = ssa_key2[['stg_name', 'ctrl']]
ssdc_df.columns = ['idx', 'ctrl']
ssdc_df = ssd.gen[['idx']].merge(right=ssdc_df, on='idx', how='left')
ssd.gen['ctrl'] = ssdc_df['ctrl'].values

# --- benchmark ssd with ssp using DCOPF ---
pp.rundcopp(ssp)
dc_comp = ssp.res_gen.copy()
gb_res = ssd.solve(disable_sfr=True, disable_ramp=True, info=False)
dc_comp['p_mw(GB)'] = ssd.mva * gb_res['pg'].values
dc_comp['gen'] = ssd.res['gen']
dcp_col = ['gen', 'p_mw', 'p_mw(GB)', 'q_mvar', 'va_degree', 'vm_pu']
dc_comp = dc_comp[dcp_col]

# --- device idx ---
ssa_dg_idx = ssa_key2.dg_idx.dropna().tolist()
ssa_syg_idx = ssa_key2.syg_idx.dropna().tolist()
ssa_gov_idx = ssa_key2.gov_idx.dropna().tolist()
ssa_stg_idx = ssa_key2.stg_idx.dropna().tolist()

# --- online and controllable device idx ---
ctrl_cond = ssa_key2.ctrl & ssa_key2.stg_u.astype(bool)
ssa_dg_idx_ctrl = ssa_key2.dg_idx[ctrl_cond].dropna().tolist()
ssa_syg_idx_ctrl = ssa_key2.syg_idx[ctrl_cond].dropna().tolist()
ssa_gov_idx_ctrl = ssa_key2.gov_idx[ctrl_cond].dropna().tolist()
ssa_stg_idx_ctrl = ssa_key2.stg_idx[ctrl_cond].dropna().tolist()

# fill NaN with False
ssa_key2.fillna(value=False, inplace=True)

# --- setup constants ---
# --- def. functions ---


def get_pe(ssa, gov_idx, dg_idx, ssa_key2):
    """Get the active power (TurbineGov/DG) after TDS, a DataFrame"""
    # TODO: may need to sum the power of same StaticGen
    # --- TurbineGov ---
    pe_syg = ssa.TurbineGov.get(src='pout', idx=gov_idx, attr='v')
    # --- DG ---
    Ip_dg = ssa.DG.get(src='Ipout_y', idx=dg_idx, attr='v')
    v_dg = ssa.DG.get(src='v', idx=dg_idx, attr='v')
    pe_dg = Ip_dg*v_dg
    # --- out ---
    pe = pd.DataFrame()
    pe['idx'] = gov_idx + dg_idx
    pe['pe'] = np.concatenate((pe_syg, pe_dg))
    ldf = pd.merge(left=ssa_key2.rename(columns={'dg_idx': 'idx'}),
                   right=pe, how='right', on='idx')
    rdf = pd.merge(left=ssa_key2.rename(columns={'gov_idx': 'idx'}),
                   right=pe, how='right', on='idx')
    pe['stg_idx'] = ldf['stg_idx'].fillna('') + rdf['stg_idx'].fillna('')
    return pe


def dp_calc(d_syn, idx_ed, intv_ed, ratio=0.1):
    """Calc SFR requirements, scalars, ``dpd_u``and ``dpd_d``, and load forecasted value ``load_exp``"""
    load = d_syn['sload'].iloc[idx_ed*intv_ed:(idx_ed*intv_ed + intv_ed)]
    load_exp = load.mean()
    load_d = load_exp * ratio
    load_u = load_exp * ratio
    return load_u, load_d, load_exp

# --- co-sim constants ---


# length of each interval
intv_ed = 300  # RTED interval, 300s
intv_agc = 4    # AGC interval, 4s
intv_pq = 1     # PQ interval, 1s; alter load and AGC
intv_step = 50  # step change interval; smooth the setpoitns
if case_name == 'npcc':
    intv_step = 100

# number of each interval
n_ed = int(t_total/intv_ed)
n_agc = int(intv_ed/intv_agc)
n_pq = int(intv_agc/intv_pq)
n_step = np.floor(intv_step/intv_agc)

# --- vars ---
# AGC table
agc_table = ssa_key2[['stg_idx', 'dg_idx', 'rg_idx', 'syg_idx',
                      'exc_idx', 'gov_idx', 'gammap', 'ctrl']]
dc_res = ssd.res
agc_table = agc_table.merge(dc_res[['gen', 'bu', 'bd']].rename(columns={'gen': 'stg_idx'}),
                            on='stg_idx', how='left')
agc_table['paux'] = 0

# AGC power of each unit
agc_in = pd.DataFrame(columns=['stg_idx'] + list(np.arange(0, t_total, 4)))
agc_in['stg_idx'] = agc_table['stg_idx']
agc_out = agc_in.copy()

# ACE vars for PI controller
ACE_integral = 0
ACE_raw = 0
Kp = 0.1  # 0.05
Ki = 0.1
# SFR boundary and total AGC input
sfr_res_data = -1 * np.ones((int(np.ceil(t_total / intv_agc)), 5))

# initial load value
ssa_p0 = ssa.PQ.p0.v.copy()
ssa_q0 = ssa.PQ.q0.v.copy()
ssa_pq_idx = ssa.PQ.idx.v
ssa_p0_sum = ssa_p0.sum()

# idx
ac_res = runopp_map(ssp, ssa_key)

cond_sch_gov = ac_res.gov_idx.fillna(False).astype(bool)
cond_sch_dg = ac_res.dg_idx.fillna(False).astype(bool)
cond_agc_gov = agc_table.ctrl * agc_table.gov_idx.fillna(False).astype(bool)
cond_agc_dg = agc_table.ctrl * agc_table.dg_idx.fillna(False).astype(bool)

sch_gov_idx = ac_res.gov_idx[cond_sch_gov].tolist()
sch_dg_idx = ac_res.dg_idx[cond_sch_dg].tolist()
agc_gov_idx = agc_table.gov_idx[cond_agc_gov].tolist()
agc_dg_idx = agc_table.dg_idx[cond_agc_dg].tolist()

# EV results, random select 3 EVs of each SOC level to record SoC
ridx = sse.ev.sample(frac=1, random_state=sse.config["seed"]).groupby('sx', sort=False).head(100).index
ev_soc_data = -1 * np.ones((t_total, len(ridx)))
ev_agc_data = -1 * np.ones((t_total, len(ridx)))
ev_na_data = -1 * np.ones((t_total, len(ridx)))
ev_soc_data[0] = sse.ev["soc"].iloc[ridx]
ev_agc_data[0] = sse.ev["agc"].iloc[ridx]
ev_na_data[0] = sse.ev["na"].iloc[ridx]

# dispatch results
rted_res = {}
