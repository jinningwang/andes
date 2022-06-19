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
    ldf = pd.merge(left=ssa_key2.rename(columns={'dg_idx':'idx'}),
         right=pe, how='right', on='idx')
    rdf = pd.merge(left=ssa_key2.rename(columns={'gov_idx':'idx'}),
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
intv_step = 50 # step change interval; smooth the setpoitns

# number of each interval
n_ed = int(t_total/intv_ed)
n_agc = int(intv_ed/intv_agc)
n_pq = int(intv_agc/intv_pq)
n_step = np.floor(intv_step/intv_agc)

# --- vars ---
# AGC table
agc_table = ssa_key2[['stg_idx', 'dg_idx', 'rg_idx', 'syg_idx', 
                      'exc_idx', 'gov_idx', 'gammap', 'ctrl']]
dcres = ssd.res
agc_table = agc_table.merge(dcres[['gen', 'bu', 'bd']].rename(columns={'gen': 'stg_idx'}),
                            on='stg_idx', how='left')
agc_table['paux'] = 0

# AGC power of each unit
agc_in = pd.DataFrame(columns=['stg_idx'] + list(np.arange(0, t_total, 4)))
agc_in['stg_idx'] = agc_table['stg_idx']
agc_out = agc_in.copy()

# ACE vars for PI controller
ACE_integral = 0
ACE_raw = 0
Kp = 0.1 # 0.05
Ki = 0.1
# SFR boundary and total AGC input
sfr_res = -1 * np.ones((int(np.ceil(t_total / intv_agc)), 5))

# initial load value
ssa_p0 = ssa.PQ.p0.v.copy()
ssa_q0 = ssa.PQ.q0.v.copy()
ssa_pq_idx = ssa.PQ.idx.v
ssa_p0_sum = ssa_p0.sum()

# EV results
ev_soc = -1 * np.ones((t_total, sse.ev.shape[0]))
ev_agc = -1 * np.ones((t_total, sse.ev.shape[0]))
ev_soc[0] = sse.ev.soc
ev_agc[0] = sse.ev.agc

# idx
ssp_res = runopp_map(ssp, ssa_key)

cond_sch_gov = ssp_res.gov_idx.fillna(False).astype(bool)
cond_sch_dg = ssp_res.dg_idx.fillna(False).astype(bool)
cond_agc_gov = agc_table.ctrl * agc_table.gov_idx.fillna(False).astype(bool)
cond_agc_dg = agc_table.ctrl * agc_table.dg_idx.fillna(False).astype(bool)

sch_gov_idx = ssp_res.gov_idx[cond_sch_gov].tolist()
sch_dg_idx = ssp_res.dg_idx[cond_sch_dg].tolist()
agc_gov_idx = agc_table.gov_idx[cond_agc_gov].tolist()
agc_dg_idx = agc_table.dg_idx[cond_agc_dg].tolist()
