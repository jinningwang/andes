# Reserve some capacity to avoid TDS crush
ssp.gen.max_p_mw = ssp.gen.max_p_mw

# store original generator data
ssp_gen0 = ssp.gen.copy()

ssp.gen['p_mw'][ssp.gen.name==ev_idx] = sse.Ptc
ssd.gen['p0'][ssd.gen.idx == ev_idx] = sse.Ptc / ssd.mva
ssa.StaticGen.set(src='p0', attr='v', idx=ev_idx, value=sse.Ptc / ssa.config.mva)

mdl_status = []

for end_time in range(3600):  # t_total
    # --- interval RTED ---
    if end_time % intv_ed == 0:
        idx_ed = end_time // intv_ed
        # --- update load ---
        sfrur, sfrdr, load_exp = dp_calc(d_syn, idx_ed, intv_ed, rsfr)
        ssp.load['scaling'] = load_exp
        ssp.gen['p_mw'][ssp.gen.name==ev_idx] = sse.Ptc
        ssd.gen['p0'][ssd.gen.idx == ev_idx] = sse.Ptc
        ssd.load['sf'] = load_exp
        ssd.update_dict()

        # --- RTED, update gen limits after SFR ---
        # set previous setpoints with `pe` from TDS
        if end_time > 0:
            vl = pd.merge(left=ssd.gen, how='left', on='idx',
                          right=ssp_res[['stg_idx', 'p']].rename(columns={'stg_idx':'idx'}))['p']
            ssd.gen.p_pre = vl.values

        # def. SFR requirements and calc. EV SFR capacities
        [prumax, prdmax] = sse.g_frc()
        # def. percentage of EV SFR capacities
        ssd.def_type2([ev_idx], [prumax*rru], [prdmax*rrd])
        ssd.def_sfr(sfrur=sfrur*ssa_p0_sum, sfrdr=sfrdr*ssa_p0_sum)

        # solve RTED
        if end_time == 0:
            dcres = ssd.solve(disable_ramp=True, info=False)
        else:
            dcres = ssd.solve(info=False)

        # reserve SFR and ramp from Generator limits in ``ssp``
        ssp_gen = pd.merge(left=ssp.gen.rename(columns={'name': 'stg_idx'}),
                           right=dcres.rename(columns={'gen': 'stg_idx'}),
                           on='stg_idx', how='left')
        # SFR limits
        ssp_gen['max_sfr'] = ssp_gen.max_p_mw - ssp_gen.pru * ssp.sn_mva
        ssp_gen['min_sfr'] = ssp_gen.min_p_mw + ssp_gen.prd * ssp.sn_mva
        # ramp limits
        if end_time > 0:
            p_pre_pp = pd.merge(left=ssp.gen.rename(columns={'name': 'stg_idx'}),
                                right=ssp_res[['stg_idx', 'p']],
                                on='stg_idx', how='left')['p']
            ssp_gen['max_ramp'] = ssp.sn_mva * (np.array(p_pre_pp) + np.array(ssd.gen['ramp_5']))
            ssp_gen['min_ramp'] = ssp.sn_mva * (np.array(p_pre_pp) - np.array(ssd.gen['ramp_5']))
            # alter generator limits
            ssp.gen.max_p_mw = ssp_gen[['max_sfr', 'max_ramp']].min(axis=1)
            ssp.gen.min_p_mw = ssp_gen[['min_sfr', 'min_ramp']].max(axis=1)
        else:
            # alter generator limits
            ssp.gen.max_p_mw = ssp_gen['max_sfr']
            ssp.gen.min_p_mw = ssp_gen['min_sfr']

        # --- ACOPF, modify setpoints ---
        # store setpoints
        if end_time > 0:
            p0 = ssp_res['p'].values  # store setpoints
        else:
            p0 = [0] * ssa_key2.shape[0]

        # solve ACOPF
        ssp_res = runopp_map(ssp, ssa_key)  # ACOPF resutls
        ssp_res['p0'] = p0                  # last setpoints
        ssp_res.fillna(False, inplace=True)  # Fill NA wil False

        # reset Generator limtis to normal limits
        ssp.gen.max_p_mw = ssp_gen0.max_p_mw
        ssp.gen.min_p_mw = ssp_gen0.min_p_mw
        mdl_status.append(ssd.mdl.status)

not_solved = []
error_status = []
for i, s in enumerate(mdl_status):
    if s != 2:
        not_solved.append(i)
        error_status.append(s)

print('RTED not solved:{}'.format(not_solved))
print('They run into: {}'.format(error_status))
