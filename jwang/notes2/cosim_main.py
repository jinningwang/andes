# Reserve some capacity to avoid TDS crush
ssp.gen.max_p_mw = ssp.gen.max_p_mw

# store original generator data
ssp_gen0 = ssp.gen.copy()

ssp.gen['p_mw'][ssp.gen.name==ev_idx] = sse.Ptc
ssd.gen['p0'][ssd.gen.idx == ev_idx] = sse.Ptc / ssd.mva

for end_time in range(t_total):  # t_total
    # --- interval RTED ---
    if end_time % intv_ed == 0:
        idx_ed = end_time // intv_ed
        # --- update load ---
        sfrur, sfrdr, load_exp = dp_calc(d_syn, idx_ed, intv_ed, rsfr)
        ssp.load['scaling'] = load_exp
        ssp.gen['p_mw'][ssp.gen.name==ev_idx] = sse.Ptc
        ssd.gen['p0'][ssd.gen.idx == ev_idx] = sse.Ptc / ssd.mva
        ssd.load['sf'] = load_exp
        ssd.update_dict()

        # --- RTED, update gen limits after SFR ---
        # set previous setpoints with `pe` from TDS
        if end_time > 0:
            pe_tds = get_pe(ssa, ssa_gov_idx, ssa_dg_idx, ssa_key2)
            pe_tds = pe_tds.merge(ssa_key2, on='stg_idx',
                                  how='right').groupby('stg_idx', as_index=False).sum()
            vl = pd.merge(left=ssd.gen, how='left', on='idx',
                          right=pe_tds[['stg_idx', 'pe']].rename(columns={'stg_idx': 'idx'}))['pe']
            ssd.gen.p_pre = vl.values

        # def. SFR requirements and calc. EV SFR capacities
        [prumax, prdmax] = sse.g_frc()
        # def. percentage of EV SFR capacities
        ssd.def_type2([ev_idx], [prumax*rru/ssd.mva], [prdmax*rrd/ssd.mva])
        ssd.def_sfr(sfrur=sfrur*ssa_p0_sum, sfrdr=sfrdr*ssa_p0_sum)

        # solve RTED
        if end_time == 0:
            dcres = ssd.solve(disable_ramp=True)
        else:
            dcres = ssd.solve()

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
                                right=pe_tds[['stg_idx', 'pe']],
                                on='stg_idx', how='left')['pe']
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

    # --- interval AGC ---
    if end_time % intv_agc == 0:
        idx_agc = end_time // intv_agc - idx_ed * n_agc
        # --- allocate AGC ---
        # assign participation factor `bu`, `bd`
        agc_table.drop(['bu', 'bd'], axis=1, inplace=True)
        agc_table = agc_table.merge(dcres[['gen', 'bu', 'bd']].rename(columns={'gen': 'stg_idx'}),
                                    on='stg_idx', how='left')
        # calc. AGC ---
        ACE_input = min(ACE_raw, dcres.pru.sum())
        if ACE_raw >= 0:
            ACE_input = min(ACE_raw, dcres.pru.sum())
            agc_table['paux'] = ACE_input * agc_table['bu'] * agc_table['gammap']
        else:
            ACE_input = max(ACE_raw, -1 * dcres.prd.sum())
            agc_table['paux'] = ACE_input * agc_table['bd'] * agc_table['gammap']
        agc_in[end_time] = agc_table['paux']
        sfr_res[end_time] = [end_time, ACE_raw, dcres.pru.sum(),
                             -1*dcres.prd.sum(), ACE_input]

        # --- record AGC ---
        if end_time > 0:
            gref = ssa.TurbineGov.get(src='pref', attr='v', idx=agc_gov_idx)
            gout = ssa.TurbineGov.get(src='pout', attr='v', idx=agc_gov_idx)
            g_aux = gout - gref
            agc_out[end_time] = np.append(g_aux, [sse.Prc/ssa.config.mva])

        # --- assign AGC ---
        # a.SynGen
        # Note: now the condition is controllable & has governor
        ssa.TurbineGov.set(src='paux0', idx=agc_gov_idx, attr='v',
                           value=agc_table.paux.values)
        # b.DG;
        ssa.DG.set(src='pext0', idx=agc_dg_idx, attr='v',
                   value=agc_table.paux.values)
        # c.EV;
        # Note: EV is in group DG, remove EV set if not used for EV SSM
        sse_agc = ssa.config.mva * agc_table[agc_table.stg_idx == ev_idx].paux.values
        # TODO: RenGen

        # --- smooth setpoints ---
        if idx_ed == 0:
            ssp_res['pref'] = ssp_res['p']
        else:
            if idx_agc == 0:
                # record `pe` from TDS in the first AGC interval
                copy = ssp_res.merge(right=pe_tds[['pe', 'stg_idx']], on='stg_idx', how='left')
                ssp_res['pe_tds'] = copy.pe
            idx_step = min((end_time - idx_ed * intv_ed) // intv_agc + 1, n_step)
            ssp_res['pref_step'] = ssp_res.p - ssp_res.p0
            # smooth change threshold: 0.01
            large_index = ssp_res['pref_step'][abs(ssp_res['pref_step']) > 0.01].index
            ssp_res['pref_delta'] = ssp_res['pref_step']
            ssp_res['pref_delta'].iloc[large_index] = ssp_res['pref_step'].iloc[large_index] / n_step * idx_step
            ssp_res['pref'] = ssp_res.p0 + ssp_res.pref_delta

            # a.SynGen
            ssa.TurbineGov.set(src='pref0', idx=sch_gov_idx,
                               attr='v', value=ssp_res.pref[cond_sch_gov].values)
            # b.DG
            ssa.DG.set(src='pref0', idx=sch_dg_idx,
                       attr='v', value=ssp_res.pref[cond_sch_dg].values)

    # --- intv_pq: alter load, run TDS ---
    # Initially, alter StaticGen: p0 and q0, run PFlow
    # Otherwise, alter Ppf and Qpf
    if end_time == 0:
        stg_opf_idx = ssp_res.stg_idx[ssp_res.controllable].tolist()
        stg_opf_val = ssp_res.p[ssp_res.controllable].tolist()
        stg_opf_v = ssp_res.vm_pu[ssp_res.controllable].tolist()
        ssa.StaticGen.set(src='p0', idx=stg_opf_idx, attr='v', value=stg_opf_val)
        ssa.StaticGen.set(src='v0', idx=stg_opf_idx, attr='v', value=stg_opf_v)
        # initial load point set as the dispatch point
        ssa.PQ.set(src='p0', idx=ssa_pq_idx, attr='v',
                   value=ssa_p0 * load_exp)
        ssa.PQ.set(src='q0', idx=ssa_pq_idx, attr='v',
                   value=ssa_q0 * load_exp)
        ssa.PFlow.run()
    else:
        ssa.PQ.set(src='Ppf', idx=ssa_pq_idx, attr='v',
                   value=ssa_p0 * d_syn['sload'].iloc[end_time])
        ssa.PQ.set(src='Qpf', idx=ssa_pq_idx, attr='v',
                   value=ssa_q0 * d_syn['sload'].iloc[end_time])
        # run `sse`
        sse.run(tf=caseH+end_time/3600, Pi=sse_agc[0],
                is_updateA=False, is_rstate=True,
                is_test=False, disable=True)
        ev_soc[end_time] = sse.ev.soc
        ev_agc[end_time] = sse.ev.agc
        sse.report(is_report=False)
        ssa.EV2.set(src='pref0', idx=ssp_res.dg_idx[ssp_res.stg_idx == ev_idx].values[0],
                    attr='v', value=sse.Ptc / ssa.config.mva)

    # run TDS
    ssa.TDS.config.tf = end_time
    if end_time == 0:
        ssa.TDS.init()
    ssa.TDS.run()
    # update AGC PI Controller
    ACE_integral = ACE_integral + ssa.ACEc.ace.v.sum()
    ACE_raw = -(Kp*ssa.ACEc.ace.v.sum() + Ki*ACE_integral)

    # ACE_raw = 0  # delete when run TDS
    # break loop if TDS run into error
    if ssa.exit_code != 0:
        raise ValueError(f"TDS error! Exit with {ssa.exit_code}, end at {end_time}s.")
