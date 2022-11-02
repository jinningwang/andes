# Reserve some capacity to avoid TDS crush
ssp.gen.max_p_mw = ssp.gen.max_p_mw

# store original generator data
ssp_gen0 = ssp.gen.copy()

ssp.gen['p_mw'][ssp.gen.name == ev_gen] = sse.data["Ptc"]
ssd.gen['p0'][ssd.gen.idx == ev_gen] = sse.data["Ptc"] / ssd.mva
ssa.StaticGen.set(src='p0', attr='v', idx=ev_gen, value=sse.data["Ptc"] / ssa.config.mva)

for end_time in tqdm(range(t_total)):  # t_total
    # --- interval RTED ---
    if end_time % intv_ed == 0:
        idx_ed = end_time // intv_ed
        # --- update load ---
        sfrur, sfrdr, load_exp = dp_calc(d_syn, idx_ed, intv_ed, rsfr)
        ssp.load['scaling'] = load_exp
        ssp.gen['p_mw'][ssp.gen.name == ev_gen] = sse.data["Ptc"]
        ssd.gen['p0'][ssd.gen.idx == ev_gen] = sse.data["Ptc"] / ssd.mva
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
        k1 = ev_num['ne'][(ev_num['time'] >= sse.data["ts"]) & (ev_num['time'] <= sse.data["ts"]+1/12)].mean()
        k0 = ev_num['ne'][ev_num['time'] >= sse.data["ts"]].iloc[0]
        k = k1 / k0
        # estiamte FRC
        [prumax, prdmax] = sse.g_frc()  # check later
        # def. percentage of EV SFR capacities
        ssd.def_type2([ev_gen], [prumax*rru/ssd.mva], [prdmax*rrd/ssd.mva])
        if case_name == 'npcc':
            ssd.gen.loc[ev_gen, 'prumax'] = prumax * rru / ssd.mva
            ssd.gen.loc[ev_gen, 'prdmax'] = prdmax * rrd / ssd.mva
        ssd.def_sfr(sfrur=sfrur*ssa_p0_sum, sfrdr=sfrdr*ssa_p0_sum)

        # solve RTED
        if end_time == 0:
            dc_res = ssd.solve(disable_ramp=True)
        else:
            dc_res = ssd.solve()

        # reserve SFR and ramp from Generator limits in ``ssp``
        ssp_gen = pd.merge(left=ssp.gen.rename(columns={'name': 'stg_idx'}),
                           right=dc_res.rename(columns={'gen': 'stg_idx'}),
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
            p0 = ac_res['p'].values  # store setpoints
        else:
            p0 = [0] * ssa_key2.shape[0]

        # solve ACOPF
        ac_res = runopp_map(ssp, ssa_key, alt=opf_alt)  # ACOPF resutls
        ac_res['p0'] = p0                  # last setpoints
        ac_res.fillna(False, inplace=True)  # Fill NA wil False

        # reset Generator limtis to normal limits
        ssp.gen.max_p_mw = ssp_gen0.max_p_mw
        ssp.gen.min_p_mw = ssp_gen0.min_p_mw

        # store dispatch results
        dc_res_copy = dc_res.copy()
        dc_res_copy['pg'] = ac_res.p
        rted_res[idx_ed] = dc_res_copy.T.to_dict()

    # --- interval AGC ---
    if end_time % intv_agc == 0:
        idx_agc = end_time // intv_agc - idx_ed * n_agc
        # --- allocate AGC ---
        # assign participation factor `bu`, `bd`
        agc_table.drop(['bu', 'bd'], axis=1, inplace=True)
        agc_table = agc_table.merge(dc_res[['gen', 'bu', 'bd']].rename(columns={'gen': 'stg_idx'}),
                                    on='stg_idx', how='left')
        # calc. AGC ---
        ACE_input = min(ACE_raw, dc_res.pru.sum())
        if ACE_raw >= 0:
            ACE_input = min(ACE_raw, dc_res.pru.sum())
            agc_table['paux'] = ACE_input * agc_table['bu'] * agc_table['gammap']
        else:
            ACE_input = max(ACE_raw, -1 * dc_res.prd.sum())
            agc_table['paux'] = ACE_input * agc_table['bd'] * agc_table['gammap']
        agc_in[end_time] = agc_table['paux']
        sfr_res_data[end_time // intv_agc] = [end_time, ACE_raw, dc_res.pru.sum(),
                                              -1*dc_res.prd.sum(), ACE_input]

        # --- record AGC ---
        if end_time > 0:
            gref = ssa.TurbineGov.get(src='pref', attr='v', idx=agc_gov_idx)
            gout = ssa.TurbineGov.get(src='pout', attr='v', idx=agc_gov_idx)
            tmp_df = pd.DataFrame()
            tmp_df['gov_idx'] = agc_gov_idx
            tmp_df[end_time] = gout - gref
            tmp_df2 = tmp_df.merge(right=agc_table,
                                on='gov_idx', how='right')
            tmp_df3 = tmp_df2[['stg_idx', end_time]].merge(right=ssa_key2,
                                    on='stg_idx', how='right')
            agc_out[end_time] = tmp_df3[end_time]

        # --- assign AGC ---
        # a.SynGen
        # Note: now the condition is controllable & has governor
        ssa.TurbineGov.set(src='paux0', idx=agc_gov_idx, attr='v',
                           value=agc_table["paux"][cond_agc_gov].values)
        # b.DG;
        ssa.DG.set(src='Pext0', idx=agc_dg_idx, attr='v',
                   value=agc_table["paux"][cond_agc_dg].values)
        # c.EV;
        # Note: EV is in group DG, remove EV set if not used for EV SSM
        sse_agc = ssa.config.mva * agc_table[agc_table.stg_idx == ev_gen].paux.values
        # TODO: RenGen

        # --- smooth setpoints ---
        if idx_ed == 0:
            ac_res['pref'] = ac_res['p']
        else:
            if idx_agc == 0:
                # record `pe` from TDS in the first AGC interval
                copy = ac_res.merge(right=pe_tds[['pe', 'stg_idx']], on='stg_idx', how='left')
                ac_res['pe_tds'] = copy.pe
            idx_step = min((end_time - idx_ed * intv_ed) // intv_agc + 1, n_step)
            ac_res['pref_step'] = ac_res.p - ac_res.p0
            # smooth change threshold: 0.01
            large_index = ac_res['pref_step'][abs(ac_res['pref_step']) > 0.01].index
            ac_res['pref_delta'] = ac_res['pref_step']
            ac_res['pref_delta'].iloc[large_index] = ac_res['pref_step'].iloc[large_index] / n_step * idx_step
            ac_res['pref'] = ac_res.p0 + ac_res.pref_delta

            # a.SynGen
            ssa.TurbineGov.set(src='pref0', idx=sch_gov_idx,
                               attr='v', value=ac_res.pref[cond_sch_gov].values)
            # b.DG
            ssa.DG.set(src='pref0', idx=sch_dg_idx,
                       attr='v', value=ac_res.pref[cond_sch_dg].values)

    # --- intv_pq: alter load, run TDS ---
    # Initially, alter StaticGen: p0 and q0, run PFlow
    # Otherwise, alter Ppf and Qpf
    if end_time == 0:
        stg_opf_idx = ac_res.stg_idx[ac_res.controllable].tolist()
        stg_opf_val = ac_res.p[ac_res.controllable].tolist()
        stg_opf_v = ac_res.vm_pu[ac_res.controllable].tolist()
        v0 = ssa.StaticGen.get(src='v0', idx=stg_opf_idx, attr='v')
        if not np.mean(stg_opf_v) == 1:
            ssa.StaticGen.set(src='p0', idx=stg_opf_idx, attr='v', value=stg_opf_val)
            ssa.StaticGen.set(src='v0', idx=stg_opf_idx, attr='v', value=stg_opf_v)
        else:
            ssa.StaticGen.set(src='p0', idx=stg_opf_idx, attr='v', value=1.1 * np.array(stg_opf_val))
            ssa.StaticGen.set(src='v0', idx=stg_opf_idx, attr='v', value=v0)
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
        ev_soc_data[end_time] = sse.ev["soc"].iloc[ridx]
        ev_agc_data[end_time] = sse.ev["agc"].iloc[ridx]
        ev_na_data[end_time] = sse.ev["na"].iloc[ridx]

        sse.report(is_report=False)
        ssa.EV2.set(src='pref0', idx=ac_res.dg_idx[ac_res.stg_idx == ev_gen].values[0],
                    attr='v', value=sse.data["Ptc"] / ssa.config.mva)

    # run TDS
    if end_time == 0:
        ssa.TDS.init()
    else:
        ssa.TDS.config.tf = end_time
        ssa.TDS.run()
    # update AGC PI Controller
    ACE_integral = ACE_integral + ssa.ACEc.ace.v.sum()
    ACE_raw = -(Kp*ssa.ACEc.ace.v.sum() + Ki*ACE_integral)

    # ACE_raw = 0  # delete when run TDS
    # break loop if TDS run into error
    if ssa.exit_code != 0:
        raise ValueError(f"TDS error! Exit with {ssa.exit_code}, end at {end_time}s.")

# save Time Series output
ssa.TDS.save_output()
