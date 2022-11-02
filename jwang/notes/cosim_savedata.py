if save_data:
    # --- save data ---
    cosim_out = pd.DataFrame()
    cosim_out['Time'] = ssa.dae.ts.t
    cosim_out['ACE'] = aced
    coif_idx = ssa.TDS.plt._process_yidx(ssa.COI.omega, None)
    cosim_out['freq'] = ssa.TDS.plt.get_values(coif_idx).reshape(-1) * ssa.config.freq

    agc_out.fillna(0, inplace=True)
    agc_out_sort = pd.merge(left=agc_out.rename(columns={'stg_idx':'idx'}),
                            right=ssd.cost[['idx']], on='idx', how='right')
    # --- agc mileage ---
    agc_mile = pd.DataFrame(columns=list(np.arange(intv_agc, t_total, intv_agc)))
    for col_id in np.arange(intv_agc, t_total, intv_agc):
        agc_mile[col_id] = np.abs(agc_out_sort[col_id] - agc_out_sort[int(col_id-intv_agc)])
    agc_mile[agc_mile.columns] *= 100

    cosim_out.to_csv(file_beging + sim_name + '_out.csv', index=False)
    sse_out.to_csv(file_beging + sim_name + '_sse.csv', index=False)
    bu_df.to_csv(file_beging + sim_name + '_bu.csv', index=False)
    bd_df.to_csv(file_beging + sim_name + '_bd.csv', index=False)
    pg_df.to_csv(file_beging + sim_name + '_pg.csv', index=False)
    agc_mile.to_csv(file_beging + sim_name + '_agcm.csv', index=False)
    sfr_res.to_csv(file_beging + sim_name + '_sfr.csv', index=False)
    ev_agc.to_csv(file_beging + sim_name + '_evagc.csv', index=False)
    ev_soc.iloc[0:800].to_csv(file_beging + sim_name + '_evsoc1.csv', index=False)
    ev_soc.iloc[800:].to_csv(file_beging + sim_name + '_evsoc2.csv', index=False)
    sse.ev.to_csv(file_beging + sim_name + '_evdata.csv', index=False)

    import csv
    new_path = open(file_beging + sim_name + '_rted.csv', 'w')
    z = csv.writer(new_path)
    for new_k, new_v in rted_res.items():
        z.writerow([new_k, new_v])
    new_path.close()
    print('Successfully save data.')
