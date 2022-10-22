if save_data:
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
