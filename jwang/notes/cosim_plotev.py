import matplotlib.font_manager

sse.plot_agc()

plt.style.use('default')

agc_ev_idx = ev_agc[ev_agc.abs().sum(axis=1) >= 1].index

fig_soc, ax_soc = plt.subplots(figsize=(5, 4))
for i in agc_ev_idx:
    ax_soc.plot(range(end_time+1), ev_soc.iloc[i])
ax_soc.set_xlabel('Time [s]')
ax_soc.set_ylabel('SOC [%]')
ax_soc.set_title('AGC EVs SOC')
ax_soc.set_xlim([0, end_time])

sse.plot(style='ieee')
