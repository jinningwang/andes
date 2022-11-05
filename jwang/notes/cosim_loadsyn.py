# --- read load from ssa ---
ptotal = ssa.PQ.as_df()['p0'].sum()
print('ANDES total load is:',ptotal.round(4),'p.u.')

# --- scale load curve ---
ddata_path = '/case/dsyn.csv'
ddata = dir_path + ddata_path
d_syn = pd.read_csv(ddata)

caseH = 18

np.random.seed(2022)
col = ['h10', 'h18', 'a10', 'a18']
col1 = ['h10', 'h18']
col2 = ['a10', 'a18']
# d_syn['a10'].iloc[200:650] *= 0.5
d_syn[col1] = (d_syn[col1] - d_syn[col1].min()) / d_syn[col1].min() + 0.8
d_syn[col2] = (d_syn[col2] - d_syn[col2].mean()) / (d_syn[col2].max() - d_syn[col2].min())

if case_name == 'ieee39':
    k = 0.1  # the coefficient can be adjusted to fit the case
    d_syn['s10'] = d_syn['h10'] + k * d_syn['a10']
    d_syn['s18'] = d_syn['h18'] + k * d_syn['a18']
    d_syn['sload'] = d_syn['s18']
    d_syn['sload'].iloc[0:300] -= 0.2 * k
    d_syn.loc[500:1200, 'sload'] = d_syn['sload'].iloc[500:1200].rolling(10).mean()
elif case_name == 'npcc':
    k = 0.005
    k2 = 0.0
    kr = 10
    kt = 0.1
    d_syn['s10'] = d_syn['h10'] + k * d_syn['a10']
    d_syn['s18'] = d_syn['h18'] + k * d_syn['a18']
    d_syn['sload'] = d_syn['s18']
    d_syn['sload'] = d_syn['sload'].rolling(kr).mean()
    d_syn['sload'] *= kt
    d_syn['sload'] += 1 - kt

# calculate expected load
step = 300
d_exp = d_syn.groupby(d_syn.index // step).mean().copy()
d_exp['time'] = range(0,3600,300)

# # align starting point of load with starting point of dispatch results
d_syn['sload'].iloc[0] = d_exp['sload'].iloc[0]
if case_name == 'ieee39':
    d_syn['sload'].iloc[1:50] = None
    d_syn['sload'].interpolate(method='linear', inplace=True)
elif case_name == 'npcc':
    d_syn['sload'].iloc[1:10] = None
    d_syn['sload'].interpolate(method='linear', inplace=True)
    d_syn.loc[1:10, 'sload'] = d_syn['sload'].iloc[0]

# --- plot load curve ---
fig_load, ax_load = plt.subplots(figsize=(5, 4))
ax_load.plot(d_syn['time'], d_syn['sload'], color='tab:blue', linestyle='-')
ystep = list(d_exp['sload'])
ystep.insert(0, d_exp['sload'].iloc[0])
ax_load.step(range(0,3900,300), ystep, color='tab:blue', linestyle='--')
ax_load.set_xlim([0, 3600])
ax_load.legend(['Actual load', 'Forecasted load'])
ax_load.set_title(f'Load profile at {caseH}H')
ax_load.set_ylabel('ratio')
ax_load.set_xlabel('Time [s]')
