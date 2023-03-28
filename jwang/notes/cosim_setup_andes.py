# --- EV Aggregator ---
if case_name == 'npcc':
    sse = ev_ssm(ts=caseH, N=10000, step=1, tp=40,
             lr=0.1, lp=60, seed=2022, name="EV1",
             n_pref=1, is_report=True,
             tt_mean=0.2, tt_var=0.05, tt_lb=0, tt_ub=0.4,
             ict=ict, ecc=ecc, agc=agc)
    sse.load_A("Aest.csv")
    case_path = '/case/npcc_ev2.xlsx'
elif case_name == 'ieee39':
    sse = ev_ssm(ts=caseH, N=50000, step=1, tp=40,
                lr=0.1, lp=60, seed=2022, name="EVA",
                n_pref=1, is_report=True,
                tt_mean=0.2, tt_var=0.05, tt_lb=0, tt_ub=0.4,
                ict=ict, ecc=ecc, agc=agc)
    sse.load_A("Aest.csv")
    case_path = '/case/ieee39_ev2.xlsx'
# historical data
ev_num = pd.read_csv("ev_num.csv")

# --- ANDES case ---
dir_path = os.path.abspath('..')

case = dir_path + case_path
ssa = andes.load(case,
                 setup=False,
                 no_output=False,
                 default_config=False)

ssa.add("Output", dict(model='ACEc', varname='ace'))
ssa.add("Output", dict(model='COI', varname='omega'))
ssa.add("Output", dict(model='TGOV1N', varname='pout'))
ssa.add("Output", dict(model='TGOV1N', varname='pref'))
ssa.add("Output", dict(model='TGOV1N', varname='paux'))

ssa.setup()

# Set output mode as 'manual', turn off TDS progress bar
ssa.TDS.config.save_mode = 'manual'
ssa.TDS.config.no_tqdm = 1
ssa.TDS.limit_store = 1
ssa.TDS.save_every = 0

# Set load as constant load.
ssa.PQ.config.p2p = 1
ssa.PQ.config.q2q = 1
ssa.PQ.config.p2z = 0
ssa.PQ.config.q2z = 0
ssa.PQ.pq2z = 0

# Turn on ``numba`` can accelerate TDS.
ssa.config.numba = 1
