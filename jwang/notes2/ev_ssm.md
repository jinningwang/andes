# EV SSM

## Intro

Reference:

1. M. Wang et al., "State Space Model of Aggregated Electric Vehicles for Frequency Regulation," in IEEE Transactions on Smart Grid, vol. 11, no. 2, pp. 981-994, March 2020, doi: 10.1109/TSG.2019.2929052.

An EV connected to the charging station has three connected status: 1) Charging state (CS); 2) Idle state (IS); 3)Discharging state (DS).

Ancillary services from EVs can be achieved by altering the conencted status. In detail, regualtion up and down can be achieved by

(a) CS $\rightarrow$ IS (b) IS $\rightarrow$ DS (c) DS $\rightarrow$ IS (d) IS $\rightarrow$ CS

It should be noted that the other swiching modes can be realized by combination of the four modes. For example, CS $\rightarrow$ DS is equivalent to (a) plus (b), and DS $\rightarrow$ CS is equivalent to (c) plus (d).

Considering the SoC as multiple levels $SoC_1$, $SoC_2$, ..., $SoC_{Ns}$, an EV can be classified in to a state that are defined by the SoC levels and connected status.

For a population of EVs, the EVs can be described by the state vector $\mathbf{x}(k)$, where each element stands for the propotion of a state.

Test markdown

## State transition

$\mathbf{x}(k+1) = \mathbf{A}\mathbf{x}(k)+\mathbf{B}\mathbf{u}(k)+\mathbf{C}\mathbf{v}(k)$

$\mathbf{y}(k) = \mathbf{D}\mathbf{x}(k)$

$\mathbf{A}$ is the state transition matrix, $3N_s \times3N_s$, and can be estimated by (1) historical data (2) analytical method.

$\mathbf{x}(k)$, $3N_s \times1$; portion of each states out of online EVs.

$\mathbf{B} = [-\mathbf{I}_{1\times N_s}, \mathbf{I}_{1\times N_s}, \mathbf{0}_{1\times N_s}] $,

$\mathbf{u}(k)$, $N_s \times1$; control matrix for modes (a) and (d), >0 means (a) and <0 means (d).

$\mathbf{C} = [\mathbf{0}_{1\times N_s}, -\mathbf{I}_{1\times N_s}, \mathbf{I}_{1\times N_s}] $,

$\mathbf{v}(k)$, $N_s \times1$; control matrix for modes (b) and (c), >0 means (b) and <0 means (c).

$\mathbf{D} = P_{ave} N_{e} [-\mathbf{1}_{1\times N_s}, \mathbf{0}_{1\times N_s}, \mathbf{1}_{1\times N_s}] $

In this form, control action $\mathbf{B}\mathbf{u}(k)+\mathbf{C}\mathbf{v}(k)$ will impact on $\mathbf{x}(k+1)$. However, the control action should be effective ***before*** the effect of $\mathbf{A}$. As a result, the SSM should be revised as:

$\mathbf{x}(k+1) = \mathbf{A}\left(\mathbf{x}(k)+\mathbf{B}\mathbf{u}(k)+\mathbf{C}\mathbf{v}(k)\right)$

$\mathbf{y}(k) = \mathbf{D}\left(\mathbf{x}(k)+\mathbf{B}\mathbf{u}(k)+\mathbf{C}\mathbf{v}(k)\right)$

NOTE:

1. in this format, the rated charging and discharging power are mixed together, can be inaccurate.

TODO: correction for random traveling behavior?

## Output

### Responding modes

(a) CS $\rightarrow$ IS, $ P_a = (\mathbf{D_a}-\mathbf{D})\mathbf{x}(k)$

(b) IS $\rightarrow$ DS, $ P_b = (\mathbf{D_b}-\mathbf{D_a})\mathbf{x}(k)$

(c) DS $\rightarrow$ IS, $ P_c = (\mathbf{D_c}-\mathbf{D})\mathbf{x}(k)$

(d) IS $\rightarrow$ CS, $ P_d = (\mathbf{D_c}-\mathbf{D_d})\mathbf{x}(k)$

$\mathbf{D_a} = P_{ave} N_{e} [\mathbf{0}_{1\times N_s}, \mathbf{0}_{1\times N_s}, \mathbf{1}_{1\times N_s}] $

$\mathbf{D_b} = P_{ave} N_{e} [\mathbf{1}_{1\times N_s}, \mathbf{1}_{1\times N_s}, \mathbf{1}_{1\times N_s}] $

$\mathbf{D_c} = P_{ave} N_{e} [-\mathbf{1}_{1\times N_s}, \mathbf{0}_{1\times N_s}, \mathbf{0}_{1\times N_s}] $

$\mathbf{D_d} = P_{ave} N_{e} [-\mathbf{1}_{1\times N_s}, -\mathbf{1}_{1\times N_s}, -\mathbf{1}_{1\times N_s}] $

NOTE:

1.$\mathbf{D_c}$ is different from the paper

TODO: The limits may need revision: ***overcharge*** and ***overdischarge*** are out of service

TODO: low charged EVs are out of service and forced to charge

### Output range

$ P_u = \mathbf{D_b}\mathbf{x}(k)$

$ P_l = \mathbf{D_d}\mathbf{x}(k)$

TODO: $\overline{PR^U}$, $\overline{PR^D}$

## Control signal

### Step1

Input signals: $\mathbf{u}$, $\mathbf{v}$, ($N_s \times1$)

If $P_i > 0$ (RegUp):

First (a) CS $\rightarrow$ IS, then (b) IS $\rightarrow$ DS. The EVs having higher SoC levels are more likely to be altered.

$ r_u = min(P_i, P_a)/(P_{ave}N_e)$

$ u_j = min(r_u -\sum_{h=j+1}^{N_s} x_{h}, x_{j}), j=[1, ..., N_s] $

$ r_v = max(P_i - P_a, 0)/(P_{ave}N_e)$

$ v_j = min(r_v -\sum_{h=j+1}^{N_s}(x_{h+N_s}+ u_h), x_{j+N_s}+ u_j), j=[1, ..., N_s] $

If $P_i < 0$ (RegDn):

First (c) DS $\rightarrow$ IS, then (d) IS $\rightarrow$ CS. The EVs having lower SoC levels are more likely to be altered.

$ r_v = max(P_i, P_c)/(P_{ave}N_e)$

$ v_j = max(r_v +\sum_{h=1}^{j-1}(x_{h-1+2N_s}), -x_{j+2N_s}), j=[1, ..., N_s] $

$ r_u = min(P_i - P_c, 0)/(P_{ave}N_e)$

$ u_j = max(r_u -\sum_{h=1}^{j} v_{h}-\sum_{h=1}^{j-1} u_{h}, -x_{j+N_s}), j=[1, ..., N_s] $

If $P_i = 0$:

$\mathbf{u} = \mathbf{0}$

$\mathbf{v} = \mathbf{0}$

NOTE:

1. Subscripts are different from the paper

### Step2

Control signals: $\mathbf{u}_{s}$, $\mathbf{v}_{s}$, ($N_s \times1$), convert the input signals to a probability.

If $P_i > 0$ (RegUp):

$ u_{s,j} = min(u_j/x_j, 1)$

$ v_{s,j} = min(v_j/(x_{j + N_s}+ u_j), 1)$

If $P_i = 0$:

$\mathbf{u}_{s} = \mathbf{0}$

$\mathbf{v}_{s} = \mathbf{0}$

If $P_i < 0$ (RegDn):

$ v_{s,j} = min(-v_j/(x_{j +2N_s}+ u_j), 1)$

$ u_{s,j} = min(-u_j/x_{j+N_s}- v_j, 1)$

### Step3

Add a direction element into $\mathbf{u}_{s}$, $\mathbf{v}_{s}$.

If $P_i > 0$ (RegUp):

$ u_{s,N_s+1} = 1$, $ v_{s,N_s+1} = 1$

If $P_i = 0$ (RegDn):

$ u_{s,N_s+1} = 0$, $ v_{s,N_s+1} = 0$

If $P_i < 0$ (RegDn):

$ u_{s,N_s+1} = -1,$ $ v_{s,N_s+1} = -1$

### Step4

fcs: C -> I

fis: I -> C

fds: D -> I

TODO: The results do not make much sense?

## Sensitivity

Impacts of EV penetration rate on the performance of AGC response

## Case Study

### Case

1. EVs do not provide SFR
2. EVs minaly provide RegUp; EVs are in high SoC level, demanded charged EVs are switched from CS to IS or DS
3. EVs minaly provide RegDn; EVs are in low SoC level, demanded charged EVs are switched from DS to IS or CS

Two scenarios: 10 H and 18 H.

Expected results: When EV providing SFR, better AGC performance; Pr

## Misc

flow_chart:

```{python}


prep grid data:

ADNES: topology,  gen. limits, ramp. limits, line limits,

Outside: gen. cost, ramp. cost,


for $t_{OPF}$ in T (interval: 5min; total: 1h; [n=12]):

    aggregate EV data (from SSM), generate $PR_{e,i,u,t}$

    Do OPF, generate $PG_{i, t}$, $PR_{g, i, u, t}$, $PR_{g, i, d, t}$


    for t in $t_{OPF}$ (interval: 4s; total: 5min; [n=75]):

        Update data into dynamic system:

            # Note, constant power model should be used in TDS.

            # Use TimeSeries as the load. 

            power change: TGOV1.paux0

            load change: 


        Run TDS: generate SFR mileage

```

Co-Sim list:

```{python}

for $t_{OPF}$ in T:

    EVA report $pru_{max}$ $prd_{max}$, eqn xxx

    TCC do OPF, eqn xxx

    Assign dispatch signal to generation units


    for t in $t_{AGC}$:

        Assign AGC signal to AGC units

        Run TDS

```
