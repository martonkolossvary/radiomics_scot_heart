# %%
import pandas as pd
import numpy as np
from sklearn import set_config
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import (
    concordance_index_ipcw,
    cumulative_dynamic_auc,
    concordance_index_censored,
)
from sklearn.model_selection import KFold, cross_val_predict
import matplotlib.pyplot as plt

# %%
# INITIALIZE -----
## Load data -----
p = pd.read_csv("")

type = "volume"  # "volume" or "burden"

p.loc[:, "Total.NCP.vol..mm3."] = p.loc[:, "Total.NCP.vol..mm3."] / 100
p.loc[:, "Total.CP.vol..mm3."] = p.loc[:, "Total.CP.vol..mm3."] / 100
p.loc[:, "Total.LD.NCP.volume..mm3."] = p.loc[:, "Total.LD.NCP.volume..mm3."] / 100
# p = p[p.loc[:, "time_MI_cardiac_death_observed_last"] <= 8.6]
## Select columns -----
eigen = [
    "brown_norm_ALL",
    "black_norm_ALL",
    "greenyellow_norm_NCP",
    "pink_norm_NCP",
    "magenta_norm_NCP",
    "brown_norm_NCP",
    "yellow_norm_NCP",
    "brown_norm_CP",
]
if type == "volume":  # Plaque volume 8y
    covars = ["risk", "cac", "obs", "Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3."]
else:  # Plaque burden 8y
    covars = ["risk", "cac", "obs", "ncppb", "cppb", "lancppb"]

p["MI_cardiac_death_observed_last"] = p["MI_cardiac_death_observed_last"].astype(bool)
vars = covars + eigen

## Create dataframes -----
d_y = np.rec.fromarrays(
    [p["MI_cardiac_death_observed_last"], p["time_MI_cardiac_death_observed_last"]],
    names=["Status", "Survival_in_years"],
)
d_x = p[vars]

# %%
# FIT MODELS -----
set_config(display="text")
times = np.linspace(0.5, 8.5, 17)  # Time points for time dependent AUC

"""
The function calculates the cumulative dynamic AUC and mean AUC.
It then plots the time-dependent AUC against the years from enrollment.
The function also adds a horizontal line for the mean AUC and sets the x and y limits for the plot.
"""
def plot_cumulative_dynamic_auc(risk_score, label, color=None):
    auc, mean_auc = cumulative_dynamic_auc(d_y, d_y, risk_score, times)
    if auc[0] < 0.5:
        auc = 1 - auc
        mean_auc = 1 - mean_auc

    plt.plot(times, auc, marker="o", color=color, label=label)
    plt.xlabel("Years from enrollment")
    plt.ylabel("Time-dependent C-index")
    # plt.axhline(mean_auc, color=color, linestyle="--")
    plt.xlim([0, 8.8])
    plt.ylim([0.65, 0.85])  # 0.45 for individual, 0.6 for multivariate
    plt.legend()

# %%
## Fit Cox models using 3-fold CV -----
var_sets = [vars[0:3], vars[0:6], vars]
model = CoxPHSurvivalAnalysis()
cv_splits = 3
CV = KFold(n_splits=cv_splits, shuffle=True, random_state=42)

models = {i: model for i in range(1, len(var_sets) + 1)}
scores = {i: [] for i in range(1, len(var_sets) + 1)}
predictions = {i: [] for i in range(1, len(var_sets) + 1)}

for i in range(1, len(var_sets) + 1):
    predictions[i] = cross_val_predict(model, d_x.loc[:, var_sets[i - 1]], d_y, cv=CV)
    print(cumulative_dynamic_auc(d_y, d_y, predictions[i], times))
    print(concordance_index_ipcw(d_y, d_y, predictions[i], tau=times[-1]))

# %%
### Compare models ----
d_xmodels = pd.DataFrame(
    {
        "Clinical": predictions[1],
        "Clinical + Plaque volume": predictions[2],
        "Clinical + Plaque volume + Eigen radiomics": predictions[3],
    }
)

for i, col in enumerate(d_xmodels.columns):  # Individual time dependent AUC values
    plot_cumulative_dynamic_auc(d_xmodels.loc[:, col], col, color=f"C{i}")
    ipwc = concordance_index_ipcw(d_y, d_y, d_xmodels.loc[:, col], tau=times[-1])
    cens = concordance_index_censored(d_y["Status"], d_y["Survival_in_years"], d_xmodels.loc[:, col])
    print(f"{col}: {ipwc}")
    print(f"{col}: {cens}")
plt.savefig("Combined_volume.pdf", dpi=600)
