# %%
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, WeibullAFTFitter

# %%
# INITIALIZE -----
## Load data -----
p = pd.read_csv("")
type = "burden"  # "volume" or "burden"

## Select columns -----
if type == "volume":  # Plaque volume 8y
    eigen = [
        "black_norm_ALL",
        "greenyellow_norm_NCP",
        "pink_norm_NCP",
        "magenta_norm_NCP",
        "brown_norm_NCP",
        "yellow_norm_NCP",
        "brown_norm_CP",
        "blue_norm_CP",
        "turquoise_norm_CP",
        "yellow_norm_CP",
    ]
    covars = ["risk", "cac", "obs", "Total.NCP.vol..mm3.", "Total.CP.vol..mm3.", "Total.LD.NCP.volume..mm3."]
else:  # Plaque burden 8y
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
    covars = ["risk", "cac", "obs", "ncppb", "cppb", "lancppb"]

p["MI_cardiac_death_observed_last"] = p["MI_cardiac_death_observed_last"].astype(bool)
vars = covars + eigen

## Create dataframes -----
d_y = np.rec.fromarrays(
    [p["MI_cardiac_death_observed_last"], p["time_MI_cardiac_death_observed_last"]],
    names=["Status", "Survival_in_years"],
)
d_x = p[vars]
d = p[vars + ["MI_cardiac_death_observed_last", "time_MI_cardiac_death_observed_last"]]

# %%
# FIT Cox MODEL -----
cph = CoxPHFitter().fit(
    d, duration_col="time_MI_cardiac_death_observed_last", event_col="MI_cardiac_death_observed_last"
)
cph.print_summary()
cph.plot_partial_effects_on_outcome(covariates="magenta_norm_NCP", values=[-2, -1, 0, 1, 2], cmap="coolwarm")

# %%
# FIT Weibull MODEL -----
wf = WeibullAFTFitter().fit(
    d, duration_col="time_MI_cardiac_death_observed_last", event_col="MI_cardiac_death_observed_last"
)
wf.print_summary()
wf.plot_partial_effects_on_outcome(covariates="magenta_norm_NCP", values=[-2, -1, 0, 1, 2], cmap="coolwarm")


# %%
