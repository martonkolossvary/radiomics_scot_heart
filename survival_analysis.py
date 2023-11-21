# %%
import pandas as pd
import numpy as np
from sklearn import set_config
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest, ComponentwiseGradientBoostingSurvivalAnalysis
from sksurv.metrics import concordance_index_ipcw, cumulative_dynamic_auc, integrated_brier_score
from sklearn.inspection import permutation_importance
from sklearn.model_selection import RepeatedStratifiedKFold
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
The function calculates the Harrell's concordance index for each feature in X using CoxPHSurvivalAnalysis model.
The function returns an array of scores, where each score corresponds to a feature in X.
"""


def fit_and_score_features(X, y):
    n_features = X.shape[1]
    scores = np.empty(n_features)
    m = CoxPHSurvivalAnalysis()
    for j in range(n_features):
        Xj = X[:, j : j + 1]
        m.fit(Xj, y)
        scores[j] = m.score(Xj, y)
    return scores


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
    plt.axhline(mean_auc, color=color, linestyle="--")
    plt.xlim([0, 8.8])
    plt.ylim([0.67, 0.85])  # 0.45 for individual, 0.6 for multivariate
    plt.legend()


"""
This code plots the coefficients of each feature for each alpha value of penalized Cox models. The top 10 coefficients are highlighted.
The top 10 features with the largest coefficients are printed, along with their corresponding coefficient value.
"""


def plot_coefficients(coefs, n_highlight):
    _, ax = plt.subplots(figsize=(9, 6))
    n_features = coefs.shape[0]
    alphas = coefs.columns
    for row in coefs.itertuples():
        ax.semilogx(alphas, row[1:], ".-", label=row.Index)

    alpha_min = alphas.min()
    top_coefs = coefs.loc[:, alpha_min].map(abs).sort_values().tail(n_highlight)
    for name in top_coefs.index:
        coef = coefs.loc[name, alpha_min]
        plt.text(alpha_min, coef, name + "   ", horizontalalignment="right", verticalalignment="center")

    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.grid(True)
    ax.set_xlabel("alpha")
    ax.set_ylabel("coefficient")


# %%
## Fit Cox model -----
### Clinical model
estimator0 = CoxPHSurvivalAnalysis()
estimator0.fit(d_x.loc[:, covars[0:3]], d_y)
pd.Series(np.exp(estimator0.coef_), index=d_x.loc[:, covars[0:3]].columns)

estimator0.score(d_x.loc[:, covars[0:3]], d_y)  # Harrell’s concordance index
# concordance_index_censored(d_y['Status'], d_y['Survival_in_years'], estimator0.predict(d_x.loc[:, covars[0:3]]))
concordance_index_ipcw(d_y, d_y, estimator0.predict(d_x.loc[:, covars[0:3]]))
scores = fit_and_score_features(d_x.loc[:, covars[0:3]].values, d_y)  # C-index for each feature
pd.Series(scores, index=d_x.loc[:, covars[0:3]].columns)

plt.figure()
for i, col in enumerate(covars[0:3]):  # Individual time dependent AUC values
    plot_cumulative_dynamic_auc(d_x.loc[:, col], col, color=f"C{i}")
    ret = concordance_index_ipcw(d_y, d_y, d_x.loc[:, col], tau=times[-1])
plt.savefig("Clinical.pdf", dpi=600)
### Clinical + Volume model
estimator1 = CoxPHSurvivalAnalysis()
estimator1.fit(d_x.loc[:, covars], d_y)
pd.Series(np.exp(estimator1.coef_), index=d_x.loc[:, covars].columns)

estimator1.score(d_x.loc[:, covars], d_y)  # Harrell’s concordance index
concordance_index_ipcw(d_y, d_y, estimator1.predict(d_x.loc[:, covars]))
scores = fit_and_score_features(d_x.loc[:, covars].values, d_y)  # C-index for each feature
pd.Series(scores, index=d_x.loc[:, covars].columns)

plt.figure()
for i, col in enumerate(covars[3:6]):  # Individual time dependent AUC values
    plot_cumulative_dynamic_auc(d_x.loc[:, col], col, color=f"C{i}")
    ret = concordance_index_ipcw(d_y, d_y, d_x.loc[:, col], tau=times[-1])
plt.savefig("Volume.pdf", dpi=600)

plt.figure()
for i, col in enumerate(covars):  # Individual time dependent AUC values
    plot_cumulative_dynamic_auc(d_x.loc[:, col], col, color=f"C{i}")
    ret = concordance_index_ipcw(d_y, d_y, d_x.loc[:, col], tau=times[-1])

### Clinical + Volume + Radiomics model
estimator2 = CoxPHSurvivalAnalysis()
estimator2.fit(d_x, d_y)
pd.Series(np.exp(estimator2.coef_), index=d_x.columns)

estimator2.score(d_x, d_y)  # Harrell’s concordance index
concordance_index_ipcw(d_y, d_y, estimator2.predict(d_x))
scores = fit_and_score_features(d_x.values, d_y)  # C-index for each feature
pd.Series(scores, index=d_x.columns)

plt.figure()
for i, col in enumerate(eigen):  # Individual time dependent AUC values
    plot_cumulative_dynamic_auc(d_x.loc[:, col], col, color=f"C{i}")
    ret = concordance_index_ipcw(d_y, d_y, d_x.loc[:, col], tau=times[-1])
plt.savefig("Radiomics.pdf", dpi=600)
# https://medium.com/the-researchers-guide/survival-analysis-in-python-km-estimate-cox-ph-and-aft-model-5533843c5d5d

# %%
### Compare models ----
d_xmodels = pd.DataFrame(
    {
        "Clinical": estimator0.predict(d_x.loc[:, covars[0:3]]),
        "Clinical + Plaque volume": estimator1.predict(d_x.loc[:, covars]),
        "Clinical + Plaque volume + Eigen radiomics": estimator2.predict(d_x),
    }
)

for i, col in enumerate(d_xmodels.columns):  # Individual time dependent AUC values
    plot_cumulative_dynamic_auc(d_xmodels.loc[:, col], col, color=f"C{i}")
    ret = concordance_index_ipcw(d_y, d_y, d_xmodels.loc[:, col], tau=times[-1])
plt.savefig("Combined_volume.pdf", dpi=600)

# %%
## Fit RF model -----
# https://scikit-survival.readthedocs.io/en/stable/user_guide/random-survival-forest.html
CV = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)
estimator_RF = RandomSurvivalForest(n_jobs=-1, random_state=42)

### RF model -----
importances = {}
for i in vars:
    importances[i] = []

for train, test in CV.split(d_x, d_y["Status"]):
    estimator_RF.fit(d_x.iloc[train], d_y[train])
    result_RF = permutation_importance(
        estimator_RF, d_x.iloc[test], d_y[test], n_repeats=100, n_jobs=-1, random_state=42
    )
    for i in range(len(result_RF["importances"])):
        importances[vars[i]] = np.append(importances[vars[i]], result_RF["importances"][i])

mean_importance = {k: np.mean(v) for k, v in importances.items()}
sd_importance = {k: np.std(v) for k, v in importances.items()}
out_mean = pd.Series(mean_importance)
out_sd = pd.Series(sd_importance)

pd.DataFrame(out_mean).to_csv("mean_importance.csv")
pd.DataFrame(out_sd).to_csv("sd_importance.csv")

# %%
## Fit LASSO model -----
cox_lasso = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=np.logspace(-4.5, 1, 10), fit_baseline_model=True)
cox_lasso.fit(d_x, d_y)

coefficients_lasso = pd.DataFrame(cox_lasso.coef_, index=d_x.columns, columns=cox_lasso.alphas_)
plot_coefficients(coefficients_lasso, n_highlight=10)

# %%
# Fit AFT model -----
est_aft_ls = ComponentwiseGradientBoostingSurvivalAnalysis(
    loss="coxph", n_estimators=300, learning_rate=1.0, random_state=0
).fit(d_x, d_y)
cindex = est_aft_ls.score(d_x, d_y)
print(round(cindex, 3))

est_aft_ls.feature_importances_
chf_funcs = est_aft_ls.predict_cumulative_hazard_function(d_x.iloc[:10])
for fn in chf_funcs:
    plt.step(fn.x, fn(fn.x), where="post")

plt.ylim(0, 1)
plt.show()

surv_funcs = est_aft_ls.predict_survival_function(d_x.iloc[:10])
for fn in surv_funcs:
    plt.step(fn.x, fn(fn.x), where="post")

plt.ylim(0, 1)
plt.show()
# %%
