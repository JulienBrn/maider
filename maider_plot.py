import pandas as pd #pandas is for tables with row names and column names
import numpy as np #numpy is to handle arrays of arbitrary dimensions
import scipy #scipy provides functions for statistical analysis (and more)
import seaborn as sns, matplotlib.pyplot as plt #These are for plotting, seaborn does nice easy plots but matplotlib is the basis
import pathlib #This is to handle files
import tqdm #tqdm handles progress bars
import pickle #pickle handles saving and loading of a python variable
import toolbox

base = pathlib.Path("/home/julien/Maider/")
parameter_cols = ["Input", "Sex", "AAV", "t", "DurationsGroup", "Session"]
result_path = "all_results.pkl"

all_results: pd.DataFrame = pickle.load((base/result_path).open("rb")).reset_index()
all_results["df/f z-scored"] = -all_results["df/f z-scored"] #TO REMOVE LATER
all_results["Input"] = np.where(all_results["Input"]=="Input1", "Distribution",
                       np.where(all_results["Input"]=="Input2", "Cue", 
                       np.where(all_results["Input"]=="Input3", "Lick", 
                       None)))


all_results= all_results.groupby(parameter_cols)["df/f z-scored"].mean().reset_index()

print(all_results)
figure_cols = list(set(parameter_cols) - set({"AAV", "Sex", "t", "DurationsGroup"}))


f = toolbox.FigurePlot(all_results, figures=figure_cols, row="AAV", col="Sex", fig_title="{" + "}{".join(figure_cols) + "}", margin_titles=True)
f.pcolormesh(x="t", y="DurationsGroup", value="df/f z-scored")


averaged_result = all_results.groupby(list(set(parameter_cols) - {"DurationsGroup"}))["df/f z-scored"].mean().reset_index()
print(averaged_result)

figure_avg_cols = list(set(figure_cols) - {"Input"})
f_avg = toolbox.FigurePlot(data=averaged_result, figures=figure_avg_cols, row="Input", col="Sex", margin_titles=True, fig_title="{" + "}{".join(figure_avg_cols) + "}" if len(figure_avg_cols) > 0 else "")
f_avg.map(sns.lineplot,  hue="AAV", x="t", y="df/f z-scored", hue_order=["a53t", "empty"])
f_avg = f_avg.add_legend()
# sns.relplot(data=averaged_result, x="t", y="df/f z-scored", hue="Input", row="AAV", col="Sex", facet_kws=dict(margin_titles=True), kind="line")
plt.show()