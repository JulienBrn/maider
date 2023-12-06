
import pandas as pd #pandas is for tables with row names and column names
import numpy as np #numpy is to handle arrays of arbitrary dimensions
import scipy #scipy provides functions for statistical analysis (and more)
import seaborn as sns, matplotlib.pyplot as plt #These are for plotting, seaborn does nice easy plots but matplotlib is the basis
import pathlib #This is to handle files
import tqdm #tqdm handles progress bars
import pickle #pickle handles saving and loading of a python variable

base = pathlib.Path("/home/julien/Maider/")
group_cols = ["Input", "Sex", "Session", "AAV", "DurationsGroup"] #Available columns: "Input", "Sex", "Session", "AAV","Side", "Subject", "Date", "TrialNum", "DurationsGroup"
result_path = "all_results.pkl"
duration_csv_path = "quantiles.csv"
zscore_method = "baseline"
files = list(base.glob("**/*Fluorescence.csv"))
plot_event=False










if not zscore_method in ["baseline", "file"]:
    raise Exception("zscore method not handled")

def get_metadata(filename): #computes session, sex, subject, date, side, aav from filename
    levels = str(filename.relative_to(base)).split("/") #splits the filename into a list of substrings separated by /
    levels = levels[:-2] #selection
    try:#try to do this, if it fails "goto" the except part
        session = int(levels[0][6:]) #selection
        sex = levels[2][0] #selection
        split = levels[2].split("_") #splits the third level into a list of substrings separated by _
        subject = split[0] #selection
        side = split[1] #selection
        date = split[2] #selection
        subject_number = int(subject[1:]) #selection
        aav = "empty" if (subject_number <6 and sex=="M") or (subject_number >=6 and sex=="F") else "a53t"
    except:
        print("Ignoring", levels)
        return [None] * 6
    return [session, sex, subject, date, side, aav]



files = pd.Series(files).to_frame()
files.columns=["File"]
tqdm.tqdm.pandas(desc="Parsing metadata")
files[["Session","Sex", "Subject", "Date", "Side", "AAV"]] = files.progress_apply(lambda row: get_metadata(row["File"]), result_type="expand", axis=1) #applies get_metadata to each filename
files= files.loc[~pd.isna(files["Date"])] #select where date is not None
# files=files.sample(10)
def get_file_data(filename):
    df = pd.read_csv(filename, header=1)
    return df

tqdm.tqdm.pandas(desc="Reading files")
data = files.groupby("File").progress_apply(lambda x: get_file_data(x.iloc[0, :]["File"]))
data = data.merge(files, how="left", on="File") #putting two tables together
print(data)
df=data


events = df["Events"].str.split(";", expand=True)
events = events.stack()
events.index.names=["line_num", "ev_num"]
events.name="Event"
events = events.reset_index()
events[["Input", "_", "End"]] = events["Event"].str.split("*", expand=True)

events = events.merge(df, how="left", right_index=True, left_on="line_num")
events= events.loc[events["End"]=="0", :]
events["is_input2"] = events["Input"] == "Input2"

events["TrialNum"]=events.sort_values("TimeStamp").groupby("File")["is_input2"].cumsum()
events = events.drop(columns="is_input2")
timestamps = events.set_index(["File", "TrialNum", "Input"])["TimeStamp"].unstack("Input")
timestamps["DurationLick"] = timestamps["Input3"] - timestamps["Input1"]
events = events.merge(timestamps, on = ["File", "TrialNum"])
events = events.drop(columns=["Input1", "Input2", "Input3"])

duration_quantiles = events["DurationLick"].quantile([0.02*i for i in range(51)])
duration_quantiles.to_csv(duration_csv_path)
events["DurationsGroup"] = np.where(~pd.isna(events["DurationLick"]), np.searchsorted(duration_quantiles.to_numpy(), events["DurationLick"].to_numpy()), None)


df["df/f"] = (df["CH1-470"]-df["CH1-410"])/df["CH1-410"]
df["df/f_zscoredfile"] = df.groupby("File")["df/f"].transform(lambda x: (x-x.mean())/x.std())



grouped_ev = events.groupby(group_cols)
all_results= {}

for group, gd in tqdm.tqdm(grouped_ev, desc="Event_groups"):
    results=[]
    ev_num=0
    for i, row in tqdm.tqdm(gd.iterrows(), desc="event", total=len(gd.index)):
        try:
            time=row["TimeStamp"]
            fdf = df.loc[df["File"] == row['File'], :].copy()
            fdf["relative_timestamp"] = fdf["TimeStamp"] - time
            selected = fdf.loc[(fdf["relative_timestamp"] >=-2000) &  (fdf["relative_timestamp"] <=4000), :].copy()
            selected["df/f z_scored"] = scipy.stats.zscore(selected["df/f"]) if zscore_method=="baseline" else selected["df/f_zscoredfile"] # Baseline zscore
            selected = selected.sort_values("relative_timestamp")
            result = np.interp(np.linspace(-2000, 4000, 500), selected["relative_timestamp"].to_numpy(), selected["df/f z_scored"].to_numpy())
            results.append(result)
        except:
            print(f"{group}, {row['File']}, {df.loc[df['File'] == row['File'], :]}")
            raise
        if plot_event:
            plt.plot(selected["relative_timestamp"], selected["df/f z_scored"])
            plt.plot(np.linspace(-2000, 4000, 500), result)
            plt.title(f"{group}, {row['File']}, ev_num={ev_num}")
            plt.show()
        ev_num+=1

    results = np.stack(results)
    results = pd.DataFrame(results.T, index = np.linspace(-2000, 4000, 500))
    results = results.mean(axis=1)
    all_results[group] = results

all_results = pd.concat(all_results)
print(all_results)
# pickle.dump(all_results, (base/"all_results.pkl").open("wb"))
all_results.index.names = group_cols +["t"]
all_results.name = "df/f z-scored"
all_results = all_results.to_frame()
print(all_results)
pickle.dump(all_results, (base/result_path).open("wb"))







# df["df/f"] = (df["CH1-410"] - df["CH1-470"])/df["CH1-410"]
# # events.columns = [f"EventNumber_{i}" for i in range(len(events.columns))]
# print(events)
# print(df)



# df[["Input", "_0", "End", "_1", "_2"]] = df["Events"].str.split("*", expand=True)
# # res["has_ev"] = res.apply(lambda x: ~pd.isna(x.to_numpy()).all(), axis=1)
# # print(res[res["has_ev"]])
# # print(res)
# df = df[~pd.isna(df["Events"])]

# print(df)
