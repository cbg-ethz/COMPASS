import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import seaborn as sns
import pandas as pd

df = pd.read_csv("/home/e840r/Documents/COMPASS_project/out_simulations/results_BITSC2sim.csv")
d={"SNVs":[],"CNAs":[],"seed":[],"rhosigma":[],"method":[],"metric":[],"similarity":[]}
for i in df.index:
    for metric in ["MP3","MP3-SNV","MP3-CNA","assignmentsAccuracy"]: # 
        d["SNVs"].append(df.loc[i,"nSNVs"])
        d["CNAs"].append(int(df.loc[i,"nCNAs"]))
        d["seed"].append(df.loc[i,"seed"])
        d["method"].append(df.loc[i,"method"])
        d["metric"].append(metric)
        d["similarity"].append(df.loc[i,metric])
        d["rhosigma"].append(df.loc[i,'rhosigma'])
df = pd.DataFrame(d)
sns.set(font="Helvetica",font_scale=2.0)
g = sns.catplot(x="rhosigma", y="similarity",hue="method",row="metric",\
    data=df, kind="box",hue_order=["COMPASS","BITSC2"],
        palette=["#E00105","#FDC508"],
    showmeans=True,meanprops={"markersize":10,"marker":"X"},
    margin_titles=False,legend=False,fliersize=2,height=3,aspect=1.6) # , height=4, aspect=.7
g.set_axis_labels("Variance in region coverage","MP3 similarity")
for AX in g.axes.flatten():
    box_patches = [patch for patch in AX.patches if type(patch) == matplotlib.patches.PathPatch]
    for i,patch in enumerate(box_patches):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        col = patch.get_facecolor()
        print(col)
        patch.set_edgecolor(col)
        patch.set_alpha(0.5)
        #artist.set_facecolor('None')
        # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        for j in range(i*7,i*7+7):
            line = AX.lines[j]
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)

for AX in g.axes[3]:
    t_labels = [x.get_text() for x in AX.get_xticklabels()]
    print(t_labels)
    def format_func(value, tick_number):
        if not value in range(len(t_labels)):
            return value
        label = float(t_labels[value])
        if label==0.5:
            return r"$\bf{"+ "{:.1f}".format(label)+r"}$"
        else:
            return "{:.1f}".format(label)
    AX.xaxis.set_major_formatter(plt.FuncFormatter(format_func))

handles, labels = g.axes[0][0].get_legend_handles_labels()

#labels, handles = zip(*sorted(zip(labels, handles), key=sort_legends))
g.axes[0][0].legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.50),
          ncol=3, fancybox=True, shadow=False)

#plt.ylim(0.33,1)

g.axes[0][0].set_ylabel("Full MP3")
g.axes[1][0].set_ylabel("SNV MP3")
g.axes[2][0].set_ylabel("CNA MP3")
g.axes[3][0].set_ylabel("Cell assignments")

#g.axes[2][0].set_xlabel("")
#g.axes[2][2].set_xlabel("")

for j in [0,1,2,3]:
    g.axes[j,0].set_title("")
#g.axes[0,0].set_title("7 SNVs")
#g.axes[0,1].set_title("20 SNVs")
#g.axes[0,2].set_title("50 SNVs")
g.axes[0,0].set_title("20 SNVs, 3 CNAs")
plt.subplots_adjust(hspace=0.1, wspace=0.06)
#plt.show()

plt.savefig("/home/e840r/Documents/COMPASS_project/Figures/simulations/BITSC2sim.svg",dpi=500,bbox_inches="tight",pad_inches=0.1)