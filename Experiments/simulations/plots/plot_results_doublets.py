import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import seaborn as sns
import pandas as pd

df = pd.read_csv("/home/e840r/Documents/COMPASS_project/out_simulations/doublets/results.csv")
d={"SNVs":[],"CNAs":[],"seed":[],"method":[],"metric":[],"similarity":[],"doublet rate":[]}
for i in df.index:
    for metric in ["MP3"]:
        if metric!="MP3-CNA" or int(df.loc[i,"nCNAs"])>0:
            d["SNVs"].append(df.loc[i,"nSNVs"])
            d["CNAs"].append(int(df.loc[i,"nCNAs"]))
            d["seed"].append(df.loc[i,"seed"])
            d["method"].append(df.loc[i,"method"])
            d["metric"].append(metric)
            d["similarity"].append(df.loc[i,metric])
            d["doublet rate"].append(df.loc[i,"doubletrate"])
df = pd.DataFrame(d)
sns.set(font="Helvetica",font_scale=2.0)
g = sns.catplot(x="doublet rate", y="similarity",hue="method",\
    data=df, kind="box",hue_order=["COMPASS","COMPASS-doublet","BITSC2","infSCITE","infSCITE-doublet"],
        palette=["#E00105","#96272D","#FDC508","#215CAF","#08407E"],
    showmeans=True,meanprops={"markersize":10,"marker":"X"},
    margin_titles=False,legend=False,fliersize=2,height=7,aspect=1.5) # , height=4, aspect=.7
g.set_axis_labels("Doublet rate","MP3 similarity to the true tree")
for AX in g.axes.flatten():
    box_patches = [patch for patch in AX.patches if type(patch) == matplotlib.patches.PathPatch]
    for i,patch in enumerate(box_patches):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        col = patch.get_facecolor()
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

for AX in g.axes[0]:
    t_labels = [x.get_text() for x in AX.get_xticklabels()]
    print(t_labels)
    def format_func(value,tick_number):
        return str(int(value))
    #AX.xaxis.set_major_formatter(plt.FuncFormatter(format_func))

handles, labels = g.axes[0][0].get_legend_handles_labels()
"""def sort_legends(key):
    name=key[0]
    print(name)
    if name=="COMPASS_0CNVs":
        return 0
    elif name=="COMPASS_3CNVs":
        return 2
    elif name=="COMPASS_8CNVs":
        return 4
    elif name=="BITSC2_0CNVs":
        return 1
    elif name=="BITSC2_3CNVs":
        return 3
    elif name=="BITSC2_8CNVs":
        return 5
    else:
        raise Exception("Invalid legend name")"""
#labels, handles = zip(*sorted(zip(labels, handles), key=sort_legends))
g.axes[0][0].legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.18),
          ncol=5, fancybox=True, shadow=False)

plt.ylim(0.0,1.03)

g.axes[0][0].set_ylabel("MP3")




#plt.show()

plt.savefig("/home/e840r/Documents/COMPASS_project/Figures/simulations/doublets_res.pdf",dpi=500,bbox_inches="tight",pad_inches=0.1)