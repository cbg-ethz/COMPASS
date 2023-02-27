import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import seaborn as sns
import pandas as pd

d={"method":[],"time":[]}
indir="/home/e840r/Documents/COMPASS_project/out_simulations/doublets/benchmarks"
for f in sorted(os.listdir(indir)):
    with open(os.path.join(indir,f),"r") as infile:
        tmp = infile.readline()
        t = float(infile.readline().split("\t")[0])
        d["time"].append(t)
        method= f.split("_")[-1][:-4]
        d["method"].append(method)

df = pd.DataFrame(d)

sns.set(font="Helvetica",font_scale=1.8)
g = sns.catplot( x="method",y="time",hue="method",dodge=False,\
    data=df, kind="box",hue_order=["COMPASS","COMPASS-doublet","BITSC2","infSCITE","infSCITE-doublet"],order=["COMPASS","COMPASS-doublet","BITSC2","infSCITE","infSCITE-doublet"],
        palette=["#E00105","#96272D","#FDC508","#215CAF","#08407E"],
    showmeans=True,meanprops={"markersize":10,"marker":"X"},
    margin_titles=False,legend=False,fliersize=2,height=6,aspect=1.0) # , height=4, aspect=.7
g.set_axis_labels("","time (s)")
g.axes[0][0].set_yscale("log")
g.axes[0][0].set(xticklabels=[])
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

for AX in g.axes[0]:
    t_labels = [x.get_text() for x in AX.get_xticklabels()]
    print(t_labels)
    def format_func(value,tick_number):
        return str(int(value))
    #AX.xaxis.set_major_formatter(plt.FuncFormatter(format_func))

handles, labels = g.axes[0][0].get_legend_handles_labels()
#labels, handles = zip(*sorted(zip(labels, handles), key=sort_legends))
g.axes[0][0].legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.35),
          ncol=5, fancybox=True, shadow=False)
g.axes[0][0].set_title("6 nodes, 7 SNVs, 0 CNA")
#plt.ylim(0.33,1)


#plt.show()




plt.savefig("/home/e840r/Documents/COMPASS_project/Figures/simulations/runtimes_doublets.svg",dpi=500,bbox_inches="tight",pad_inches=0.1)