import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import seaborn as sns
import pandas as pd

d={"method":[],"n_cells":[],"time":[]}
indir="/home/e840r/Documents/COMPASS_project/out_simulations/ncells/benchmarks"
for f in sorted(os.listdir(indir)):
    with open(os.path.join(indir,f),"r") as infile:
        tmp = infile.readline()
        t = float(infile.readline().split("\t")[0])
        d["time"].append(t)
        ncells = int(f.split("_")[0].rstrip("ncells"))
        d["n_cells"].append(ncells)
        method= f.split("_")[-1][:-4]
        d["method"].append(method)

df = pd.DataFrame(d)

sns.set(font="Helvetica",font_scale=1.6)
g = sns.catplot(x="n_cells", y="time",hue="method",dodge=False,\
    data=df, kind="box",hue_order=["COMPASS","BITSC2","infSCITE"],
        palette=["#E00105","#FDC508","#215CAF"],
    showmeans=True,meanprops={"markersize":10,"marker":"X"},
    margin_titles=False,legend=False,fliersize=2,height=5,aspect=1.05) # , height=4, aspect=.7
g.set_axis_labels("Number of cells","time (s)")
g.axes[0][0].set_yscale("log")
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
#g.axes[0][0].legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.35),
#          ncol=5, fancybox=True, shadow=False)
g.axes[0][0].set_title("6 nodes, 20 SNVs")
#plt.ylim(0.33,1)


#plt.show()




plt.savefig("/home/e840r/Documents/COMPASS_project/Figures/simulations/runtimes.svg",dpi=500,bbox_inches="tight",pad_inches=0.1)