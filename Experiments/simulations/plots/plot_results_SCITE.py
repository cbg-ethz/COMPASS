import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import seaborn as sns
import pandas as pd

colors = ["#A51C30","#EF0107","#FA8072","#FFAD2E","#FFEC00","#F9FF20","#0039a6","#4169E1","#87CEFA"]

colors = ["#96272D","#B7352D","#D48681","#FFAD2E","#FFEC00","#F9FF20","#08407E","#215CAF","#7A9DCF"] #ETH

df = pd.read_csv("/home/e840r/Documents/COMPASS_project/out_simulations/results_SCITEsim.csv")
sns.set(font="Helvetica",font_scale=1.8)
g = sns.catplot(x="nodeprob", y="MP3",hue="method", col="dropoutrate", data=df, kind="box",palette=["#EF0107","#4169E1"],showmeans=True,meanprops={"markersize":10,"marker":"X"},
        margin_titles=False,legend=False,fliersize=2,height=3,aspect=1.2)
g.set_axis_labels("Dirichlet concentration for the node probabilities","MP3 Similarity to the true tree")
g.set(ylim=(0.25,1.03))
sns.set_style("whitegrid")
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
    def format_func(value, tick_number):
        if not value in range(len(t_labels)):
            return value
        label = float(t_labels[value])
        if label==0.3:
            return r"$\bf{"+ "{:.1f}".format(label)+r"}$"
        else:
            return "{:.1f}".format(label)
    AX.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
g.axes[0][1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.80),
          ncol=2, fancybox=True, shadow=False)
plt.ylim(0,1)

g.axes[0][0].set_ylabel("MP3")

g.axes[0][0].set_xlabel("")
g.axes[0][2].set_xlabel("")

g.axes[0][0].set_title("dropout rate = 0.04")
g.axes[0][1].set_title("dropout rate = 0.07")
g.axes[0][2].set_title("dropout rate = 0.1")

#plt.show()

plt.savefig("/home/e840r/Documents/COMPASS_project/Figures/simulations/SCITE_simulation.svg",dpi=500,bbox_inches="tight",pad_inches=0.1)