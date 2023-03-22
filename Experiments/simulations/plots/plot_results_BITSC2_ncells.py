import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd

l_ncells=[100,200,500,1000]
MP3_results=[]

df = pd.read_csv("/home/e840r/Documents/COMPASS_project/out_simulations/ncells/results_BITSC2.csv",sep=",")

sns.set(font="Helvetica",font_scale=1.8)

g = sns.catplot(x="ncells", y="MP3",data=df, kind="box",hue="method",palette=["#FDC508"], #palette="Set1"
    showmeans=True,meanprops={"markersize":10,"marker":"X","markerfacecolor":"black", "markeredgecolor":"black"},
    margin_titles=False,legend=False,fliersize=2,height=5,aspect=1.2) # , height=4, aspect=.7
g.set_axis_labels("Number of cells","MP3 similarity")

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
plt.ylim(0.5,1.1)
handles, labels = g.axes[0][0].get_legend_handles_labels()
#labels, handles = zip(*sorted(zip(labels, handles), key=sort_legends))
g.axes[0][0].legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.45),
          ncol=5, fancybox=True, shadow=False)
plt.savefig("/home/e840r/Documents/COMPASS_project/Figures/simulations/BITSC2_ncells.svg",bbox_inches="tight")
plt.show()






df = pd.read_csv("Experiments/workflow_BITSC2/results_BITSC2_merged.csv")
df = df[df["nodes"]==6]
df = df[df["SNVs"]==20]
df = df[df["CNVs"]==3]
sns.set(font="Helvetica",font_scale=1.4)
g = sns.catplot(x="rho", y="MP3",hue="method", col="nodes",row="SNVs",\
    data=df, kind="box",\
    palette="Set1",\
    showmeans=True,meanprops={"markersize":10,"marker":"X"},
    margin_titles=False,legend=False,fliersize=2) # , height=4, aspect=.7
g.set_axis_labels("Variance in coverage between regions","MP3 similarity to the true tree")
for AX in g.axes.flatten():
    for i,artist in enumerate(AX.artists):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        col = artist.get_facecolor()
        print(col)
        artist.set_edgecolor(col)
        artist.set_alpha(0.5)
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
        if label==0.5:
            return r"$\bf{"+ "{:.1f}".format(label)+r"}$"
        else:
            return "{:.1f}".format(label)
    AX.xaxis.set_major_formatter(plt.FuncFormatter(format_func))


g.axes[0][0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.45),
          ncol=3, fancybox=True, shadow=False)

#plt.ylim(0.33,1)
#plt.show()
plt.setp(g.axes[0][0].get_legend().get_texts(), fontsize='12')

plt.savefig("Comparison_BITSC2_500cells.png",dpi=500,bbox_inches="tight",pad_inches=0.1)
