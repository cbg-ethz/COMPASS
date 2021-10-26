import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_csv("COMPASS/Experiments/workflow_BITSC2/results_BITSC2_CNVleaves.csv")
sns.set(font="Helvetica",font_scale=2.0)
g = sns.catplot(x="rho", y="BD",hue="method", col="propCNVleaves",row="nodes", data=df, kind="box",palette="Set1",showmeans=True,meanprops={"markersize":10,"marker":"X"},
        margin_titles=False,legend=False,fliersize=2,height=3,aspect=1.7) # , height=4, aspect=.7
g.set_axis_labels("                                                Variance in coverage between regions","Bourque Distance to the true tree")
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

for AX in g.axes[2]:
    t_labels = [x.get_text() for x in AX.get_xticklabels()]
    def format_func(value, tick_number):
        if not value in range(len(t_labels)):
            return value
        label = float(t_labels[value])
        if label==0.5:
            return r"$\bf{"+ "{:.1f}".format(label)+r"}$"
        else:
            return "{:.1f}".format(label)
    AX.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
g.axes[0][2].legend(loc='upper center', bbox_to_anchor=(-0.05, 1.60),
          ncol=2, fancybox=True, shadow=False)

g.axes[0][0].set_ylabel("")
g.axes[2][0].set_ylabel("")

g.axes[2][0].set_xlabel("")
g.axes[2][2].set_xlabel("")
g.axes[2][3].set_xlabel("")

for x in g.axes:
    for y in x:
        y.title.set_fontsize(18)
plt.ylim(0,12)
plt.show()




plt.savefig("Comparison_BITSC2_CNVleaves.png",dpi=500,bbox_inches="tight",pad_inches=0.1)