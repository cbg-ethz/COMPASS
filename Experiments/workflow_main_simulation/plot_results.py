import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_csv("Experiments/workflow_main_simulation/results.csv")
#df=df[df["CNLOHs"]==0].reset_index()
methodCNVs=[]
for i in range(df.shape[0]):
    methodCNVs.append(df.loc[i,"method"]+"_"+str(df.loc[i,"CNVs"]) + "CNVs")
df["methodsCNV"] = methodCNVs
sns.set(font="Helvetica",font_scale=2.0)
# hue_order=["COMPASS_0CNVs","COMPASS_1CNVs","COMPASS_3CNVs","BITSC2_0CNVs","BITSC2_1CNVs","BITSC2_3CNVs"],
#  palette=["#ff9998","#f37676","#e72129","#aec6e6","#4c8fc1","#2276b1"],\
g = sns.catplot(x="rho", y="MP3",hue="methodsCNV", col="nodes",row="SNVs",\
    data=df, kind="box",hue_order=["COMPASS_0CNVs","COMPASS_3CNVs","COMPASS_8CNVs","BITSC2_0CNVs","BITSC2_3CNVs","BITSC2_8CNVs","SCITE_0CNVs","SCITE_3CNVs","SCITE_8CNVs"],
        palette=["#ee3d6a","#c12949","#932938","#eebd4f","#c1a236","#916e2f","#aec6e6","#4c8fc1","#2276b1"],
    showmeans=True,meanprops={"markersize":10,"marker":"X"},
    margin_titles=False,legend=False,fliersize=2,height=3,aspect=1.7) # , height=4, aspect=.7
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

for AX in g.axes[2]:
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
def sort_legends(key):
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
        raise Exception("Invalid legend name")
#labels, handles = zip(*sorted(zip(labels, handles), key=sort_legends))
g.axes[0][1].legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 2.15),
          ncol=3, fancybox=True, shadow=False)

#plt.ylim(0.33,1)

g.axes[0][0].set_ylabel("")
g.axes[2][0].set_ylabel("")

g.axes[2][0].set_xlabel("")
g.axes[2][2].set_xlabel("")
plt.show()

plt.savefig("Evaluation.png",dpi=500,bbox_inches="tight",pad_inches=0.1)