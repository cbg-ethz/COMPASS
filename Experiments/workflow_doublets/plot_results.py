import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_csv("Experiments/workflow_doublets/results.csv")
#df=df[df["CNLOHs"]==0].reset_index()
#methodCNVs=[]
#for i in range(df.shape[0]):
#    methodCNVs.append(df.loc[i,"method"]+"_"+str(df.loc[i,"CNVs"]) + "CNVs")
#df["methodsCNV"] = methodCNVs
sns.set(font="Helvetica",font_scale=1.0)
# hue_order=["COMPASS_0CNVs","COMPASS_1CNVs","COMPASS_3CNVs","BITSC2_0CNVs","BITSC2_1CNVs","BITSC2_3CNVs"],
#  palette=["#ff9998","#f37676","#e72129","#aec6e6","#4c8fc1","#2276b1"],\
g = sns.catplot(x="doublet_rate", y="MP3",hue="method", col="nodes",row="SNVs",\
    data=df, kind="box", hue_order = ["COMPASS","SCITE","BITSC2","COMPASS-nodoublet","SCITE-nodoublet"],palette=["#932938","#2276b1","#eebd4f","#ee3d6a","#aec6e6"],
    showmeans=True,meanprops={"markersize":10,"marker":"X"},
    margin_titles=False,legend=False,fliersize=2,height=3,aspect=1.7) # , height=4, aspect=.7
g.set_axis_labels("Doublet rate","MP3 similarity to the true tree")
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

handles, labels = g.axes[0][0].get_legend_handles_labels()
g.axes[0][0].legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.55),
          ncol=5, fancybox=True, shadow=False)

#plt.ylim(0.33,1)

#plt.show()

plt.savefig("Evaluation_doublets.png",dpi=500,bbox_inches="tight",pad_inches=0.1)