from glob import glob
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1,num=1,figsize=(6,4))

colors = {1.25:"red",1.5:"green",2.0:"black",3.0:"blue",4.0:"orange"}

for fn in glob("test-*.out"):
    d = []
    r = float(fn[:-4][5:])
    with open(fn,"r",encoding="ascii") as fh:
        for line in fh.readlines()[11:]:
            d.append([float(x) for x in line.strip().split()])
    d = np.array(d)

    ax.loglog(d[:,0],d[:,1],color=colors[r],
              label=f"r={r}",lw=0.25,marker=".",markersize=1.5)

ax.set_xlabel("$t_D$")
ax.set_ylabel("$p_D$")
ax.set_ylim([1.0E-2,1.5])
ax.grid(True)
ax.legend(loc=0,fontsize="small")

fig.savefig("powerlaw-prediction.png",dpi=200)
plt.close(1)
