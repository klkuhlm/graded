from glob import glob
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1,num=1,figsize=(6,4))

colors = {1.25:"red",1.5:"green",2.0:"black",3.0:"blue",4.0:"orange"}

Lc = 1.0
Tc = 0.05 * (0.05 * 5.0E-10 + 1.0E-6) * Lc**2 * 1.2E-3 / 2.0E-17
Pc = 1.5E+6

print(Lc,Tc,Pc)

for fn in glob("test-*.out"):
    d = []
    r = float(fn[:-4][5:])
    with open(fn,"r",encoding="ascii") as fh:
        for line in fh.readlines()[11:]:
            d.append([float(x) for x in line.strip().split()])
    d = np.array(d)

    fn = f"*-{r:.3g}r-drawdown.csv"
    fns = glob(fn)
    print(fn,fns)
    dd = np.loadtxt(fns[0]) # drawdown data (sec,MPa)
    
    ax.loglog(d[:,0]*Tc,d[:,1]*Pc * 1.0E-6,color=colors[r],
              label=f"r={r} mod",lw=0.5)
    ax.loglog(dd[:,0], dd[:,1] , color=colors[r],label=f"r={r} obs",marker=".")

ax.set_xlabel("$t$ [sec]")
ax.set_ylabel("$p$ [MPa]")
ax.set_ylim([1.0E-3 * Pc * 1.0E-6,2.5 * Pc * 1.0E-6])
ax.grid(True)
ax.legend(loc=0,fontsize="small")
ax.set_title(f"$T_c=${Tc:.3E} $P_c=${Pc:.3E}")

fig.savefig("powerlaw-prediction.png",dpi=200)
plt.close(1)
