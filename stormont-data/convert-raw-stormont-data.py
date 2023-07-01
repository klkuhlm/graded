from glob import glob
import numpy as np
import matplotlib.pyplot as plt

DAY2SEC = 24.0 * 60.0 * 60.0

# digitized "mineby" line
mb = []
with open("mineby.csv", "r", encoding="ascii") as fh:
    for line in fh.readlines():
        mb.append([float(x) for x in line.strip().split(",")])

mb = np.array(mb)
mineby_time = np.mean(mb[:, 0]) * DAY2SEC
print(f"one day in seconds {DAY2SEC:.4e}")
print(f"mineby time (sec): {mineby_time:.4e}")

d = {}
delta = []
bores = {}

fig, ax = plt.subplots(1, 1, num=1, figsize=(6, 4))

for fn in glob("??-*r.csv"):
    bore_no = fn.split("-")[0]
    bore_r = float(fn.split("-")[1][:-5])  # drop "r.csv" from end
    bores[bore_no] = bore_r
    d[bore_no] = []
    with open(fn, "r", encoding="ascii") as fh:
        for line in fh.readlines():
            d[bore_no].append([float(x) for x in line.strip().split(",")])
    d[bore_no] = np.array(d[bore_no])
    d[bore_no][:, 0] *= DAY2SEC
    print(bore_no, bore_r, d[bore_no].shape)
    print(d[bore_no][:5, :])

    delta = np.zeros_like(d[bore_no])
    delta[:, 0] = d[bore_no][:, 0] - mineby_time
    idx = np.argmin(np.abs(delta[:, 0]))
    delta[:, 1] = d[bore_no][idx, 1] - d[bore_no][:, 1]

    sidx = np.argsort(delta[:, 0])
    sdelta = delta[sidx, :]
    pos_mask = sdelta[:, 0] > 0.0

    ax.loglog(
        sdelta[pos_mask, 0],
        sdelta[pos_mask, 1],
        label=f"{bore_no} {bore_r}r",
        lw=0.25,
        marker=".",
    )

    np.savetxt(fn.replace(".csv", "-drawdown.csv"), sdelta[pos_mask], fmt="%.3e")

ax.set_xlabel("time since mineby [sec]")
ax.set_ylabel("$\\Delta P$ since mineby [MPa]")
ax.grid(True)
ax.legend(loc=0, fontsize="small")

fig.savefig("drawdown.png", dpi=200)
plt.close(1)
