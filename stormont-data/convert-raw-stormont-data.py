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

fig2, axes = plt.subplots(1, 2, num=2, figsize=(8, 4), constrained_layout=True)

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
    delta[:, 0] = d[bore_no][:, 0] - mineby_time  # delta time

    if 0:
        idx = np.argmin(np.abs(delta[:, 0]))  # value closest to zero time (pos or neg)
        norm_t = d[bore_no][idx, 0]
        norm_p = d[bore_no][idx, 1]
    if 0:
        idx = np.argmax(d[bore_no][:, 1])
        norm_t = d[bore_no][idx, 0]
        norm_p = d[bore_no][idx, 1]
    if 1:
        avg_mask = np.logical_and(
            d[bore_no][:, 0] >= 140 * DAY2SEC, d[bore_no][:, 0] <= mineby_time + 14 * DAY2SEC
        )
        norm_t = np.median(d[bore_no][avg_mask, 0])
        #norm_p = np.median(d[bore_no][avg_mask, 1])
        norm_p = np.max(d[bore_no][avg_mask, 1])

    print(norm_p, norm_t / DAY2SEC)
    delta[:, 1] = (norm_p - d[bore_no][:, 1]) / norm_p  # delta relative pressure

    sidx = np.argsort(delta[:, 0])
    sdelta = delta[sidx, :]
    pos_mask = sdelta[:, 0] > 0.0

    # original data
    axes[0].plot(
        d[bore_no][sidx, 0] / DAY2SEC, d[bore_no][sidx, 1], lw=0.25, marker="."
    )
    axes[0].plot(mb[:, 0], mb[:, 1], "k--")

    # drop in pressure (not normalized)
    axes[1].plot(d[bore_no][sidx, 0] / DAY2SEC, delta[sidx, 1] * norm_p)
    axes[1].plot(mb[:, 0], mb[:, 1], "k--")

    ax.loglog(
        sdelta[pos_mask, 0] / DAY2SEC,
        sdelta[pos_mask, 1],
        label=f"{bore_no} {bore_r}r",
        lw=0.25,
        marker=".",
    )

    np.savetxt(fn.replace(".csv", "-drawdown.csv"), sdelta[pos_mask], fmt="%.3e")

axes[0].set_xlabel("days")
axes[0].set_ylabel("pressure [MPa]")
axes[0].grid(True)
axes[0].set_ylim([mb[1, 1], mb[0, 1]])
axes[1].set_xlabel("days")
axes[1].set_ylabel("$\\Delta$ pressure [MPa]")
axes[1].grid(True)
axes[1].set_ylim([-0.5, 3.5])

ax.set_xlabel("time since mineby [day]")
ax.set_ylabel("relative $\\Delta P$ since mineby [MPa]")
ax.set_title(f"assuming far-field pressure is {norm_p:.1f} MPa")
ax.grid(True)
ax.legend(loc=0, fontsize="small")

fig2.savefig("pressure.png", dpi=200)
plt.close(2)

fig.savefig("scaled_drawdown.png", dpi=200)
plt.close(1)
