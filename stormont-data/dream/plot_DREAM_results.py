import scipy.io as sio
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import corner
import subprocess


def acf(x, length=30):
    return np.array([1] + [np.corrcoef(x[:-i], x[i:])[0, 1] for i in range(1, length)])

# read in MATLAB .mat file
base = f"powerlaw_stormont_uniform"
print(base)

READ_FINAL = False
burnin = 10000

if READ_FINAL:
    m = sio.loadmat(f"{base}_results.mat")
    niter, nvars, nchains = m["Sequences"].shape
    print(f"reading final results {niter} {nvars-2} {nchains}")
else:
    m = sio.loadmat(f"DREAM_ZS.mat")
    niter, nvars, nchains = m["Sequences"].shape
    # intermediate results padded with zeros at end, cut these off
    niter = niter - np.argmax(np.abs(m["Sequences"][::-1, 0, 0]) > 0)
    print(f"reading intermediate results {niter} {nvars-2} {nchains}")

npar = nvars - 2

norm = matplotlib.colors.Normalize(
    vmin=-(m["Sequences"][burnin:, npar, :].max()),
    vmax=-(m["Sequences"][burnin:, npar, :].min()),
)
gray = matplotlib.colormaps["gray"]

fn = f"{base}"

svn = ["k0", "n0", "eta", "tau", "cm", "m"]
varnames = [
    "log$_{10} k_0$",
    "log$_{10} n_0$",
    "$\\eta$",
    "$\\tau$",
    "log$_{10} c_m$",
    "$m$",
]
colors = ["red", "green", "blue", "cyan", "magenta", "black"]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot of all chains
if 1:

    INCLUDE_LIKELIHOOD = False
    nlag = 1000
    lags = 1 + np.arange(nlag)

    fig, ax = plt.subplots(nrows=npar, ncols=2, sharex="col", figsize=(6, 7), num=1)
    print("ax.shape", ax.shape)
    for i in range(npar):
        for j in range(nchains):
            ax[i, 0].plot(
                m["Sequences"][:niter, i, j],
                color=colors[j],
                linestyle="-",
                linewidth=0.5,
            )
            if INCLUDE_LIKELIHOOD:
                if j == 0 and i == 0:
                    axl = ax[i, 0].twinx()
                    axl.semilogy(
                        -m["Sequences"][:niter, npar, j],
                        color="black",
                        linestyle="-",
                        linewidth=0.75,
                    )

            ac = acf(m["Sequences"][burnin:niter, i, j], length=nlag)
            ax[i, 1].plot(lags, ac, color=colors[j], alpha=0.5)

        ax[i, 0].set_ylabel(varnames[i])
        ax[i, 0].set_xlim([0, niter])
        ax[i, 0].grid(True)
        ax[i, 0].axvline(burnin, color="black", lw=1, linestyle="--")

        ax[i, 1].set_xlim([0, nlag])
        ax[i, 1].set_ylim([-0.5, 1])
        ax[i, 1].grid(True)

    ax[-1, 0].set_xlabel("Iterations per chain")
    ax[-1, 1].set_xlabel("lag")
    fig.suptitle(fn)
    fig.tight_layout()
    fig.savefig(f"{fn}-all-chains.png", dpi=200)
    plt.close(1)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plot of posteriors as histograms
if 1:

    ncols = 2
    nrows = int(np.ceil(npar / ncols))

    def gaussian(x, m, v):
        return np.exp(-((x - m) ** 2) / (2.0 * v)) / np.sqrt(2.0 * np.pi * v)

    bins = 30
    fig, ax = plt.subplots(
        nrows=nrows, ncols=ncols, sharex=False, sharey=False, figsize=(6, 7), num=1
    )
    fh = open(f"{fn}-stats.txt", "w")
    fh.write("i,variable,mean,norm,variance,Gelman-Rubin\n")
    L = niter - burnin

    grand_mean = []

    for i in range(npar):

        chain_mean = np.mean(m["Sequences"][burnin:niter, i, :], axis=0)
        xx = m["Sequences"][burnin:niter, i, :].flatten()
        grand_mean.append(np.mean(xx))
        btw_chain_var = L / 2.0 * np.sum((chain_mean - grand_mean[-1]) ** 2)
        win_chain_var = np.var(m["Sequences"][burnin:niter, i, :], axis=0)
        W = np.mean(win_chain_var)
        gelman_rubin = (float(L - 1) / L * W + 1.0 / L * btw_chain_var) / W

        mn = grand_mean[-1]
        vr = np.var(xx)
        xout = np.linspace(xx.min(), xx.max())
        r = i // ncols
        c = i % ncols

        hn, hb, _ = ax[r, c].hist(
            xx, bins=bins, density=True, histtype="step", linewidth=2.0
        )
        ax[r, c].plot(xout, gaussian(xout, mn, vr), "r-")
        ax[r, c].set_xlabel(varnames[i])
        # ax[r,c].set_ylabel("probability")
        ax[r, c].grid(True)
        md = np.argmax(hn)  # tallest bin
        bin_centers = (hb[1:] + hb[:-1]) / 2.0
        nr = bin_centers[md]  # bin center (hist returns edges)
        fh.write(f"{i+1},{svn[i]},{mn:.7E},{nr:.7E},{vr:.7E},{gelman_rubin:.7E}\n")
        print(i + 1, svn[i], mn, nr, vr, gelman_rubin)

    if npar % 2 > 0:
        ax[nrows - 1, ncols - 1].set_axis_off()

    fh.close()
    fig.suptitle(fn)
    fig.tight_layout()
    fig.savefig(f"{fn}-posterior-histograms.png", dpi=200)
    plt.close(1)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# corner plot of joint posterior distributions
if 1:

    xx = np.concatenate(
        (
            m["Sequences"][burnin:niter, 0:npar, 0],
            m["Sequences"][burnin:niter, 0:npar, 1],
            m["Sequences"][burnin:niter, 0:npar, 2],
        ),
        axis=0,
    )

    fig = plt.figure(1, figsize=(8, 8))
    fig = corner.corner(
        data=xx, labels=varnames, quantiles=[0.16, 0.5, 0.84], nbins=bins, fig=fig
    )

    fig.suptitle(fn)
    fig.savefig(f"{fn}-corner.png", dpi=300)
    plt.close(1)

# ****************************************
# horsetail plot
if 1:
    Lc = 0.4825  # borehole radius (m)
    mu = 1.75e-03  # brine viscosity (Pa*sec)
    cf = 3.1e-10  # fluid compressibility (1/Pa)

    colors = ["red","green","blue","black"]

    nhorsetail = 50

    # data tables are time (seconds) and change in pressure (MPa)
    data53 = np.loadtxt('53-1.25r-drawdown.txt')
    data54 = np.loadtxt('54-1.5r-drawdown.txt')
    data56 = np.loadtxt('56-2r-drawdown.txt')
    data58 = np.loadtxt('58-3r-drawdown.txt')
    
    obs = [data53, data54, data56, data58] # don't concatenate
    rvec = [1.25, 1.5, 2.0, 3.0] # rD for each observation

    for ob in obs:
        print(ob.shape)
        print(ob[:5,:])
    
    # rough non-dimensionalizing for ballparking time range from data
    c = 10.0 ** grand_mean[1] * cf + 10.0 ** grand_mean[4]  # (n*cf+cm)
    # n0*c*Lc^2*mu/k0
    Tc = (10.0 ** grand_mean[1] * c * Lc ** 2 * mu / 10.0 ** grand_mean[0])  
    tD = np.logspace(
        np.log10(0.75 * obs[0][0, 0] / Tc),
        np.log10(10.0 * obs[0][-1, 0] / Tc), 100
    )
    np.savetxt("tD.in", tD)
    print("tD.shape",tD.shape)
    print(f"c={c:.3e} Tc={Tc:.3e} tD0={tD[0]} tDN={tD[-1]}")

    # random sample from chain
    np.random.seed(1234)
    sidx = np.random.randint(low=burnin - 1, high=niter - 1, size=nhorsetail)

    norm = matplotlib.colors.Normalize(
        vmin=-(m["Sequences"][burnin:, npar, :].max()),
        vmax=-(m["Sequences"][burnin:, npar, :].min()),
    )

    fig, axes = plt.subplots(1, 2, num=1, figsize=(6, 4))
    lw = 0.25

    for i in range(nhorsetail):
        if i % 10 == 0:
            print(i)
        for j in range(nchains):
            
            params = m["Sequences"][sidx[i], :, j]
            p = [10.0**params[0], 10.0**params[1], params[2],
                 params[3], 10.0**params[4], params[5]]
            pout = [params[2], params[2]*params[3], params[5], -999.0] # eta, kappa, m, rD

            #print("params",params)
            #print("p",p)
            #print("pout",pout)
            
            #cg = gray(norm(-params[npar]))
            c = p[1] * cf + p[4]  # (n*cf + cm)
            Tc = p[1] * c * Lc ** 2 * mu / p[0]  # n0*c*Lc^2*mu/k0
            P_scale = 3.0e+6 

            print(f"c={c:.3e} Tc={Tc:.3e} Pc={P_scale:.3e}")
            
            tt = tD * Tc
            
            for ii,rv in enumerate(rvec):
                pout[3] = rv
                
                fh = open("parameters.in", "w")
                fh.write(" ".join([f"{x:.12E}" for x in pout]) + "\n")
                fh.close()
                
                subprocess.run("./drive-powerlaw.sh", shell=True)

                sim = np.loadtxt("powerlaw.out", usecols=(1,))
                #print(sim)
                
                SIM = sim * P_scale

                # re-dimensionalize model output
                axes[0].semilogx(tt, SIM, linestyle="-", color=colors[ii], lw=lw)
                axes[1].loglog(tt, SIM, linestyle="-", color=colors[ii], lw=lw)

    for ii in range(4):
        axes[0].semilogx(
            obs[ii][:, 0],
            obs[ii][:, 1]*1.0E+6,
            color=colors[ii],
            marker=".",
            markerfacecolor=colors[ii],
            linestyle="none",
            markersize=4,
        )
        axes[1].loglog(
            obs[ii][:, 0],
            obs[ii][:, 1]*1.0E+6,
            color=colors[ii],
            marker=".",
            markerfacecolor=colors[ii],
            linestyle="none",
            markersize=4,
        )

    axes[0].set_xlabel("time (s)")
    axes[0].set_ylabel("p (Pa)")
    axes[1].set_xlabel("time (s)")
    axes[1].set_ylabel("p (Pa)")

    axes[0].grid(True)
    axes[1].grid(True)

    ##fig.suptitle(bore)
    fig.tight_layout()
    fig.savefig(f"{fn}-horsetail.png", dpi=200)
    plt.close(1)
