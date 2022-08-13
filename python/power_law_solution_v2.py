import mpmath as mp

method = "dehoog"
mp.mp.dps = 16  # decimal digits of precision

FRACTURED = False

# sigma=0 is plain Type-II solution, non-zero is wellbore storage
sigma = mp.mpmathify("1.0E-3")  

lamda = mp.mpmathify("1.0E-6")  # inter-porosity exchange coeff
omega = mp.mpmathify("1.0E-3")  # fracture porosity ratio
invomega = 1 - omega

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

print(
    f"FRACTURED:{FRACTURED}, dps:{mp.mp.dps}, method:{method}"
)

def sol_p(p):
    # prediction of gradient of pressure at wellbore
    # for specified pressure source
    # alpha = (kappa - m + 1) / 2
    # gamma = (kappa - eta + 2) / 2
    if FRACTURED:
        beta = mp.sqrt(
            (omega * p + invomega * p * lamda / (invomega * p + lamda)) / gammasq
        )
    else:
        beta = mp.sqrt(p / gammasq)
    # nu = alpha / gamma
    return mp.sqrt(p) * mp.besselk(nu - 1, beta) / (mp.besselk(nu, beta) * p)

def sol_q(p):
    # prediction of pressure at wellbore
    # for specified flowrate source
    # alpha = (kappa - m + 1) / 2
    # gamma = (kappa - eta + 2) / 2
    if FRACTURED:
        beta = mp.sqrt(
            (omega * p + invomega * p * lamda / (invomega * p + lamda)) / gammasq
        )
    else:
        beta = mp.sqrt(p / gammasq)
    # nu = alpha / gamma
    bkn = mp.besselk(nu, beta)
    return bkn / (
        mp.power(p, 1.5) * mp.besselk(nu - 1, beta) - sigma * p ** 2 * bkn
    )

def sol_qr(p):
    # prediction of pressure at general r
    # for specified flowrate source
    # alpha = (kappa - m + 1) / 2
    # gamma = (kappa - eta + 2) / 2
    if FRACTURED:
        beta = mp.sqrt(
            (omega * p + invomega * p * lamda / (invomega * p + lamda)) / gammasq
        )
    else:
        beta = mp.sqrt(p / gammasq)
    # nu = alpha / gamma
    return (
        rD ** alpha
        * mp.besselk(nu, beta * rD ** gamma)
        / (
            mp.power(p, 1.5) * mp.besselk(nu - 1, beta)
            + sigma * p ** 2 * mp.besselk(nu, beta)
        )
    )

def sol_pr(p):
    # prediction of pressure at general r
    # for specified pressure source
    # alpha = (kappa - m + 1) / 2
    # gamma = (kappa - eta + 2) / 2
    if FRACTURED:
        beta = mp.sqrt(
            (omega * p + invomega * p * lamda / (invomega * p + lamda)) / gammasq
        )
    else:
        beta = mp.sqrt(p / gammasq)
    # nu = alpha / gamma
    return (
        rD ** (alpha + gamma - 1)
        * mp.besselk(nu - 1, beta * rD ** gamma)
        / (mp.besselk(nu, beta) * mp.sqrt(p))
    )

tDvec = [10.0 ** x for x in mp.linspace(-4.0, 8.0, 20)]

colors = ["red", "green", "blue", "cyan", "magenta", "black", "orange"]
line_styles = ["-", ":", "--", "-."]

color_lines = []
line_lines = []

tau_vec = [mp.mpmathify(x) for x in ["1", "1.5", "2", "2.5", "3"]]
eta_vec = [mp.mpmathify(x) for x in ["-0.5", "0", "0.5", "1", "2"]]
m_vec = [mp.mpmathify(x) for x in ["0", "1", "2", "3"]]

fhp = open(f"power_law_tD_typeI_python_{FRACTURED}.txt", "w")
fhq = open(f"power_law_tD_typeII_python_{FRACTURED}.txt", "w")

header = "kappa_eta_m " + " ".join([f"{float(x):.8E}" for x in tDvec]) + "\n"
fhp.write(header)
fhq.write(header)

for tau in tau_vec:
    print(f"tau:{tau}")

    for eta in eta_vec:
        figL = plt.figure(1, figsize=(5, 5))
        figR = plt.figure(2, figsize=(5, 5))
        axp = figL.add_subplot(111)
        axq = figR.add_subplot(111)

        for ii, (color, m) in enumerate(zip(colors, m_vec)):

            kappa = eta * tau  # physically motivated, ensures kappa >= eta
            alpha = (1 + kappa - m) / 2
            gamma = (2 + kappa - eta) / 2
            gammasq = gamma * gamma
            nu = alpha / gamma  # this simplification requires gamma > 0

            key = f"\u03c4{float(tau):.2g}_\u03BA{float(kappa):.2g}_\u03B7{float(eta):.2g}_m{float(m):.2g} "

            print(
                f"{key} \u03B1:{alpha}, \u03B3:{gamma}, \u03BD:{float(nu):.3g}"
            )

            pvec = []
            qvec = []
            for tD in tDvec:
                pvec.append(mp.invertlaplace(sol_p, tD, method=method))
                qvec.append(mp.invertlaplace(sol_q, tD, method=method))

            label = f"$m={float(m):.2g}$ $\\eta={float(eta):.2g}$"

            fhq.write(key + " ".join([f"{float(x):.8E}" for x in qvec]) + "\n")
            fhp.write(key + " ".join([f"{float(x):.8E}" for x in pvec]) + "\n")

            axp.loglog(tDvec, pvec, ls='-', color=color, label=label)
            axq.loglog(tDvec, qvec, ls='-', color=color, label=label)

        axp.set_xlabel("$t_D$")
        axq.set_xlabel("$t_D$")
        axp.set_ylabel("$\\partial p/\\partial x$")
        # axp.set_title("specified pressure")
        axp.legend(
            loc="lower left",
            fontsize="xx-small",
            ncol=2,
            handlelength=3,
        )
        axp.set_ylim(ymin=1.0e-6)
        axp.grid()

        axq.set_ylabel("$p$")
        # axq.set_title("specified flowrate")
        axq.set_ylim(ymin=0.01)
        axq.legend(
            loc="lower right",
            fontsize="xx-small",
            ncol=2,
            handlelength=3,
        )
        axq.grid()

        figL.tight_layout()
        figR.tight_layout()
        figL.savefig(
            f"power_law_p_tD_solution_tau{float(tau):.2g}_eta{float(eta):.2g}_{FRACTURED}.png", dpi=200
        )
        figR.savefig(
            f"power_law_q_tD_solution_tau{float(tau):.2g}_eta{float(eta):.2g}_{FRACTURED}.png", dpi=200
        )
        plt.close(1)
        plt.close(2)

fhp.close()
fhq.close()
print("")
