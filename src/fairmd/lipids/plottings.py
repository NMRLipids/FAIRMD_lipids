"""
@DRAFT
Network communication. Downloading files. Checking links etc.
"""

import json
import os

import matplotlib.pyplot as plt
import numpy as np

from fairmd.lipids import FMDL_EXP_PATH, FMDL_SIMU_PATH


def plotFormFactor(  # noqa: N802
    exp_form_factor,
    k,
    legend,
    plot_color,
):
    """:meta private:"""
    x_vals = []
    y_vals = []
    for i in exp_form_factor:
        x_vals.append(i[0])
        y_vals.append(k * i[1])
    plt.plot(x_vals, y_vals, label=legend, color=plot_color, linewidth=4.0)
    plt.xlabel(r"$q_{z} [Å^{-1}]$", size=20)
    plt.ylabel(r"$|F(q_{z})|$", size=20)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlim([0, 0.69])
    plt.ylim([-10, 250])
    plt.legend(loc="upper right")
    plt.savefig("FormFactor.pdf")


def plotOrderParameters(OPsim, OPexp):  # noqa
    """:meta private:"""
    xValuesHG = []  # noqa: N806
    xValuesSN1 = []  # noqa: N806
    xValuesSN2 = []  # noqa: N806

    yValuesHGsim = []  # noqa: N806
    yValuesSN1sim = []  # noqa: N806
    yValuesSN2sim = []  # noqa: N806
    yValuesHGsimERR = []  # noqa: N806
    yValuesSN1simERR = []  # noqa: N806
    yValuesSN2simERR = []  # noqa: N806
    yValuesHGexp = []  # noqa: N806
    yValuesSN1exp = []  # noqa: N806
    yValuesSN2exp = []  # noqa: N806
    xValuesHGexp = []  # noqa: N806
    xValuesSN1exp = []  # noqa: N806
    xValuesSN2exp = []  # noqa: N806

    sn1carbons = {
        "M_G1C3_M M_G1C3H1_M": 2,
        "M_G1C3_M M_G1C3H2_M": 2,
        "M_G1C4_M M_G1C4H1_M": 3,
        "M_G1C4_M M_G1C4H2_M": 3,
        "M_G1C5_M M_G1C5H1_M": 4,
        "M_G1C5_M M_G1C5H2_M": 4,
        "M_G1C6_M M_G1C6H1_M": 5,
        "M_G1C6_M M_G1C6H2_M": 5,
        "M_G1C7_M M_G1C7H1_M": 6,
        "M_G1C7_M M_G1C7H2_M": 6,
        "M_G1C8_M M_G1C8H1_M": 7,
        "M_G1C8_M M_G1C8H2_M": 7,
        "M_G1C9_M M_G1C9H1_M": 8,
        "M_G1C9_M M_G1C9H2_M": 8,
        "M_G1C10_M M_G1C10H1_M": 9,
        "M_G1C10_M M_G1C10H2_M": 9,
        "M_G1C11_M M_G1C11H1_M": 10,
        "M_G1C11_M M_G1C11H2_M": 10,
        "M_G1C12_M M_G1C12H1_M": 11,
        "M_G1C12_M M_G1C12H2_M": 11,
        "M_G1C13_M M_G1C13H1_M": 12,
        "M_G1C13_M M_G1C13H2_M": 12,
        "M_G1C14_M M_G1C14H1_M": 13,
        "M_G1C14_M M_G1C14H2_M": 13,
        "M_G1C15_M M_G1C15H1_M": 14,
        "M_G1C15_M M_G1C15H2_M": 14,
        "M_G1C16_M M_G1C16H1_M": 15,
        "M_G1C16_M M_G1C16H2_M": 15,
        "M_G1C17_M M_G1C17H1_M": 16,
        "M_G1C17_M M_G1C17H2_M": 16,
        "M_G1C17_M M_G1C17H3_M": 16,
    }

    sn2carbons = {
        "M_G2C3_M M_G2C3H1_M": 2,
        "M_G2C3_M M_G2C3H2_M": 2,
        "M_G2C4_M M_G2C4H1_M": 3,
        "M_G2C4_M M_G2C4H2_M": 3,
        "M_G2C5_M M_G2C5H1_M": 4,
        "M_G2C5_M M_G2C5H2_M": 4,
        "M_G2C6_M M_G2C6H1_M": 5,
        "M_G2C6_M M_G2C6H2_M": 5,
        "M_G2C7_M M_G2C7H1_M": 6,
        "M_G2C7_M M_G2C7H2_M": 6,
        "M_G2C8_M M_G2C8H1_M": 7,
        "M_G2C8_M M_G2C8H2_M": 7,
        "M_G2C9_M M_G2C9H1_M": 8,
        "M_G2C9_M M_G2C9H2_M": 8,
        "M_G2C10_M M_G2C10H1_M": 9,
        "M_G2C10_M M_G2C10H2_M": 9,
        "M_G2C11_M M_G2C11H1_M": 10,
        "M_G2C11_M M_G2C11H2_M": 10,
        "M_G2C12_M M_G2C12H1_M": 11,
        "M_G2C12_M M_G2C12H2_M": 11,
        "M_G2C13_M M_G2C13H1_M": 12,
        "M_G2C13_M M_G2C13H2_M": 12,
        "M_G2C14_M M_G2C14H1_M": 13,
        "M_G2C14_M M_G2C14H2_M": 13,
        "M_G2C15_M M_G2C15H1_M": 14,
        "M_G2C15_M M_G2C15H2_M": 14,
        "M_G2C16_M M_G2C16H1_M": 15,
        "M_G2C16_M M_G2C16H2_M": 15,
        "M_G2C17_M M_G2C17H1_M": 16,
        "M_G2C17_M M_G2C17H2_M": 16,
        "M_G2C17_M M_G2C17H3_M": 16,
        "M_G2C18_M M_G2C18H1_M": 17,
        "M_G2C18_M M_G2C18H2_M": 17,
        "M_G2C18_M M_G2C18H3_M": 17,
        "M_G2C19_M M_G2C19H1_M": 18,
        "M_G2C19_M M_G2C19H2_M": 18,
        "M_G2C19_M M_G2C19H3_M": 18,
    }

    HGcarbons = {  # noqa: N806
        "M_G3N6C1_M M_G3N6C1H1_M": 1,
        "M_G3N6C1_M M_G3N6C1H2_M": 1,
        "M_G3N6C1_M M_G3N6C1H3_M": 1,
        "M_G3N6C2_M M_G3N6C2H1_M": 1,
        "M_G3N6C2_M M_G3N6C2H2_M": 1,
        "M_G3N6C2_M M_G3N6C2H3_M": 1,
        "M_G3N6C3_M M_G3N6C3H1_M": 1,
        "M_G3N6C3_M M_G3N6C3H2_M": 1,
        "M_G3N6C3_M M_G3N6C3H3_M": 1,
        "M_G3C5_M M_G3C5H1_M": 2,
        "M_G3C5_M M_G3C5H2_M": 2,
        "M_G3C4_M M_G3C4H1_M": 3,
        "M_G3C4_M M_G3C4H2_M": 3,
        "M_G3_M M_G3H1_M": 4,
        "M_G3_M M_G3H2_M": 4,
        "M_G2_M M_G2H1_M": 5,
        "M_G1_M M_G1H1_M": 6,
        "M_G1_M M_G1H2_M": 6,
    }

    for key in OPsim:
        if "M_G1C" in key:
            try:
                xValuesSN1.append(sn1carbons[key])
                yValuesSN1sim.append(float(OPsim[key][0][0]))
                yValuesSN1simERR.append(float(OPsim[key][0][2]))
                yValuesSN1exp.append(OPexp[key][0][0])
                xValuesSN1exp.append(sn1carbons[key])
            except Exception:
                pass
        elif "M_G2C" in key:
            try:
                xValuesSN2.append(sn2carbons[key])
                yValuesSN2sim.append(float(OPsim[key][0][0]))
                yValuesSN2simERR.append(float(OPsim[key][0][2]))
                yValuesSN2exp.append(OPexp[key][0][0])
                xValuesSN2exp.append(sn2carbons[key])
            except Exception:
                pass
        elif "M_G3" in key or "M_G2_M" in key or "M_G1_M" in key:
            try:
                xValuesHG.append(HGcarbons[key])
                yValuesHGsim.append(float(OPsim[key][0][0]))
                yValuesHGsimERR.append(float(OPsim[key][0][2]))
                yValuesHGexp.append(OPexp[key][0][0])
                xValuesHGexp.append(HGcarbons[key])
            except Exception:
                pass
    plt.rc("font", size=15)
    plt.errorbar(
        xValuesHGexp,
        yValuesHGexp,
        yerr=0.02,
        fmt=".",
        color="black",
        markersize=25,
    )
    plt.errorbar(
        xValuesHG,
        yValuesHGsim,
        yerr=yValuesHGsimERR,
        fmt=".",
        color="red",
        markersize=20,
    )
    my_xticks = ["\u03b3", "\u03b2", "\u03b1", "$g_{1}$", "$g_{2}$", "$g_{3}$"]
    plt.xticks([1, 2, 3, 4, 5, 6], my_xticks, size=20)
    plt.yticks(size=20)
    plt.ylabel(r"$S_{CH}$", size=25)
    plt.savefig("HG.pdf")
    plt.show()

    plt.text(2, -0.04, "sn-1", fontsize=25)
    plt.xticks(np.arange(min(xValuesSN1), max(xValuesSN1) + 1, 2.0))
    plt.plot(xValuesSN1, yValuesSN1sim, color="red")
    plt.plot(xValuesSN1exp, yValuesSN1exp, color="black")
    plt.errorbar(
        xValuesSN1,
        yValuesSN1sim,
        yerr=yValuesSN1simERR,
        fmt=".",
        color="red",
        markersize=25,
    )
    plt.errorbar(
        xValuesSN1exp,
        yValuesSN1exp,
        yerr=0.02,
        fmt=".",
        color="black",
        markersize=20,
    )
    plt.ylabel(r"$S_{CH}$", size=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.savefig("sn-1.pdf")
    plt.show()

    plt.text(2, -0.04, "sn-2", fontsize=25)
    plt.xticks(np.arange(min(xValuesSN2), max(xValuesSN2) + 1, 2.0))
    plt.plot(xValuesSN2, yValuesSN2sim, color="red")
    plt.plot(xValuesSN2exp, yValuesSN2exp, color="black")
    plt.errorbar(
        xValuesSN2,
        yValuesSN2sim,
        yValuesSN2simERR,
        fmt=".",
        color="red",
        markersize=25,
    )
    plt.errorbar(
        xValuesSN2exp,
        yValuesSN2exp,
        yerr=0.02,
        fmt=".",
        color="black",
        markersize=20,
    )
    plt.xlabel("Carbon", size=25)
    plt.ylabel(r"$S_{CH}$", size=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.savefig("sn-2.pdf")
    plt.show()


def plotSimulation(system, lipid: str):  # noqa: N802
    """
    Creates plots of form factor and C-H bond order parameters for the selected
    ``lipid`` from a simulation given by system.

    :param system: FAIRMD Lipids ID number of the simulation
    :param lipid: universal molecul name of the lipid

    """
    path = os.path.join(FMDL_SIMU_PATH, system["path"])
    ff_path_sim = os.path.join(path, "FormFactor.json")
    op_path_sim = os.path.join(path, lipid + "OrderParameters.json")
    ffqual_fpath = os.path.join(path, "FormFactorQuality.json")

    print("DOI: ", system["DOI"])

    try:
        with open(ffqual_fpath) as json_file:
            ff_quality = json.load(json_file)
        print("Form factor quality: ", ff_quality[0])
        ffdir = os.path.join(FMDL_EXP_PATH, "FormFactors", system["EXPERIMENT"]["FORMFACTOR"])
        for subdir, _, files in os.walk(ffdir):
            for filename in files:
                if filename.endswith("_FormFactor.json"):
                    ff_path_exp = subdir + "/" + filename
        with open(ff_path_exp) as json_file:
            ff_exp = json.load(json_file)
    except Exception:
        print("Force field quality not found")

    with open(op_path_sim) as json_file:
        op_sim = json.load(json_file)

    op_exp = {}
    for exp_op_folder in list(system["EXPERIMENT"]["ORDERPARAMETER"][lipid].values()):
        op_path_exp = os.path.join(FMDL_EXP_PATH, "OrderParameters", exp_op_folder, lipid + "_Order_Parameters.json")
        with open(op_path_exp) as json_file:
            op_exp.update(json.load(json_file))

    try:
        with open(ff_path_sim) as json_file:
            ff_sim = json.load(json_file)
        plotFormFactor(ff_sim, 1, "Simulation", "red")
        plotFormFactor(ff_exp, ff_quality[1], "Experiment", "black")
        plt.show()
    except Exception:
        plt.show()
        print("Form factor plotting failed")

    plotOrderParameters(op_sim, op_exp)
