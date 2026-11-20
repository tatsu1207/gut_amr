#!/usr/bin/env python3
"""Write every figure panel as its own PNG (300 dpi) to ~/vault/kjm/figures/panels_png/.

Figs 3-5 already exist as separate panel PNGs and are copied. Figs 1, 2 and 6 are drawn as
one matplotlib figure each; their plotting scripts are re-run with saving redirected to a
temporary folder (the submitted images are not touched), and each panel is cropped to the
tight bounding box of its axes.
"""
import os, sys, runpy, shutil, tempfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mf
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

E = os.path.expanduser
OUT = E("~/vault/kjm/figures/panels_png/")
shutil.rmtree(OUT, ignore_errors=True)  # stale panel names must not survive a layout change
os.makedirs(OUT, exist_ok=True)
TMP = tempfile.mkdtemp()
D, H = E("~/vault/kjm/miseq/dada2_output/figures/"), E("~/vault/kjm/hiseq/figures/")

COPY = {
    "fig3A_genus_trajectories": D + "genus_trajectories.png",
    "fig3B_genus_heatmap": D + "genus_heatmap.png",
    "fig4A_total_arg_by_diet": H + "fig4a_total_arg.png",
    "fig4B_arg_class_composition": H + "fig4b_arg_class_composition.png",
    "fig4C_arg_foldchange_heatmap": H + "fig4c_arg_foldchange_heatmap.png",
    "figS1_rarefaction": D + "figS1_rarefaction.png",
    "figS2_scfa": D + "figS2_scfa.png",
    "figS3_taxonomy_stacked_bars": D + "taxonomy_stacked_bars.png",
    "figS4_top_args_heatmap": H + "fig4d_top_args_heatmap.png",
    "fig5A_assembly_arg_classes": H + "fig5a_assembly_arg_classes.png",
    "fig5B_top_arg_genera": H + "fig5b_top_arg_genera.png",
    "fig5C_lactococcus_vs_mls": H + "fig5c_lactococcus_vs_mls.png",
    "fig5D_genomic_context": H + "fig5d_genomic_context.png",
    "fig5E_assembly_args_by_diet": H + "fig5e_assembly_args_by_diet.png",
}
# script, name of the saved PNG to capture, panels as (name, axes indices)
SPLIT = [
    (D + "plot_fig1_physiology.py", "fig1_physiology.png",
     [("fig1A_body_weight", [0]), ("fig1B_LDL", [1]), ("fig1C_HDL", [2]),
      ("fig1D_acetate", [3]), ("fig1E_propionate", [4]), ("fig1F_butyrate", [5])]),
    (D + "plot_fig2_diversity.py", "fig2_diversity.png",
     [("fig2A_shannon", [0]), ("fig2B_richness", [1]), ("fig2C_PCoA", [2]), ("fig2D_displacement", [3])]),
    (E("~/vault/kjm/hiseq/butyrate_analysis/plot_fig6_butyrate.py"), "fig6_butyrate_capacity.png",
     [("fig6A_but_trajectories", [0, 1, 2, 3]), ("fig6B_but_vs_cecal_butyrate", [4]), ("fig6C_but_by_genus", [5])]),
]

for name, src in COPY.items():
    shutil.copyfile(src, OUT + name + ".png")

orig = mf.Figure.savefig
for script, target, panels in SPLIT:
    caught = {}
    def hooked(self, fname, *a, **k):
        caught[os.path.basename(str(fname))] = self
        return orig(self, os.path.join(TMP, os.path.basename(str(fname))), *a, **k)
    mf.Figure.savefig = hooked
    _close = plt.close
    plt.close = lambda *a, **k: None
    cwd = os.getcwd(); os.chdir(os.path.dirname(script))
    try:
        runpy.run_path(script, run_name="__main__")
    finally:
        os.chdir(cwd); mf.Figure.savefig = orig; plt.close = _close
    fig = caught[target]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    axes = [ax for ax in fig.axes if ax.get_visible()]
    for pname, idx in panels:
        bb = Bbox.union([axes[i].get_tightbbox(r) for i in idx])
        bb = bb.transformed(fig.dpi_scale_trans.inverted()).padded(0.08)
        fig.savefig(OUT + pname + ".png", dpi=300, bbox_inches=bb)
    plt.close(fig)

shutil.rmtree(TMP)
print("\n".join(sorted(os.listdir(OUT))))
