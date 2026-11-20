"""Display labels for sampling timepoints (2026-11-20).

Data files and sample IDs keep the original timepoint codes (T0h, T24h, T72h, S72h, SD7, SD10,
SD14), which match the SRA submission. Figures show the manuscript's labels instead:
T = days since the first antibiotic dose, R = days since the last dose.

install() wraps Figure.savefig so that, just before saving, every text element of the figure
(tick labels, titles, legends, annotations, heatmap labels) has old codes replaced by new ones.
"""
import re
import matplotlib.figure as _mf
import matplotlib.text as _mt

MAP = {"T0h": "T0", "T24h": "T1", "T72h": "T3", "S72h": "R3", "SD7": "R7", "SD10": "R10", "SD14": "R14"}
_PAT = re.compile(r"(?<![A-Za-z0-9])(T0h|T24h|T72h|S72h|SD7|SD10|SD14)(?![A-Za-z0-9])")


def label(s):
    return _PAT.sub(lambda m: MAP[m.group(1)], str(s))


def relabel(fig):
    # Tick labels are regenerated from the axis formatter at draw time, so fix them on the axis.
    fig.canvas.draw()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            labs = [tl.get_text() for tl in axis.get_ticklabels()]
            if any(_PAT.search(s) for s in labs):
                axis.set_ticks(axis.get_ticklocs())
                axis.set_ticklabels([label(s) for s in labs])
    for t in fig.findobj(_mt.Text):
        s = t.get_text()
        if s and _PAT.search(s):
            t.set_text(label(s))


def install():
    if getattr(_mf.Figure.savefig, "_tplabels", False):
        return
    orig = _mf.Figure.savefig

    def savefig(self, *a, **k):
        relabel(self)
        return orig(self, *a, **k)
    savefig._tplabels = True
    _mf.Figure.savefig = savefig
