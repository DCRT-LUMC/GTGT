from collections.abc import Sequence

import svg
from exonviz.draw import draw_exons
from exonviz.mutalyzer import build_exons
from mutalyzer.description import Description

config = {
    "width": 768,
    "height": 20,
    "scale": 1.0,
    "noncoding": True,
    "gap": 0,
    "color": "#4C72B7",
    "exonnumber": True,
    "firstexon": 1,
    "lastexon": 9999,
    "variantcolors": ["#BA1C30", "#DB6917", "#EBCE2B", "#702C8C", "#C0BD7F"],
    "variantshape": "pin",
}


def draw(d: Description) -> str:
    exons, dropped_variants = build_exons(
        hgvs=d.input_description,
        mutalyzer=d.output()["selector_short"],
        config=config,
    )

    # Determine the minimum scale at which we can draw every exon
    min_scale = max([e.min_scale() for e in exons])
    config["scale"] = max(min_scale, 1)

    fig = draw_exons(exons, config=config)

    # Make the figure scalable with CSS by resetting the width and height
    # and setting a viewBox
    width = fig.width
    height = fig.height
    fig.viewBox = svg.ViewBoxSpec(0, 0, width, height)

    fig.width = None
    fig.height = None

    return str(fig)
