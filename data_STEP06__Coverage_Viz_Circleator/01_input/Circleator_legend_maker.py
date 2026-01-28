#!/usr/bin/env python3
"""
Create a standalone SVG legend for a Circleator figure.

Output: circleator_legend.svg

Requires: svgwrite
    pip install svgwrite
"""

import svgwrite

# Colors (hex) - updated for the current Circleator config
COLORS = {
    "Genes": "#4daf4a",                 # green
    "rRNAs": "#e41a1c",                 # red
    "GC content": "#ff8c1a",            # orange
    "Illumina Coverage": "#3aa0d5",     # blue
    "ONT Coverage": "#984ea3",          # violet
}

OUT_SVG = "circleator_legend.svg"

# Layout parameters
padding = 12
swatch = 14            # square size (px)
row_gap = 10           # vertical gap between rows (px)
text_dx = 10           # gap between swatch and text (px)
font_size = 14
font_family = "Arial"

# Compute canvas size
rows = list(COLORS.keys())
width = 320
height = padding * 2 + len(rows) * swatch + (len(rows) - 1) * row_gap

dwg = svgwrite.Drawing(OUT_SVG, size=(f"{width}px", f"{height}px"))
dwg.viewbox(0, 0, width, height)

# White background
dwg.add(dwg.rect(insert=(0, 0), size=(width, height), fill="white", opacity=1.0))

y = padding
x = padding

for label in rows:
    color = COLORS[label]

    # Color swatch
    dwg.add(
        dwg.rect(
            insert=(x, y),
            size=(swatch, swatch),
            fill=color,
            stroke="black",
            stroke_width=0.8,
        )
    )

    # Text label
    dwg.add(
        dwg.text(
            label,
            insert=(x + swatch + text_dx, y + swatch * 0.78),
            font_size=font_size,
            font_family=font_family,
            fill="black",
        )
    )

    y += swatch + row_gap

dwg.save()
print(f"Wrote {OUT_SVG}")
