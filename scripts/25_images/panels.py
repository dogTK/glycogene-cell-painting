"""Genes of the image panels, read from select_genes.py output (results/25_images/gene_selection.csv).
panel 'oe' = overexpression, 'ko' = knockout (both shown in Figure 3c); order of the slots: AGP up, AGP down, ER up.
"""
from pathlib import Path

import pandas as pd

R = Path.home() / 'research/glyco_cell_painting/results/25_images'
SOURCE = {'ORF': 'source_4', 'CRISPR': 'source_13'}


def panel_genes(panel):
    """[(gene, modality, channel, direction), ...] in slot order."""
    s = pd.read_csv(R / 'gene_selection.csv')
    s = s[(s.panel == panel) & s.chosen]
    order = {('AGP', 'up'): 0, ('AGP', 'down'): 1, ('ER', 'up'): 2}
    s = s.assign(o=[order[(c, d)] for c, d in zip(s.channel, s.direction)]).sort_values('o')
    return list(s[['gene', 'modality', 'channel', 'direction']].itertuples(index=False, name=None))
