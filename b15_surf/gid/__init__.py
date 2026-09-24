"""B15-GID — blind grammar identification from pole / zero / split fingerprints.

The pass reads ONLY the label-stripped dataset (``b15_surf.zero.dataset.label_stripped``)
and never imports ``b15_surf.zero`` evaluators (WP3.1). Detection of poles,
zeros and split ranks is SOUND (exact rational arithmetic through the frozen
``pir.symbolic.linear`` bridge and the exact RANK_TEST of ``b3_electroweak.rank``).
"""
