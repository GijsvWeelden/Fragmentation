This directory contains plotting code for data. It also contains a number of scripts to simplify the interaction with the plotting code.

# plotCutVariation.C
Plots purity and efficiency of V0 sample when cutting on one or more of their variables.

# v0cutvarsingle
Runs `plotCutVariation.C`
Specify input, data set, hadron, variable to investigate, pt, fit settings and rebinning.

# v0cutvar
Same as `v0cutvarsingle`, but runs all variables.

# plotJet.C
# plotPurity.C
# plotV0.C
# plotV0MassinJet.C
# plotV0inJet.C

# plotTrackQA.C
QA plots for V0 tracks.

# merge
Merges plots outputted by `plotCutVariation.C`
Mode 1: fits across pt
Mode 2: efficiency
Mode 3: signal-over-background
Mode 4: significance (S/sqrt(S+B))
