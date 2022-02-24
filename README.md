# Hybrid_Model

IMPORTANT: ALWAYS ACTIVATE THE ALICE ENVIRONMENT WITH `ali`, where `.zshrc` should contain:

```
source /cvmfs/alice.cern.ch/etc/login.sh
ali=alienv enter VO_ALICE@AliPhysics::vAN-20200108_ROOT6-1,VO_ALICE@HepMC::HEPMC_02_06_10-1
```

Compile `<code>.cxx` with:

```
make <code>
```

# Jewel simulation
Simulate events with jewel with:

```
./analyze_hepmc_jet_shapes_constsub_eventwise_treeout /user/marcovl/jewel/run_pp_5tev02/1/example.hepmc jet_tree_pp5tev
```
This will output `jet_tree_pp5tev<settings>.root` (with `<settings>` "`_nobkg`" etc).

This code can be ran on the cluster with the `submit_files.sh` and `analyze_files_batch.sh` scripts.

```
./submit_files.sh /dcache/alice/marcovl/jewel/<set> <start_index> <end_index> <dirs_per_job>
```

# Directory Structure

## Plotting macros
`trees_to_hists.C`: Reads in trees from root file and makes histograms for each observable (for the complete set, leading jets, and awayside jets separately). The histograms are 2d (pt:obs) to allow for pt cuts/binning in later analysis.
`plot_jetprops.C`: Reads in histograms from root file and makes plots for various pt ranges, distinguishing jewel settings (`pp`, `AA`) and leading/awayside jets. The plots are saved in the `plots` directory.

## Plots
Contains plots for various jewel settings, e.g. `2tev76_full_nobkg` meaning full jets, without background subtraction, for energy 2.76 TeV. Distinguishes two main categories of plots: `Away` and `ppAAnrAAr`.
### Away
These plots compare the leading and awayside jets for each respective jewel setting, e.g. `pp` leading vs `pp` awayside.
### ppAAnrAAr
These plots compare jets across jewel settings, i.e. `pp`, `AAnr`, `AAr`. There are separate plots for comparing the complete jet sample, as well as only the leading or awayside jets.