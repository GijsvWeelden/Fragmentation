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

# Accessing TH1 from a TList
In files with many histograms, `TList` objects are used to create structure. When in a root file, they can be accessed by:

```
TFile *f = TFile::Open("sample.root");
TList *myList = (TList*)f->Get("list");
TH1F *hist = (TH1F*)myList->FindObject("hist");
```

# ppClass and macros
The `ppClass.C` will extract the important information from the simulation `.root` file and save it in histograms. It uses the `MakeClass` method.

To make a class:
```
root -l "myfile.root"
.ls
T->MakeClass("MyClass")
```

To use the class, modify the `MyClass::Loop()` function and run it with:
```
root –l MyClass.C
MyClass m
m.Loop()
```

Other macros can be used similarly: `root -l <macro>.C`


# Directory Structure

## Plots
`ppAAnr`, `ppAAr`, `ppAAnrAAr`: compares different settings
`awayside/pp`, `AA_recoil`, `AA_norecoil`: compares leading/awayside jets in a single setting
`awayside/Awayside`, `Leading`: compares leading or awayside jets across different settings

## Macros
`plotting_macros/jetprops_2dhists.C`:
Reads in trees from `run*/jet_shapes*.root` and saves the different variables in 2D histograms of type (pt,observable) in `2dhists_*.root`. When looking into single variables across pt bins and/or settings (pp, AA), this is much faster than reading in the whole tree.

`plotting_macros/awayside_jetprops_hists.C`:
Reads in the trees from `run*/jet_shapes*.root` and saves the different variables in 1D histograms, for leading and awayside jets separately. The histograms are saved in `awayside_*.root`.

`plotting_macros/jetprops_hists.C`:
Reads in trees from `run*/jet_shapes*.root` for various settings and saves the different variables in various pt bins in 1D histograms in the `compare_*.root` files.

`plotting_macros/plot_jetprops.C`:
Reads histograms from `compare_*.root` and plots them for pp, AA. It also plots the ratios AA/pp. The plots are saved in the `plots` directory. Should be absorbed in other plotting macro.

`plotting_macros/plot_awayside_jetprops_jetwise.C`:
Reads histograms from `awayside_*.root` and plots variables in various pt bins, comparing pp, AAnr and AAr.

`plotting_macros/plot_awayside_jetprops_settingwise.C`:
Reads histograms from `awayside_*.root` and plots variables in various pt bins, comparing leading, awayside and full jets. The plots are saved in `plots/<setting>`, distinguishing plots containing the full sample (ALF) and plots only containing leading and awayside jets (AL).

## Marco's macros
`jewel_plot_utils.C`
`legend_utils.C`
`zplot_shapes.C`
`plot_zg.C`
`style.C`

## Old/junk macros
`plot_jet_variables.C`: Deprecated and terrible, don't use! Should be deleted.
`test.C`: For quick tests only.
