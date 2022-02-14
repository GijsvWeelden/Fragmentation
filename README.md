# Hybrid_Model

IMPORTANT: ALWAYS ACTIVATE THE ALICE ENVIRONMENT WITH `ali=alienv enter VO_ALICE@AliPhysics::vAN-20200108_ROOT6-1,VO_ALICE@HepMC::HEPMC_02_06_10-1`

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
