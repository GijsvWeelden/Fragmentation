# Hybrid_Model

IMPORTANT: ALWAYS ACTIVATE THE ALICE ENVIRONMENT WITH `ali`

Compile `<code>.cxx` with:

```
make <code>
```

# Jewel simulation
Simulate events with jewel with:

```
./analyze_hepmc_jet_shapes_constsub_eventwise_treeout user/marcovl/jewel/run_pp_5tev02/1/example.hepmc jet_tree_pp5tev
```
This will output `jet_tree_pp5tev<settings>.root` (with `<settings>` "`_nobkg`" etc).

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
