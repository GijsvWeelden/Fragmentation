# PYTHIA

The main simulation code in this repo is `/pythia/pythia8_simple.cpp`. It takes the following arguments:

```
pythia8_simple <nEvents> <outFile> <ptHatMin> <ptHatMax>
```

inside the environment `alias o2="alienv enter VO_ALICE@O2Physics::nightly-20220928-1, VO_ALICE@GCC-Toolchain::v10.2.0-alice2-14"`.

compile `<code>.cxx` with:

```
make <code>
```

# Submitting jobs
To achieve higher statistics, we can run the pythia code as a job on stoomboot. For this, use `/pythia/batch_pythia8_simple.sh`.

# Plotting macros
