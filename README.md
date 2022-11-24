# PYTHIA

The main simulation code in this repo is `/pythia/pythia8_simple.cpp`. It takes the following arguments:

```
pythia8_simple <nEvents> <outFile> <ptHatMin> <ptHatMax>
```

Run inside the environment `o2`

```
alias o2="alienv enter VO_ALICE@O2Physics::nightly-20220928-1, VO_ALICE@GCC-Toolchain::v10.2.0-alice2-14"
```

And compile `<code>.cxx` with:

```
make <code>
```

# Submitting jobs
To achieve higher statistics, we can run the pythia code as a job on stoomboot. For this, use `/pythia/submit.sh`. In particular, it is useful to run a large number of jobs with `/pythia/batch.sh`:

```
./batch <nJobs> <jobSize>
```

# Plotting macros
