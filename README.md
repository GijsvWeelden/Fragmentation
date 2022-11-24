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

# KKP Fragmentation Functions
In the `brick` directory, we have access to the KKP fragmentation functions. First, compile the code in the `kkp` dir, with `./make`. Then, plot the FFs with the `plot_kkp_bm.C`. Note that `plot_kkp_bm_sub.C` does not run standalone, but is run through `plot_kkp_bm.C`.

More useful is to access the fragmentation functions themselves, which can be done as follows:

```
TF1 *frag_f = new TF1("frag_f",kkp_func,0,1,3);
frag_f->SetParameter(0,E); // Q for fragmentation
frag_f->SetParameter(1,1); // parton flavour 1=u
frag_f->SetParameter(2,1); // hadron type: 1=pi+pi-; 4 = p+pbar
```

where the various settings are found in `kkp/kkp.h`.

# Plotting macros
