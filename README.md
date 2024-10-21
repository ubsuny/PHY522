# Tutorial for `fastjet`, `fastjet-contrib`, `pythia`, and `root`


This is a tutorial for PHY522, Advanced Topics in Particle Physics at the University at Buffalo, State University of New York.
It will download two programs that demonstrate the use of
`fastjet` ([software](https://fastjet.fr), [Phys.Lett.B641:57-61,2006](https://arxiv.org/abs/hep-ph/0512210), [CERN-PH-TH/2011-297](https://arxiv.org/abs/1111.6097)),
`fastjet-contrib` ([software](https://fastjet.hepforge.org/contrib/)), specifically the N-subjettiness [JHEP 1103:015,2011](https://arxiv.org/abs/1011.2268) and soft drop [JHEP 1405 (2014) 146](https://arxiv.org/abs/1402.2657) packages, 
`pythia` ([software](https://pythia.org), [LU-TP 22-16, MCNET-22-04](https://arxiv.org/abs/2203.11601), [LU TP 14-36, MCNET-14-22, CERN-PH-TH-2014-190, FERMILAB-PUB-14-316-CD, DESY 14-178, SLAC-PUB-16122](https://arxiv.org/abs/1410.3012))
and `root` [link](https://root.cern) together.

   * `pythia2fastjet` : use `pythia` to generate events with a given configuration and cluster them with various algorithms in `fastjet`. 
   * `pythia2root`: a more industrial-scale analysis that will generate events, cluster stable particles, and write the events and jets to a `.root` file for later analysis.

## Details of physics selection

   * Events will be simulated with `pythia` as per the user-supplied configuration. An example `qcd_multijets.cfg` is provided. 
   * Jets will be constructed with the anti-$k_T$ jet algorithm with R=0.8 [JHEP 0804:063,2008](https://arxiv.org/abs/0802.1189). A minimum jet $p_T$ cut is applied of 50 GeV. 
   * The clustered AK8 jets will then be reclustered with the soft drop algorithm [JHEP 1405 (2014) 146](https://arxiv.org/abs/1402.2657) with various $\beta$ parameters and $z_{cut}=0.1$. 
   * To avoid confusion, jets that contain 90\% of their energy from electrons or muons are not considered as jets.
   * Various n-subjettiness variables are also calculated [JHEP 1103:015,2011](https://arxiv.org/abs/1011.2268) . 

## Get the software

Execute these commands to download the correct software environment and docker image

```
git clone https://github.com/ubsuny/PHY522.git
cd PHY522
docker pull srappoccio/rivet-pythia-uproot:latest
docker run --rm -it -v ${PWD}:${PWD} -w ${PWD} --entrypoint "/bin/bash" -p 8888:8888 srappoccio/rivet-pythia-uproot:latest
```

## Compile the software

Then when you are inside the docker image, compile the executables: 

```
make pythia2fastjet
make pythia2root
```

## Running the software

After you are done compiling, you can use these two programs to do studies.


```
$ ./pythia2fastjet                     
usage: ./pythia2fastjet config_file n_events (verbose?)
```

```
$ ./pythia2fastjet
usage: ./pythia2fastjet config_file n_events (verbose?)
```


An example is provided in `qcd_multijets.cfg`, containing this text:

```
HardQCD:all = on
PhaseSpace:pTHatMin = 100
```



