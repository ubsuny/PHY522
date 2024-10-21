# Tutorial for `fastjet`, `pythia`, and `root`




## Get the software

Execute these commands to download the correct software environment and docker image

```
git clone https://github.com/ubsuny/PHY522.git
cd PHY522
docker pull srappoccio/rivet-pythia-uproot:latest
docker run --rm -it -v ${PWD}:${PWD} -w ${PWD} --entrypoint "/bin/bash" -p 8888:8888 srappoccio/rivet-pythia-uproot:latest
```

Then when you are inside the docker image, compile the executables: 

```
make pythia2fastjet
make pythia2root
```

