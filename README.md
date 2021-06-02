# IPFP-for-HPC

An IPFP implentation based on OpenMPI and OpenMP to speed up the process of computing the aggregated matrices used in the paper: Estimating, monitoring, and forecasting COVID-19 epidemics (2020) by Albani et al.

#### Prerequisites

you need to have `OpenMPI` and `OpenMP` on your machine (new versions)

in the branch 'older_version' you can find a versions that can compile with older mpicc versions.

#### Usage

you can try the code in this way:

```
$ git clone https://github.com/thomasreolon/IPFP-for-HPC
$ cd IPFP-for-HPC
$ make run
```

_obviously this application should be run on a cluster with the real dataset and not on your local machine_
