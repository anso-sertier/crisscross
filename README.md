# crisscross

crisscross is Somatic structural variants (SV) tool that uses WGS data and two
complementary signals from the read alignments to detect structural variants:
(a) discordant pair mapping (wrong read orientation or incorrect insert-size) and
(b) soft-clipping (unmapped first or last bases of read) that allows resolving
SV breakpoints at the base pair. A cluster of discordant pairs and one or two
clusters of soft-clipped reads defined an SV candidate: the discordant pairs
cluster defined two associated regions, possibly on different chromosomes and
the soft-clipped reads cluster(s), located in these regions, pinpointed the
potential SV breakpoint positions. We further checked that the soft-clipped
bases at each SV breakpoint were correctly aligned in the neighborhood of the
associated region. SV events were then classified as germline or somatic
depending on their presence in the matched normal sample.


## Prerequisites

### Required software for compiling and running

* gcc compiler version 4.9.2
* python version 2.7.x
* [bwa version 0.7.5a](https://sourceforge.net/projects/bio-bwa/files/)
* c++ library :
** [boost version 1.55.0](https://www.boost.org/users/history/version_1_55_0.html)
** [bamtools API version 2.3.0](https://github.com/pezmaster31/bamtools)
* [bedtools version 2.25.0](https://github.com/arq5x/bedtools2/releases/tag/v2.25.0)
* [samtools](https://github.com/samtools/samtools) 

#### Compilation

Two C++ scripts needs to be compiled (no makefile available yet ...)

```sh
g++ -o bam_extract_ab_sc_check bam_extract_ab_sc_check.cpp -I /PATH/bamtools/include -L /PATH/bamtools/lib -I /PATH/boost_1_55_0/include/ -L /PATH/boost_1_55_0/lib -lbamtools -lz -lboost_regex -lboost_program_options
g++ -o ab_clustering           ab_clustering.cpp           -I /PATH/bamtools/include -L /PATH/bamtools/lib -I /PATH/boost_1_55_0/include/ -L /PATH/boost_1_55_0/lib -lbamtools -lz -lboost_program_options
``` 

### Running

This tool need paired tumor/normal bams to run (bwa alignment if possible).


#### Parallel running

#### Linear running

Change variables name or path in ``make_config_linear.sh```, then run

```sh
make_config_linear.sh
launch_linear_pipeline.sh
```


## Authors

* Anne-Sophie Sertier (anne-sophie.sertier[at]lyon.unicancer.fr)

## Licence

This project is licensed under the GPL License - see the [LICENSE](LICENSE)
file for details

