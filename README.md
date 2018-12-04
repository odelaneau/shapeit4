# Segmented HAPlotype Estimation and Imputation Tools version 4 (SHAEPIT4)

SHAPEIT4 is a refactored version of SHAPEIT version 1 to 3 with multiple additional features.
SHAPEIT4:
- Quiclky selects conditioning haplotypes using a PBWT of the haplotypes (sub-linear scaling with sample size).
- Leverages phase sets when specified in the input VCF/BCF (uses phase informative reads).
- Accounts for sets of pre-phased variants (uses haplotype scaffold). 
- Reads/Writes files using HTSlib (Better I/O performance).
- Re-implements genotype graph and HMM routines (much faster implementation).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

SHAPEIT4 is built on two libraries that need to be installed on your system:

* [HTSlib] (http://www.htslib.org/) - A C library for reading/writing high-throughput sequencing data
* [BOOST] (https://www.boost.org/) - Peer-reviewed portable C++ source libraries

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

## History of development


This document details the development history of the SHAPEIT software and web package. We refer to the following people using their initials as follows - Olivier Delaneau (OD), Jean-Francois Zagury (JFZ), Jonathan Marchini (JM), Jared O’Connell (JO).

Jan 2007-Dec 2008 : The ground work for the development of SHAPEIT was carried out during this time while OD was graduate student under the supervision of JFZ, leading to an initial ShapeIT publication [3]. During this time OD was paid by CNAM. 

Jan 2009 – Dec 2010 - After his PhD, OD has pursued the development of SHAPEIT as a post-doc, under the supervision of JFZ. This work has lead to the first package SHAPEIT v1 released on the website on 20/11/2011. During these years OD was paid on some of JFZ’s grants and by Peptinov. The main methodological developments in SHAPEIT v1 are (i) using a compact representation of the inferred haplotypes using a graphical model, (ii) representing the space of possible haplotypes using a markov chain, which has the downstream benefit of making the core calculations linear in the number of conditioning states. The SHAPEIT v1 paper [5] was written by OD, JM and JFZ. The paper was written between Oct 2010 and March 2011, when it was submitted. The paper was published on 4 Dec 2011.

Jan 2011 – Dec 2011 – OD moved to Oxford University to work also in collaboration with Jonathan Marchini, while still being paid by Peptinov. During the period the SHAPEIT v1 method continued to be refined and tested ahead of release and publication. In addition, the development of SHAPEIT v2 occured, which involved fusing together the ideas of SHAPEIT v1 and IMPUTE v2 (local selection of conditionning haplotypes using Hamming distance minimization). SHAPEIT v2 was released on 25/09/2013.

Jan 2012 - August 2013 - OD has worked under the supervision of JM and was paid from JM’s MRC grants and Leverhulme award . The SHAPEIT v2 method continued to be developed during this period. The SHAPEIT v2 paper was written by OD, JM and JFZ. The paper was submitted on 16 June 2012. The paper [6] was published on 27 Nov 2011. In addition, OD and JM worked on (i) adapting the method to allow for phasing genotype likelihoods derived from sequence data, which resulted in a publication [9], (ii) a method of phasing using sequencing reads, which resulted in a publication [7]. During this same period JO worked with JM and OD to develop a way of using SHAPEIT v2 for phasing related individuals [8], and for scaling the method to handle many hundreds of thousands of samples.

The published scientific articles corresponding to these works are indicated  on the next page.

Publications related to SHAPEIT

-The background work has been published in these publications :

[1] Stephens M, Smith NJ, Donnelly P (2001) A new statistical method for haplotype reconstruction from population data. American journal of human genetics 68(4) 978-989

[2] Stephens M, Donnelly P (2003) A comparison of bayesian methods for haplotype reconstruction from population genotype data. American journal of human genetics 73(5):1162-1169.

[3] Delaneau O, Coulonges C, Zagury JF (2008). SHAPE-IT : a new rapid and accurate haplotyping software. BMC Bioinformatics, 9:540-53.

[4] B. Howie, P. Donnelly, J. Marchini (2009) A Flexible and Accurate Genotype Imputation Method for the Next Generation of Genome-Wide Association Studies. PLoS Genetics  5(6): e1000529

-The core of machinery of the ShapeIT algorithm is described in the two following publications :

[5] O. Delaneau, J. Marchini, JF. Zagury (2012) A linear complexity phasing method for thousands of genomes. Nat Methods. 9(2):179-81. doi: 10.1038/nmeth.1785 
(The paper was written between Oct 2010 and March 2011, when it was submitted. The paper was published on 4 Dec 2011).

[6] O. Delaneau, JF. Zagury, J. Marchini (2013) Improved whole chromosome phasing for disease and population genetic studies. Nat Methods. 10(1):5-6. doi: 10.1038/nmeth.2307 
(The paper was submitted on 16 June 2012. The paper was published on 27 Nov 2012).


-The applications and extensions of the core algorithm to various popular types of data is described in the 3 following papers:

[7] O. Delaneau, B. Howie, A. Cox, J-F. Zagury, J. Marchini (2013) Haplotype estimation using sequence reads. American Journal of Human Genetics 93 (4) 787-696

[8] J. O'Connell, D. Gurdasani, O. Delaneau, et al. (2014) A general approach for haplotype phasing across the full spectrum of relatedness. PLoS Genetics 

[9] O. Delaneau, J. Marchini, The 1000 Genomes Project Consortium (2014) Integrating sequence and array data to create an improved 1000 Genomes Project haplotype reference panel. (in review)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc