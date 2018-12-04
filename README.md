# Segmented HAPlotype Estimation and Imputation Tools version 4 (SHAEPIT4)

SHAPEIT4 is a fast and accurate method for estimation of haplotypes (aka phasing) for SNP array and sequencing data. The version 4 is a refactored and improved version of the SHAPEIT algorithm with multiple key additional features:
- It includes a Positional Burrow Wheeler Transform (PBWT) based approach to quickly select a small set of informative conditioning haplotypes to be used when updating the phase of an individual.
- It leverages phase information contained in sequencing reads. It is designed to be run after a first step of haplotype assembly that regroups heterozygous genotypes into phase sets.
- It accounts for sets of pre-phased genotypes (i.e. haplotype scaffold). The scaffold can be derived either from family data or large reference panels.
- It reads and writes files using HTSlib for better I/O performance in either VCF or BCF formats.
- It re-implements genotype graph and HMM routines or better hardware usage and performance.
- It provides the code in an open source format (licence MIT) on github.

If you use the SHAPEIT4 in your research work, please cite the following paper:

O. Delaneau, M. Robinson, JF. Zagury, J. Marchini, ET. Dermitzakis. Integrative haplotype estimation with sub-linear complexity. BioRxiv 2018.

## Documentation

GO ON THIS [WEBSITE] (https://odelaneau.github.io/shapeit4/).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
