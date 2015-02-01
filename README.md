ICGC-TCGA-PCAWG SMC Validation
==============================


Synapse project:  [syn2875157](https://www.synapse.org/#!Synapse:syn2875157/wiki/)
Evaluation ID:    3060780


Need:

* Synapse entities for VCF files (to be kept in GNOS)
* May have to compute tumor-only VCF files
* "Truth" VCF files
* Provenance linking back to TCGA and ICGC data sources


* Submit VCFs to evaluation queue
* Score tumor-only VCF files against "Truth" files
* Compute confusion matrix for each VCF and store (in a table?)





Plan:
* write script to create table of:
  - synapse id
  - tumor sample id
  - normal sample id
  - variant type
  - source (DKFZ/EMBL/SANGER/UCSC)
  - alorithm / workflow


Questions:
* Looks like all VCFs under Annai are germline?
* Any more info needed?

