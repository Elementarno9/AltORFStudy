# AltORFStudy

**Students**

* Akhmetgaliev Eduard (F)

* Kershinskaya Ekaterina (S)

**Supervisor**

Konstantin Senkevich, McGill

### Aim

Analyze data from OpenVar to make a list of SNPs that are nonpathogenic in the canonical ORF along with altering the alternative ORF with moderate/high impact and affecting expression in brain tissues as candidates for future functionality studies.

**About OpenVar**

OpenVar is the first tool for genomic variant annotation and functional effect prediction supporting deep open reading frame (ORF) annotation and polycistronic annotation of Human, Mouse, Rat and Fruit fly transcripts. 
OpenVar builds on the well-known and extensively used SNPEff tool ([Cingolani et al., 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3679285/)), but also offers the possibility to predict variant effect in alternative ORFs as defined in [OpenProt](https://openprot.org) ([Brunet et al., 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5932603/)).


**Data**

The available data at the start of the project were: OpenVar results for 16 VCF with GWAS SNPs associated with 7 neurological and 9 psychiatric disorders.
OpenVar results for every disorder contained:
* the submitted input vcf (input_vcf.vcf);
* a text file of all analysis warnings (warnings.txt). This file contain all SNPs that were not included in the analysis and the reason why (e.g. none of the alleles match the allele of the reference genome at the given position);
* a .tsv file listing all consequences for each variant (study_name_annOnePerLine.tsv). This file thus contain as many lines as there are consequences for all of the submitted variants;
* a .tsv file listing only the maximal impact on canonical and alternative ORFs for each variant (study_name_max_impact.tsv). This file thus contain as many lines as the number of submitted variants; 

*OpenVar results publication was not allowed*.

### Workflow

<img width="459" alt="image" src="https://user-images.githubusercontent.com/90474946/169116816-b1801eb1-d814-425d-b68f-42a45b8c06c8.png">

Aggregation and filtering was made with python3.9 and pandas. GTEx API can be found [here](https://gtexportal.org/home/api-docs/index.html).

### Results summary

**Neurological disorders**

Disorder | Total number of variants | Number of SNP-candidates |
------ | --- | --- |
Alzheimer's disease | 2369 | 72 |
Dementia with Lewy Bodies | 190 | 6 |
Headache | 3346 | 183 |
Neuromyelitis | 653 | 43 |
Parkinsons's disease | 3465 | 163 |
Prion | 41 | 3 |
Stroke | 360 | 5 |

**Psychiatric disorders**

Disorder | Total number of variants | Number of SNP-candidates |
------ | --- | --- |
Attention deficit hyperactivity disorder | 317 | 4 |
Alcohol Dependence | 6 | 2 |
Anorexia Nervosa | 294 | 28 |
Autism spectrum disorder | 92 | 3 |
Bipolar Disorder | 240 | 7 |
Major Depressive Disorder | 49 | 0 |
Posttraumatic Stress | 5 | 0 |
Schizoaffective Disorder | 22337 | 791 |
Tourette syndrome | 1 | 0 |


**Graphs**

For each disorder there is a Sankey diagram showing impact change between refORF and altORF for SNP-candidates. For example, the diagram for Parkinsons's disease, only for max impacted SNPs. 

<img width="937" alt="image" src="https://user-images.githubusercontent.com/90474946/169359077-49e8d5a0-f805-4446-8f7f-29a33c398cc0.png">

### Discussion

There is an opinion that the alternative ORFs represent an overlooked layer of complexity in the coding potential of genomes and are transforming our current vision of the nature of coding genes. Our investigation shows that there are a lot of SNPs that potentially could cause a disorder due to their effect on the alternative ORFs. These SNPs are planned to become candidates for future functionality studies.

### Program

Our research can be repeated with the provided program.

#### Installation and example run

The program requires **python 3.6+** and all libraries from `requirements.txt`.

**Important:** Input data required to be in a format of OpenVar.

There is an example parallel run for Parkinson's disease data in folder 'Neurol_PD', where output will also be saved:

```
python main.py --input Neurol_PD --parallel
```

**Parameters**

**`--input`**

Path to the folder with OpenVar data.

**`--output`**

Path to the output folder. Default is same as input path.

**`--name`**

The prefix for .tsv file, if it is not the same as input path.

**`--parallel`**

Using a few cores for faster run.

**`--only-sankey`**

Only replot sankey for existing output. Can be used when plotting parameters need to be changed.


### References

1. Brunet MA, Leblanc S, Roucou X. Reconsidering proteomic diversity with functional investigation of small ORFs and alternative ORFs. Exp Cell Res. 2020 Aug 1;393(1):112057. doi: 10.1016/j.yexcr.2020.112057. Epub 2020 May 6. PMID: 32387289.
2. Marie A Brunet, Sébastien Leblanc, Xavier Roucou et al. Openvar: Functional Annotation Of Variants In Non-Canonical Open Reading Frames, 22 March 2022, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-1412249/v1]
3. Zhang D, Guelfi S, Garcia-Ruiz S, Costa B, Reynolds RH, D'Sa K, Liu W, Courtin T, Peterson A, Jaffe AE, Hardy J, Botía JA, Collado-Torres L, Ryten M. Incomplete annotation has a disproportionate impact on our understanding of Mendelian and complex neurogenetic disorders. Sci Adv. 2020 Jun 10;6(24):eaay8299. doi: 10.1126/sciadv.aay8299. PMID: 32917675; PMCID: PMC7286675.
4. Sheshukova, E.V., Shindyapina, A.V., Komarova, T.V. et al. “Matreshka” genes with alternative reading frames. Russ J Genet 52, 125–140 (2016). https://doi.org/10.1134/S1022795416020149
