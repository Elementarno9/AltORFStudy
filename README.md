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

For each disorder there is a Sankey diagram showing impact change between refORF and altORF for SNP-candidates. For example, diagram for Parkinsons's disease.

<img width="937" alt="image" src="https://user-images.githubusercontent.com/90474946/169359077-49e8d5a0-f805-4446-8f7f-29a33c398cc0.png">


### Installation and example run

The program requires **python 3.6+** and all libraries from `requirements.txt`.

**Important: ** Input data required to be in a format of OpenVar.

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
