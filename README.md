# ExpEvoAnalyzer
A snakemake pipeline to analyse experimental evolution data.

## Dependencies:

* python>=3.9
* bakta
* ska2
* spades
* pandas
* biopython
* snakemake>=9.6
* snpeff
* shovill
* pyvcf

## To install:

Install the required packages using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[mamba](https://github.com/mamba-org/mamba):

```
git clone https://github.com/samhorsfield96/ExpEvoAnalyzer.git
cd ExpEvoAnalyzer
mamba env create -f environment.yml
mamba activate expevoanalyzer
```

## Running:

Update `config.yaml` to specify workflow and directory paths.
- `output_dir`: path to output directory. Does not need to exist prior to running.
- `input_dir`: path to directory containing paired-end FASTQ files (gzipped or uncompressed), one per line.
- `reference_reads_R1`: path to FASTQ of read 1 for the reference isolate.
- `reference_reads_R2`: path to FASTQ of read 2 for the reference isolate.
- `shovill_params`: additional parameters for [shovill](https://github.com/tseemann/shovill) e.g. `--trim`
- `ksize`: integer, k-mer size used for [SKA2](https://github.com/bacpop/ska.rust)
- `snpEFF_config`: path to snpEff config file, usually `scripts/snpEff.config`
- `bakta_db`: path to [bakta](https://github.com/oschwengers/bakta) database. This needs to be downloaded and updated prior to running (see [here](https://github.com/oschwengers/bakta#database-download))

Run snakemake:

```
snakemake --cores <cores>
```