# CGI\_supernova\_Pipeline


This pipeline is for 10X's supernova software for use with stLFR data.
This could definitely be extended to do more.
However, for samples we aren't trying to compare to Sentieon this is about what needs to be done.


## Directory Contents

- Snakefile
    - This snakefile contains all the rules for running supernova with stLFR reads
    - This pipeline is comparatively more barebones and just outputs a supernovas assembly
- config.yaml
    - The config file contains all the parameters necessary for modular assembly and evaluation
    - More info can be found in the comments and below


## Running the pipeline

Most fastq files are stored at `/research/rv-02/Projects/Project_stLFR/`.
Analyses are run at `/research/rv-02/Projects/Project_stLFR/pipeline_analysis/`.
Within `pipeline_analysis/` directories are organized by slide and lanes pertinent to any given samples.
For example, four samples are sequenced on two slides, such that LaneN of both slides contain the same sample.
These will be found under `pipeline_analysis/` as `slide1_slide2/`, within which will be subdirectories `slide1_lane1_slide2_lane1/`, `slide1_lane2_slide2_lane2/`, etc.

Within any given sample analysis directory will be a `fastq/` directory and one or more analysis directories.
The analysis directory for Supernova analysis is typically `Supernova_Analysis`.
Within the fastq directory, the directory containing the appropriate lanes fastq information are symlinked as `slideN_laneN/`.
For example, the  directory `/research/rv-02/Projects/Project_stLFR/pipeline_analysis/slide1_slide2/slide1_lane1_slide2_lane1/fastq/` would contain two subdirectories, `slide1_lane1` and `slide2_lane1`.
These would point to `/research/rv-02/Projects/Project_stLFR/slide1/lane1` and `/research/rv-02/Projects/Project_stLFR/slide2/lane1` respectively.

On the same level as `fastq/` you'll potentially find `Supernova_Analysis/`.
The `config.yaml` file should be copied to within `Supernova_Analysis`.
At this point, all that's needed is to modify `config.yaml` and run snakemake appropriately.
An example of running snakemake is below.

```
# -j specifies the number of threads to use
# -k specifies keep going with independent jobs after an error
# -s specifies the location of the snakefile
# you can also copy the snakefile to the current directory and omit -s
snakemake -j 20 -k -s /research/rv-02/home/eanderson/CGI_supernova_Pipeline/Snakefile 2>&1 | tee snakemake.err.txt
```

## Modifying config.yaml

- samples
    - lib_id:
        - sample or library name, definitely update
    - fq_path: `"../fastq"`
        - shouldn't have to be changed
        - can be changed if you'd like the pipeline to use reads in a different directory structure
    - lanes:
        - Lanes within fastq path to use
- params
    - read_len:
        - modify appropriately
    - min_reads:
        - modify if you have too high of coverage and want to subset a fewer selection of reads
    - max_reads:
        - modify if you have too high of coverage and don't want to guess and check with min_reads
- threads
    - just the threads for supernova, modify at will
