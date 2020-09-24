# import config file
configfile: "./config.yaml"

# These are the targets we want to make
rule run_all:
    input:
        "SuperNova/{}/outs/report.txt".format(config['samples']['lib_id']),
        expand("{id}_pseudohap2.{hap}.fasta.gz", id=config['samples']['lib_id'], hap = [1, 2])

# As with the some other pipelines we can get weird file names
# and this pipeline starts from scratch
def get_fastqs_one(wildcards):
    from pathlib import Path
    fq_path = Path(config['samples']['fq_path'])
    lanes = config['samples']['lanes']
    fq_files = []
    for lane in lanes:
        loc_path = fq_path / lane
        for fp in loc_path.iterdir():
            # search fastq path lanes for read 1
            if fp.is_file and str(fp).endswith("_1.fq.gz"):
                print(fp)
                fq_files.append(str(fp))
    return fq_files


# As with the some other pipelines we can get weird file names
# and this pipeline starts from scratch
def get_fastqs_two(wildcards):
    from pathlib import Path
    fq_path = Path(config['samples']['fq_path'])
    lanes = config['samples']['lanes']
    fq_files = []
    for lane in lanes:
        loc_path = fq_path / lane
        for fp in loc_path.iterdir():
            # search fastq path lanes for read 2
            if fp.is_file and str(fp).endswith("_2.fq.gz"):
                print(fp)
                fq_files.append(str(fp))
    return fq_files


# aggregate all fq1 files in data/
rule cat_read_one:
    input:
        get_fastqs_one
    output:
        "data/read_1.fq.gz"
    run:
        if len(config['samples']['lanes']) > 1:
            shell("cat {input} > {output}")
        else:
            shell("ln -s ../{input} {output}")


# aggregate all fq2 files in data/
rule cat_read_two:
    input:
        get_fastqs_two
    output:
        "data/read_2.fq.gz"
    run:
        if len(config['samples']['lanes']) > 1:
            shell("cat {input} > {output}")
        else:
            shell("ln -s ../{input} {output}")


# barcodes can also get messy
# sometimes they're from the barcode list and sometimes they're from the RC list
# This function samples the barcodes and attempts to determine the appropriate list
def determine_barcode_list(n_samples):
    import gzip
    import sys

    # set potential barcode files based on the config file, also set the fastq path
    barcode_file = config['params']['toolsdir'] + config['params']['barcode']
    barcode_rc_file = config['params']['toolsdir'] + config['params']['barcode_RC']
    fastq_path = "data/read_2.fq.gz"

    # Read in barcodes fro, the barcodes file
    # We allow for one mismatch for every barcode
    # so we add those as well
    def get_barcodes(barcodes_file):
        print(f"reading in {barcodes_file}", file=sys.stderr)
        nucs = ['A', 'C', 'G', 'T']
        barcodes = []
        with open(barcodes_file, "r") as bcs_file:
            for line in bcs_file:
                bc = line.strip().split()[0]
                # iterate through each nucleotide in the barcode
                for i in range(0, len(bc)):
                    # iterate through potential mismatch nucleotides
                    for nuc in nucs:
                        bc_alt = list(bc)
                        bc_alt[i] = nuc
                        # add mismatched barcode to the barcodes list
                        barcodes.append("".join(bc_alt))

            # return the list as a set, their shouldn't be any duplicates
            # this is mostly to make the search faster
            barcodes_set = set(barcodes)
            return barcodes_set

    # scrape the fastq file for barcodes
    def get_fq_barcodes(fastq_path, n_samples):
        with gzip.open(fastq_path, "rb") as fastq:
            fq_bcs = []
            # keep track of number of barcodes we've looked at
            counter = 1
            print("Getting fastq barcodes", file=sys.stderr)
            for line in fastq:
                counter += 1
                # break after we've sampled enough barcodes
                if counter == n_samples + 1:
                    break
                if counter % 400000 == 0:
                    print(f"Read in {counter/4} lines", file=sys.stderr)
                # append barcodes for each read
                if counter % 4 == 3:
                    seq = line.decode('utf-8').strip()
                    fq_bcs.append(seq[-10:])
                    fq_bcs.append(seq[-26:-16])
                    fq_bcs.append(seq[-42:-32])

            return fq_bcs

    # get barcodes, rc_barcodes, and fq_barcodes
    barcodes = get_barcodes(barcode_file)
    barcodes_rc = get_barcodes(barcode_rc_file)
    fq_bcs = get_fq_barcodes(fastq_path, n_samples)

    bcs_found = 0
    rc_bcs_found = 0
    # check fq barcodes against the barcodes and rc_barcodes
    for bc in fq_bcs:
        if bc in barcodes:
            bcs_found += 1
        if bc in barcodes_rc:
            rc_bcs_found += 1

    print(f"Barcodes: {bcs_found}\nRC_Barcodes: {rc_bcs_found}\n",file=sys.stderr)
    # if more barcodes were found than RC barcodes return the barcodes path
    if bcs_found > rc_bcs_found:
        print("Barcode_list: {}".format(config['params']['barcode']), file=sys.stderr)
        return config['params']['barcode']
    # if more RC barcodes found retunr RC barcodes path
    elif rc_bcs_found > bcs_found:
        print("Barcode_list: {}".format(config['params']['barcode_RC']), file=sys.stderr)
        return config['params']['barcode_RC']
    # if the same number of barcodes were found try again with more samples (unlikely)
    # if we sample more than 10 mil reads, exit - something is wrong
    else:
        if n_samples > 40000000:
            print("Can't determine correct barcode file for sample", file=sys.stderr)
            sys.exit(1)
        else:
            return determine_barcode_list(n_samples * 2)


# We can also get barcodes in two formats
# BGI started omitting the 6bp spacers in between barcodes
# We want to check ti see which format we're getting and run the correct command accordingly
def get_read_diff(read1, read2):
    import sys
    import gzip
    with gzip.open(read1, "r") as r1, gzip.open(read2, "r") as r2:
        r1.readline()
        r2.readline()
        seq1 = r1.readline().strip()
        seq2 = r2.readline().strip()
        # check to see that our supplied read length is the same as the determined barcode
        # exit if they're not
        if len(seq1) != config['params']['read_len']:
            print(f"interpreted read length ({len(seq1)}) does not match supplied read length ({config['params']['read_len']})", file=sys.stderr)
            sys.exit(1)
        # check to see if there are 42 or 30 more bases in fq2 than fq1
        if len(seq2) - len(seq1) == 42:
            print(f"42BP barcode detected", file=sys.stderr)
            return True
        elif len(seq2) - len(seq1) == 30:
            print(f"30BP barcode detected", file=sys.stderr)
            return False
        # if we get something other than 42 or 30 something is wrong
        else:
            print(f"Unknown barcode length detected", file=sys.stderr)
            sys.exit(1)


# rule for splitting reads
rule split_reads:
    input:
        expand("data/read_{i}.fq.gz", i=range(1,3))
    output:
        expand("data/split_read.{i}.fq.gz", i=range(1,3))
    params:
        len = config['params']['read_len'],
        toolsdir = config['params']['toolsdir']
    run:
        # determine which barcode list to use
        params.barcode = determine_barcode_list(400000)
        # check the difference in length then use the appropriate command
        if get_read_diff(input[0], input[1]):
            # this uses the 42 bp split script
            shell("perl {params.toolsdir}/tools/split_barcode_PEXXX_42_reads.pl "
                "{params.toolsdir}{params.barcode} "
                "{input} {params.len} data/split_read "
                "2> data/split_stat_read.err")
        else:
            # this uses the 30 bp split script
            shell("perl {params.toolsdir}/tools/split_barcode_PEXXX_30_reads.pl "
                "{params.toolsdir}{params.barcode} "
                "{input} {params.len} data/split_read "
                "2> data/split_stat_read.err")


# parse the log to keep only barcodes with > min_reads
# (the log has read pairs, not reads so we divide by 2)
rule get_barcode_whitelist:
    input:
        "split_stat_read1.log"
    output:
        "Convert_Reads/barcodes_to_keep.txt"
    params:
        min_reads = config['params']['min_reads']
    run:
        min_reads = params.min_reads / 2
        shell("awk '(NF==3 && $3~/[0-9]+_[0-9]+_[0-9]+/ && $2>={min_reads}){{print}}' {input} > {output}")


# pair the barcodes we'll keep with the supernova whitelist
rule make_tmp_barcodes:
    input:
        whitelist = "Convert_Reads/barcodes_to_keep.txt",
        supernova_bcs = config['params']['supernova_whitelist']
    output:
        "Convert_Reads/tmp_barcode.txt"
    run:
        import subprocess as sbp
        import sys

        def get_lines(infile):
            # get a count of barcodes in a file
            lines = sbp.Popen(['wc', '-l', infile], stdout=sbp.PIPE,
                              stderr=sbp.STDOUT, universal_newlines=True)
            stdout,stderr = lines.communicate()
            line_count = int(stdout.strip().split()[0])
            return line_count

        # get lines in our whitelist
        bc_whitelist_count = get_lines(input.whitelist)
        print(f"Usable barcodes: {bc_whitelist_count}", file=sys.stderr)
        # get lines in 10X whitelist
        supernova_bc_count = get_lines(input.supernova_bcs)
        print(f"SuperNova barcodes: {supernova_bc_count}", file=sys.stderr)
        # get number of times we'll have to repeat the 10X whitelist
        copy_number = (bc_whitelist_count // supernova_bc_count) + 1
        print(f"Copy number: {copy_number}", file=sys.stderr)

        # concatenate the 10X whitelist file that many times
        with open(str(output), "w") as f:
            for i in range(copy_number):
                cat = sbp.Popen(['cat', input.supernova_bcs], stdout=f, universal_newlines=True)
                cat.wait()

        outlines = get_lines(str(output))
        print(f"{output} lines: {outlines}", file=sys.stderr)


# create a linked file of the 10X barcodes and ours
rule barcode_id_link:
    input:
        "Convert_Reads/barcodes_to_keep.txt",
        "Convert_Reads/tmp_barcode.txt"
    output:
        "Convert_Reads/barcode_id_link.txt"
    shell:
        "paste {input} | "
        "awk '(NF==4 && length($1)>0){{print $3,$4}}' OFS='\t' > {output}"


# convert the reads to 10X format
rule convert_reads:
    input:
        expand("data/split_read.{i}.fq.gz", i=range(1,3)),
        "Convert_Reads/barcode_id_link.txt"
    output:
        "Convert_Reads/read-I1_si-TTCACGCG_lane-001-chunk-001.fastq.gz",
        "Convert_Reads/read-R1_si-TTCACGCG_lane-001-chunk-001.fastq.gz",
        "Convert_Reads/read-R2_si-TTCACGCG_lane-001-chunk-001.fastq.gz"
    params:
        perl_script = config['params']['convert_reads_script']
    run:
        inputs = []
        for entry in input:
            inputs.append("../" + entry)
        shell("""cd Convert_Reads ;
                 perl {params.perl_script} {inputs}""")

# symlink the 10X reads to where supernova wants them
rule link_reads:
    input:
        "Convert_Reads/read-I1_si-TTCACGCG_lane-001-chunk-001.fastq.gz",
        "Convert_Reads/read-R1_si-TTCACGCG_lane-001-chunk-001.fastq.gz",
        "Convert_Reads/read-R2_si-TTCACGCG_lane-001-chunk-001.fastq.gz"
    output:
        "SuperNova/Fastqs/sample_S1_L001_I1_001.fastq.gz",
        "SuperNova/Fastqs/sample_S1_L001_R1_001.fastq.gz",
        "SuperNova/Fastqs/sample_S1_L001_R2_001.fastq.gz"
    run:
        for fastq, link in zip(input, output):
            shell("ln -s $(realpath {fastq}) {link}")


# run supernova
rule run_supernova:
    input:
        "SuperNova/Fastqs/sample_S1_L001_I1_001.fastq.gz",
        "SuperNova/Fastqs/sample_S1_L001_R1_001.fastq.gz",
        "SuperNova/Fastqs/sample_S1_L001_R2_001.fastq.gz"
    output:
        "SuperNova/{}/outs/report.txt".format(config['samples']['lib_id'])
    params:
        supernova = config['params']['supernova_install'],
        lib_id = config['samples']['lib_id'],
        max_reads = config['params']['max_reads'],
        max_mem = config['params']['max_mem']
    threads:
        config['threads']['supernova']
    shell:
        "cd SuperNova; "
        "rm -r {params.lib_id}; "
        "{params.supernova} run --id={params.lib_id} "
            "--fastqs=Fastqs "
            "--maxreads={params.max_reads} "
            "--localmem={params.max_mem} "
            "--localcores={threads} "
            "--nopreflight"


# output supernova's pseudohap2 style output
rule generate_psuedohap2:
    input:
        "SuperNova/{}/outs/report.txt".format(config['samples']['lib_id'])
    output:
        expand("{id}_psuedohap2.{hap}.fasta.gz", id=config['samples']['lib_id'], hap = [1, 2])
    params:
        supernova = config['params']['supernova_install'],
        lib_id = config['samples']['lib_id']
    shell:
        "{params.supernova} mkoutput "
            "--style=pseudohap2 "
            "--asmdir=SuperNova/{params.lib_id}/outs/assembly "
            "--outprefix={params.lib_id}_pseudohap2"
