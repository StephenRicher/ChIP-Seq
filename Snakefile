import pandas as pd

BASE = workflow.basedir

ENVS = f'{BASE}/workflow/envs'
SCRIPTS = f'{BASE}/workflow/scripts'

configfile : f'{BASE}/config/config.yaml'
GENOME = config['genome']
BUILD = config['build']
FASTQ_SCREEN_CONFIG = config['fastq_screen_config']

try:
    workdir : config['outdir']
except KeyError:
    pass

DATA_TABLE = config['data']
DATA = pd.read_table(config['data'], sep = ',', dtype = {'rep' : str})

# Validate read file input with wildcard definitions
if not DATA['group'].str.match(r'[^\/\s.-]+').all():
    sys.exit(f'Invalid group definition in {DATA_TABLE}.')
if not DATA['rep'].str.match(r'\d+').all():
    sys.exit(f'Invalid replicate definition in {DATA_TABLE}.')
if not DATA['read'].str.match(r'R[12]').all():
    sys.exit(f'Invalid read definition in {DATA_TABLE}.')
if not DATA['type'].str.match(r'input|bound').all():
    sys.exit(f'Invalid read definition in {DATA_TABLE}.')

seq_type = 'se'

if seq_type == "pe":
    single_regex = '[^\/\s.-]+-\d+-(input|bound)-R[12]'
    DATA['single'] = (DATA[['group', 'rep', 'type', 'read']]
        .apply(lambda x: '-'.join(x), axis = 1))
    trimmed_out = ['fastq/trimmed/{sample_type}-R1.trim.fastq.gz',
                   'fastq/trimmed/{sample_type}-R2.trim.fastq.gz']
    cutadapt_cmd = ('cutadapt '
            '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA '
            '-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT '
            '{params.others} --cores {THREADS} '
            '-o {output.trimmed[0]} -p {output.trimmed[1]} '
            '{input} > {output.qc} '
        '2> {log}')
    hisat2_cmd = ('(hisat2 '
            '-x {params.index} -p 6 '
            '-1 {input.reads_in[0]} -2 {input.reads_in[1]} '
            '--summary-file {output.qc} '
        '| samtools view -b > {output.bam}) '
        '2> {log}')
else:
    single_regex = '[^\/\s.-]+-\d+-(input|bound)'
    DATA['single'] = (DATA[['group', 'rep', 'type']]
        .apply(lambda x: '-'.join(x), axis = 1))
    trimmed_out = ['fastq/trimmed/{sample_type}.trim.fastq.gz']
    cutadapt_cmd = ('cutadapt '
            '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA '
            '{params.others} --cores {THREADS} '
            '-o {output.trimmed[0]} {input} > {output.qc} '
        '2> {log}')
    hisat2_cmd = ('(hisat2 '
            '-x {params.index} -p 6 '
            '-U {input.reads_in[0]} '
            '--summary-file {output.qc} '
        '| samtools view -b > {output.bam}) '
        '2> {log}')

# Define single and sample names from definitions.
DATA['sample_type'] = (DATA[['group', 'rep', 'type']]
    .apply(lambda x: '-'.join(x), axis = 1))
DATA['sample'] = (DATA[['group', 'rep']]
    .apply(lambda x: '-'.join(x), axis = 1))
DATA = DATA.set_index(['group', 'sample', 'sample_type', 'single'], drop = False)

# Extract groups and replicates.
GROUPS = {}
for group in DATA['group']:
    GROUPS[group] = list(DATA.loc[group]['rep'].unique())
# Extract sample names
SAMPLES_TYPE = list(DATA['sample_type'].unique())

THREADS = 1

wildcard_constraints:
    group = '[^\/\s.-]+',
    sample = '[^\/\s.-]+-\d+',
    sample_type = '[^\/\s.-]+-\d+-(input|bound)',
    single = single_regex,
    rep = '\d+',
    read = 'R[12]',
    type = 'input|bound'

rule all:
    input:
        ['qc/multiqc', 'qc/multibamqc', f'genome/{BUILD}.fa.gz.fai']

rule bgzip_genome:
    input:
        GENOME
    output:
        f'genome/{BUILD}.fa.gz'
    log:
        f'logs/bgzip_genome/{BUILD}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        '(zcat -f {input} | bgzip --stdout > {output}) 2> {log}'

rule faidx:
    input:
        rules.bgzip_genome.output
    output:
         multiext(f'{rules.bgzip_genome.output}', '.fai', '.gzi')
    log:
        'logs/index_genome/index_genome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'

rule fastqc:
    input:
        lambda wc: DATA.xs(wc.single, level = 3)['path']
    output:
        html = 'qc/fastqc/{single}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'

# Modify FASTQ filename to match {sample}-{read} for multiQC
rule modify_fastqc:
    input:
        'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    output:
        'qc/fastqc/{single}.raw_fastqc.zip'
    params:
        name = lambda wc: f'{wc.single}'
    log:
        'logs/modify_fastqc/{single}.raw.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        '{SCRIPTS}/modify_fastqc.sh {input} {output} {params.name} &> {log}'

rule cutadapt:
    input:
        lambda wc: DATA.xs(wc.sample_type, level = 3)['path']
    output:
        trimmed = trimmed_out,
        qc = 'qc/cutadapt/unmod/{sample_type}.cutadapt.txt'
    group:
        'cutadapt'
    params:
        others = '--minimum-length 20 --quality-cutoff 20 '
                 '--gc-content 46 --overlap 6 --error-rate 0.1'
    log:
        'logs/cutadapt/{sample_type}.log'
    conda:
        f'{ENVS}/cutadapt.yaml'
    threads:
        THREADS
    shell:
        cutadapt_cmd

rule modify_cutadapt:
    input:
        rules.cutadapt.output.qc
    output:
        'qc/cutadapt/{sample_type}.cutadapt.txt'
    group:
        'cutadapt'
    log:
        'logs/modify_cutadapt/{sample_type}.log'
    conda:
        f'{ENVS}/coreutils.yaml'
    shell:
        'awk -v sample={wildcards.sample_type} '
            '-f {SCRIPTS}/modify_cutadapt.awk {input} > {output} 2> {log}'

rule fastq_screen:
    input:
        'fastq/trimmed/{single}.trim.fastq.gz'
    output:
        txt = 'qc/fastq_screen/{single}.fastq_screen.txt',
        png = 'qc/fastq_screen/{single}.fastq_screen.png'
    params:
        fastq_screen_config = FASTQ_SCREEN_CONFIG,
        subset = 100000,
        aligner = 'bowtie2'
    log:
        'logs/fastq_screen/{single}.log'
    threads:
        8
    wrapper:
        "0.49.0/bio/fastq_screen"

rule fastqc_trimmed:
    input:
        'fastq/trimmed/{single}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{single}.trim_fastqc.html',
        zip = 'qc/fastqc/{single}.trim_fastqc.zip'
    log:
        'logs/fastqc_trimmed/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'

rule unzip_genome:
    input:
        rules.bgzip_genome.output
    output:
        f'genome/{BUILD}.fa'
    log:
        f'logs/unzip_genome/{config["build"]}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        'bgzip -cd {input} > {output}'

rule hisat2_build:
    input:
        rules.unzip_genome.output
    output:
        expand('genome/index/{build}.{n}.ht2',
            build = config["build"],
            n = ['1', '2', '3', '4', '5', '6', '7', '8'])
    params:
        basename = f'genome/index/{config["build"]}'
    log:
        f'logs/hisat2_build/{config["build"]}.log'
    conda:
        f'{ENVS}/hisat2.yaml'
    threads:
        6
    shell:
        'hisat2-build --threads {threads} {input} {params.basename} '
            '&> {log}'

rule hisat2_map:
    input:
        reads_in = trimmed_out,
        index = rules.hisat2_build.output
    output:
        bam = 'mapped/{sample_type}.bam',
        qc = 'qc/hisat2/{sample_type}.hisat2.txt'
    params:
        index = f'genome/index/{config["build"]}',
    log:
        'logs/hisat2_map/{sample_type}.log'
    conda:
        f'{ENVS}/hisat2.yaml'
    threads:
        12
    shell:
        hisat2_cmd

rule samtools_sort:
    input:
        rules.hisat2_map.output.bam
    output:
        'mapped/{sample_type}.sort.bam',
    log:
        'logs/samtools_sort/{sample_type}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        12
    shell:
        'samtools sort -@ 6 {input} > {output} 2> {log}'

rule index_bam:
    input:
        rules.samtools_sort.output
    output:
        f'{rules.samtools_sort.output}.bai'
    log:
        'logs/index_bam/{sample_type}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        12
    shell:
        'samtools index -@ 6 {input} &> {log}'

rule samtools_stats:
    input:
        rules.samtools_sort.output
    output:
        'qc/samtools/stats/{sample_type}.stats.txt'
    group:
        'SAM_QC'
    log:
        'logs/samtools_stats/{sample_type}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'

rule samtools_idxstats:
    input:
        bam = rules.samtools_sort.output,
        index = rules.index_bam.output
    output:
        'qc/samtools/idxstats/{sample_type}.idxstats.txt'
    group:
        'SAM_QC'
    log:
        'logs/samtools_idxstats/{sample_type}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools idxstats {input.bam} > {output} 2> {log}'

rule samtools_flagstat:
    input:
        rules.samtools_sort.output
    output:
        'qc/samtools/flagstat/{sample_type}.flagstat.txt'
    group:
        'SAM_QC'
    log:
        'logs/samtools_flagstat/{sample_type}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools flagstat {input} > {output} 2> {log}'

rule bamqc:
    input:
        bam = rules.samtools_sort.output
    output:
        directory('qc/bamqc/{sample_type}')
    resources:
        mem_mb = 3000
    log:
        'logs/bamqc/{sample_type}.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    threads:
        12
    shell:
        'qualimap bamqc '
            '-bam {input.bam} -outdir {output} '
            '-nt 6 --java-mem-size={resources.mem_mb}M '
            '--paint-chromosome-limits '
        '&> {log}'

rule multibamqc_config:
    input:
        expand('qc/bamqc/{sample_type}', sample_type = SAMPLES_TYPE)
    output:
        'qc/bamqc/multibamqc.config'
    group:
        'multiBAM_QC'
    log:
        'logs/multibamqc_config/multibamqc_config.log'
    shell:
        '{SCRIPTS}/multibamqc_config.py {input} > {output} 2> {log}'

rule multibamqc:
    input:
        rules.multibamqc_config.output
    output:
        directory('qc/multibamqc')
    group:
        'multiBAM_QC'
    log:
        'logs/multibamqc/multibamqc.log'
    conda:
        f'{ENVS}/qualimap.yaml'
    shell:
        'qualimap multi-bamqc '
            '--data {input} -outdir {output} '
        '&> {log}'

rule multiqc:
    input:
        expand('qc/fastqc/{single}.raw_fastqc.zip',
            single =  DATA['single']),
        expand('qc/fastq_screen/{single}.fastq_screen.txt',
            single =  DATA['single']),
        expand('qc/cutadapt/{single}.cutadapt.txt', single =  DATA['single']),
        expand('qc/fastqc/{single}.trim_fastqc.zip',
            single =  DATA['single']),
        expand('qc/hisat2/{sample_type}.hisat2.txt',
            sample_type = SAMPLES_TYPE),
        expand('qc/samtools/stats/{sample_type}.stats.txt',
            sample_type = SAMPLES_TYPE),
        expand('qc/samtools/idxstats/{sample_type}.idxstats.txt',
            sample_type = SAMPLES_TYPE),
        expand('qc/samtools/flagstat/{sample_type}.flagstat.txt',
            sample_type = SAMPLES_TYPE),
        expand('qc/bamqc/{sample_type}', sample_type = SAMPLES_TYPE),
    output:
        directory('qc/multiqc')
    log:
        'logs/multiqc/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc --outdir {output} '
            '--force --config {BASE}/config/multiqc_config.yaml {input} '
        '&> {log}'
