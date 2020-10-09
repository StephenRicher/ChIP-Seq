import math
import tempfile
import itertools
import pandas as pd
from snake_setup import set_config, load_samples

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

if not config:
    configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
default_config = {
    'workdir':           workflow.basedir,
    'tmpdir':            tempfile.gettempdir(),
    'threads':           workflow.cores       ,
    'data':              ''          ,
    'paired':            ''          ,
    'computeMatrixScale':
        {'exon':         True        ,},
    'genome':
        {'build':        'genome'    ,
         'index':        ''          ,
         'annotation':   ''          ,
         'sequence':     ''          ,
         'genes':        ''          ,
         'blacklist':    None        ,},
    'cutadapt':
        {'forwardAdapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
         'reverseAdapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
         'overlap':         3                                 ,
         'errorRate':       0.1                               ,
         'minimumLength':   0                                 ,
         'qualityCutoff':  '0,0'                              ,
         'GCcontent':       50                                ,},
    'fastq_screen':      None,
}

config = set_config(config, default_config)

workdir: config['workdir']
BUILD = config['genome']['build']

# Read path to samples in pandas
samples = load_samples(config['data'])

# Extract groups and replicates.
GROUPS = {}
for group in DATA['group']:
    GROUPS[group] = list(DATA.loc[group]['rep'].unique())

# Get group-rep-type for all samples
SAMPLES = list(samples['sample'].unique())
# Get group-rep for input samples
INPUTS = [sample[:-6] for sample in SAMPLES if sample.endswith('input')]
# Get group-rep for bound sample
BOUNDS = [sample[:-6] for bound in SAMPLES if sample.endswith('bound')]
# Generate list of group comparisons - this avoids self comparison
COMPARES = [f'{i[0]}-{i[1]}' for i in itertools.combinations(list(GROUPS), 2)]

wildcard_constraints:
    group = '|'.join(GROUPS),
    bound = '|'.join(BOUNDS),
    sample = '|'.join(SAMPLES),
    single = '|'.join(config['single']),
    rep = '\d+',
    read = 'R[12]',
    type = 'input|bound',
    stage = 'markdup|filtered'

rule all:
    input:
        ['qc/multiqc', f'genome/{BUILD}.fa.gz.fai',
          expand('qc/deeptools/compare{compare}/plot{type}{compare}-{mode}.png',
            mode=['scaled', 'reference'], type=['Heatmap', 'Profile'],
            compare=['Input', 'Group'])]


rule fastQC:
    input:
        lambda wc: samples.xs(wc.single, level=3)['path']
    output:
        html = 'qc/fastqc/{single}.raw_fastqc.html',
        zip = 'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    log:
        'logs/fastqc/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule modifyFastQC:
    input:
        'qc/fastqc/unmod/{single}.raw.fastqc.zip'
    output:
        'qc/fastqc/{single}.raw_fastqc.zip'
    params:
        name = lambda wc: f'{wc.single}'
    log:
        'logs/modifyFastQC/{single}.raw.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/modifyFastQC.py {input} {output} {params.name} &> {log}'


def cutadaptOutput():
    if config['paired']:
        return ['fastq/trimmed/{bound}-R1.trim.fastq.gz',
                'fastq/trimmed/{bound}-R2.trim.fastq.gz']
    else:
        return ['fastq/trimmed/{bound}.trim.fastq.gz']


def cutadaptCmd():
    if config['paired']:
        cmd = ('cutadapt -a {params.forwardAdapter} -A {params.reverseAdapter} '
            '-o {output.trimmed[0]} -p {output.trimmed[1]} ')
    else:
        cmd = 'cutadapt -a {params.forwardAdapter} -o {output.trimmed[0]} '

    cmd += ('--overlap {params.overlap} --error-rate {params.errorRate} '
        '--minimum-length {params.minimumLength} '
        '--quality-cutoff {params.qualityCutoff} '
        '--gc-content {params.GCcontent} '
        '--cores {threads} {input} > {output.qc} 2> {log}')
    return cmd


rule cutadapt:
    input:
        lambda wc: samples.xs(wc.sample, level=1)['path']
    output:
        trimmed = cutadaptOutput(),
        qc = 'qc/cutadapt/unmod/{bound}.cutadapt.txt'
    params:
        forwardAdapter = config['cutadapt']['forwardAdapter'],
        reverseAdapter = config['cutadapt']['reverseAdapter'],
        overlap = config['cutadapt']['overlap'],
        errorRate = config['cutadapt']['errorRate'],
        minimumLength = config['cutadapt']['minimumLength'],
        qualityCutoff = config['cutadapt']['qualityCutoff'],
        GCcontent = config['cutadapt']['GCcontent']
    log:
        'logs/cutadapt/{bound}.log'
    conda:
        f'{ENVS}/cutadapt.yaml'
    threads:
        config['threads']
    shell:
        cutadaptCmd()


rule modifyCutadapt:
    input:
        rules.cutadapt.output.qc
    output:
        'qc/cutadapt/{bound}.cutadapt.txt'
    group:
        'cutadapt'
    log:
        'logs/modifyCutadapt/{bound}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/modifyCutadapt.py {wildcards.sample} {input} '
        '> {output} 2> {log}'


if config['fastq_screen']:
    rule fastqScreen:
        input:
            'fastq/trimmed/{single}.trim.fastq.gz'
        output:
            txt = 'qc/fastq_screen/{single}.fastq_screen.txt',
            png = 'qc/fastq_screen/{single}.fastq_screen.png'
        params:
            fastq_screen_config = config['fastq_screen'],
            subset = 100000,
            aligner = 'bowtie2'
        log:
            'logs/fastq_screen/{single}.log'
        threads:
            config['threads']
        wrapper:
            '0.49.0/bio/fastq_screen'


rule fastQCtrimmed:
    input:
        'fastq/trimmed/{single}.trim.fastq.gz'
    output:
        html = 'qc/fastqc/{single}.trim_fastqc.html',
        zip = 'qc/fastqc/{single}.trim_fastqc.zip'
    log:
        'logs/fastQCtrimmed/{single}.log'
    wrapper:
        '0.49.0/bio/fastqc'


rule bowtie2Build:
    input:
        config['genome']['sequence']
    output:
        expand('genome/index/{build}.{n}.bt2',
               n=['1', '2', '3', '4', 'rev.1', 'rev.2'], build=BUILD)
    params:
        basename = f'genome/index/{BUILD}'
    log:
        'logs/bowtie2Build.log'
    conda:
        f'{ENVS}/bowtie2.yaml'
    threads:
        config['threads']
    shell:
        'bowtie2-build --threads {threads} {input} {params.basename} &> {log}'


def bowtie2Cmd():
    if config['paired']:
        return ('bowtie2 -x {params.index} -1 {input.reads_in[0]} '
            '-2 {input.reads_in[1]} --threads {threads} '
            '> {output.sam} 2> {log}; cp {log} {output.qc}')
    else:
        return ('bowtie2 -x {params.index} -U {input.reads[0]} '
            '--threads {threads} > {output.sam} 2> {log}; '
            'cp {log} {output.qc}')


rule bowtie2Map:
    input:
        reads = rules.cutadapt.output.trimmed,
        index = rules.bowtie2Build.output
    output:
        sam = pipe('mapped/{bound}.sam'),
        qc = 'qc/bowtie2/{bound}.bowtie2.txt'
    params:
        index = f'genome/index/{BUILD}'
    group:
        'map'
    log:
        'logs/bowtie2Map/{bound}.log'
    conda:
        f'{ENVS}/bowtie2.yaml'
    threads:
        min(1, (config['threads'] / 2) - 1)
    shell:
        bowtie2Cmd()


rule fixBAM:
    input:
        rules.hisat2.output.sam
    output:
        pipe('mapped/{bound}.fixmate.bam')
    log:
        'logs/fixBAM/{bound}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools fixmate -O bam,level=0 -m {input} {output} &> {log}'


rule sortBAM:
    input:
        rules.fixBAM.output
    output:
        'mapped/{bound}.sort.bam'
    log:
        'logs/sortBAM/{bound}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        min(1, (config['threads'] / 2) - 1)
    shell:
        'samtools sort -@ {threads} {input} > {output} 2> {log}'


rule markdupBAM:
    input:
        rules.sortBAM.output
    output:
        bam = 'mapped/{bound}.markdup.bam',
        qc = 'qc/deduplicate/{bound}.txt'
    log:
        'logs/markdupBAM/{bound}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        config['threads']
    shell:
        'samtools markdup -@ {threads} '
        '-sf {output.qc} {input} {output.bam} &> {log}'


rule indexBAM:
    input:
        'mapped/{bound}.{stage}.bam'
    output:
        'mapped/{bound}.{stage}.bam.bai'
    log:
        'logs/indexBAM/{bound}-{stage}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        min(1, (config['threads'] / 2) - 1,)
    shell:
        'samtools index -@ {threads} {input} &> {log}'

rule samtoolsStats:
    input:
        rules.markdupBAM.output.bam
    output:
        'qc/samtools/stats/{bound}.stats.txt'
    group:
        'samQC'
    log:
        'logs/samtools_stats/{bound}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools stats {input} > {output} 2> {log}'


rule samtoolsIdxstats:
    input:
        bam = 'mapped/{bound}.markdup.bam',
        index = 'mapped/{bound}.markdup.bam.bai'
    output:
        'qc/samtools/idxstats/{bound}.idxstats.txt'
    group:
        'samQC'
    log:
        'logs/samtools_idxstats/{bound}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools idxstats {input.bam} > {output} 2> {log}'


rule samtoolsFlagstat:
    input:
        rules.markdupBAM.output.bam
    output:
        'qc/samtools/flagstat/{bound}.flagstat.txt'
    group:
        'samQC'
    log:
        'logs/samtoolsFlagstat/{bound}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools flagstat {input} > {output} 2> {log}'



def setBlacklistCommand():
    if config['genome']['blacklist']:
        cmd = ('bedtools merge -d {params.distance} -i {input} '
               '> {output} 2> {log}')
    else:
        cmd = 'touch {output} &> {log}'
    return cmd


rule processBlacklist:
    input:
        config['genome']['blacklist'] if config['genome']['blacklist'] else []
    output:
        f'genome/{BUILD}-blacklist.bed'
    params:
        # Max distance between features allowed for features to be merged.
        distance = 1000
    log:
        f'logs/processBlacklist/{BUILD}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        setBlacklistCommand()


rule estimateReadFiltering:
    input:
        bam = rules.markdupBAM.output.bam,
        index = 'mapped/{bound}.markdup.bam.bai',
        blacklist = rules.processBlacklist.output
    output:
        'qc/deeptools/estimateReadFiltering/{bound}.txt'
    params:
        minMapQ = 15,
        binSize = 10000,
        distanceBetweenBins = 0,
        properPair = '--samFlagInclude 2' if config['paired'] else '',
    log:
        'logs/estimateReadFiltering/{bound}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'estimateReadFiltering --bamfiles {input.bam} --outFile {output} '
        '--binSize {params.binSize} --blackListFileName {input.blacklist} '
        '--distanceBetweenBins {params.distanceBetweenBins} '
        '--minMappingQuality {params.minMapQ} --ignoreDuplicates '
        '--samFlagExclude 260 {params.properPair} '
        '--numberOfProcessors {threads} &> {log}'


rule alignmentSieve:
    input:
        bam = rules.markdupBAM.output.bam,
        index = 'mapped/{sample_type}.markdup.bam.bai',
        blacklist = rules.processBlacklist.output
    output:
        bam = 'mapped/{sample_type}.filtered.bam',
        qc = 'qc/deeptools/{sample_type}-filter-metrics.txt'
    params:
        minMapQ = 15,
        maxFragmentLength = 2000,
        properPair = '--samFlagInclude 2' if config['paired'] else '',
    log:
        'logs/alignmentSieve/{sample_type}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'alignmentSieve --bam {input.bam} --outFile {output.bam} '
        '--minMappingQuality {params.minMapQ} --ignoreDuplicates '
        '--samFlagExclude 260 {params.properPair} '
        '--maxFragmentLength {params.maxFragmentLength} '
        '--blackListFileName {input.blacklist} '
        '--numberOfProcessors {threads} --filterMetrics {output.qc} &> {log}'


rule indexBAM:
    input:
        'mapped/{sample_type}.{stage}.bam'
    output:
        'mapped/{sample_type}.{stage}.bam.bai'
    log:
        'logs/indexBAM/{sample_type}-{stage}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        config['threads']
    shell:
        'samtools index -@ {threads} {input} &> {log}'


rule multiBamSummary:
    input:
        bams = expand('mapped/{sample}.filtered.bam', sample=SAMPLES),
        indexes = expand('mapped/{sample}.filtered.bam.bai', sample=SAMPLES)
    output:
        'qc/deeptools/multiBamSummary.npz'
    params:
        binSize = 10000,
        distanceBetweenBins = 0,
        labels = ' '.join(SAMPLES),
        extendReads = 150
    log:
        'logs/multiBamSummary.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'multiBamSummary bins --bamfiles {input.bams} --outFileName {output} '
        '--binSize {params.binSize} --labels {params.labels} '
        '--distanceBetweenBins {params.distanceBetweenBins} '
        '--extendReads {params.extendReads} --numberOfProcessors {threads} &> {log}'


rule plotCorrelation:
    input:
        rules.multiBamSummary.output
    output:
        plot = 'qc/deeptools/plotCorrelation.png',
        matrix = 'qc/deeptools/plotCorrelation.tsv'
    params:
        corMethod = 'pearson',
        colorMap = 'viridis',
        labels = ' '.join(SAMPLES)
    log:
        'logs/plotCorrelation.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    shell:
        'plotCorrelation --corData {input} --corMethod {params.corMethod} '
        '--whatToPlot heatmap --labels {params.labels} --skipZeros '
        '--colorMap {params.colorMap} --plotNumbers '
        '--plotFile {output.plot} --outFileCorMatrix {output.matrix} &> {log}'


def setColours(samples):
    """ Find group (and input) and assign to specific colour. """
    colours = ''
    colourPool = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
                  '#0072B2', '#D55E00', '#CC79A7']
    # Add double quotes for compatibility with shell
    colourPool = [f'"{colour}"' for colour in colourPool]
    usedColours = {}
    for sample in samples:
        if sample.endswith('-input'):
            # Inputs are always black
            colours += '"#000000" '
        else:
            group = sample.split('-')[0]
            if group not in usedColours:
                usedColours[group] = colourPool[0]
                # Remove from pool
                colourPool = colourPool[1:]
            colours += f'{usedColours[group]} '
    return f'{colours}'


rule plotPCA:
    input:
        rules.multiBamSummary.output
    output:
        plot = 'qc/deeptools/plotPCA.png',
        data = 'qc/deeptools/plotPCA.tab'
    params:
        labels = ' '.join(SAMPLES),
        colours = setColours(SAMPLES)
    log:
        'logs/plotPCA.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    shell:
        'plotPCA --corData {input} --colors {params.colours} '
        '--labels {params.labels} --plotFile {output.plot} '
        '--outFileNameData {output.data} &> {log}'


rule plotCoverage:
    input:
        bams = expand('mapped/{sample}.filtered.bam', sample=SAMPLES),
        indexes = expand('mapped/{sample}.filtered.bam.bai', sample=SAMPLES)
    output:
        plot = 'qc/deeptools/plotCoverage.png',
        data = 'qc/deeptools/plotCoverage.tab',
        info = 'qc/deeptools/plotCoverage.info'
    params:
        nSamples = 1000000,
        labels = ' '.join(SAMPLES)
    log:
        'logs/plotCoverage.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'plotCoverage --bamfiles {input.bams} --labels {params.labels} '
        '--plotFile {output.plot} --outRawCounts {output.data} '
        '--numberOfSamples {params.nSamples} --numberOfProcessors {threads} '
        '> {output.info} 2> {log}'


rule plotFingerprint:
    input:
        bams = expand('mapped/{sample}.filtered.bam', sample=SAMPLES),
        indexes = expand('mapped/{sample}.filtered.bam.bai', sample=SAMPLES)
    output:
        plot = 'qc/deeptools/plotFingerprint.png',
        data = 'qc/deeptools/plotFingerprint.tab'
    params:
        nSamples = 500000,
        labels = ' '.join(SAMPLES)
    log:
        'logs/plotFingerprint.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'plotFingerprint --bamfiles {input.bams} --labels {params.labels} '
        '--plotFile {output.plot} --outRawCounts {output.data} '
        '--numberOfSamples {params.nSamples} --skipZeros '
        '--numberOfProcessors {threads} &> {log}'


if INPUTS:
    rule mergeInput:
        input:
            expand('mapped/{input}-input.markdup.bam', input=INPUTS)
        output:
            'mapped/input/all-input.sort.bam'
        log:
            'logs/mergeInput.log'
        conda:
            f'{ENVS}/samtools.yaml'
        shell:
            'samtools merge {output} {input} &> {log}'


    rule indexInputBAM:
        input:
            rules.mergeInput.output
        output:
            f'{rules.mergeInput.output}.bai'
        log:
            'logs/indexInputBAM.log'
        conda:
            f'{ENVS}/samtools.yaml'
        threads:
            config['threads']
        shell:
            'samtools index -@ {threads} {input} &> {log}'


    rule bamCompare:
        input:
            treatment ='mapped/{bound}-bound.filtered.bam',
            treatmentIndex = 'mapped/{bound}-bound.filtered.bam.bai',
            control = rules.mergeInput.output,
            controlIndex = rules.indexInputBAM.output,
        output:
            'bigwig/compareInput/{bound}-Input.bigwig',
        params:
            binSize = 10,
            scale = 'SES',
            extendReads = 150,
            operation = 'log2',
            genomeSize = 2652783500,
        log:
            'logs/bamCompare/{bound}.log'
        conda:
            f'{ENVS}/deeptools.yaml'
        threads:
            config['threads']
        shell:
            'bamCompare --bamfile1 {input.treatment} --bamfile2 {input.control} '
            '--outFileName {output} --binSize {params.binSize} '
            '--extendReads {params.extendReads} --operation {params.operation} '
            '--effectiveGenomeSize {params.genomeSize} '
            '--numberOfProcessors {threads} &> {log}'


    rule computeMatrixScaled:
        input:
            expand('bigwig/compareInput/{bound}-Input.bigwig', bound=BOUNDS)
        output:
            scaledGZ = 'deeptools/computeMatrix/matrix-scaled.gz',
            scaled = 'deeptools/computeMatrix/matrix-scaled.tab',
            sortedRegions = 'deeptools/computeMatrix/genes-scaled.bed'
        params:
            binSize = 10,
            regionBodyLength = 5000,
            upstream = 2000,
            downstream = 2000,
            metagene = '--metagene' if config['computeMatrixScale']['exon'] else '',
            samplesLabel = ' '.join(BOUNDS),
            genes = config['genome']['genes'],
            averageType = 'mean'
        log:
            'logs/computeMatrixScaled.log'
        conda:
            f'{ENVS}/deeptools.yaml'
        threads:
            config['threads']
        shell:
            'computeMatrix scale-regions --scoreFileName {input} '
            '--regionsFileName {params.genes} --outFileName {output.scaledGZ} '
            '--outFileNameMatrix {output.scaled} --skipZeros '
            '--regionBodyLength {params.regionBodyLength} '
            '--samplesLabel {params.samplesLabel} --binSize {params.binSize} '
            '--averageTypeBins {params.averageType} {params.metagene} '
            '--upstream {params.upstream} --downstream {params.downstream} '
            '--outFileSortedRegions {output.sortedRegions} '
            '--numberOfProcessors {threads} &> {log}'


    rule computeMatrixReference:
        input:
            expand('bigwig/compareInput/{bound}-Input.bigwig', bound=BOUNDS)
        output:
            referenceGZ = 'deeptools/computeMatrix/matrix-reference.gz',
            reference = 'deeptools/computeMatrix/matrix-reference.tab',
            sortedRegions = 'deeptools/computeMatrix/genes-reference.bed'
        params:
            binSize = 10,
            upstream = 5000,
            downstream = 5000,
            samplesLabel = ' '.join(BOUNDS),
            genes = config['genome']['genes'],
            averageType = 'mean',
            referencePoint = 'TSS'
        log:
            'logs/computeMatrix.log'
        conda:
            f'{ENVS}/deeptools.yaml'
        threads:
            config['threads']
        shell:
            'computeMatrix reference-point --scoreFileName {input} '
            '--regionsFileName {params.genes} --outFileName {output.referenceGZ} '
            '--outFileNameMatrix {output.reference} --skipZeros '
            '--samplesLabel {params.samplesLabel} --binSize {params.binSize} '
            '--averageTypeBins {params.averageType} '
            '--referencePoint {params.referencePoint} '
            '--upstream {params.upstream} --downstream {params.downstream} '
            '--outFileSortedRegions {output.sortedRegions} '
            '--numberOfProcessors {threads} &> {log}'


    rule plotProfile:
        input:
            'deeptools/computeMatrix/matrix-{mode}.gz'
        output:
            plot = 'qc/deeptools/compareInput/plotProfileInput-{mode}.png',
            data = 'qc/deeptools/compareInput/plotProfileInput-data-{mode}.tab',
            bed = 'qc/deeptools/compareInput/plotProfileInput-regions-{mode}.bed'
        params:
            dpi = 300,
            plotsPerRow = 2,
            averageType = 'mean',
            referencePoint = 'TSS'
        log:
            'logs/plotProfile-{mode}.log'
        conda:
            f'{ENVS}/deeptools.yaml'
        shell:
            'plotProfile --matrixFile {input} --outFileName {output.plot} '
            '--outFileSortedRegions {output.bed} --outFileNameData {output.data} '
            '--dpi {params.dpi} --averageType {params.averageType} '
            '--refPointLabel {params.referencePoint} '
            '--numPlotsPerRow {params.plotsPerRow} &> {log}'


    rule plotHeatmap:
        input:
            'deeptools/computeMatrix/matrix-{mode}.gz'
        output:
            plot = 'qc/deeptools/compareInput/plotHeatmapInput-{mode}.png',
            data = 'qc/deeptools/compareInput/plotHeatmapInput-data-{mode}.tab',
            bed = 'qc/deeptools/compareInput/plotHeatmapInput-regions-{mode}.bed'
        params:
            dpi = 300,
            zMax = 3,
            zMin = -3,
            kmeans = 3,
            width = max(4, len(BOUNDS) * 2),
            colorMap = 'RdBu_r',
            averageType = 'mean',
            referencePoint = 'TSS',
            interpolationMethod = 'auto'
        log:
            'logs/plotHeatmap-{mode}.log'
        conda:
            f'{ENVS}/deeptools.yaml'
        threads:
            config['threads']
        shell:
            'plotHeatmap --matrixFile {input} --outFileName {output.plot} '
            '--outFileSortedRegions {output.bed} --outFileNameMatrix {output.data} '
            '--interpolationMethod {params.interpolationMethod} '
            '--dpi {params.dpi} --kmeans {params.kmeans} '
            '--zMin {params.zMin} --zMax {params.zMax} '
            '--colorMap {params.colorMap} --refPointLabel {params.referencePoint} '
            '--averageTypeSummaryPlot {params.averageType} '
            '--heatmapWidth {params.width} '
            '--refPointLabel {params.referencePoint} &> {log}'


rule mergeReplicates:
    input:
        lambda wc: expand('mapped/{group}-{rep}-bound.filtered.bam',
            group=wc.group, rep=GROUPS[wc.group])
    output:
        'mapped/{group}-bound.filtered.bam'
    log:
        'logs/mergeReplicates/{group}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools merge {output} {input} &> {log}'


rule indexMergedReplicates:
    input:
        rules.mergeReplicates.output
    output:
        f'{rules.mergeReplicates.output}.bai'
    log:
        'logs/indexMergedReplicates/{group}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        config['threads']
    shell:
        'samtools index -@ {threads} {input} &> {log}'


rule bamCompareGroups:
    input:
        group1 = 'mapped/{group1}-bound.filtered.bam',
        group1Index = 'mapped/{group1}-bound.filtered.bam.bai',
        group2 = 'mapped/{group2}-bound.filtered.bam',
        group2Index = 'mapped/{group2}-bound.filtered.bam.bai',
    output:
        'bigwig/compareGroup/{group1}-{group2}.bigwig',
    params:
        binSize = 10,
        scale = 'SES',
        extendReads = 150,
        operation = 'log2',
        genomeSize = 2652783500,
    log:
        'logs/bamCompareGroups/{group1}-{group2}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'bamCompare --bamfile1 {input.group1} --bamfile2 {input.group2} '
        '--outFileName {output} --binSize {params.binSize} '
        '--extendReads {params.extendReads} --operation {params.operation} '
        '--effectiveGenomeSize {params.genomeSize} '
        '--numberOfProcessors {threads} &> {log}'


rule computeMatrixScaledGroups:
    input:
        expand('bigwig/compareGroup/{compare}.bigwig', compare=COMPARES)
    output:
        scaledGZ = 'deeptools/computeMatrixGroups/matrix-scaled.gz',
        scaled = 'deeptools/computeMatrixGroups/matrix-scaled.tab',
        sortedRegions = 'deeptools/computeMatrixGroups/genes-scaled.bed'
    params:
        binSize = 10,
        regionBodyLength = 5000,
        upstream = 2000,
        downstream = 2000,
        metagene = '--metagene' if config['computeMatrixScale']['exon'] else '',
        samplesLabel = ' '.join(COMPARES),
        genes = config['genome']['genes'],
        averageType = 'mean'
    log:
        'logs/computeMatrixScaledGroups.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'computeMatrix scale-regions --scoreFileName {input} '
        '--regionsFileName {params.genes} --outFileName {output.scaledGZ} '
        '--outFileNameMatrix {output.scaled} --skipZeros '
        '--regionBodyLength {params.regionBodyLength} '
        '--samplesLabel {params.samplesLabel} --binSize {params.binSize} '
        '--averageTypeBins {params.averageType} {params.metagene} '
        '--upstream {params.upstream} --downstream {params.downstream} '
        '--outFileSortedRegions {output.sortedRegions} '
        '--numberOfProcessors {threads} &> {log}'


rule computeMatrixReferenceGroups:
    input:
        expand('bigwig/compareGroup/{compare}.bigwig', compare=COMPARES)
    output:
        referenceGZ = 'deeptools/computeMatrixGroups/matrix-reference.gz',
        reference = 'deeptools/computeMatrixGroups/matrix-reference.tab',
        sortedRegions = 'deeptools/computeMatrixGroups/genes-reference.bed'
    params:
        binSize = 10,
        upstream = 5000,
        downstream = 5000,
        samplesLabel = ' '.join(COMPARES),
        genes = config['genome']['genes'],
        averageType = 'mean',
        referencePoint = 'TSS'
    log:
        'logs/computeMatrixReferenceGroups.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'computeMatrix reference-point --scoreFileName {input} '
        '--regionsFileName {params.genes} --outFileName {output.referenceGZ} '
        '--outFileNameMatrix {output.reference} --skipZeros '
        '--samplesLabel {params.samplesLabel} --binSize {params.binSize} '
        '--averageTypeBins {params.averageType} '
        '--referencePoint {params.referencePoint} '
        '--upstream {params.upstream} --downstream {params.downstream} '
        '--outFileSortedRegions {output.sortedRegions} '
        '--numberOfProcessors {threads} &> {log}'

rule plotProfileGroups:
    input:
        'deeptools/computeMatrixGroups/matrix-{mode}.gz'
    output:
        plot = 'qc/deeptools/compareGroup/plotProfileGroup-{mode}.png',
        data = 'qc/deeptools/compareGroup/plotProfileGroup-data-{mode}.tab',
        bed = 'qc/deeptools/compareGroup/plotProfileGroup-regions-{mode}.bed'
    params:
        dpi = 300,
        plotsPerRow = 2,
        averageType = 'mean',
        referencePoint = 'TSS'
    log:
        'logs/plotProfileGroups-{mode}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    shell:
        'plotProfile --matrixFile {input} --outFileName {output.plot} '
        '--outFileSortedRegions {output.bed} --outFileNameData {output.data} '
        '--dpi {params.dpi} --averageType {params.averageType} '
        '--refPointLabel {params.referencePoint} '
        '--numPlotsPerRow {params.plotsPerRow} &> {log}'


rule plotHeatmapGroups:
    input:
        'deeptools/computeMatrixGroups/matrix-{mode}.gz'
    output:
        plot = 'qc/deeptools/compareGroup/plotHeatmapGroup-{mode}.png',
        data = 'qc/deeptools/compareGroup/plotHeatmapGroup-data-{mode}.tab',
        bed = 'qc/deeptools/compareGroup/plotHeatmapGroup-regions-{mode}.bed'
    params:
        dpi = 300,
        zMax = 3,
        zMin = -3,
        kmeans = 3,
        width = max(4, len(BOUNDS) * 2),
        colorMap = 'RdBu_r',
        averageType = 'mean',
        referencePoint = 'TSS',
        interpolationMethod = 'auto'
    log:
        'logs/plotHeatmapGroups-{mode}.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'plotHeatmap --matrixFile {input} --outFileName {output.plot} '
        '--outFileSortedRegions {output.bed} --outFileNameMatrix {output.data} '
        '--interpolationMethod {params.interpolationMethod} '
        '--dpi {params.dpi} --kmeans {params.kmeans} '
        '--zMin {params.zMin} --zMax {params.zMax} '
        '--colorMap {params.colorMap} --refPointLabel {params.referencePoint} '
        '--averageTypeSummaryPlot {params.averageType} '
        '--heatmapWidth {params.width} '
        '--refPointLabel {params.referencePoint} &> {log}'


def macs2Control(wc):
    if INPUTS:
        return rules.mergeInput.output
    else:
        return []

def macs2Command():
    if INPUTS:
        return ('macs2 callpeak '
            '--treatment {input.bound} '
        	'--control {input.input} '
         	'--format BAM --gsize {params.genomeSize} '
        	'--name {wildcards.sample} '
        	'--outdir {params.dir} &> {log}')
    else:
        return ('macs2 callpeak '
            '--treatment {input.bound} '
         	'--format BAM --gsize {params.genomeSize} '
        	'--name {wildcards.sample} '
        	'--outdir {params.dir} &> {log}')

rule macs2:
    input:
        input = macs2Control,
        bound = 'mapped/{bound}-bound.filtered.bam'
    output:
        summits = 'macs2/{bound}/{bound}_summits.bed',
        narrowPeak = 'macs2/{bound}/{bound}_peaks.narrowPeak',
        xlsPeak = 'macs2/{bound}/{bound}_peaks.xls',
        model = 'macs2/{bound}/{bound}_model.r'
    params:
        dir = directory('macs2/{bound}'),
        genomeSize = 2652783500
    log:
        'logs/macs/{bound}.log'
    conda:
        f'{ENVS}/macs2.yaml'
    shell:
        macs2Command()


rule consensusPeaks:
    input:
        lambda wc: expand(
            'macs2/{group}-{rep}/{group}-{rep}_peaks.narrowPeak',
            group=wc.group, rep=GROUPS[wc.group])
    output:
        'macs2/{group}-consensus_peaks.narrowPeak'
    params:
        # Mininium replicates required to retain peak
        min = lambda wc: math.ceil(len(GROUPS[wc.group]) / 2)
    log:
        'logs/consensusPeaks/{group}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools multiinter -i {input} '
        "| awk -v n={params.min} '{{if ($4 >= n) {{print}}}}' "
        "> {output} 2> {log}"


rule intersectConsensus:
    input:
        expand('macs2/{group}-consensus_peaks.narrowPeak',
            group=list(GROUPS)),
    output:
        'macs2/consensus/all-consensus_peaks.narrowPeak'
    log:
        f'logs/intersectConsensus.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools multiinter -i {input} > {output} 2> {log}'


rule plotEnrichment:
    input:
        bams = expand('mapped/{sample}.filtered.bam',
            sample=SAMPLES),
        indexes = expand('mapped/{sample}.filtered.bam.bai',
            sample=SAMPLES),
        peaks = rules.intersectConsensus.output
    output:
        plot = 'qc/deeptools/plotEnrichment.png',
        data = 'qc/deeptools/plotEnrichment.tab'
    params:
        plotsPerRow = 2,
        labels = ' '.join(SAMPLES),
        genes = config['genome']['genes'],
        regionLabel = '"MACS2 consensus peaks" "Gene features"'
    log:
        'logs/plotEnrichment.log'
    conda:
        f'{ENVS}/deeptools.yaml'
    threads:
        config['threads']
    shell:
        'plotEnrichment --bamfiles {input.bams} '
        '--BED {input.peaks} {params.genes} '
        '--plotFile {output.plot} --outRawCounts {output.data} '
        '--labels {params.labels} --numPlotsPerRow {params.plotsPerRow} '
        '--regionLabels {params.regionLabel} --numberOfProcessors {threads} '
        '&> {log}'


rule multiQC:
    input:
        expand('qc/fastqc/{single}.raw_fastqc.zip', single= samples['single']),
        expand('qc/fastq_screen/{single}.fastq_screen.txt',
            single=samples['single']) if config['fastq_screen'] else [],
        expand('qc/cutadapt/{single}.cutadapt.txt', single=samples['single']),
        expand('qc/fastqc/{single}.trim_fastqc.zip', single=samples['single']),
        expand('qc/bowtie2/{sample}.bowtie2.txt', sample=SAMPLESs),
        expand('macs2/{bound}/{bound}_peaks.xls', sample=BOUNDS),
        expand('qc/deeptools/estimateReadFiltering/{sample}.txt',
            sample=SAMPLES_),
        'qc/deeptools/plotCorrelation.tsv',
        'qc/deeptools/plotPCA.tab',
        'qc/deeptools/plotCoverage.info',
        'qc/deeptools/plotCoverage.tab',
        'qc/deeptools/plotFingerprint.png',
        'qc/deeptools/plotFingerprint.tab',
        'qc/deeptools/compareInput/plotProfileInput-data-scaled.tab',
        'qc/deeptools/compareGroup/plotProfileGroup-data-scaled.tab',
        'qc/deeptools/plotEnrichment.tab',
        expand('qc/samtools/stats/{sample}.stats.txt', sample=SAMPLES),
        expand('qc/samtools/idxstats/{sample}.idxstats.txt', sample=SAMPLES),
        expand('qc/samtools/flagstat/{sample}.flagstat.txt', sample=SAMPLES)
    output:
        directory('qc/multiqc')
    log:
        'logs/multiqc.log'
    conda:
        f'{ENVS}/multiqc.yaml'
    shell:
        'multiqc --outdir {output} '
            '--force --config {BASE}/config/multiqc_config.yaml {input} '
        '&> {log}'
