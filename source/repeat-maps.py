#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Identifying RBP enrichment at repeat families and other multi-copy elements.

Designed as an orthogonal approach to peak calling, this pipeline maps trimmed fastq files to a set of manually
curated sequences belonging to distinct repeat families, including ribosomal RNAs (e.g. 18S, 28S, 5S, 5.8S, and
the 45S precursor), snRNAs (e.g. U1, U2), snoRNAs, tRNAs, YRNAs, as well as repetitive elements (e.g. Alu, LINE,
ERV). It then merges these repeat-mapped reads with genome-mapped reads and performs its own PCR-deduplication
step, resulting in a table of enriched binding sites within repeat families or unique genomic elements.
"""

import os
import argparse
import glob
import os
import gzip
import sys
import time
import subprocess

import datetime
import tempfile
import shutil
import itertools
import math

import pysam
import cmder
from seqflow import Flow, task, logger


def file_path(p):
    ps = p.split(':') if ':' in p else [p]
    for path in ps:
        if not os.path.isfile(path):
            logger.error(f'FASTQ file "{path}" may not be a file or does not exist.')
            sys.exit(1)
    return p


parser = argparse.ArgumentParser(description=__doc__.splitlines()[1], prog='repeat-maps')
parser.add_argument('--fastq', required=True, nargs='+', type=file_path,
                    help='Path to adapters trimmed fastq file(s). For PE dataset, fastq files for R1 and R2 need to '
                         'be separated via ":", e.g., fastq1:fastq2.')
parser.add_argument('--bam', required=True, nargs='+', type=file_path,
                    help='path to genomic mapping bam file(s) corresponding to fastq file(s).')
parser.add_argument('--dataset', nargs='+', help='Name of the each dataset.')
parser.add_argument('--species', default='hg19', choices=('hg19', 'hg38'),
                    help="Short name for species, e.g. hg19, hg38, mm10, ..., default: hg19 .")
parser.add_argument('--database', help='The basename of the bowtie index files, default: auto populate for species.')
parser.add_argument('--outdir',
                    help="Path to the output directory. Default to the current work directory and if the specified "
                         "path does not exist, it will try to create it first.",
                    default=os.getcwd())
parser.add_argument('--cpus', type=int, help='Maximum number of CPUs can be used, default: 8.', default=8)
parser.add_argument('--time', type=int, help='Time (in integer hours) for running your job.', default=4)
parser.add_argument('--memory', type=int, help='Amount of memory (in GB) for all CPU cores.', default=32)
parser.add_argument('--email', help='Email address for notifying you the start, end, and abort of you job.')
parser.add_argument('--scheduler', choices=('pbs', 'qsub', 'slurm', 'sbatch', 'local'), default='local',
                    help='Name (case insensitive) of the scheduler for job, default: local.')
parser.add_argument('--job', help="Name of your job", default='repeat.map')
parser.add_argument('--hold_submit', action='store_true',
                    help='Generate the submit script but hold it without submitting to the job scheduler.')
parser.add_argument('--debug', action='store_true', help='Invoke debug mode (for development use only).')
parser.add_argument('--dry_run', action='store_true',
                    help='Print out steps and files involved in each step without actually running the pipeline.')


options = parser.parse_args()
if not os.path.isdir(options.outdir):
    try:
        os.mkdir(options.outdir)
    except OSError as e:
        logger.error(e)
        sys.exit(1)
os.chdir(options.outdir)

if len(options.fastq) != len(options.bam):
    logger.error(f'Number of FASTQ files does not match the number of BAM files.')
    sys.exit(1)
    
F1_TO_FASTQ = {fastq.split(':')[0]: fastq if ':' in fastq else f'{fastq}:' for fastq in options.fastq}
datasets = options.dataset if options.dataset else [os.path.basename(fastq).replace('.fastq.gz', '')
                                                    for fastq in F1_TO_FASTQ]
if len(datasets) != len(options.fastq):
    logger.error(f'Number of FASTQ files does not match the number of dataset names.')
    sys.exit(1)

FASTQ_TO_SAM = {fastq: f'{ds}.repetitive.elements.sam' for fastq, ds in zip(F1_TO_FASTQ, datasets)}
DATASET_TO_SAM = {dataset: (f'{dataset}.repetitive.elements.sam', bam)
                  for dataset, fastq, bam in zip(datasets, F1_TO_FASTQ, options.bam)}
DATA_TYPES = {dataset: 'PE' if ':' in fastq else 'SE' for dataset, fastq in zip(datasets, options.fastq)}
SAM_TO_BAM = {f'{dataset}.repetitive.elements.{x}{y}.sam': os.path.basename(bam).replace('.bam', f'.genome.{x}{y}.sam')
              for dataset in datasets for x in 'ATGCN' for y in 'ATGCN' for bam in options.bam}
BAM_TO_SAM = {v: k for k, v in SAM_TO_BAM.items()}
GENOME_BAM = {os.path.basename(bam).replace('.bam', f'.genome.{x}{y}.sam'): bam
              for bam in options.bam for x in 'ATGCN' for y in 'ATGCN'}
dbs = {'hg19': '/storage/vannostrand/reference_data/repeat-map-hg19/bowtie2_index/repeat.elements',
       'hg38': '/storage/vannostrand/reference_data/repeat-family-mapping-grch38/bowtie2_index/repeat.elements'}
if options.species not in dbs:
    logger.error(f'Invalid species or species has not been implemented yest.')
    sys.exit(1)
DB = options.database if options.database else dbs[options.species]
setattr(options, 'database', DB)


@task(inputs=options.fastq, outputs=lambda i: FASTQ_TO_SAM[i])
def bowtie_family_map(fastq, sam):
    fastq1, fastq2 = F1_TO_FASTQ[fastq].split(':')
    fastq = f'{fastq1}:{fastq2}' if fastq2 else fastq1
    cmd = f'bowtie_family_map.pl {fastq} {DB} {sam} {options.species} {options.cpus}'
    cmder.run(cmd, msg='Mapping and annotating repeat elements ...')
    
    
def split_bam(bam, flag='SE', out=''):
    exes = {'SE': 'split_bam_se.pl', 'PE': 'split_bam_se.pl'}
    out = f' {out}' if out else ''
    cmder.run(f'{exes[flag]} {bam}{out}', msg=f'Splitting {bam} ...')
    
    
@task(inputs=bowtie_family_map, outputs=lambda i: i.replace('.sam', '.NN.sam'))
def split_repeat_bam(sam, out):
    split_bam(sam, flag=DATA_TYPES[out.replace('.repetitive.elements.NN.sam', '')])


@task(inputs=[], outputs=list(GENOME_BAM.keys())[:1], parent=split_repeat_bam)
def split_genome_bam(inputs, out):
    bam = GENOME_BAM[out]
    flag = cmder.run(f'samtools view -c -f 1 {bam}').stdout.read().strip()
    flag = 'SE' if flag == '0' else 'PE'
    split_bam(bam, flag=flag, out=os.path.basename(bam).replace('.bam', '.genome'))


@task(inputs=[], parent=split_genome_bam, cpus=options.cpus,
      outputs=[k.replace('.repetitive.elements', '.repetitive.elements.combine.with.unique.map.dedup')
               for k in SAM_TO_BAM])
def dedup(bam, out):
    sam = out.replace('.combine.with.unique.map.dedup', '')
    flag = DATA_TYPES[sam.split('.repetitive.elements.')[0]]
    bam = SAM_TO_BAM[sam]
    cmder.run(f'duplicate_removal.pl {sam} {bam} {flag} {options.species}',
              msg=f'Deduplicating {sam} ...')


@task(inputs=[], parent=dedup,
      outputs=[f'{dataset}.repetitive.elements.combine.with.unique.map.count.tsv' for dataset in datasets])
def merge_result(inputs, out):
    pattern = out.replace('.count.tsv', '.dedup*sam')
    tsv = out.replace('.count.tsv', '.tsv.gz')
    cmder.run(f'cat {pattern} | pigz -p {options.cpus} > {tsv}')
    cmder.run(f'merge_sam.pl {out} *.count.txt')


def cleanup():
    cmder.run(f'rm *.repetitive.elements.sam')
    cmder.run('rm *.[ATGCN][ATGCN].sam *.[ATGCN][ATGCN].count.txt* failed_jobs*')


def schedule():
    sbatch = """#!/usr/bin/env bash

#SBATCH -n {cpus}                       # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t {runtime}                  # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=32G                   # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name={job}            # Short name for the job
#SBATCH --output=%j.{job}.log
"""
    sbatch_email = """
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL
"""
    pbs = """ #!/usr/bin/env bash

#PBS -l nodes=1:ppn={cores}
#PBS -l walltime={runtime}
#PBS -l vmem={memory}gb
#PBS -j oe
#PBS -N {jobname}
"""
    pbs_email = """
#PBS -M {email}
#PBS -m abe
"""
    
    if options.scheduler in ('pbs', 'qsub'):
        runtime, directive, exe, mail = f'{options.time}:00:00', pbs, 'qsub', pbs_email
        project = '/project/vannostrand'
    elif options.scheduler in ('slurm', 'sbatch'):
        days, hours = divmod(options.time, 24)
        runtime, directive, exe, mail = f'{days}-{hours:02}:00', sbatch, 'sbatch', sbatch_email
        project = '/storage/vannostrand'
    else:
        raise ValueError(f'Unsupported scheduler: {options.scheduler}, see help for supported schedulers.')
    setattr(options, 'scheduler', 'local')
    if options.debug:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    else:
        root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    program = os.path.join(root, 'repeat-maps')
    setattr(options, 'runtime', runtime)
    setattr(options, 'project', project)
    
    arguments = ('fastq', 'bam', 'dataset', 'database', 'species', 'outdir', 'cpus', 'scheduler')
    d = vars(options)
    args = []
    for k in arguments:
        v = d[k]
        if k == 'dataset':
            v = ' '.join(v)
        elif isinstance(v, (list, tuple)):
            separator = ' \\\n  ' + ' ' * (len(k) + 3)
            v = separator.join(v)
        args.append(f'--{k} {v}')
    arguments = ' \\\n  '.join(args)
    # arguments = ' \\\n  '.join([f'--{k} {" ".join(d[k]) if isinstance(d[k], (list, tuple)) else d[k]}'
    #                             for k in arguments if d.get(k, None)])
    setattr(options, 'arguments', arguments)
    
    code = fr"""
export TMPDIR={project}/tmp
export TEMP={project}/tmp
export TMP={project}/tmp

{program} \
  {arguments}
"""
    
    text = [directive, mail, code] if options.email else [directive, code]
    text = ''.join(text).format(**vars(options))
    
    submitter = os.path.join(options.outdir, 'repeat.map.submit.sh')
    with open(submitter, 'w') as o:
        o.write(text)
    
    print(f'Job submit script was saved to:\n    {submitter}')
    if options.hold_submit:
        print(f'Job {options.job} was not submitted yet, submit it after carefully review the submit script using:')
        print(f'    {exe} {submitter}')
    else:
        subprocess.run([exe, submitter])
        print(f'Job {options.job} was successfully submitted with the following resources:')
        data = {'Job name:': options.job, 'Output directory:': options.outdir,
                'Number of cores:': options.cpus, 'Job memory:': options.memory,
                'Job runtime:': f'{runtime} (D-HH:MM)'}
        for k, v in data.items():
            print(f'{k:>20} {v}')


@logger.catch()
def main():
    if options.scheduler == 'local':
        flow = Flow('repeat-maps', description=__doc__.strip())
        flow.run(dry_run=options.dry_run, cpus=options.cpus)
        logger.info('Cleaning up and finalizing ...')
        cleanup()
        logger.info('Mission accomplished!')
    else:
        schedule()


if __name__ == '__main__':
    main()
