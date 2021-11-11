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

import pandas as pd
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
                    help="Short name for species, e.g. hg19, hg38, mm10, ..., default: hg19.")
parser.add_argument('--database', help='The basename of the bowtie index files, default: auto populate for species.')
parser.add_argument('--outdir',
                    help="Path to the output directory. Default to the current work directory and if the specified "
                         "path does not exist, it will try to create it first.",
                    default=os.getcwd())
parser.add_argument('--summary_html', help='Path to html file that summary will be saved.',
                    default='repetitive.elements.enrichment.summary.html')
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
    
R1 = [fastq.split(':')[0] for fastq in options.fastq]
R1_TO_FASTQ = {fastq.split(':')[0]: fastq if ':' in fastq else f'{fastq}:' for fastq in options.fastq}
DATASETS = options.dataset if options.dataset else [os.path.basename(fastq).replace('.fastq.gz', '')
                                                    for fastq in R1_TO_FASTQ]
if len(DATASETS) != len(options.fastq):
    logger.error(f'Number of FASTQ files does not match the number of dataset names.')
    sys.exit(1)
    
FASTQ_TO_LINK = {r1: f'{dataset}.r1.fastq.gz' for r1, dataset in zip(R1_TO_FASTQ, DATASETS)}
FASTQ_TO_BAM = {r1: bam for r1, bam in zip(R1_TO_FASTQ, options.bam)}
DATASET_TO_TYPES = {dataset: 'PE' if ':' in fastq else 'SE' for dataset, fastq in zip(DATASETS, options.fastq)}
DBS = {'hg19': '/storage/vannostrand/reference_data/repeat-map-hg19/bowtie2_index/repeat.elements',
       'hg38': '/storage/vannostrand/reference_data/repeat-family-mapping-grch38/bowtie2_index/repeat.elements'}
DB = options.database if options.database else DBS[options.species]
setattr(options, 'database', DB)


def _link(src, dst):
    if src:
        if os.path.exists(src):
            if os.path.exists(dst):
                if os.path.islink(dst):
                    logger.debug(f'Link {dst} for source {src} already exists.')
                else:
                    logger.warning(f'Both link {dst} and source {src} point to the same file.')
            else:
                os.symlink(src, dst)
        else:
            logger.error(f'File {src} does not exist.')
            sys.exit(1)
        

@task(inputs=list(FASTQ_TO_LINK.keys()), outputs=lambda i: FASTQ_TO_LINK[i])
def soft_link(fastq, link):
    (fastq1, fastq2), bam = R1_TO_FASTQ[fastq].split(':'), FASTQ_TO_BAM[fastq]
    dataset = link.replace('.r1.fastq.gz', '')
    _link(fastq1, os.path.join(options.outdir, link))
    _link(fastq2, os.path.join(options.outdir, f'{dataset}.r2.fastq.gz'))
    _link(bam, os.path.join(options.outdir, f'{dataset}.genome.bam'))
    return link


@task(inputs=soft_link, outputs=lambda i: i.replace('.r1.fastq.gz', '.repetitive.elements.sam'))
def bowtie_family_map(fastq, sam):
    fastq2 = fastq.replace('.r1.fastq.gz', '.r2.fastq.gz')
    if os.path.exists(fastq2):
        fastq = f'{fastq}:{fastq2}'
    cmd = f'bowtie_family_map.pl {fastq} {DB} {sam} {options.species} {options.cpus}'
    cmder.run(cmd)
    
    
def split_bam(bam, flag='SE', out=''):
    exes = {'SE': 'split_bam_se.pl', 'PE': 'split_bam_se.pl'}
    out = f' {out}' if out else ''
    cmder.run(f'{exes[flag]} {bam}{out}')
    
    
@task(inputs=bowtie_family_map, outputs=lambda i: i.replace('.sam', '.NN.sam'))
def split_repeat_bam(sam, out):
    split_bam(sam, flag=DATASET_TO_TYPES[out.replace('.repetitive.elements.NN.sam', '')])


@task(inputs=[], outputs=[f'{dataset}.genome.NN.sam' for dataset in DATASETS],
      parent=split_repeat_bam, cpus=options.cpus)
def split_genome_bam(inputs, out):
    bam = out.replace('.genome.NN.sam', '.genome.bam')
    flag = cmder.run(f'samtools view -c -f 1 {bam}').stdout.read().strip()
    flag = 'SE' if flag == '0' else 'PE'
    split_bam(bam, flag=flag, out=out.replace('.NN.sam', ''))


@task(inputs=[], parent=split_genome_bam, cpus=options.cpus,
      outputs=[f'{dataset}.repetitive.elements.combine.with.unique.map.dedup.{x}{y}.sam'
               for dataset in DATASETS for x in 'ATGCN' for y in 'ATGCN'])
def dedup(bam, out):
    sam = out.replace('.repetitive.elements.combine.with.unique.map.dedup.', '.repetitive.elements.')
    bam = out.replace('.repetitive.elements.combine.with.unique.map.dedup.', '.genome.')
    flag = DATASET_TO_TYPES[out.split('.repetitive.elements.')[0]]
    cmder.run(f'duplicate_removal.pl {sam} {bam} {flag} {options.species}')


@task(inputs=[], parent=dedup,
      outputs=[f'{dataset}.repetitive.elements.combine.with.unique.map.count.tsv.gz' for dataset in DATASETS])
def merge_result(inputs, out):
    pattern = out.replace('.count.tsv.gz', '.dedup*sam')
    map_out = out.replace('.count.tsv.gz', '.tsv.gz')
    cmder.run(f'cat {pattern} | pigz -p {options.cpus} > {map_out}')
    tsv = out.replace('.tsv.gz', '.tsv')
    cmder.run(f'merge_count.pl {tsv} *.[ATGCN][ATGCN].count.txt')
    cmder.run(f'pigz -p {options.cpus} {tsv}')
    

TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.1/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https:////cdn.datatables.net/1.11.1/css/jquery.dataTables.min.css" rel="stylesheet">

    <title>Repetitive Elements Enrichment Summary</title>
</head>

<body>
    <div class="container">
        <h1 class="py-3">RBP enrichment at repeat families and other multi-copy elements</h1>
        {table}
    </div>
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.1/dist/js/bootstrap.bundle.min.js"></script>
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
    <script src="https://cdn.datatables.net/1.11.1/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/fixedheader/3.2.0/js/dataTables.fixedHeader.min.js"></script>
    <script src="https://cdn.datatables.net/fixedcolumns/4.0.0/js/dataTables.fixedColumns.min.js"></script>
    <script>
        $(document).ready( function () {{
            $('table.dataframe').DataTable({{
                paging: false,
                ordering: false,
                fixedHeader: true,
                scrollX: true,
                scrollY: 400,
                scrollCollapse: true,
                fixedColumns: {{left: 1}}
            }});
        }} );
    </script>
</body>
</html>
"""


@task(inputs=[], outputs=[options.summary_html], parent=merge_result)
def summarize(inputs, out):
    outputs = [f'{dataset}.repetitive.elements.combine.with.unique.map.count.tsv.gz' for dataset in DATASETS]
    df = pd.DataFrame()
    for path in outputs:
        sample = os.path.basename(path).split('.trim.repetitive.')[0]
        dd = pd.read_csv(path, sep='\t', comment='#', header=None, usecols=[0, 1, 3],
                         names=['group', 'element', 'percentage'])
        dd = dd.rename(columns={'percentage': sample})
        dd = dd[(dd['group'] == 'TOTAL') & (~dd['element'].str.contains('\|'))]
        dd = dd.drop(columns=['group'])
        df = dd if df.empty else pd.merge(df, dd, on='element', how='outer')
    df = df.fillna(0)
    df[df.select_dtypes(include=['number']).columns] *= 3
    df.to_csv(out.replace('.html', '.tsv'), index=False, sep='\t', float_format='%.6f')
    with open(out, 'w') as o:
        table = df.to_html(index=False, escape=False, float_format='%.4f', justify='center',
                           classes='dataframe table table-light table-bordered text-center w-100')
        o.write(TEMPLATE.format(table=table))


@task(inputs=[], parent=summarize, outputs=['repetitive.elements.multi.map.delete.tsv.gz'])
def cleanup(tsv, out):
    cmder.run(f'rm *.repetitive.elements.sam || true')
    cmder.run('rm *.[ATGCN][ATGCN].sam *.[ATGCN][ATGCN].count.txt* failed_jobs* || true')
    _ = [os.unlink(p) for p in glob.iglob('*.fastq.gz') if os.path.islink(p)]
    _ = [os.unlink(p) for p in glob.iglob('*.bam') if os.path.islink(p)]
    cmder.run(f'tail -n +1 *repetitive.elements.bowtie.log > repetitive.elements.bowtie.map.log')
    cmder.run(f'rm *repetitive.elements.bowtie.log')
    cmder.run(f'tail -n +1 *.map.delete.tsv | pigz -p {options.cpus} > {out}')
    cmder.run(f'rm *.map.delete.tsv')


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
    else:
        schedule()


if __name__ == '__main__':
    main()
