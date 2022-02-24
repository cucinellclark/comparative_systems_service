#!/usr/bin/env python

import copy
import gzip
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request as request
from contextlib import closing
from multiprocessing import Process

import requests

from bvbrc_api import authenticateByEnv, getHostManifest



def gzipMove(source, dest):
    with open(source, "rb") as f_in:
        with gzip.open(dest, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def moveRead(filepath):
    directory, base = os.path.split(filepath)
    new_base = base.replace(")", "_").replace("(", "_")
    if new_base != base:
        newpath = os.path.join(directory, new_base)
        os.rename(filepath, newpath)
        return newpath
    else:
        return filepath


def run_compare_systems(job_data, output_dir, tool_params={}):
    # arguments:
    # list of genomes [{"genome":somefile,"annotation":somefile}]
    # dictionary of library dictionaries structured as {libraryname:{library:libraryname, replicates:[{read1:read1file, read2:read2file}]}}
    # parametrs_file is a json keyed parameters list.
    # Example tool_params: '{"fastqc":{"-p":"2"},"trim_galore":{"-p":"2"},"bowtie2":{"-p":"2"},"hisat2":{"-p":"2"},"samtools_view":{"-p":"2"},"samtools_index":{"-p":"2"}}'
    output_dir = os.path.abspath(output_dir)
    subprocess.call(["mkdir", "-p", output_dir])
    genome_list, read_list, recipe = setup(job_data, output_dir, tool_params)
    # print("genome_list: {}\nread_list: {}\nrecipe: {}".format(genome_list, read_list, recipe), file=sys.stdout)
    # sys.stdout.flush()
    for step in recipe:
        step = step.upper()
        if step == "PATHWAYS":
            run_trim(genome_list, output_dir, job_data, tool_params)
        elif step == "SUBSYSTEMS":
            run_subsystems(genome_list, output_dir, job_data, tool_params)
        elif step == "FAMILIES":
            run_families(read_list, output_dir, job_data, tool_params)
        else:
            print("Skipping step. Not found: {}".format(step), file=sys.stderr)
