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

from bvbrc_api import authenticateByEnv,getGenomeGroupIds 

def run_pathways(genome_ids,output_file,output_dir, session):
    
    ### queries
    base_query = "https://www.patricbrc.org/api/pathway/?in(genome_id,("
    end_query = "))&eq(annotation,PATRIC)&limit(200000000)&http_accept=text/tsv" 
    query = base_query + ",".join(genome_ids) + end_query
    print("Pathways Query:\n{0}".format(query))
    req = session.get(query)
    pathways_file = os.path.join(output_dir,output_file+"_pathways.tsv")
    with open(pathways_file,"w") as o:
        o.write(req.text)

def get_genome_group_ids(group_list,s):
    genome_group_ids = []
    for genome_group in group_list:
        list_text = getGenomeGroupIds(genome_group,s,genomeGroupPath=True).strip('][').split(',')
        # ['{"genome_id":"562.79202"}', '{"genome_id":"562.80445"}', '{"genome_id":"562.80446"}']
        for genome_entry in list_text:
            genome_entry = genome_entry.strip('}{').split(":")
            genome_group_ids.append(genome_entry[1].strip('\"'))
    return genome_group_ids

def run_compare_systems(job_data, output_dir):

    ###Setup session
    s = requests.Session()
    authenticateByEnv(s)

    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        subprocess.call(["mkdir", "-p", output_dir])
    output_file = job_data["output_file"]


    print("run_systems: job_data = {0}".format(job_data)) 
    print("run_systems: output_dir = {0}".format(output_dir)) 
    
    genome_ids = job_data["genome_ids"]
    if len(job_data["genome_groups"]) > 0:
        genome_group_ids = get_genome_group_ids(job_data["genome_groups"],s)
        # make ids unique 
        genome_ids = genome_ids + genome_group_ids
        genome_ids = list(set(genome_ids))

    run_pathways(genome_ids,output_file,output_dir,s)
    #run_subsystems(genome_ids,output_dir)
    #run_families(genome_ids,output_dir)
