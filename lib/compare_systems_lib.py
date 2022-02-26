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
import pandas as pd
import numpy as np

from bvbrc_api import authenticateByEnv,getGenomeGroupIds 

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def build_key(component_list):
    mod_list = []
    for x in component_list:
        if not isinstance(x,str):
            x = 'nan'
        mod_list.append(x)
    #key_parts = [x.lower() for x in mod_list]
    return ":".join(mod_list)

def run_families(genome_ids, output_file, output_dir, session):
    base_query = "https://www.patricbrc.org/api/genome_feature/?in(genome_id,("
    end_query = "))&limit(20000000)&http_accept=text/tsv"
    query = base_query + ",".join(genome_ids) + end_query
    print("GenomeFeatures Query:\n{0}".format(query))
    feature_df = pd.read_csv(query,sep="\t")
    base_query2 = "https://www.patricbrc.org/api/genome_feature/?in(feature_id,("
    end_query2 = "))&limit(20000000)&http_accept=text/tsv"
    print_query = True
    proteinfams_file = os.path.join(output_dir,output_file+"_proteinfams.tsv")
    if os.path.exists(proteinfams_file):
        os.remove(proteinfams_file)
    proteinfams_list = []
    for fids in chunker(feature_df.feature_id.tolist(), 10):
        query2 = base_query2 + ",".join(fids) + end_query2
        if print_query: #print first query
            print("ProteinFamilies Query:\n{0}".format(query2))
            print_query = False
        tmp_df = pd.read_csv(query2,sep="\t")
        proteinfams_list.append(tmp_df)
    # TODO: remove feature_df to save memory???
    # TODO: check the results from concat are correct
    proteinfams_df = pd.concat(proteinfams_list)
    proteinfams_df.to_csv(proteinfams_file, header=True, sep="\t")
    

def run_subsystems(genome_ids, output_file, output_dir, session):
    
    # json(facet,{"stat":{"type":"field","field":"superclass","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","class":{"type":"field","field":"class","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","gene_count":"unique(feature_id)","subclass":{"type":"field","field":"subclass","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","gene_count":"unique(feature_id)"}}}}}}}):  

    # query 
    base_query = "https://www.patricbrc.org/api/subsystem/?in(genome_id,("
    end_query = "))&limit(200000000)&http_accept=text/tsv"
    query = base_query + ",".join(genome_ids) + end_query
    print("Subsystems Query:\n{0}".format(query))
    #subsystems_df = pd.read_csv(query,sep="\t")
    subsystems_file = os.path.join(output_dir,output_file+"_subsystems.tsv")
    #subsystems_df.to_csv(subsystems_file, header=True, sep="\t")

    # TODO: remove, used for testing
    subsystems_df = pd.read_csv(subsystems_file, sep="\t", index_col=0)

    # TODO: combine fact dictionaries later

    # faceting for subsystems overview
    overview_dict = {}
    key_set = set()
    for genome_id in subsystems_df['genome_id'].unique():
        genome_df = subsystems_df.loc[subsystems_df['genome_id'] == genome_id]
        overview_dict[genome_id] = {}
        overview_dict[genome_id]["superclass_counts"] = len(genome_df['superclass'].unique())
        overview_dict[genome_id]["gene_counts"] = len(genome_df['gene'].unique())
        for superclass in genome_df['superclass'].unique():
            superclass_df = subsystems_df.loc[subsystems_df['superclass'] == superclass]
            overview_dict[genome_id][superclass] = {}
            overview_dict[genome_id][superclass]["class_counts"] = len(superclass_df['class'].unique())
            overview_dict[genome_id][superclass]["gene_counts"] = len(superclass_df['gene'].unique())
            for clss in subsystems_df['class'].unique():
                class_df = superclass_df.loc[superclass_df['class'] == clss]
                overview_dict[genome_id][superclass][clss] = {}
                overview_dict[genome_id][superclass][clss]['subclass_counts'] = len(class_df['subclass'].unique())
                overview_dict[genome_id][superclass][clss]['gene_counts'] = len(class_df['gene'].unique())
                for subclass in class_df['subclass'].unique():
                    subclass_df = class_df.loc[class_df['subclass'] == subclass]
                    overview_dict[genome_id][superclass][clss][subclass] = {}
                    overview_dict[genome_id][superclass][clss][subclass]['subsystem_name_counts'] = len(subclass_df['subsystem_name'].unique()) 
                    overview_dict[genome_id][superclass][clss][subclass]['gene_counts'] = len(subclass_df['gene'].unique())
                    for name in subsystems_df['subsystem_name'].unique(): 
                        key = build_key([superclass,clss,subclass,name]) #TODO: might not need 
                        key_set.add(key)

    # faceting for subsystems table

    st_list = [] #subsystem table list
    # get stats on a per-genome basis
    # use unique subsystems_ids to create subsystem table for each genome
    for genome_id in subsystems_df['genome_id'].unique(): 
        print('---{0}'.format(genome_id))
        genome_df = subsystems_df.loc[subsystems_df['genome_id'] == genome_id]
        genome_table = genome_df.drop(['feature_id','public','refseq_locus_tag','role_id','genome_id','taxon_id','gene','role_name','owner','product','patric_id','genome_name','id','_version_','date_inserted','date_modified'], axis=1)
        genome_table = genome_table.drop_duplicates()
        
        # TODO: what is genome_count? Seems like it is always 1
        genome_table['genome_count'] = [1]*genome_table.shape[0]
        genome_table['gene_count'] = [0]*genome_table.shape[0]
        genome_table['role_count'] = [0]*genome_table.shape[0]
    
        genome_table['genome_id'] = [genome_id]*genome_table.shape[0]
        
        for sub_id in genome_table['subsystem_id']:
            tmp_df = genome_df.loc[genome_df['subsystem_id'] == sub_id] 
            genome_table.loc[genome_table['subsystem_id'] == sub_id,'gene_count'] = len(tmp_df['gene'])
            genome_table.loc[genome_table['subsystem_id'] == sub_id,'role_count'] = len(tmp_df['role_id'])

        st_list.append(genome_table)    

    subsystems_table = pd.concat(st_list)
    subsystems_table.to_csv("test.txt",sep="\t")

def run_pathways(genome_ids,output_file,output_dir, session):
    
    ### queries
    base_query = "https://www.patricbrc.org/api/pathway/?in(genome_id,("
    end_query = "))&eq(annotation,PATRIC)&limit(200000000)&http_accept=text/tsv" 
    query = base_query + ",".join(genome_ids) + end_query
    print("Pathways Query:\n{0}".format(query))
    pathway_df = pd.read_csv(query,sep="\t")
    pathways_file = os.path.join(output_dir,output_file+"_pathways.tsv")
    pathway_df.to_csv(pathways_file, header=True, sep="\t")

def get_genome_group_ids(group_list,s):
    genome_group_ids = []
    for genome_group in group_list:
        list_text = getGenomeGroupIds(genome_group,s,genomeGroupPath=True).strip('][').split(',')
        # ['{"genome_id":"562.79202"}', '{"genome_id":"562.80445"}', '{"genome_id":"562.80446"}']
        for genome_entry in list_text:
            genome_entry = genome_entry.strip('}{').split(":")
            genome_group_ids.append(genome_entry[1].strip('\"'))
    return genome_group_ids

# TODO:
#   - Write out pandas dfs without indices

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

    # TODO: add chunking
    #run_pathways(genome_ids,output_file,output_dir,s)
    run_subsystems(genome_ids,output_file,output_dir,s)
    #run_families(genome_ids,output_file,output_dir,s)
