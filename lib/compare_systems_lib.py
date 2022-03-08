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

# Returns the entry with the maximum occurences in a pandas dataframe column
# df: pandas dataframe
# col: column name
def get_maximum_value(df, col):
    # value counts returns an ordered list, descending
    col_counts = df[col].value_counts() 
    max_label = col_counts[0].index
    return max_label

def run_families(genome_ids, output_file, output_dir, session):
    base_query = "https://www.patricbrc.org/api/genome_feature/?in(genome_id,("
    end_query = "))&limit(20000000)&http_accept=text/tsv"
    query = base_query + ",".join(genome_ids) + end_query
    print("GenomeFeatures Query:\n{0}".format(query))
    #feature_df = pd.read_csv(query,sep="\t")
    base_query2 = "https://www.patricbrc.org/api/genome_feature/?in(feature_id,("
    end_query2 = "))&limit(20000000)&http_accept=text/tsv"
    print_query = True
    proteinfams_file = os.path.join(output_dir,output_file+"_proteinfams.tsv")
    '''
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
    '''
    # TODO: remove, used for testing
    proteinfams_df = pd.read_csv(proteinfams_file,sep="\t",index_col=0)
    
    plfam_list = [] 
    pgfam_list = []
    figfam_list = []
    for genome_id in proteinfams_df['genome_id'].unique():
        print("---{0}".format(genome_id))    
    
        genome_df = proteinfams_df.loc[proteinfams_df['genome_id'] == genome_id]

        plfam_table = genome_df.drop(['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','figfam_id','pgfam_id','protein_id',
                                    'aa_length','na_length','gene','go'], axis=1) 
        pgfam_table = genome_df.drop(['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','figfam_id','plfam_id','protein_id',
                                    'aa_length','na_length','gene','go'], axis=1)
        '''
        figfam_table = genome_df.drop(['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','plfam_id','pgfam_id','protein_id',
                                    'aa_length','na_length','gene','go'], axis=1)
        '''
        import pdb
        pdb.set_trace()
        # TODO: different product field for same pl/pg_fam id
        plfam_table = plfam_table.drop(['product'],axis=1)
        pgfam_table = pgfam_table.drop(['product'],axis=1)
        ###

        # TODO: floor or ceiling or decimal for mean and stddev stats
        # TODO: need to check the stats calculationgs: min and max are having issues
        # plfam_stats 
        plfam_table['genome_count'] = [1]*plfam_table.shape[0] 
        plfam_table['min_aa_length'] = [0]*plfam_table.shape[0] 
        plfam_table['max_aa_length'] = [0]*plfam_table.shape[0] 
        plfam_table['mean_aa_length'] = [0]*plfam_table.shape[0] 
        plfam_table['stddev_aa_length'] = [0]*plfam_table.shape[0] 
        # TODO: fix product
        plfam_table['product'] = ['PRODUCT DESCRIPTION']*plfam_table.shape[0]

        for plfam_id in plfam_table['plfam_id']:
            tmp_df = genome_df.loc[genome_df['plfam_id'] == plfam_id]
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'min_aa_length'] = np.min(tmp_df['aa_length'])
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'max_aa_length'] = np.max(tmp_df['aa_length'])
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'mean_aa_length'] = np.mean(tmp_df['aa_length'])
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'max_aa_length'] = np.std(tmp_df['aa_length'])
        # pgfam_stats 
        pgfam_table['genome_count'] = [1]*pgfam_table.shape[0] 
        pgfam_table['min_aa_length'] = [0]*pgfam_table.shape[0] 
        pgfam_table['max_aa_length'] = [0]*pgfam_table.shape[0] 
        pgfam_table['mean_aa_length'] = [0]*pgfam_table.shape[0] 
        pgfam_table['stddev_aa_length'] = [0]*pgfam_table.shape[0] 
        # TODO: fix product
        pgfam_table['product'] = ['PRODUCT DESCRIPTION']*pgfam_table.shape[0]
        '''
        for pgfam_id in pgfam_table['pgfam_id']:
            tmp_df = genome_df.loc[genome_df['pgfam_id'] == pgfam_id]
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'min_aa_length'] = np.min(tmp_df['aa_length'])
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'max_aa_length'] = np.max(tmp_df['aa_length'])
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'mean_aa_length'] = np.mean(tmp_df['aa_length'])
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'max_aa_length'] = np.std(tmp_df['aa_length'])
        '''
        # figfam_stats 
        '''
        figfam_table['genome_count'] = [1]*figfam_table.shape[0] 
        figfam_table['min_aa_length'] = [0]*figfam_table.shape[0] 
        figfam_table['max_aa_length'] = [0]*figfam_table.shape[0] 
        figfam_table['mean_aa_length'] = [0]*figfam_table.shape[0] 
        figfam_table['stddev_aa_length'] = [0]*figfam_table.shape[0] 

        for figfam_id in figfam_table['figfam_id']:
            tmp_df = figfam_table.loc[figfam_table['figfam_id'] == figfam_id]
            figfam_table.loc[figfam_table['figfam_id'] == figfam_id,'min_aa_length'] = np.min(tmp_df['aa_length'])
            figfam_table.loc[figfam_table['figfam_id'] == figfam_id,'max_aa_length'] = np.max(tmp_df['aa_length'])
            figfam_table.loc[figfam_table['figfam_id'] == figfam_id,'mean_aa_length'] = np.mean(tmp_df['aa_length'])
            figfam_table.loc[figfam_table['figfam_id'] == figfam_id,'max_aa_length'] = np.std(tmp_df['aa_length'])
        '''


        plfam_list.append(plfam_table)
        pgfam_list.append(pgfam_table)
        #figfam_list.append(figfam_table)
    # TODO: across genome stats adjustments

    # write out tables
    plfam_output = pd.concat(plfam_list)
    pgfam_output = pd.concat(pgfam_list)
    #figfam_output = pd.concat(figfam_list)
    plfam_summary_file = proteinfams_file.replace(".tsv","_plfam_summary.tsv")
    pgfam_summary_file = proteinfams_file.replace(".tsv","_pgfam_summary.tsv")
    #figfam_summary_file = proteinfams_file.replace(".tsv","_figfam_summary.tsv")

    plfam_output.to_csv(plfam_summary_file,sep="\t",index=False)
    pgfam_output.to_csv(pgfam_summary_file,sep="\t",index=False)
    #figfam_output.to_csv(figfam_summary_file,sep="\t",index=False)

    print("ProteinFamilies Complete")

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
                    #for name in subsystems_df['subsystem_name'].unique(): 
                    #    key = build_key([superclass,clss,subclass,name]) #TODO: might not need 
                    #    key_set.add(key)

    subsystems_overview_file = subsystems_file.replace('.tsv','_overview.json')
    with open(subsystems_overview_file,'w') as o:
        json.dump(overview_dict,o)

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
        
        # add stats columns 
        for sub_id in genome_table['subsystem_id']:
            tmp_df = genome_df.loc[genome_df['subsystem_id'] == sub_id] 
            genome_table.loc[genome_table['subsystem_id'] == sub_id,'gene_count'] = len(tmp_df['gene'])
            genome_table.loc[genome_table['subsystem_id'] == sub_id,'role_count'] = len(tmp_df['role_id'])
            # TODO: genome count calculation
        
        # TODO: genes tab table

        st_list.append(genome_table)    

    # TODO: genes table
    subsystems_table = pd.concat(st_list)
    subsystems_table_output_file = subsystems_file.replace('.tsv','_summary.tsv')
    subsystems_table.to_csv(subsystems_table_output_file,sep="\t",index=False)

def run_pathways(genome_ids,output_file,output_dir, session):
    
    ### queries
    base_query = "https://www.patricbrc.org/api/pathway/?in(genome_id,("
    end_query = "))&eq(annotation,PATRIC)&limit(200000000)&http_accept=text/tsv" 
    query = base_query + ",".join(genome_ids) + end_query
    print("Pathways Query:\n{0}".format(query))
    pathway_df = pd.read_csv(query,sep="\t")
    pathways_file = os.path.join(output_dir,output_file+"_pathways.tsv")
    pathway_df.to_csv(pathways_file, header=True, sep="\t", index=False)

    # TODO: create alt_locus_tag query

    #TODO: remove, reading in file for testing
    #pathway_df = pd.read_csv(pathways_file,sep="\t")
    pathways_list = []
    ecnum_list = []
    genes_list = []
    for genome_id in genome_ids: 
        genome_id = float(genome_id)
        print('---Faceting GenomeId: {0}---'.format(genome_id))
        genome_df = pathway_df.loc[pathway_df['genome_id'] == genome_id]

        pathway_table = genome_df.drop(['pathway_ec','genome_name','accession','genome_ec',
                                        'product','gene','public','patric_id','sequence_id','ec_number',
                                        'feature_id','taxon_id','ec_description','refseq_locus_tag','owner',
                                        'id','_version_','date_inserted','date_modified'], axis=1)

        # TODO: add index column
        ec_table = genome_df.drop(['pathway_ec','genome_name','accession','genome_ec','product','feature_id',
                                    'gene','public','patric_id','sequence_id','taxon_id','refseq_locus_tag',
                                    'owner','id','_version_','date_inserted','date_modified'], axis=1)
        # TODO: add alt_locus_tag column
        genes_table = genome_df.drop(['pathway_ec','genome_ec','public','sequence_id','feature_id',
                                        'taxon_id','owner','id','_version_','date_inserted','date_modified'], axis=1)

        pathway_table = pathway_table.drop_duplicates()
        ec_table = ec_table.drop_duplicates()
        genes_table = genes_table.drop_duplicates()
        # TODO: other fillna?
        genes_table.gene = genes_table.gene.fillna('')

        # add pathway stats columns 
        pathway_table['genome_count'] = [1]*pathway_table.shape[0]
        pathway_table['gene_count'] = [0]*pathway_table.shape[0]
        pathway_table['ec_count'] = [0]*pathway_table.shape[0]
        pathway_table['genome_ec'] = [0]*pathway_table.shape[0]
        for pathway_id in pathway_table['pathway_id']:
            tmp_df = genome_df.loc[genome_df['pathway_id'] == pathway_id]
            pathway_table.loc[pathway_table['pathway_id'] == pathway_id,'gene_count'] = len(tmp_df['feature_id'].unique())
            pathway_table.loc[pathway_table['pathway_id'] == pathway_id,'ec_count'] = len(tmp_df['ec_number'].unique())
            pathway_table.loc[pathway_table['pathway_id'] == pathway_id,'genome_ec'] = len(tmp_df['ec_number'].unique())

        # get first ec_number entry data for ec_number duplicates
        # pathway_id is not duplicated for each ec_number
        keep_rows = []
        ec_list = []
        for i in range(0,ec_table.shape[0]):
            if not ec_table.iloc[i]['ec_number'] in ec_list:
                ec_list.append(ec_table.iloc[i]['ec_number'])
                keep_rows.append(i)
        ec_table = ec_table.iloc[keep_rows]

        # add ec_number stats columns
        ec_table['genome_count'] = [1]*ec_table.shape[0]
        ec_table['gene_count'] = [0]*ec_table.shape[0]
        ec_table['ec_count'] = [0]*ec_table.shape[0]
        ec_table['genome_ec'] = [0]*ec_table.shape[0]
        for ec_number in ec_table['ec_number']:
            tmp_df = genome_df.loc[genome_df['ec_number'] == ec_number]
            ec_table.loc[ec_table['ec_number'] == ec_number,'gene_count'] = len(tmp_df['feature_id'].unique()) 
            ec_table.loc[ec_table['ec_number'] == ec_number,'ec_count'] = len(tmp_df['ec_number'].unique()) 
            ec_table.loc[ec_table['ec_number'] == ec_number,'genome_ec'] = len(tmp_df['ec_number'].unique())

        # get first gene entry data for gene duplicates
        # pathway_id is not duplicated for each gene
        keep_rows = []
        g_list = []
        for i in range(0,genes_table.shape[0]):
            if not genes_table.iloc[i]['gene'] in g_list:
                g_list.append(genes_table.iloc[i]['gene'])
                keep_rows.append(i)
        genes_table = genes_table.iloc[keep_rows]

        # genes table stats
        genes_table['genome_count'] = [1]*genes_table.shape[0]
        genes_table['gene_count'] = [0]*genes_table.shape[0]
        genes_table['ec_count'] = [0]*genes_table.shape[0]
        genes_table['genome_ec'] = [0]*genes_table.shape[0]
        genes_table['alt_locus_tag'] = ['TMP_Alt_LOCUS']*genes_table.shape[0]
        for gene in genes_table['gene']:
            tmp_df = genome_df.loc[genome_df['gene'] == gene]
            genes_table.loc[genes_table['gene'] == gene,'gene_count'] = len(tmp_df['feature_id'].unique())
            genes_table.loc[genes_table['gene'] == gene,'ec_count'] = len(tmp_df['ec_number'].unique())
            genes_table.loc[genes_table['gene'] == gene,'genome_ec'] = len(tmp_df['ec_number'].unique())
        
        # append to lists
        pathways_list.append(pathway_table)
        ecnum_list.append(ec_table)
        genes_list.append(genes_table)

    # concatenate tables and do final calculations: genome_count, genome_ec
    # counting is done per-genome, multi-genome calculation adjustments are done on the front end
    pathway_output = pd.concat(pathways_list)
    ec_output = pd.concat(ecnum_list)
    genes_output = pd.concat(genes_list)

    output_json = {}
    output_json['pathway'] = pathway_output.to_csv(index=False,sep='\t')
    output_json['ecnumber'] = ec_output.to_csv(index=False,sep='\t')
    output_json['genes'] = genes_output.to_csv(index=False,sep='\t')
    output_json['genome_ids'] = genome_ids

    output_json_file = pathways_file.replace('.tsv','_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    # TODO: change to json file
    #return
    # write out tables
    #pathway_summary_file = pathways_file.replace(".tsv","_pathway_summary.tsv")
    #ec_summary_file = pathways_file.replace(".tsv","_ec_summary.tsv")
    # TODO: genes_output

    #pathway_output.to_csv(pathway_summary_file,sep="\t",index=False)
    #ec_output.to_csv(ec_summary_file,sep="\t",index=False)
    # TODO: write genes summary file
    print("Pathways Complete")

def get_genome_group_ids(group_list,session):
    genome_group_ids = []
    for genome_group in group_list:
        list_text = getGenomeGroupIds(genome_group,session,genomeGroupPath=True).strip('][').split(',')
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


    print("Run ComparativeSystems:\njob_data = {0}".format(job_data)) 
    print("output_dir = {0}".format(output_dir)) 
    
    genome_ids = job_data["genome_ids"]
    if len(job_data["genome_groups"]) > 0:
        genome_group_ids = get_genome_group_ids(job_data["genome_groups"],s)
        # make ids unique 
        genome_ids = genome_ids + genome_group_ids
        genome_ids = list(set(genome_ids))

    # TODO: add chunking
    # TODO: add recipe
    run_pathways(genome_ids,output_file,output_dir,s)
    #run_subsystems(genome_ids,output_file,output_dir,s)
    #run_families(genome_ids,output_file,output_dir,s)
