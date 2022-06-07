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

from bvbrc_api import authenticateByEnv,getGenomeIdsByGenomeGroup,getFeatureDf,getSubsystemsDf,getPathwayDf

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
    #base_query = "https://www.patricbrc.org/api/genome_feature/?in(genome_id,("
    #end_query = "))&limit(-1)&http_accept=text/tsv"
    #query = base_query + ",".join(genome_ids) + end_query
    #print("GenomeFeatures Query:\n{0}".format(query))
    #feature_df = pd.read_csv(query,sep="\t")
    feature_list = []
    # proteinfams_df = getFeatureDf(genome_ids,session, limit=2500000)
    # change column names
    column_map = {
        'Genome': 'genome_name',
        'Genome ID': 'genome_id',
        'Accession': 'accession',
        'BRC ID': 'patric_id',
        'RefSeq Locus Tag': 'refseq_locus_tag',
        'Alt Locus Tag': 'alt_locus_tag',
        'Feature ID': 'feature_id',
        'Annotation': 'annotation',
        'Feature Type': 'feature_type',
        'Start': 'start',
        'End': 'end',
        'Length': 'length',
        'Strand': 'strand',
        'FIGfam ID': 'figfam_id',
        'PATRIC genus-specific families (PLfams)': 'plfam_id',
        'PATRIC cross-genus families (PGfams)': 'pgfam_id',
        'Protein ID': 'protein_id',
        'AA Length': 'aa_length',
        'Gene Symbol': 'gene',
        'Product': 'product',
        'GO': 'go'
    }
    proteinfams_df.rename(columns=column_map, inplace=True)

    proteinfams_file = os.path.join(output_dir,output_file+"_proteinfams.tsv")
    proteinfams_df.to_csv(proteinfams_file, index=False, header=True, sep="\t")
    # TODO: remove, used for testing
    #proteinfams_df = pd.read_csv(proteinfams_file,sep="\t")
    
    plfam_list = [] 
    pgfam_list = []
    for genome_id in genome_ids:
        print("---{0}".format(genome_id))    
        genome_df = proteinfams_df.loc[proteinfams_df['genome_id'] == genome_id]

        plfam_table = genome_df.drop(['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','figfam_id','pgfam_id','protein_id',
                                    'aa_length','gene','go'], axis=1) 
        pgfam_table = genome_df.drop(['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','figfam_id','plfam_id','protein_id',
                                    'aa_length','gene','go'], axis=1)

        # Get unique family ids, first row for information 
        keep_rows = []
        plfam_id_list = []
        for i in range(0,plfam_table.shape[0]):
            if not plfam_table.iloc[i]['plfam_id'] in plfam_id_list:
                plfam_id_list.append(plfam_table.iloc[i]['plfam_id'])
                keep_rows.append(i)
        plfam_table = plfam_table.iloc[keep_rows]

        keep_rows = []
        pgfam_id_list = []
        for i in range(0,pgfam_table.shape[0]):
            if not pgfam_table.iloc[i]['pgfam_id'] in pgfam_id_list:
                pgfam_id_list.append(pgfam_table.iloc[i]['pgfam_id'])
                keep_rows.append(i)
        pgfam_table = pgfam_table.iloc[keep_rows]

        # plfam_stats 
        plfam_table['feature_count'] = [0]*plfam_table.shape[0]
        plfam_table['genome_count'] = [1]*plfam_table.shape[0] 
        plfam_table['genomes'] = [0]*plfam_table.shape[0]
        plfam_table['aa_length_min'] = [0]*plfam_table.shape[0] 
        plfam_table['aa_length_max'] = [0]*plfam_table.shape[0] 
        plfam_table['aa_length_mean'] = [0]*plfam_table.shape[0] 
        plfam_table['aa_length_std'] = [0]*plfam_table.shape[0] 
        for plfam_id in plfam_table['plfam_id']:
            tmp_df = genome_df.loc[genome_df['plfam_id'] == plfam_id]
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'aa_length_min'] = np.min(tmp_df['aa_length'])
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'aa_length_max'] = np.max(tmp_df['aa_length'])
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'aa_length_mean'] = np.mean(tmp_df['aa_length'])
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'aa_length_std'] = np.std(tmp_df['aa_length'])
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'feature_count'] = len(tmp_df['feature_id'])
            # genomes used in Heatmap viewer
            plfam_table.loc[plfam_table['plfam_id'] == plfam_id,'genomes'] = format(len(tmp_df['feature_id']),'#04x').replace('0x','')

        # pgfam_stats 
        pgfam_table['feature_count'] = [0]*pgfam_table.shape[0] 
        pgfam_table['genome_count'] = [1]*pgfam_table.shape[0] 
        pgfam_table['genomes'] = [0]*pgfam_table.shape[0]
        pgfam_table['aa_length_min'] = [0]*pgfam_table.shape[0] 
        pgfam_table['aa_length_max'] = [0]*pgfam_table.shape[0] 
        pgfam_table['aa_length_mean'] = [0]*pgfam_table.shape[0] 
        pgfam_table['aa_length_std'] = [0]*pgfam_table.shape[0] 
        for pgfam_id in pgfam_table['pgfam_id']:
            tmp_df = genome_df.loc[genome_df['pgfam_id'] == pgfam_id]
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'aa_length_min'] = np.min(tmp_df['aa_length'])
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'aa_length_max'] = np.max(tmp_df['aa_length'])
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'aa_length_mean'] = np.mean(tmp_df['aa_length'])
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'aa_length_std'] = np.std(tmp_df['aa_length'])
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'feature_count'] = len(tmp_df['feature_id'])
            pgfam_table.loc[pgfam_table['pgfam_id'] == pgfam_id,'genomes'] = format(len(tmp_df['feature_id']),'#04x').replace('0x','')


        plfam_list.append(plfam_table)
        pgfam_list.append(pgfam_table)
        #figfam_list.append(figfam_table)

    # write out tables
    # counting is done per-genome, multi-genome calculation adjustments are done on the front end
    plfam_output = pd.concat(plfam_list)
    pgfam_output = pd.concat(pgfam_list)
    
    output_json = {}
    output_json['plfam'] = plfam_output.to_csv(index=False,sep='\t')
    output_json['pgfam'] = pgfam_output.to_csv(index=False,sep='\t')
    output_json['genome_ids'] = genome_ids

    output_json_file = proteinfams_file.replace('.tsv','_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    print("ProteinFamilies Complete")

def run_subsystems(genome_ids, output_file, output_dir, session):
    
    # json(facet,{"stat":{"type":"field","field":"superclass","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","class":{"type":"field","field":"class","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","gene_count":"unique(feature_id)","subclass":{"type":"field","field":"subclass","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","gene_count":"unique(feature_id)"}}}}}}}):  

    # subsystems_df = getSubsystemsDf(genome_ids,session) 

    # Superclass, class, and subclass can be different cases: convert all to lower case
    subsystems_df['superclass'] = subsystems_df['superclass'].str.lower()
    subsystems_df['class'] = subsystems_df['class'].str.lower()
    subsystems_df['subclass'] = subsystems_df['subclass'].str.lower()

    # query 
    #base_query = "https://www.patricbrc.org/api/subsystem/?in(genome_id,("
    #end_query = "))&limit(200000000)&http_accept=text/tsv"
    #query = base_query + ",".join(genome_ids) + end_query
    #print("Subsystems Query:\n{0}".format(query))
    #subsystems_df = pd.read_csv(query,sep="\t")
    subsystems_file = os.path.join(output_dir,output_file+"_subsystems.tsv")
    subsystems_df.to_csv(subsystems_file, header=True, sep="\t")

    # TODO: remove, used for testing
    #subsystems_df = pd.read_csv(subsystems_file, sep="\t", index_col=0)

    # faceting for subsystems overview
    overview_dict = {}
    key_set = set()
    for genome_id in subsystems_df['genome_id'].unique():
        genome_df = subsystems_df.loc[subsystems_df['genome_id'] == genome_id]
        overview_dict[genome_id] = {}
        overview_dict[genome_id]["superclass_counts"] = len(genome_df['superclass'].unique())
        # TODO: check that this is correct for each level
        overview_dict[genome_id]["gene_counts"] = len(genome_df['gene'].unique())
        for superclass in genome_df['superclass'].unique():
            superclass_df = subsystems_df.loc[subsystems_df['superclass'] == superclass]
            overview_dict[genome_id][superclass] = {}
            overview_dict[genome_id][superclass]["class_counts"] = len(superclass_df['class'].unique())
            overview_dict[genome_id][superclass]["gene_counts"] = len(superclass_df['gene'].unique())
            for clss in superclass_df['class'].unique():
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

    #subsystems_overview_file = subsystems_file.replace('.tsv','_overview.json')
    #with open(subsystems_overview_file,'w') as o:
    #    json.dump(overview_dict,o)

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
    #subsystems_table_output_file = subsystems_file.replace('.tsv','_summary.tsv')
    #subsystems_table.to_csv(subsystems_table_output_file,sep="\t",index=False)

    output_json_file = subsystems_file.replace('.tsv','_tables.json')
    
    output_json = {}
    output_json['genome_ids'] = genome_ids
    output_json['overview'] = overview_dict 
    output_json['subsystems'] = subsystems_table.to_csv(index=False,sep='\t')
    with open(output_json_file,'w') as o:
        o.write(json.dumps(output_json))

def run_pathways(genome_ids,output_file,output_dir, session):
    
    pathways_file = os.path.join(output_dir,output_file+'_pathways.tsv')
    # TODO: create alt_locus_tag query
    # pathway_df = getPathwayDf(genome_ids,session, limit=2500000)
    # TODO:
    # - move this to p3_core/lib/bvbrc_api.py
    # convert pathway_id to string and pad with leading zeros
    pathway_df['pathway_id'] = pathway_df['pathway_id'].apply(lambda x: '{0:0>5}'.format(x)) 
    pathway_df.to_csv(pathways_file,sep='\t',index=False)

    #TODO: remove, reading in file for testing
    #pathway_df = pd.read_csv(pathways_file,sep="\t")
    pathways_list = []
    ecnum_list = []
    genes_info_dict = {} # key is patric_id
    # TODO:
    # - make sure to check genome_id type issue 
    for genome_id in genome_ids: 
        print('---Faceting GenomeId: {0}---'.format(genome_id))
        genome_df = pathway_df.loc[pathway_df['genome_id'] == genome_id]
        '''
        pathway_table = genome_df.drop(['pathway_ec','genome_name','accession','genome_ec',
                                        'product','gene','public','patric_id','sequence_id','ec_number',
                                        'feature_id','taxon_id','ec_description','refseq_locus_tag','owner',
                                        'id','_version_','date_inserted','date_modified'], axis=1)
        # TODO: add index column
        ec_table = genome_df.drop(['pathway_ec','genome_name','accession','genome_ec','product','feature_id',
                                    'gene','public','patric_id','sequence_id','taxon_id','refseq_locus_tag',
                                    'owner','id','_version_','date_inserted','date_modified'], axis=1)
        # TODO: add alt_locus_tag column
        '''
        pathway_table = genome_df[['genome_id','annotation','pathway_class','pathway_name','pathway_id']]
        ec_table = genome_df[['genome_id','annotation','pathway_class','pathway_name','pathway_id','ec_number','ec_description']]

        pathway_table = pathway_table.drop_duplicates()
        ec_table = ec_table.drop_duplicates()

        # add pathway stats columns 
        pathway_table['genome_count'] = [1]*pathway_table.shape[0]
        pathway_table['gene_count'] = [0]*pathway_table.shape[0]
        pathway_table['ec_count'] = [0]*pathway_table.shape[0]
        pathway_table['genome_ec'] = [0]*pathway_table.shape[0]
        for pathway_id in pathway_table['pathway_id']:
            tmp_df = genome_df.loc[genome_df['pathway_id'] == pathway_id]
            pathway_table.loc[pathway_table['pathway_id'] == pathway_id,'gene_count'] = len(tmp_df['feature_id'].unique())
            pathway_table.loc[pathway_table['pathway_id'] == pathway_id,'ec_count'] = len(tmp_df['ec_number'].unique())
            pathway_table.loc[pathway_table['pathway_id'] == pathway_id,'genome_ec'] = len(tmp_df['genome_ec'].unique())
            # for genes info: take first record
            p_id = tmp_df.iloc[0]['patric_id']
            if p_id not in genes_info_dict:
                genes_info_dict[p_id] = tmp_df.iloc[0] 

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

        # append to lists
        pathways_list.append(pathway_table)
        ecnum_list.append(ec_table)

    # concatenate tables and do final calculations: genome_count, genome_ec
    # counting is done per-genome, multi-genome calculation adjustments are done on the front end
    pathway_output = pd.concat(pathways_list)
    ec_output = pd.concat(ecnum_list)

    # Get gene data
    feature_list = []
    gene_df = getFeatureDf(genome_ids,session, limit=2500000)
    # change column names
    column_map = {
        'Genome': 'genome_name',
        'Genome ID': 'genome_id',
        'Accession': 'accession',
        'BRC ID': 'patric_id',
        'RefSeq Locus Tag': 'refseq_locus_tag',
        'Alt Locus Tag': 'alt_locus_tag',
        'Feature ID': 'feature_id',
        'Annotation': 'annotation',
        'Feature Type': 'feature_type',
        'Start': 'start',
        'End': 'end',
        'Length': 'length',
        'Strand': 'strand',
        'FIGfam ID': 'figfam_id',
        'PATRIC genus-specific families (PLfams)': 'plfam_id',
        'PATRIC cross-genus families (PGfams)': 'pgfam_id',
        'Protein ID': 'protein_id',
        'AA Length': 'aa_length',
        'Gene Symbol': 'gene',
        'Product': 'product',
        'GO': 'go'
    }
    gene_df.rename(columns=column_map, inplace=True)
    # Parse gene data
    # TODO:
    # - make sure to check genome_id type issue 
    genes_list = []
    for genome_id in genome_ids:
        print('---Faceting GenomeId Genes Table: {0}---'.format(genome_id))
        genome_df = gene_df.loc[gene_df['genome_id'] == genome_id]
        
        genes_table = genome_df.drop_duplicates()
        genes_table.gene = genes_table.gene.fillna('')
    
        # get first gene entry data for gene duplicates
        # pathway_id is not duplicated for each gene
        keep_rows = []
        g_list = []
        for i in range(0,genes_table.shape[0]):
            if not genes_table.iloc[i]['gene'] in g_list:
                g_list.append(genes_table.iloc[i]['gene'])
                keep_rows.append(i)
            if genes_table.iloc[i]['gene'] == '':
                keep_rows.append(i)
        genes_table = genes_table.iloc[keep_rows]

        # fill in ec_description, ec_number, index, pathwayid, pathway name 
        genes_table['ec_description'] = ['']*genes_table.shape[0]
        genes_table['ec_number'] = ['']*genes_table.shape[0]
        genes_table['index'] = ['']*genes_table.shape[0]
        genes_table['pathway_id'] = ['']*genes_table.shape[0]
        genes_table['pathway_name'] = ['']*genes_table.shape[0]
        # genes table stats
        genes_table['genome_count'] = [1]*genes_table.shape[0]
        genes_table['gene_count'] = [0]*genes_table.shape[0]
        genes_table['ec_count'] = [0]*genes_table.shape[0]
        genes_table['genome_ec'] = [0]*genes_table.shape[0]
        #genes_table['alt_locus_tag'] = ['']*genes_table.shape[0]
        # use patric_id to account for '' genes
        for p_id in genes_table['patric_id']:
            tmp_df = pathway_df.loc[pathway_df['genome_id'] == genome_id]
            tmp_df = tmp_df[tmp_df['patric_id'] == p_id]
            genes_table.loc[genes_table['patric_id'] == p_id,'gene_count'] = len(tmp_df['feature_id'].unique())
            genes_table.loc[genes_table['patric_id'] == p_id,'ec_count'] = len(tmp_df['ec_number'].unique())
            genes_table.loc[genes_table['patric_id'] == p_id,'genome_ec'] = len(tmp_df['ec_number'].unique())
            if p_id in genes_info_dict:
                genes_table.loc[genes_table['patric_id'] == p_id,'ec_description'] = genes_info_dict[p_id]['ec_description']
                genes_table.loc[genes_table['patric_id'] == p_id,'ec_number'] = genes_info_dict[p_id]['ec_number']
                genes_table.loc[genes_table['patric_id'] == p_id,'index'] = genes_info_dict[p_id]['id']
                genes_table.loc[genes_table['patric_id'] == p_id,'pathway_id'] = genes_info_dict[p_id]['pathway_id']
                genes_table.loc[genes_table['patric_id'] == p_id,'pathway_name'] = genes_info_dict[p_id]['pathway_name']
 
        genes_list.append(genes_table)

    genes_output = pd.concat(genes_list)

    output_json = {}
    output_json['pathway'] = pathway_output.to_csv(index=False,sep='\t')
    output_json['ecnumber'] = ec_output.to_csv(index=False,sep='\t')
    output_json['genes'] = genes_output.to_csv(index=False,sep='\t')
    output_json['genome_ids'] = genome_ids
    

    output_json_file = pathways_file.replace('.tsv','_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    print("Pathways Complete")

# Store pathways, subsystems, and protein families queries in a dictionary
def run_all_queries(genome_ids, session):
    query_dict = {}
    ### Run pathways query
    if True:
        print('pathways query')
        pathway_df = getPathwayDf(genome_ids,session, limit=2500000)
        if not pathway_df is None:
            query_dict['pathway'] = pathway_df
        else:
            sys.stderr.write('Pathways dataframe is None\n')
    ### Run subsystems query
    if True:
        print('subsystems query')
        subsystems_df = getSubsystemsDf(genome_ids,session) 
        if not subsystems_df is None:
            query_dict['subsystems'] = subsystems_df
        else:
            sys.stderr.write('Subsystems dataframe is None\n')
    ### Run features query
    if True:
        print('features query')
        feature_df = getFeatureDf(genome_ids,session, limit=2500000)
        if not feature_df is None:
            query_dict['feature'] = feature_df
        else:
            sys.stderr.write('Features dataframe is None\n')
    return query_dict

def get_genome_group_ids(group_list,session):
    genome_group_ids = []
    for genome_group in group_list:
        genome_id_list = getGenomeIdsByGenomeGroup(genome_group,session,genomeGroupPath=True)
        genome_group_ids = genome_group_ids + genome_id_list
    return genome_group_ids

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
        if len(genome_group_ids) == 0:
            sys.stderr.write('FAILED to get genome ids for genome groups: exiting')
            sys.exit(-1)
        # make ids unique 
        genome_ids = genome_ids + genome_group_ids
        genome_ids = list(set(genome_ids))

    query_dict = run_all_queries(genome_ids, s)
    import pdb
    pdb.set_trace()

    # TODO: add chunking
    # TODO: add recipe
    # TODO: add multithreading
    #run_pathways(genome_ids,output_file,output_dir,s)
    run_subsystems(genome_ids,output_file,output_dir,s)
    #run_families(genome_ids,output_file,output_dir,s)
