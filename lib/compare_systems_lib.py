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

from bvbrc_api import authenticateByEnv,getGenomeIdsByGenomeGroup,getFeatureDataFrame,getSubsystemsDataFrame,getPathwayDataFrame,getDataForGenomes,getQueryData

import time
import io

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

def return_columns_to_remove(system,columns):
    if system is 'subsystems_genes':
        drop_subsystem_columns = ['date_inserted','date_modified','genome_name','gene','owner','patric_id','public','product','refseq_locus_tag','taxon_id','_version_']
        table_columns = list(set.intersection(set(drop_subsystem_columns),set(columns))) 
        return table_columns
    elif system is 'subsystems_subsystems':
        drop_subsystem_columns = ['feature_id','public','role_id','genome_id','taxon_id','role_name','owner','product','patric_id','genome_name','id','_version_','date_inserted','date_modified']
        table_columns = list(set.intersection(set(drop_subsystem_columns),set(columns))) 
        return table_columns
    elif system is 'proteinfamilies_plfams':
        drop_plfams_columns = ['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','figfam_id','pgfam_id','protein_id',
                                    'aa_length','gene','go']
        table_columns = list(set.intersection(set(drop_plfams_columns),set(columns))) 
        return table_columns
    elif system is 'proteinfamilies_pgfams':
        drop_pgfams_columns = ['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','figfam_id','plfam_id','protein_id',
                                    'aa_length','gene','go']
        table_columns = list(set.intersection(set(drop_pgfams_columns),set(columns))) 
        return table_columns
    elif system is 'pathways_genes':
        drop_gene_columns = ['genome_name','accession','alt_locus_tag','refseq_locus_tag','feature_id','annotation','product']
        table_columns = list(set.intersection(set(drop_gene_columns),set(columns)))
        return table_columns
    else: # pathways does not have drop columns
        sys.stderr.write("Error, system is not a valid type\n")
        return [] 

def get_plfam_stats(row,stats_df,stats_name):
    plfam_stats = stats_df.loc[row['plfam_id']]
    is_dataframe = isinstance(plfam_stats,pd.DataFrame)
    if stats_name is 'aa_length_min':
        min_value = np.min(plfam_stats['aa_length']) if is_dataframe else plfam_stats['aa_length']
        return min_value
    elif stats_name is 'aa_length_max':
        max_value = np.max(plfam_stats['aa_length']) if is_dataframe else plfam_stats['aa_length']
        return max_value
    elif stats_name is 'aa_length_mean':
        mean_value = np.mean(plfam_stats['aa_length']) if is_dataframe else plfam_stats['aa_length']
        return mean_value
    elif stats_name is 'aa_length_std':
        std_value = np.std(plfam_stats['aa_length']) if is_dataframe else 0
        return std_value
    elif stats_name is 'feature_count': 
        count_value = len(plfam_stats) if is_dataframe else 1
        return count_value
    elif stats_name is 'genomes':
        # genomes used in Heatmap viewer
        genomes_value = format(len(plfam_stats['feature_id']),'#04x').replace('0x','') if is_dataframe else format(1,'#04x').replace('0x','')
        return genomes_value
    else:
        sys.stderr.write(f'invalid stats name: {stats_name}')
        return None
    

def run_families(genome_ids, query_dict, output_file, output_dir, genome_data, session):
    plfam_dict = {}
    pgfam_dict = {}
    genome_dict = {}
    plfam_dict['unique_set'] = set()
    pgfam_dict['unique_set'] = set()
    present_genome_ids = set()
    for gids in chunker(genome_ids, 20):
        base = "https://www.patricbrc.org/api/genome_feature/?http_download=true"
        query = f"in(genome_id,({','.join(gids)}))&limit(2500000)&sort(+feature_id)&eq(annotation,PATRIC)"
        headers = {"accept":"text/tsv", "content-type":"application/rqlquery+x-www-form-urlencoded", 'Authorization': session.headers['Authorization']}
        result_header = True
        for line in getQueryData(base,query,headers):
            if result_header:
                result_header = False
                print(line)
                continue
            line = line.strip().split('\t')
            # 20 entries in query result with pgfam and plfam data
            if len(line) < 20: 
                continue
            try:
                genome_id = line[1].replace('\"','')
                plfam_id = line[14].replace('\"','')
                pgfam_id = line[15].replace('\"','')
                aa_length = line[17].replace('\"','')
                product = line[19].replace('\"'.'')
            except Exception as e:
                sys.stderr.write(f'Error with the following line:\n{e}\n{line}\n')
                continue
            if genome_id not in genome_dict:
                genome_dict[genome_id] = {}
                genome_dict[genome_id]['plfam_set'] = set()
                genome_dict[genome_id]['pgfam_set'] = set()
            present_genome_ids.add(genome_id)
            genome_dict[genome_id]['plfam_set'].add(plfam_id) 
            genome_dict[genome_id]['pgfam_set'].add(pgfam_id) 
            if genome_id not in plfam_dict:
                plfam_dict[genome_id] = {}
            if genome_id not in pgfam_dict:
                pgfam_dict[genome_id] = {}
            if plfam_id and plfam_id not in plfam_dict[genome_id]:
                plfam_dict[genome_id][plfam_id] = {}
                plfam_dict[genome_id][plfam_id]['aa_length_list'] = []
                # TODO: check if I need to check for duplicate features
                plfam_dict[genome_id][plfam_id]['feature_count'] = 0
                plfam_dict[genome_id][plfam_id]['genome_count'] = 1
                plfam_dict[genome_id][plfam_id]['product'] = product
            if pgfam_id and pgfam_id not in pgfam_dict[genome_id]:
                pgfam_dict[genome_id][pgfam_id] = {}
                pgfam_dict[genome_id][pgfam_id]['aa_length_list'] = []
                # TODO: check if I need to check for duplicate features
                pgfam_dict[genome_id][pgfam_id]['feature_count'] = 0
                pgfam_dict[genome_id][pgfam_id]['genome_count'] = 1
                pgfam_dict[genome_id][pgfam_id]['product'] = product
            if plfam_id:
                plfam_dict['unique_set'].add(plfam_id)
                plfam_dict[genome_id][plfam_id]['aa_length_list'].append(int(aa_length))
                plfam_dict[genome_id][plfam_id]['feature_count'] = plfam_dict[genome_id][plfam_id]['feature_count'] + 1
            if pgfam_id:
                pgfam_dict['unique_set'].add(pgfam_id)
                pgfam_dict[genome_id][pgfam_id]['aa_length_list'].append(int(aa_length))
                pgfam_dict[genome_id][pgfam_id]['feature_count'] = pgfam_dict[genome_id][pgfam_id]['feature_count'] + 1
        
    plfam_line_list = []        
    pgfam_line_list = []
    header = 'family_id\tgenome_id\tfeature_count\tgenome_count\tproduct\taa_length_min\taa_length_max\taa_length_mean\taa_length_std'
    plfam_line_list.append(header)
    pgfam_line_list.append(header)
    for plfam_id in plfam_dict['unique_set']: 
        for gid in genome_ids: 
            if gid not in plfam_dict:
                continue
            #plfam_data['plfam_id'] = plfam_id
            #plfam_data['genome_id'] = gid
            if plfam_id in plfam_dict[gid]:
                aa_length_list = plfam_dict[gid][plfam_id]['aa_length_list'] 
                aa_length_max = max(aa_length_list)
                aa_length_min = min(aa_length_list)
                aa_length_mean = np.mean(aa_length_list)
                aa_length_std = np.std(aa_length_list)
                feature_count = plfam_dict[gid][plfam_id]['feature_count']
                genome_count = plfam_dict[gid][plfam_id]['genome_count']
                genomes = format(feature_count,'#04x').replace('0x','')
                product = plfam_dict[gid][plfam_id]['product']
                plfam_str = f'{plfam_id}\t{gid}\t{feature_count}\t{genome_count}\t{product}\t{aa_length_min}\t{aa_length_max}\t{aa_length_mean}\t{aa_length_std}'
                plfam_line_list.append(plfam_str)
    
    for pgfam_id in pgfam_dict['unique_set']: 
        for gid in genome_ids: 
            if gid not in pgfam_dict:
                continue
            pgfam_data = {}
            pgfam_data['pgfam_id'] = pgfam_id
            pgfam_data['genome_id'] = gid
            if pgfam_id in pgfam_dict[gid]:
                aa_length_list = pgfam_dict[gid][pgfam_id]['aa_length_list'] 
                aa_length_max = max(aa_length_list)
                aa_length_min = min(aa_length_list)
                aa_length_mean = np.mean(aa_length_list)
                aa_length_std = np.std(aa_length_list)
                feature_count = pgfam_dict[gid][pgfam_id]['feature_count']
                genome_count = pgfam_dict[gid][pgfam_id]['genome_count']
                genomes = format(feature_count,'#04x').replace('0x','')
                product = pgfam_dict[gid][pgfam_id]['product']
                pgfam_str = f'{pgfam_id}\t{gid}\t{feature_count}\t{genome_count}\t{product}\t{aa_length_min}\t{aa_length_max}\t{aa_length_mean}\t{aa_length_std}'
                pgfam_line_list.append(pgfam_str)
    
    output_json = {}
    output_json['plfam'] = '\n'.join(plfam_line_list) 
    output_json['pgfam'] = '\n'.join(pgfam_line_list) 
    #output_json['genome_ids'] = genome_ids
    output_json['genome_ids'] = list(set(genome_ids).intersection(present_genome_ids)) 
    output_json['job_name'] = output_file

    output_json_file = os.path.join(output_dir,output_file+'_proteinfams_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    print("ProteinFamilies Complete")

def run_subsystems(genome_ids, query_dict, output_file, output_dir, genome_data, session):
    
    # json(facet,{"stat":{"type":"field","field":"superclass","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","class":{"type":"field","field":"class","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","gene_count":"unique(feature_id)","subclass":{"type":"field","field":"subclass","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","gene_count":"unique(feature_id)"}}}}}}}):  

    # subsystems_df = getSubsystemsDataFrame(genome_ids,session) 
    subsystems_df = query_dict['subsystems']
    
    # Superclass, class, and subclass can be different cases: convert all to lower case
    #subsystems_df['superclass'] = subsystems_df['superclass'].str.lower()
    #subsystems_df['class'] = subsystems_df['class'].str.lower()
    #subsystems_df['subclass'] = subsystems_df['subclass'].str.lower()

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
        overview_dict[genome_id]["superclass_counts"] = len(genome_df['subsystem_id'].unique())
        # TODO: check that this is correct for each level
        overview_dict[genome_id]["gene_counts"] = genome_df.shape[0] 
        for superclass in genome_df['superclass'].str.lower().unique():
            superclass_df = subsystems_df.loc[subsystems_df['superclass'].str.lower() == superclass]
            overview_dict[genome_id][superclass] = {}
            overview_dict[genome_id][superclass]["class_counts"] = len(superclass_df['subsystem_id'].unique())
            overview_dict[genome_id][superclass]["gene_counts"] = superclass_df.shape[0] 
            for clss in superclass_df['class'].str.lower().unique():
                class_df = superclass_df.loc[superclass_df['class'].str.lower() == clss]
                overview_dict[genome_id][superclass][clss] = {}
                overview_dict[genome_id][superclass][clss]['subclass_counts'] = len(class_df['subsystem_id'].unique())
                overview_dict[genome_id][superclass][clss]['gene_counts'] = class_df.shape[0] 
                for subclass in class_df['subclass'].str.lower().unique():
                    subclass_df = class_df.loc[class_df['subclass'].str.lower() == subclass]
                    overview_dict[genome_id][superclass][clss][subclass] = {}
                    overview_dict[genome_id][superclass][clss][subclass]['subsystem_name_counts'] = len(subclass_df['subsystem_id'].unique()) 
                    overview_dict[genome_id][superclass][clss][subclass]['gene_counts'] = subclass_df.shape[0] 
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
        #genome_table = genome_df.drop(['feature_id','public','refseq_locus_tag','role_id','genome_id','taxon_id','gene','role_name','owner','product','patric_id','genome_name','id','_version_','date_inserted','date_modified'], axis=1)
        genome_table = genome_df.drop(return_columns_to_remove('subsystems_subsystems',genome_df.columns.tolist()), axis=1)
        genome_table = genome_table.drop_duplicates()
        
        # TODO: what is genome_count? Seems like it is always 1
        genome_table['genome_count'] = [1]*genome_table.shape[0]
        genome_table['gene_count'] = [0]*genome_table.shape[0]
        genome_table['role_count'] = [0]*genome_table.shape[0]
    
        genome_table['genome_id'] = [genome_id]*genome_table.shape[0]

        # TODO: unknown if this needs to stay or go
        # Get unique subsystem ids, first row for information
        keep_rows = []
        sub_id_list = []
        for i in range(0,genome_table.shape[0]):
            if not genome_table.iloc[i]['subsystem_id'] in sub_id_list:
                sub_id_list.append(genome_table.iloc[i]['subsystem_id'])
                keep_rows.append(i)
        genome_table = genome_table.iloc[keep_rows]

        # add stats columns 
        for sub_id in genome_table['subsystem_id']:
            tmp_df = genome_df.loc[sub_id] 
            if isinstance(tmp_df, pd.DataFrame): # if only one record, returns as Series
                if 'gene' in tmp_df.columns.tolist():
                    genome_table.loc[sub_id,'gene_count'] = len(tmp_df['gene'].unique()) # unsure if this is correct)
                else: # unsure if this is correct
                    genome_table.loc[sub_id,'gene_count'] = tmp_df.shape[0]
                genome_table.loc[sub_id,'role_count'] = len(tmp_df['role_id'].unique())
            else:
                genome_table.loc[sub_id,'gene_count'] = 1
                genome_table.loc[sub_id,'role_count'] = 1
            # TODO: genome count calculation
        
        st_list.append(genome_table)    

    subsystems_table = pd.concat(st_list)
    #subsystems_table_output_file = subsystems_file.replace('.tsv','_summary.tsv')
    #subsystems_table.to_csv(subsystems_table_output_file,sep="\t",index=False)

    gene_df = query_dict['feature']
    # change column names
    
    #gene_df.drop(['_version_'],inplace=True)
    # Add subsystems columns to genes table
    #remove_columns = ['date_inserted','date_modified','genome_name','gene','owner','patric_id','public','product','refseq_locus_tag','taxon_id','_version_']
    gene_df = pd.merge(gene_df,subsystems_df.drop(return_columns_to_remove('subsystems_genes',subsystems_df.columns.tolist()),axis=1),on=['genome_id','feature_id'],how='inner')

    output_json_file = subsystems_file.replace('.tsv','_tables.json')
    
    output_json = {}
    output_json['genome_ids'] = genome_ids
    output_json['genome_names'] = genome_data.set_index('Genome ID').loc[genome_ids]['Genome Name'].tolist() # returns a list of genome names in the same order as the genome ids 
    output_json['overview'] = overview_dict 
    output_json['job_name'] = output_file
    output_json['subsystems'] = subsystems_table.to_csv(index=False,sep='\t')
    output_json['genes'] = gene_df.to_csv(index=False,sep='\t')
    with open(output_json_file,'w') as o:
        o.write(json.dumps(output_json))

    print('Subsystems complete')

def run_pathways(genome_ids, query_dict, output_file,output_dir, genome_data, session):
    
    pathways_file = os.path.join(output_dir,output_file+'_pathways.tsv')
    # TODO: create alt_locus_tag query
    # pathway_df = getPathwayDataFrame(genome_ids,session, limit=2500000)
    pathway_df = query_dict['pathway']
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
        ec_table = genome_df[['genome_id','annotation','pathway_class','pathway_name','pathway_id','ec_number','ec_description','ec_index']]

        pathway_table = pathway_table.drop_duplicates()
        ec_table = ec_table.drop_duplicates()

        # add pathway stats columns 
        pathway_table['genome_count'] = [1]*pathway_table.shape[0]
        pathway_table['gene_count'] = [0]*pathway_table.shape[0]
        pathway_table['ec_count'] = [0]*pathway_table.shape[0]
        pathway_table['genome_ec'] = [0]*pathway_table.shape[0]
        for pathway_id in pathway_table['pathway_id']: # should be unique
            tmp_df = genome_df.loc[pathway_id]
            if isinstance(tmp_df, pd.DataFrame): # if only one entry, returns a Series 
                pathway_table.loc[pathway_id,'gene_count'] = len(tmp_df['feature_id'].unique())
                pathway_table.loc[pathway_id,'ec_count'] = len(tmp_df['ec_number'].unique())
                pathway_table.loc[pathway_id,'genome_ec'] = len(tmp_df['genome_ec'].unique())
                # for genes info: take first record
                p_id = tmp_df.iloc[0]['patric_id']
            else:
                pathway_table.loc[pathway_id,'gene_count'] = 1 
                pathway_table.loc[pathway_id,'ec_count'] = 1 
                pathway_table.loc[pathway_id,'genome_ec'] = 1 
                # for genes info: take first record
                p_id = tmp_df['patric_id']
            if p_id not in genes_info_dict:
                genes_info_dict[p_id] = tmp_df.iloc[0] 

        # TODO: maybe modify once the multiple-rows thing has been decided
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
        genome_df.set_index('ec_index',inplace=True)
        ec_table.set_index('ec_index',inplace=True)
        for ec_number in ec_table['ec_number']:
            tmp_df = genome_df.loc[ec_number]
            if isinstance(tmp_df, pd.DataFrame): # if only one entry, returns a Series
                ec_table.loc[ec_number,'gene_count'] = len(tmp_df['feature_id'].unique()) 
                ec_table.loc[ec_number,'ec_count'] = len(tmp_df['ec_number'].unique()) 
                ec_table.loc[ec_number,'genome_ec'] = len(tmp_df['ec_number'].unique())
            else:
                ec_table.loc[ec_number,'gene_count'] = 1 
                ec_table.loc[ec_number,'ec_count'] = 1 
                ec_table.loc[ec_number,'genome_ec'] = 1 
        # append to lists
        pathways_list.append(pathway_table)
        ecnum_list.append(ec_table)

    # concatenate tables and do final calculations: genome_count, genome_ec
    # counting is done per-genome, multi-genome calculation adjustments are done on the front end
    pathway_output = pd.concat(pathways_list)
    ec_output = pd.concat(ecnum_list)

    # Get gene data
    feature_list = []
    # gene_df = getFeatureDataFrame(genome_ids,session, limit=2500000)
    gene_df = query_dict['feature']
    
    genes_output = pd.merge(gene_df.drop(return_columns_to_remove('pathways_genes',gene_df.columns.tolist()), axis=1),pathway_df,on=['genome_id','patric_id'],how='inner')

    if False:
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

    #genes_output = pd.concat(genes_list)

    output_json = {}
    output_json['pathway'] = pathway_output.to_csv(index=False,sep='\t')
    output_json['ecnumber'] = ec_output.to_csv(index=False,sep='\t')
    output_json['genes'] = genes_output.to_csv(index=False,sep='\t')
    output_json['genome_ids'] = genome_ids
    output_json['job_name'] = output_file
    

    output_json_file = pathways_file.replace('.tsv','_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    print("Pathways Complete")

# Store pathways, subsystems, and features queries in a dictionary
def run_all_queries(genome_ids, session):
    query_dict = {}
    ### Run pathways query
    if True:
        print('pathways query')
        pathway_df = getPathwayDataFrame(genome_ids,session, limit=2500000)
        if not pathway_df is None:
            pathway_df['pathway_index'] = pathway_df['pathway_id']
            pathway_df['ec_index'] = pathway_df['ec_number']
            pathway_df.set_index('pathway_index', inplace=True)
            query_dict['pathway'] = pathway_df
        else:
            sys.stderr.write('Pathways dataframe is None\n')
    ### Run subsystems query
    if True:
        print('subsystems query')
        subsystems_df = getSubsystemsDataFrame(genome_ids,session) 
        if not subsystems_df is None:
            subsystems_df['subsystem_index'] = subsystems_df['subsystem_id']
            subsystems_df.set_index('subsystem_index', inplace=True)
            query_dict['subsystems'] = subsystems_df
        else:
            sys.stderr.write('Subsystems dataframe is None\n')
    ### Run features query
    if True:
        print('features query')
        feature_df = getFeatureDataFrame(genome_ids,session, limit=2500000)
        if not feature_df is None:
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
            if 'Genome ID' in feature_df.columns:
                feature_df.rename(columns=column_map, inplace=True)
            feature_df['plfam_index'] = feature_df['plfam_id']
            feature_df['pgfam_index'] = feature_df['pgfam_id']
            #feature_df.set_index('plfam_index', inplace=True)
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

    # optionally add more genome info to output 
    genome_data = getDataForGenomes(genome_ids,s) 

    query_dict = run_all_queries(genome_ids, s)

    # TODO: add chunking
    # TODO: add recipe
    # TODO: add multithreading
    run_pathways(genome_ids, query_dict, output_file, output_dir, genome_data, s)
    run_subsystems(genome_ids, query_dict, output_file, output_dir, genome_data, s)
    run_families(genome_ids, query_dict, output_file, output_dir, genome_data, s)
