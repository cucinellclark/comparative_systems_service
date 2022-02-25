import os, sys

def get_genome_group_ids(group_list):
    print("group_list = {0}".format(group_list))
    base_url = "https://patricbrc.org"
    for genome_group in group_list:
        group_name = os.path.basename(genome_group)
        group_path = os.path.dirname(genome_group)
        genomeGroupSpecifier = "GenomeGroup(" + group_name.replace("/", "%2f") + ")"
        print("group_name = {0}".format(group_name))
        #query = base_url + 

def run_pathways(job_data,output_dir):
    
    # genome_ids will either be an empty list or populated in the job_data 
    genome_ids = job_data["genome_ids"]

    ### get genome ids from genome groups
    genome_group_list = job_data["genome_groups"] 
    genome_group_ids = get_genome_group_ids(genome_group_list)    
    # TODO: change to a set then back to remove duplicates??
    #genome_ids = genome_ids + genome_group_ids

    ### queries
