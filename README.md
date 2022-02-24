# Comparative Systems
```
Mode 1 (Service)
Inputs:
#1 Job data
json file for job {"reference_genome_id": "1310806.3", \
                    "output_file": "rnaseq_baumanii_1505311", \
                    "recipe": [], "output_path": "/anwarren@patricbrc.org/home/test",\
                    }
#2 Server string to get genome from using PATRIC ID system
#3 Override parameter string to govern number of processes used
#4 Output directory

Outputs:
#1 Pathways JSON file
#2 Subsystems JSON file
#3 Protein families JSON file

This program will:
#1 Use the given genome_id , genome_groups to create the JSON structures specified

Dependencies:
Python 3.x: requests, pandas
```
