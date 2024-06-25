
# Application specification: ComparativeSystems

This is the application specification for service with identifier ComparativeSystems.

The backend script implementing the application is [App-ComparativeSystems.pl](../service-scripts/App-ComparativeSystems.pl).

The raw JSON file for this specification is [ComparativeSystems.json](ComparativeSystems.json).

This service performs the following task:   Create datastructures to decompose genomes

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |
| genome_ids | Genome Ids | list  |  | ARRAY(0x55d0a0f580c8) |
| genome_groups | Genome Groups | list  |  | ARRAY(0x55d0a0feb448) |

