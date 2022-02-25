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

from bvbrc_api import authenticateByEnv 
from pathways import run_pathways


def run_compare_systems(job_data, output_dir):
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        subprocess.call(["mkdir", "-p", output_dir])

    print("run_systems: job_data = {0}".format(job_data)) 
    print("run_systems: output_dir = {0}".format(output_dir)) 
    
    run_pathways(job_data,output_dir)
    #run_subsystems(job_data,output_dir)
    #run_families(job_data,output_dir)
