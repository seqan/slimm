#!/usr/bin/env python

import argparse
import glob
import os
import re
import datetime
import sys
import urllib
import shutil
import tarfile
import Queue
import thread
import threading

def extract_then_delete(source, destination):
    if (source.endswith("tar.gz")):
        tar = tarfile.open(source, "r:gz")
        tar.extractall(path=destination)
        tar.close()
    elif (source.endswith("tar")):
        tar = tarfile.open(source, "r:")
        tar.extractall(path=destination)
        tar.close()
    os.remove(source)



def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def taxonomy_download(taxa_type, working_dir, today_string):
    taxa_folder_path = working_dir + "/" + taxa_type + "_" + today_string
    new_file_path = taxa_folder_path + ".tar.gz"
    print "Downloading " + taxa_type + " file ..."

    # DELETE FOR TESTING ONLY
    # old_taxa_folder_path = working_dir + "/.old/" + taxa_type + "_" + today_string
    # shutil.copytree(old_taxa_folder_path, taxa_folder_path)
    urllib.urlretrieve("ftp://ftp.ncbi.nih.gov/pub/taxonomy/" + taxa_type + ".tar.gz", new_file_path)
    extract_then_delete(new_file_path, taxa_folder_path)
    return taxa_folder_path

def download_one(download_queue, genomes_dir):
    q = download_queue.get()
    try:
        taxa_id = q[0]
        url = q[1]
        destination = genomes_dir + "/" + str(taxa_id) + ".fna.gz"
        urllib.urlretrieve(url, destination)
    except Exception, e:
        download_queue.put(q)

