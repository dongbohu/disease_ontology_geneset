#!/usr/bin/env python3

import requests

def get_url_content(url):
    """Request an input URL and return text lines in the reponse."""

    resp = requests.get(url)
    text_lines = resp.text.strip('\n').split('\n')
    return text_lines


def get_obo_release():
    """Parse HumanDO.obo data file and return its release date."""

    url = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo"
    text_lines = get_url_content(url)
    release_date = ""

    # Find the line that is in the following format:
    # "data-version: doid/releases/YYYY-MM-DD/doid-non-classified.obo"
    # and extract "YYYY-MM-DD" part as the release string.
    for line in text_lines:
        if line.startswith("data-version: "):
            full_version = line.strip().split(' ')[1]
            release_date = full_version.split('/')[2]
            break

    return release_date


def get_genemap2_release():
    """Parse genemap2.txt data file and return its release date."""

    url = "https://raw.githubusercontent.com/greenelab/disease_ontology_geneset/master/data/latest/genemap2.txt"
    text_lines = get_url_content(url)
    genemap2_release = ""

    # Find the line that is in the following format:
    # "data-version: doid/releases/YYYY-MM-DD/doid-non-classified.obo"
    # and extract "YYYY-MM-DD" part as the release string.
    for line in text_lines:
        if line.startswith("# Generated: "):
            release_date = line.strip().split(': ')[1]
            break

    return release_date


def get_release(self):
    """
    Return a string that combines the release dates of both obo and
    genemap2 data files.
    """

    obo_release = get_obo_release()
    genemap2_release = get_genemap2_release()
    return "obo-" + obo_release + "_" + "genemap2-" + genemap2_release


# Test harness
if __name__ == "__main__":
    print(get_release(None))
