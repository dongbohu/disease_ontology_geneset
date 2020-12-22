#!/usr/bin/env python3


import json
import os
import re
import sys

from datetime import date
from go import go

import logging
logger = logging.getLogger(__name__)

tax_id = 9606

#do_obo_url = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo"
#genemap_url = "https://data.omim.org/downloads/z9hwkkLwTHyKrrmsmXkYiQ/genemap2.txt"

# dhu: hardcoded paths of input files
obo_filename = "data/HumanDO.obo"

# "confidence" column (#7) in "genemap.txt" is deprecated. All values in this
# column are empty.  For differences between "genemap.txt" and "genemap2.txt",
# see the end of both files.
#
# "genemap2.txt" is better because it already maps "MIM Number" to "Entrez Gene ID".
genemap_filename = "data/omim/genemap2.txt"

# Two varibles when searching MIM Disease ID from "Phenotypes" column in `genemap_filename`
FIND_MIMID = re.compile('\, [0-9]* \([1-4]\)')  # Regex pattern
PHENOTYPE_FILTER = '(3)'

# dhu: Do we still need this?
tag_mapping_filename = "data/tissue-disease_curated-associations.txt"


# This function is copied from `annotation-refinery/utils.py`
def build_tags_dictionary(
        tag_mapping_filename,
        geneset_id_column,
        geneset_name_column,
        tag_column,
        header
):
    tags_dict = {}
    tag_file_fh = open(tag_mapping_filename, 'r')

    if header:
        next(tag_file_fh)

    for line in tag_file_fh:
        toks = line.strip().split('\t')
        gs_id = toks[geneset_id_column]
        gs_name = toks[geneset_name_column]
        # Underscores may be used in files in place of spaces
        gs_name = gs_name.replace('_', ' ')
        gs_tag = toks[tag_column]

        if gs_id not in tags_dict:
            tags_dict[gs_id] = {'gs_name': gs_name, 'gs_tags': [gs_tag]}
        else:
            tags_dict[gs_id]['gs_tags'].append(gs_tag)

    tag_file_fh.close()

    return tags_dict


# dhu: This function is copied from "annotation-refinery/process_do.py"
def build_doid_omim_dict(obo_filename):
    """
    Function to read in DO OBO file and build dictionary of DO terms
    from OBO file that have OMIM cross-reference IDs

    Arguments:
    obo_filename -- A string. Location of the DO OBO file to be read in.

    Returns:
    doid_omim_dict -- A dictionary of only the DO terms in the OBO file
    that have OMIM xrefs. The keys in the dictionary are DOIDs, and the
    values are sets of OMIM xref IDs.
    """
    obo_fh = open(obo_filename, 'r')
    doid_omim_dict = {}

    # This statement builds a list of the lines in the file and reverses
    # its order. This is because the list 'pop()' method pops out list
    # elements starting from the end. This way the lines will be read in
    # the following loop in order, from top to bottom of the file.
    obo_reversed_str_array = obo_fh.readlines()[::-1]

    while obo_reversed_str_array:  # Loop adapted from Dima @ Princeton
        line = obo_reversed_str_array.pop()
        line = line.strip()
        if line == '[Term]':
            while line != '' and obo_reversed_str_array:
                line = obo_reversed_str_array.pop()

                if line.startswith('id:'):
                    doid = re.search('DOID:[0-9]+', line)
                    if doid:
                        doid = doid.group(0)

                if line.startswith('xref: OMIM:'):
                    # If term has OMIM xref, get it and add it to the
                    # doid_omim_dict. Otherwise, ignore.
                    omim = re.search('[0-9]+', line).group(0)

                    if doid not in doid_omim_dict:
                        doid_omim_dict[doid] = set()
                    if omim not in doid_omim_dict[doid]:
                        doid_omim_dict[doid].add(omim)

    return doid_omim_dict


class MIMdisease:
    def __init__(self):
        self.id = ''
        self.phenotype = ''  # Phenotype mapping method
        self.genes = []      # list of gene IDs


def build_mim_diseases_dict(genemap_filename):
    """
    Function to parse genemap file and build a dictionary of MIM
    diseases.

    Arguments:
    genemap_filename -- A string. Location of the genemap file to read in.

    Returns:
    mim_diseases -- A dictionary. The keys are MIM disease IDs, and the
    values are `MIMdisease` objects, defined by the class above.

    *N.B. MIM IDs are not all one type of object (unlike Entrez IDs,
    for example) - they can refer to phenotypes/diseases, genes, etc.
    """

    mim_diseases = {}

    genemap_fh = open(genemap_filename, 'r')
    for line in genemap_fh:  # Loop based on Dima's @ Princeton
        tokens = line.strip('\n').split('\t')

        try:
            mim_geneid = tokens[5].strip()
            entrez_id = tokens[9].strip()
            disorders = tokens[12].strip()
        except IndexError:
            continue

        # Skip line if disorders field is empty
        if disorders == '':
            continue

        # Log message for empty "Entrez Gene ID" column
        if entrez_id == '':
            logger.info(f"Empty Entrez Gene ID for MIM NUmber {mim_geneid}")
            continue

        # Split disorders and handle them one by one
        disorders_list = disorders.split(';')
        for disorder in disorders_list:
            if '[' in disorder or '?' in disorder:
                continue

            # This next line returns a re Match object:
            # It will be None if no match is found.
            mim_info = re.search(FIND_MIMID, disorder)

            if mim_info:
                split_mim_info = mim_info.group(0).split(' ')
                mim_disease_id = split_mim_info[1].strip()
                mim_phenotype = split_mim_info[2].strip()

                # Check if the mim_phenotype number is the one
                # in our filter. If not, skip and continue
                if mim_phenotype != PHENOTYPE_FILTER:
                    continue

                if mim_disease_id not in mim_diseases:
                    mim_diseases[mim_disease_id] = MIMdisease()
                    mim_diseases[mim_disease_id].id = mim_disease_id
                    mim_diseases[mim_disease_id].phenotype = mim_phenotype

                if entrez_id not in mim_diseases[mim_disease_id].genes:
                    mim_diseases[mim_disease_id].genes.append(entrez_id)

    return mim_diseases


def add_do_term_annotations(doid_omim_dict, disease_ontology, mim_diseases):
    """
    Function to add annotations to only the disease_ontology terms found in
    the doid_omim_dict (created by the build_doid_omim_dict() function).

    Arguments:
    doid_omim_dict -- Dictionary mapping DO IDs to OMIM xref IDs. Only DOIDs
    with existing OMIM xref IDs are present as keys in this dictionary.

    disease_ontology -- A Disease Ontology that has parsed the DO OBO file.
    This is actually just a go.go() object (see imports for this file) that
    has parsed a DO OBO file instead of a GO OBO file.

    mim_diseases -- Dictionary of MIM IDs as the keys and MIMdisease
    objects (defined above) as values.

    Returns:
    Nothing, only adds annotations to DO terms.

    """
    logger.debug(disease_ontology.go_terms)

    for doid in doid_omim_dict.keys():
        term = disease_ontology.get_term(doid)

        if term is None:
            continue

        logger.info("Processing %s", term)

        omim_id_list = doid_omim_dict[doid]

        for omim_id in omim_id_list:
            # If omim_id is not in mim_diseases dict, ignore it.
            if omim_id not in mim_diseases:
                continue

            mim_entry = mim_diseases[omim_id]
            for gene_id in mim_entry.genes:
                entrez = int(gene_id)
                term.add_annotation(gid=entrez, ref=None)


def create_term_title(do_term):
    """
    Small function to create the DO term title in the desired
    format: DO-<DO integer ID>:<DO term full name>
    Example: DO-9351:diabetes mellitus

    Arguments:
    do_term -- This is a go_term object from the go() class (go.go)

    Returns:
    title -- A string of the DO term's title in the desired format.
    """
    do_id = do_term.go_id
    do_num = do_id.split(':')[1]
    title = 'DO' + '-' + do_num + ':' + do_term.full_name

    return title


def create_term_abstract(do_term, doid_omim_dict):
    """
    Function to create the DO term abstract in the desired
    format.

    Arguments:
    do_term -- This is a go_term object from the go() class (go.go)

    doid_omim_dict -- A dictionary of DO terms mapping to sets of OMIM xrefs.
    This is returned by the build_doid_omim_dict() function above.

    Returns:
    abstract -- A string of the DO term's abstract in the desired format.
    """
    omim_clause = ''

    doid = do_term.go_id
    if doid in doid_omim_dict:
        # `omim_list` is sorted to make return value reproducible.
        omim_list = sorted(list(doid_omim_dict[doid]))
    else:
        omim_list = []

    if len(omim_list):
        omim_clause = ' Annotations directly to this term are provided ' + \
                        'by the OMIM disease ID'  # Is that sentence right?

        if len(omim_list) == 1:
            omim_clause = omim_clause + ' ' + omim_list[0]
        else:
            omim_clause = omim_clause + 's ' + ', '.join(omim_list[:-1]) + \
                ' and ' + omim_list[-1]
        omim_clause = omim_clause + '.'

    abstract = ''

    if do_term.description:
        abstract += do_term.description

    else:
        logger.info("No OBO description for term %s", do_term)

    abstract += (
        ' Annotations from child terms in the disease ontology are'
        + ' propagated through transitive closure.'
        + omim_clause
    )

    logger.info(abstract)
    return abstract


def process_do_terms():
    """
    Function to read in config INI file and run the other functions to
    process DO terms.
    """

    disease_ontology = go()
    loaded_obo_bool = disease_ontology.load_obo(obo_filename)

    if loaded_obo_bool is False:
        logger.error('DO OBO file could not be loaded.')

    doid_omim_dict = build_doid_omim_dict(obo_filename)

    #mim2entrez_dict = build_mim2entrez_dict(mim2gene_filename)

    mim_diseases = build_mim_diseases_dict(genemap_filename)

    add_do_term_annotations(doid_omim_dict, disease_ontology, mim_diseases)

    disease_ontology.populated = True
    disease_ontology.propagate()

    tags_dictionary = None
    if tag_mapping_filename:
        # The following 4 variables are copied from `human.ini`
        do_id_column = 2
        do_name_column = 3
        tag_column = 1
        header = True

        tags_dictionary = build_tags_dictionary(
            tag_mapping_filename, do_id_column, do_name_column, tag_column, header
        )

    do_terms = []
    for term_id, term in disease_ontology.go_terms.items():
        do_term = {}
        do_term['_id'] = create_term_title(term)
        do_term['is_public'] = True
        do_term['creator'] = 'disease_ontology_parser'
        do_term['date'] = date.today().isoformat()
        do_term['taxid'] = tax_id
        do_term['genes'] = []
        do_term['disease_ontology'] = {
            'id': term_id,
            'abstract': create_term_abstract(term, doid_omim_dict)
        }

        for annotation in term.annotations:
            if annotation.gid not in do_term['genes']:
                do_term['genes'].append(annotation.gid)

        if do_term['genes']:
            if tags_dictionary and term_id in tags_dictionary:
                do_term['disease_ontology']['tags'] = tags_dictionary[term_id]['gs_tags']

            do_term['genes'].sort()  # sort genes to make output reproducible
            do_terms.append(do_term)

    return do_terms


# Test harness
if __name__ == "__main__":
    """
    doid_omim_dict = build_doid_omim_dict(obo_filename)
    #for k, v in doid_omim_dict.items(): print(k, ":", v)

    mim2entrez = build_mim2entrez_dict(mim2gene_filename)
    #for k, v in mim2entrez.items(): print(k, ":", v)

    mim_diseases = build_mim_diseases_dict(genemap_filename, mim2entrez)
    #for k, v in mim_diseases.items(): print(k, ":", v)
    """


    do_terms = process_do_terms()
    print(json.dumps(do_terms, indent=2))

    #print("\nTotal number of gs:", len(do_terms))

    # genemap.txt:  4,192 genesets
    # genemap2.txt: 4,194 genesets
