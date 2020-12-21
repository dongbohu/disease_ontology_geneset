#!/usr/bin/env python3

import json
import os
import re
import sys

from go import go
from slugify import slugify

import logging
logger = logging.getLogger(__name__)

# Type of OMIM term filter
TYPE_FILTER = set(['gene', 'gene/phenotype'])

# Phenotype filter
PHENO_FILTER = '(3)'

# Confidence filter
# dhu: "confidence" column in `genemap.txt` is deprecated now.
#CONF_FILTER = ['C', 'P']
CONF_FILTER = None

# dhu: The following global variables were defined in `human.ini`
FILTER_BY_TAXONOMY_ID = True

organism = "Homo sapiens"
tax_id = 9606
xrdb = "Entrez"

#do_obo_url = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo"
#mim2gene_url = "https://omim.org/static/omim/data/mim2gene.txt"
#genemap_url = "https://data.omim.org/downloads/z9hwkkLwTHyKrrmsmXkYiQ/genemap.txt"

# dhu: hardcoded paths of input files
obo_filename = "data/HumanDO.obo"
mim2gene_filename = "data/omim/mim2gene.txt"

# dhu: "confidence" column (#7) in "genemap.txt" is deprecated. All values in
#      this column are empty.
#      For differences between "genemap.txt" and "genemap2.txt", see the end of both files.
#      "genemap2.txt" is probably a better choice, because it already maps "MIM Number" to
#      "Entrez Gene ID".
genemap_filename = "data/omim/genemap.txt"

# dhu: ignore "tag_mapping_filename" for now
tag_mapping_filename = None


# This function is copied from `annotation-refinery/utils.py`
def build_tags_dictionary(
        tag_mapping_file,
        geneset_id_column,
        geneset_name_column,
        tag_column,
        header
):
    tags_dict = {}
    tag_file_fh = open(tag_mapping_file, 'r')

    if header:
        tag_file_fh.next()

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


# dhu: This function is copied from "annotation-refinery/process_do.py".
def build_mim2entrez_dict(mim2gene_filename):
    """
    Function to parse mim2gene.txt file and build dictionary of MIM
    numbers to Entrez IDs. The file itself is called mim2gene, and it
    includes Entrez, Symbol and Ensembl identifiers, but we will use
    Entrez.

    Arguments:
    mim2gene_filename -- A string. Location of the mim2gene.txt file to read in.

    Returns:
    mim2entrez_dict -- A dictionary mapping MIM IDs to Entrez IDs for MIM
    Entry Types that pass the TYPE_FILTER. The keys are MIM IDs and the
    values are Entrez IDs.
    """
    mim2entrez_dict = {}

    mim2gene_fh = open(mim2gene_filename, 'r')

    for line in mim2gene_fh:  # Loop based on loop from Dima @ Princeton
        toks = line.split('\t')

        try:  # This is to catch lines that are not in the format we want.
            mim = toks[0]
            mim_type = toks[1]
            entrez_gid = toks[2]
        except IndexError:
            continue

        if mim_type in TYPE_FILTER:
            if entrez_gid == '':
                logger.warning(
                    "Gene Entrez ID was blank for MIM ID '%s' in %s"
                    " mim-to-gene mapping file",
                    mim,
                    mim2gene_filename
                )
                continue
            if mim in mim2entrez_dict:
                logger.warning("MIM already exists in mim2entrez_dict: %s", mim)
            mim2entrez_dict[mim] = entrez_gid
    return mim2entrez_dict


class MIMdisease:
    def __init__(self):
        self.mimid = ''
        self.phe_mm = ''  # Phenotype mapping method
        self.genetuples = []  # (Gene ID, Gene confidence)


def build_mim_diseases_dict(genemap_filename, mim2entrez_dict):
    """
    Function to parse genemap file and build a dictionary of MIM
    diseases.

    Arguments:
    genemap_filename -- A string. Location of the genemap file to read in.

    mim2entrez_dict -- A dictionary mapping MIM gene IDs to Entrez IDs.
    This is the output of the build_mim2entrez_dict() function above.

    Returns:
    mim_diseases -- A dictionary. The keys are MIM disease IDs, and the
    values are `MIMdisease` objects, defined by the class above.

    *N.B. MIM IDs are not all one type of object (unlike Entrez IDs,
    for example) - they can refer to phenotypes/diseases, genes, etc.
    """
    FIND_MIMID = re.compile('\, [0-9]* \([1-4]\)')
    mim_diseases = {}

    genemap_fh = open(genemap_filename, 'r')
    for line in genemap_fh:  # Loop based on Dima's @ Princeton
        # The choice of fields relies on info from the genemap.key
        # file from omim
        toks = line.split('\t')

        try:  # This is to catch lines that are not in the format we want.
            confidence = toks[6].strip()
            mim_geneid = toks[8].strip()
            disorders = toks[11].strip()
        except IndexError:
            continue

        # Skip line if disorders field is empty, or confidence is
        # something other than the one(s) in our filter.
        #
        # dhu: 'confidence' field is deprecated, so its value is always empty now
        #if disorders == '' or confidence not in CONF_FILTER:
        if disorders == '':
            continue

        # A lot of MIM IDs probably will not be in the mim2entrez_dict,
        # since it did not include MIM IDs that are labelled as 'phenotype'
        # (see TYPE_FILTER above)
        if mim_geneid not in mim2entrez_dict:
            logger.warn('Entrez ID for MIM gene ID ' + str(mim_geneid) +
                        ' was not found in mim2entrez_dict')
            continue

        entrezid = mim2entrez_dict[mim_geneid]
        tuple_gid_conf = (entrezid, confidence)

        disorders_list = disorders.split(';')
        for disorder in disorders_list:

            if '[' in disorder or '?' in disorder:
                continue

            # This next line returns a re Match object:
            # It will be None if no match is found.
            mim_info = re.search(FIND_MIMID, disorder)

            if mim_info:
                split_mim_info = mim_info.group(0).split(' ')

                # mim_dis stands for MIM disease ID
                mim_dis_id = split_mim_info[1].strip()
                mim_phetype = split_mim_info[2].strip()

                # Check if the mim_phetype number is the one
                # in our filter. If not, skip and continue
                if mim_phetype != PHENO_FILTER:
                    continue

                if mim_dis_id not in mim_diseases:
                    mim_diseases[mim_dis_id] = MIMdisease()
                    mim_diseases[mim_dis_id].mimid = mim_dis_id
                    mim_diseases[mim_dis_id].phe_mm = mim_phetype

                if tuple_gid_conf not in mim_diseases[mim_dis_id].genetuples:
                    mim_diseases[mim_dis_id].genetuples.append(tuple_gid_conf)

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
            # Ignore if omim_id is not present in mim_diseases dictionary
            if omim_id not in mim_diseases:
                continue

            mim_entry = mim_diseases[omim_id]

            for gene_tuple in mim_entry.genetuples:
                entrez = int(gene_tuple[0])  # The first item is the Entrez ID
                term.add_annotation(gid=entrez, ref=None)


def create_do_term_title(do_term):
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


def create_do_term_abstract(do_term, doid_omim_dict):
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
        omim_list = list(doid_omim_dict[doid])
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

    conf_clause = ''
    if CONF_FILTER and len(CONF_FILTER):
        conf_clause = ' Only annotations with confidence labeled '

        if len(CONF_FILTER) == 1:
            conf_clause = conf_clause + CONF_FILTER[0]
        else:
            conf_clause = conf_clause + ', '.join(CONF_FILTER[:-1]) + \
                ' or ' + CONF_FILTER[-1]

        conf_clause += ' by OMIM have been added.'

    abstract = ''

    if do_term.description:
        abstract += do_term.description

    else:
        logger.info("No OBO description for term %s", do_term)

    abstract += ' Annotations from child terms in the disease ontology are' + \
        ' propagated through transitive closure.' + omim_clause + conf_clause

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

    doid_omim_dict = build_doid_omim_dict(do_obo_filename)

    mim2entrez_dict = build_mim2entrez_dict(mim2gene_filename)

    mim_diseases = build_mim_diseases_dict(genemap_filename, mim2entrez_dict)

    add_do_term_annotations(doid_omim_dict, disease_ontology, mim_diseases)

    disease_ontology.populated = True
    disease_ontology.propagate()

    # dhu: "TAG_MAPPING_FILE" is ignored for now
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
        if len(term.annotations):
            print(f"term_id:", term_id)
            print(f"term.annotations:", term.annotations, "\n")

        do_term = {}
        do_term['title'] = create_do_term_title(term)
        do_term['abstract'] = create_do_term_abstract(term, doid_omim_dict)
        do_term['xrdb'] = xrdb
        do_term['organism'] = organism
        do_term['slug'] = slugify(term_id + '-' + organism)
        do_term['annotations'] = {}

        for annotation in term.annotations:
            if annotation.gid not in do_term['annotations']:
                do_term['annotations'][annotation.gid] = []
            else:
                do_term['annotations'][annotation.gid].append(annotation.ref)

        if do_term['annotations']:
            if tags_dictionary and term_id in tags_dictionary:
                do_term['tags'] = tags_dictionary[term_id]['gs_tags']
            do_terms.append(do_term)

    return do_terms


# Test harness
if __name__ == "__main__":
    #doid_omim_dict = build_doid_omim_dict(obo_filename)
    #for k, v in doid_omim_dict.items(): print(k, ":", v)

    mim2entrez = build_mim2entrez_dict(mim2gene_filename)
    #for k, v in mim2entrez.items(): print(k, ":", v)
    mim_diseases = build_mim_diseases_dict(genemap_filename, mim2entrez)
    for k, v in mim_diseases.items(): print(k, ":", v)


    """
    do_terms = process_do_terms()

    #for gs in do_terms:
    #    print(json.dumps(gs, indent=2))

    print("\nTotal number of gs:", len(do_terms))
    """
