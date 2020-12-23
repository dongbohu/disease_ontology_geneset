#!/bin/bash

# This script automatically updates `data/HumanDO.obo` and `data/omim/genemap2.txt`.

FORMAT_STR="############"

# Update `./data/HumanDO.obo`
echo
echo "${FORMAT_STR} Updating HumanDO.obo ... ${FORMAT_STR}"
wget "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo"

if [ "$?" -gt 0 ]; then
    echo "${FORMAT_STR} ERROR: Failed to download HumanDO.obo ${FORMAT_STR}"
elif [ "$(diff -q ./HumanDO.obo ./data/HumanDO.obo)" ]; then
    mv -f ./HumanDO.obo ./data
    git ci -m "Update HumanDO.obo" ./data/HumanDO.obo
    echo "${FORMAT_STR} HumanDO.obo updated ${FORMAT_STR}"
else
    echo "${FORMAT_STR} No need to update HumanDO.obo ${FORMAT_STR}"
    rm -f ./HumanDO.obo
fi

# Update `./data/omim/genemap2.txt`
echo
echo "${FORMAT_STR} Updating genemap2.txt ... ${FORMAT_STR}"
wget "https://data.omim.org/downloads/z9hwkkLwTHyKrrmsmXkYiQ/genemap2.txt"

if [ "$?" -gt 0 ]; then
    echo "${FORMAT_STR} ERROR: Failed to download genemap2.txt ${FORMAT_STR}"
elif [ $(diff -q ./genemap2.txt ./data/omim/genemap2.txt) ]; then
    mv -f ./genemap2.txt ./data/omim/
    git ci -m "Update genemap2.txt" ./data/omim/genemap2.txt
    echo "${FORMAT_STR} genemap2.txt updated ${FORMAT_STR}"
else
    echo "${FORMAT_STR} No need to update genemap2.txt ${FORMAT_STR}"
    rm -f ./genemap2.txt
fi
