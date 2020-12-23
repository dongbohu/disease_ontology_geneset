#!/usr/bin/env python3

def get_release(self):
    import requests

    DO_OBO_URL = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo"
    release = ""
    resp = requests.get(DO_OBO_URL)
    text_lines = resp.text.strip('\n').split('\n')

    # Find the line that is in the following format:
    # "data-version: doid/releases/YYYY-MM-DD/doid-non-classified.obo"
    # and extract "YYYY-MM-DD" part as the release string.
    for line in text_lines:
        if line.startswith("data-version: "):
            full_version = line.strip().split(' ')[1]
            release = full_version.split('/')[2]
            break

    return release


# Test harness
if __name__ == "__main__":
    print(get_release(None))
