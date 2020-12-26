## This directory includes the data files used to generate disease ontology genesets.

1. `latest/` includes the latest `HumanDO.obo` and `genemap2.txt` file.
   - `HumanDO.obo` is downloaded from [here]
(https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.obo)
   - `genemap2.txt` is downloaded from [here](https://omim.org/downloads) (registration required)

2. `snapshot-for-test/` includes a few data files for tests only.

3. Note that `confidence` column (#7) in `genemap.txt` is deprecated. All values in
   this column are empty. `genemap.txt` is now replaced by `genemap2.txt`. The latter
   also maps `MIM Number` to `Entrez Gene ID`, so `mim2gene.txt` is not needed either.

   The differences between `genemap.txt` and `genemap2.txt` are described at the end
   of both files.

4. When `https://data.omim.org/downloads/<your_download_key>/genemap2.txt` is downloaded
   more than 10 times a day, the response will become:
     > This data account exceeded its download cap, please contact us at
     > https://omim.org/contact if this is an issue

   OMIM maintainers confirmed that the download cap is 10 times per day, and
   `genemap2.txt` is updated daily. So there is no reason to download it more than
   once per day.
