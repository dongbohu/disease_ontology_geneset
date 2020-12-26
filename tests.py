#!/usr/bin/env python3

"""This script include a few simple tests for local use only."""

import json
import os
import unittest
from parser import get_genesets

class TestResult(unittest.TestCase):
    def test_snapshot(self):
        """20 seconds to run this test."""

        data_dir = "./data/snapshot-for-test"
        obo_filename = os.path.join(data_dir, "2020-12-22.HumanDO.obo")
        genemap_filename = os.path.join(data_dir, "2020-12-23.genemap2.txt")
        genesets = get_genesets(obo_filename, genemap_filename)

        # 4,222 genesets
        self.assertEqual(len(genesets), 4222)

        # Confirm information in the first geneset
        first_gs = genesets[0]
        self.assertEqual(first_gs['_id'], "DO-10124:corneal disease")
        self.assertEqual(first_gs['is_public'], True)
        self.assertEqual(first_gs['creator'], 'disease_ontology_parser')
        self.assertEqual(first_gs['taxid'], 9606)

        genes_in_first_gs = [ g['source'] for g in first_gs['genes'] ]
        self.assertEqual(len(genes_in_first_gs), 21)
        self.assertEqual(genes_in_first_gs[0], '1296')
        self.assertEqual(genes_in_first_gs[-1], '200576')
        self.assertEqual(
            first_gs['disease_ontology']['abstract'],
            "An eye disease that affects the cornea, which is the transparent surface of the eye that assists in light refraction. Annotations from child terms in the disease ontology are propagated through transitive closure."
        )

        # Confirm information in the last geneset
        last_gs = genesets[-1]
        self.assertEqual(last_gs['_id'], "DO-9955:hypoplastic left heart syndrome")
        self.assertEqual(last_gs['is_public'], True)
        self.assertEqual(last_gs['creator'], 'disease_ontology_parser')
        self.assertEqual(last_gs['taxid'], 9606)

        genes_in_last_gs = [ g['source'] for g in last_gs['genes'] ]
        self.assertEqual(genes_in_last_gs, ['1482', '2697'])

        self.assertEqual(
            last_gs['disease_ontology']['abstract'],
            'A congenital heart disease characterized by abnormal development of the left-sided structures of the heart. Annotations from child terms in the disease ontology are propagated through transitive closure. Annotations directly to this term are provided by the OMIM disease IDs 241550 and 614435.'
        )


# Test harness
if __name__ == '__main__':
    unittest.main()


# 2020-12-22:
#   - genemap.txt:  4,192 genesets
#   - genemap2.txt: 4,194 genesets

# 2020-12-23:
#   - 2,959 unique Entrez gene IDs to query
#   - 4,222 genesets
