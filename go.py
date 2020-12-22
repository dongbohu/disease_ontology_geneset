import logging
import re

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class go:
    heads = None
    go_terms = None
    alt_id2std_id = None
    populated = None
    s_orgs = None

    # populate this field if you want to mark this GO as organism specific
    go_organism_tax_id = None

    def __init__(self):
        """
        Initialize data structures for storing the tree.
        """
        self.heads = []
        self.go_terms = {}
        self.alt_id2std_id = {}
        self.populated = False
        self.s_orgs = []

    def load_obo(self, path):
        """Load obo from the defined location. """
        obo_fh = None
        try:
            obo_fh = open(path)
        except IOError:
            logger.error('Could not open %s on the local filesystem.', path)

        if obo_fh is None:
            logger.error('Could not open %s.', path)
            return False

        self.parse(obo_fh)
        return True

    def parse(self, obo_fh):
        """
        Parse the passed obo handle.
        """
        inside = False
        gterm = None
        for line in obo_fh:
            fields = line.rstrip().split()

            if len(fields) < 1:
                continue
            elif fields[0] == '[Term]':
                if gterm:
                    if gterm.head:
                        self.heads.append(gterm)
                inside = True
            elif fields[0] == '[Typedef]':
                if gterm:
                    if gterm.head:
                        self.heads.append(gterm)
                inside = False

            elif inside and fields[0] == 'id:':
                if fields[1] in self.go_terms:
                    logger.debug("Term %s exists in go()", fields[1])
                    gterm = self.go_terms[fields[1]]
                else:
                    logger.debug("Adding term %s to go()", fields[1])
                    gterm = GOTerm(fields[1])
                    self.go_terms[gterm.get_id()] = gterm
            elif inside and fields[0] == 'def:':
                desc = ' '.join(fields[1:])
                desc = desc.split('"')[1]
                gterm.description = desc
            elif inside and fields[0] == 'name:':
                fields.pop(0)
                name = '_'.join(fields)
                name = re.sub('[^\w\s_-]', '_', name).strip().lower()
                name = re.sub('[-\s_]+', '_', name)
                gterm.name = name
                gterm.full_name = ' '.join(fields)
            elif inside and fields[0] == 'namespace:':
                gterm.namespace = fields[1]
            elif inside and fields[0] == 'def:':
                gterm.desc = ' '.join(fields[1:]).split('\"')[1]
            elif inside and fields[0] == 'alt_id:':
                gterm.alt_id.append(fields[1])
                self.alt_id2std_id[fields[1]] = gterm.get_id()
            elif inside and fields[0] == 'is_a:':
                logger.debug("Making term.head for term %s = False", gterm)
                gterm.head = False
                fields.pop(0)
                pgo_id = fields.pop(0)
                if pgo_id not in self.go_terms:
                    self.go_terms[pgo_id] = GOTerm(pgo_id)

                gterm.is_a.append(self.go_terms[pgo_id])
                self.go_terms[pgo_id].parent_of.add(gterm)
                gterm.child_of.add(self.go_terms[pgo_id])
            elif inside and fields[0] == 'relationship:':
                if fields[1].find('has_part') != -1:
                    # Has part is not a parental relationship --
                    # it is actually for children.
                    continue
                logger.debug("Making term.head for term %s = False", gterm)
                gterm.head = False
                pgo_id = fields[2]
                if pgo_id not in self.go_terms:
                    self.go_terms[pgo_id] = GOTerm(pgo_id)
                # Check which relationship you are with this parent go term
                if (fields[1] == 'regulates' or
                        fields[1] == 'positively_regulates' or
                        fields[1] == 'negatively_regulates'):
                    gterm.relationship_regulates.append(self.go_terms[pgo_id])
                elif fields[1] == 'part_of':
                    gterm.relationship_part_of.append(self.go_terms[pgo_id])
                else:
                    logger.info("Unkown relationship %s",
                                self.go_terms[pgo_id].name)

                self.go_terms[pgo_id].parent_of.add(gterm)
                gterm.child_of.add(self.go_terms[pgo_id])
            elif inside and fields[0] == 'is_obsolete:':
                logger.debug("Making term.head for term %s = False", gterm)
                gterm.head = False
                del self.go_terms[gterm.get_id()]

        # This loop checks that all terms that have been marked as head=True
        # have been added to self.heads
        for term_id, term in self.go_terms.items():
            if term.head:
                if term not in self.heads:
                    logger.debug("Term %s not in self.heads, adding now", term)
                    self.heads.append(term)

        logger.debug("Terms that are heads: %s", self.heads)

    def propagate(self):
        """
        propagate all gene annotations
        """
        logger.info("Propagate gene annotations")
        logger.debug("Head term(s) = %s", self.heads)
        for head_gterm in self.heads:
            logger.info("Propagating %s", head_gterm.name)
            self.propagate_recurse(head_gterm)

    def propagate_recurse(self, gterm):
        if not len(gterm.parent_of):
            logger.debug("Base case with term %s", gterm.name)
            return

        for child_term in gterm.parent_of:
            self.propagate_recurse(child_term)
            new_annotations = set()

            regulates_relation = (gterm in child_term.relationship_regulates)
            part_of_relation = (gterm in child_term.relationship_part_of)

            for annotation in child_term.annotations:
                copied_annotation = None
                # If this relation with child is a regulates(and its sub class)
                # filter annotations
                if regulates_relation:
                    # only add annotations that didn't come from a part of or
                    # regulates relationship
                    if annotation.ready_regulates_cutoff:
                        continue
                    else:
                        copied_annotation = annotation.prop_copy(
                            ready_regulates_cutoff=True)
                elif part_of_relation:
                    copied_annotation = annotation.prop_copy(
                        ready_regulates_cutoff=True)
                else:
                    copied_annotation = annotation.prop_copy()

                new_annotations.add(copied_annotation)
            gterm.annotations = gterm.annotations | new_annotations

    def get_term(self, tid):
        logger.debug('get_term: %s', tid)
        term = None
        try:
            term = self.go_terms[tid]
        except KeyError:
            try:
                term = self.go_terms[self.alt_id2std_id[tid]]
            except KeyError:
                logger.error('Term name does not exist: %s', tid)
        return term


class Annotation(object):
    def __init__(self, xdb=None, gid=None, ref=None, evidence=None, date=None,
                 direct=False, cross_annotated=False, origin=None,
                 ortho_evidence=None, ready_regulates_cutoff=False):
        super(Annotation, self).__setattr__('xdb', xdb)
        super(Annotation, self).__setattr__('gid', gid)
        super(Annotation, self).__setattr__('ref', ref)
        super(Annotation, self).__setattr__('evidence', evidence)
        super(Annotation, self).__setattr__('date', date)
        super(Annotation, self).__setattr__('direct', direct)
        super(Annotation, self).__setattr__('cross_annotated', cross_annotated)
        super(Annotation, self).__setattr__('origin', origin)
        super(Annotation, self).__setattr__('ortho_evidence', ortho_evidence)
        super(Annotation, self).__setattr__('ready_regulates_cutoff',
                                            ready_regulates_cutoff)

    def prop_copy(self, ready_regulates_cutoff=None):
        if ready_regulates_cutoff is None:
            ready_regulates_cutoff = self.ready_regulates_cutoff

        return Annotation(xdb=self.xdb, gid=self.gid, ref=self.ref,
                          evidence=self.evidence, date=self.date, direct=False,
                          cross_annotated=False,
                          ortho_evidence=self.ortho_evidence,
                          ready_regulates_cutoff=ready_regulates_cutoff)

    def __hash__(self):
        return hash((self.xdb, self.gid, self.ref, self.evidence, self.date,
                     self.direct, self.cross_annotated, self.ortho_evidence,
                     self.ready_regulates_cutoff, self.origin))

    def __eq__(self, other):
        return (self.xdb, self.gid, self.ref, self.evidence, self.date,
                self.direct, self.cross_annotated, self.ortho_evidence,
                self.ready_regulates_cutoff, self.origin).__eq__((
                    other.xdb, other.gid, other.ref, other.evidence,
                    other.date, other.direct, other.cross_annotated,
                    other.ortho_evidence, other.ready_regulates_cutoff,
                    other.origin))

    def __setattr__(self, *args):
        raise TypeError("Attempt to modify immutable object.")
    __delattr__ = __setattr__


class GOTerm:
    go_id = ''
    is_a = None
    relationship = None
    parent_of = None
    child_of = None
    annotations = None
    alt_id = None
    namespace = ''
    included_in_all = None
    valid_go_term = None
    cross_annotated_genes = None
    head = None
    name = None
    full_name = None
    description = None
    base_counts = None
    counts = None
    summary = None
    desc = None
    votes = None

    def __init__(self, go_id):
        self.head = True
        self.go_id = go_id
        self.annotations = set([])
        self.cross_annotated_genes = set([])
        self.is_a = []
        self.relationship_regulates = []
        self.relationship_part_of = []
        self.parent_of = set()
        self.child_of = set()
        self.alt_id = []
        self.included_in_all = True
        self.valid_go_term = True
        self.name = None
        self.full_name = None
        self.description = None
        self.base_counts = None
        self.counts = None
        self.desc = None
        self.votes = set([])

    def __cmp__(self, other):
        return cmp(self.go_id, other.go_id)

    def __hash__(self):
        return(self.go_id.__hash__())

    def __repr__(self):
        return(self.go_id + ': ' + self.name)

    def get_id(self):
        return self.go_id

    def map_genes(self, id_name):
        mapped_annotations_set = set([])
        for annotation in self.annotations:
            mapped_genes = id_name.get(annotation.gid)
            if mapped_genes is None:
                logger.warning('No matching gene id: %s', annotation.gid)
                continue
            for mgene in mapped_genes:
                mapped_annotations_set.add(
                    Annotation(
                        xdb=None,
                        gid=mgene,
                        direct=annotation.direct,
                        ref=annotation.ref,
                        evidence=annotation.evidence,
                        date=annotation.date,
                        cross_annotated=annotation.cross_annotated
                    )
                )
        self.annotations = mapped_annotations_set

    def get_annotated_genes(self, include_cross_annotated=True):
        genes = []
        for annotation in self.annotations:
            if (not include_cross_annotated) and annotation.cross_annotated:
                continue
            genes.append(annotation.gid)
        return genes

    def add_annotation(
            self,
            gid,
            ref=None,
            cross_annotated=False,
            allow_duplicate_gid=True,
            origin=None,
            ortho_evidence=None
    ):
        if not allow_duplicate_gid:
            for annotated in self.annotations:
                if annotated.gid == gid:
                    return
        self.annotations.add(
            Annotation(
                gid=gid,
                ref=ref,
                cross_annotated=cross_annotated,
                origin=origin,
                ortho_evidence=ortho_evidence
            )
        )

    def get_annotation_size(self):
        return len(self.annotations)

    def get_namespace(self):
        return self.namespace
