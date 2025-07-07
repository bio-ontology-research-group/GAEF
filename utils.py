from collections import deque, Counter
import warnings
import pandas as pd
import numpy as np
from xml.etree import ElementTree as ET
import math

BIOLOGICAL_PROCESS = 'GO:0008150'
MOLECULAR_FUNCTION = 'GO:0003674'
CELLULAR_COMPONENT = 'GO:0005575'
HAS_FUNCTION = 'http://mowl.borg/has_function'

FUNC_DICT = {
    'cc': CELLULAR_COMPONENT,
    'mf': MOLECULAR_FUNCTION,
    'bp': BIOLOGICAL_PROCESS}

NAMESPACES = {
    'cc': 'cellular_component',
    'mf': 'molecular_function',
    'bp': 'biological_process'
}


class Ontology(object):

    def __init__(self, filename='data/go-basic.obo', with_rels=True):
        self.ont = self.load(filename, with_rels)
        self.ic = None
        self.ic_norm = 0.0
        self.ancestors = {}

    def has_term(self, term_id):
        return term_id in self.ont

    def get_term(self, term_id):
        if self.has_term(term_id):
            return self.ont[term_id]
        return None

    def load(self, filename, with_rels):
        ont = dict()
        obj = None
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line == '[Term]':
                    if obj is not None:
                        ont[obj['id']] = obj
                    obj = dict()
                    obj['is_a'] = list()
                    obj['part_of'] = list()
                    obj['regulates'] = list()
                    obj['alt_ids'] = list()
                    obj['is_obsolete'] = False
                    continue
                elif line == '[Typedef]':
                    if obj is not None:
                        ont[obj['id']] = obj
                    obj = None
                else:
                    if obj is None:
                        continue
                    l = line.split(": ")
                    if l[0] == 'id':
                        obj['id'] = l[1]
                    elif l[0] == 'alt_id':
                        obj['alt_ids'].append(l[1])
                    elif l[0] == 'namespace':
                        obj['namespace'] = l[1]
                    elif l[0] == 'is_a':
                        obj['is_a'].append(l[1].split(' ! ')[0])
                    elif with_rels and l[0] == 'relationship':
                        it = l[1].split()
                        # add all types of relationships
                        obj['is_a'].append(it[1])
                    elif l[0] == 'name':
                        obj['name'] = l[1]
                    elif l[0] == 'is_obsolete' and l[1] == 'true':
                        obj['is_obsolete'] = True
            if obj is not None:
                ont[obj['id']] = obj
        for term_id in list(ont.keys()):
            for t_id in ont[term_id]['alt_ids']:
                ont[t_id] = ont[term_id]
            if ont[term_id]['is_obsolete']:
                del ont[term_id]
        for term_id, val in ont.items():
            if 'children' not in val:
                val['children'] = set()
            for p_id in val['is_a']:
                if p_id in ont:
                    if 'children' not in ont[p_id]:
                        ont[p_id]['children'] = set()
                    ont[p_id]['children'].add(term_id)
     
        return ont

    def get_ancestors(self, term_id):
        if term_id not in self.ont:
            return set()
        if term_id in self.ancestors:
            return self.ancestors[term_id]
        term_set = set()
        q = deque()
        q.append(term_id)
        while(len(q) > 0):
            t_id = q.popleft()
            if t_id not in term_set:
                term_set.add(t_id)
                for parent_id in self.ont[t_id]['is_a']:
                    if parent_id in self.ont:
                        q.append(parent_id)
        self.ancestors[term_id] = term_set
        return term_set

    def get_prop_terms(self, terms):
        prop_terms = set()

        for term_id in terms:
            prop_terms |= self.get_anchestors(term_id)
        return prop_terms


    def get_parents(self, term_id):
        if term_id not in self.ont:
            return set()
        term_set = set()
        for parent_id in self.ont[term_id]['is_a']:
            if parent_id in self.ont:
                term_set.add(parent_id)
        return term_set


    def get_namespace_terms(self, namespace):
        terms = set()
        for go_id, obj in self.ont.items():
            if obj['namespace'] == namespace:
                terms.add(go_id)
        return terms

    def get_namespace(self, term_id):
        return self.ont[term_id]['namespace']
    
    def get_term_set(self, term_id):
        if term_id not in self.ont:
            return set()
        term_set = set()
        q = deque()
        q.append(term_id)
        while len(q) > 0:
            t_id = q.popleft()
            if t_id not in term_set:
                term_set.add(t_id)
                for ch_id in self.ont[t_id]['children']:
                    q.append(ch_id)
        return term_set



#### GET ANCESTORS ####

def get_ancestors(protein_go_terms, obo_file='data/go-basic.obo'):
    """
    Expands GO terms to include their ancestors for each protein.

    Parameters:
    - protein_go_terms (dict): {protein_id: set(GO_terms)}
    - obo_file (str): Path to GO ontology file

    Returns:
    - dict: {protein_id: set(GO_terms + ancestors)}
    """
    go = Ontology(obo_file, with_rels=True)
    result = {}
    for protein_id, terms in protein_go_terms.items():
        expanded = set()
        for term in terms:
            expanded.update(go.get_ancestors(term))
        result[protein_id] = expanded
    return result

#### GET SPECIFIC TERMS ####
def get_specific(protein_go_terms, obo_file='data/go-basic.obo'):
    """
    Filters GO terms to retain only the most specific (non-ancestor) terms.

    Parameters:
    - protein_go_terms (dict): {protein_id: set(GO_terms)}
    - obo_file (str): Path to GO ontology file

    Returns:
    - dict: {protein_id: set(most specific GO_terms)}
    """
    go = Ontology(obo_file, with_rels=True)
    result = {}
    for protein_id, terms in protein_go_terms.items():
        specific = set(terms)
        for term in terms:
            ancestors = go.get_ancestors(term)
            ancestors.discard(term)
            specific -= ancestors
        result[protein_id] = specific
    return result
