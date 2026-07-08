# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Pubmed to ontologies
'''

class Pubmed2OLS:
    def __init__(self):
      super().__init__()
      if not hasattr(self, "ontologies"):
         self.ontologies = {}
      self.ontologies["EFO"] = {"Term Source Name": "EFO", "Term Source File": "http://www.ebi.ac.uk/efo/efo.owl", "Term Source Version": "3.89.0"}
      self.ontologies["MeSH"] = {"Term Source Name": "MeSH", "Term Source File": "https://id.nlm.nih.gov/mesh/", "Term Source Version": "2024-08-09"}
    
    def pubstatus2efo(self, pub_status: str) -> list:
       '''
       take pubmed publication status and map to EFO publication status term. Return efo term, source ref, accession number. 
       '''
       if not pub_status:
        return [None, None, None]
       
       statuses = [status.strip() for status in pub_status.split("+")]
       primary_status = statuses[0]

       onto_terms = {
          "published": ["published", "EFO", "EFO_0001796"],
          "preprint": ["preprint", "EFO", "EFO_0010558"],
          "submitted": ["submitted", "EFO", "EFO_0001794"],
          "in preparation": ["in preparation", "EFO", "EFO_0001795"],
          "retracted": ["Retracted Publication", "MeSH", "D016441"],
       }
       
       pub_status_mapping = {
          "ppublish": onto_terms["published"], 
          "epublish" : onto_terms["published"],
          "aheadofprint" : onto_terms["published"],
          "retracted" : onto_terms["retracted"],
          "pmc" : onto_terms["published"],
          "pmcr": onto_terms["published"],
          "pubmed": onto_terms["published"],
          "medline": onto_terms["published"],
          "premedline": onto_terms["published"],
          "publisher": onto_terms["submitted"],
          "inprocess": onto_terms["published"],
          "entrez": onto_terms["submitted"],
       }

       return pub_status_mapping.get(primary_status.lower(), [primary_status, None, None])
