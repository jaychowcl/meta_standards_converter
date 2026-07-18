# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
GEO to ontologies
'''

class GEO2OLS:
    def __init__(self):
        super().__init__()
        if not hasattr(self, "ontologies"):
            self.ontologies = {}
        self.ontologies["EFO"] = {"Term Source Name": "EFO", "Term Source File": "http://www.ebi.ac.uk/efo/efo.owl", "Term Source Version": "3.90.0"}
        self.ontologies["OBI"] = {"Term Source Name": "OBI", "Term Source File": "http://purl.obolibrary.org/obo/obi.owl", "Term Source Version": "2026-03-19"}
        pass

    def geoprotocols2efo(self, protocol_type: str) -> list:
       '''
       take GEO protocol types and map to EFO protocol type term. Return efo term, source ref, accession number. 
       '''
       if not protocol_type:
        raise ValueError("data not given")
       
       onto_terms = {
          "sample treatment protocol" : ["sample treatment protocol", "EFO", "EFO_0003809"],
          "sample collection protocol" : ["sample collection protocol", "EFO", "EFO_0005518"],
          "nucleic acid extraction protocol" : ["nucleic acid extraction protocol", "EFO", "EFO_0002944"],
          "nucleic acid sequencing protocol" : ["nucleic acid sequencing protocol", "EFO", "EFO_0004170"],
          "growth protocol" : ["growth protocol", "EFO", "EFO_0003789"],
          "labelling protocol" : ["labelling protocol", "EFO", "EFO_0003808"],
          "hybridization protocol" : ["hybridization protocol", "EFO", "EFO_0003790"],
          "array scanning and feature extraction protocol" : ["array scanning and feature extraction protocol", "EFO", "EFO_0003814"],
          "normalization data transformation protocol" : ["normalization data transformation protocol", "EFO", "EFO_0003816"],
          "protocol" : ["protocol", "OBI", "OBI_0000272"],
          }

       protocol_type_mapping = {
         "Manufacture-Protocol": onto_terms["protocol"],
         "Sample-Collection-Protocol": onto_terms["sample collection protocol"],
         "Treatment-Protocol": onto_terms["sample treatment protocol"],
         "Growth-Protocol": onto_terms["growth protocol"],
         "Extract-Protocol": onto_terms["nucleic acid extraction protocol"],
         "Library-Construction-Protocol": ["nucleic acid library construction protocol", "EFO", "EFO_0004184"],
         "Label-Protocol": onto_terms["labelling protocol"],
         "Hybridization-Protocol": onto_terms["hybridization protocol"],
         "Scan-Protocol": onto_terms["array scanning and feature extraction protocol"],
         "Data-Processing": onto_terms["normalization data transformation protocol"],
         "Nucleic-Acid-Sequencing-Protocol": onto_terms["nucleic acid sequencing protocol"],
         }
       
       return protocol_type_mapping.get(protocol_type, [protocol_type, None, None])
       
