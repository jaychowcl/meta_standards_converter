# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Unified storage class for metadata
'''

class MetaStore:
    def __init__(self):
        pass

    def validate_investigation_metadata(self, investigation_metadata:dict) -> bool:
        '''
        Takes investigation_metadata dict and validates json structure, 
        '''
        assert self._validate_investigation_metadata_structure(investigation_metadata), "investigation_metadata dict does not have correct structure"


        return True

    def _validate_investigation_metadata_structure(self, investigation_metadata:dict) -> bool:
        '''
        Validates that the investigation_metadata dict has the correct structure (keys and value types).
        '''
        
        pass