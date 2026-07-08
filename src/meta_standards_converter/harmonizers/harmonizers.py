# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
All harmonizer classes
'''

from .pubmed2ols import Pubmed2OLS
from .geo2ols import GEO2OLS

class Harmonizer(Pubmed2OLS, GEO2OLS):
    def __init__(self):
        self.ontologies = {}
        super().__init__()
        pass

    
