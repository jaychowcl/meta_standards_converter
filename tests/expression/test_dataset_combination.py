# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
import unittest
import anndata
import pandas
from scipy import sparse
import pytest
from meta_standards_converter.converters.dataset_combination import DatasetCombinationPolicy, DatasetCompatibilityError


def policy():
    def forbidden(*args):
        raise AssertionError("Evidence inspection must not invoke matrix combination collaborators")
    return DatasetCombinationPolicy(scientific_modules=forbidden, attach_sample_values=forbidden,
        package_version=lambda: "6.0.0", metadata_schema_version="1.0")

@pytest.mark.parametrize("allow_unverified", [False, True])
def test_policy_refuses_matrix_combination(allow_unverified):
    with pytest.raises(DatasetCompatibilityError, match="combination is disabled"):
        policy().combine({}, allow_unverified=allow_unverified)

class TestCombinationEvidence(unittest.TestCase):
    anndata = anndata
    pandas = pandas
    sparse = sparse

    def test_all_unknown_samples_are_not_falsely_compatibility_verified(self):
        adatas = {
            sample_id: self.anndata.AnnData(
                X=self.sparse.csr_matrix([[1]]),
                obs=self.pandas.DataFrame(index=[f"{sample_id}-cell"]),
                var=self.pandas.DataFrame(index=["feature:1"]),
            )
            for sample_id in ("GSM1", "GSM2")
        }

        missing = policy().missing_combination_evidence(adatas)

        self.assertEqual(
            {
                "organism": ["GSM1", "GSM2"],
                "reference": ["GSM1", "GSM2"],
                "modality": ["GSM1", "GSM2"],
                "feature_namespace": ["GSM1", "GSM2"],
            },
            missing,
        )


    def test_feature_namespace_distinguishes_entrez_ids_from_gene_symbols(self):
        converter = policy()
        entrez = self.anndata.AnnData(
            X=self.sparse.csr_matrix([[1, 2, 3]]),
            obs=self.pandas.DataFrame(index=["cell"]),
            var=self.pandas.DataFrame(index=["7157", "1956", "7422"]),
        )
        symbols = self.anndata.AnnData(
            X=self.sparse.csr_matrix([[1, 2, 3]]),
            obs=self.pandas.DataFrame(index=["cell"]),
            var=self.pandas.DataFrame(index=["TP53", "EGFR", "VEGFA"]),
        )

        self.assertEqual("entrez", converter.feature_namespace(entrez))
        self.assertEqual("symbol", converter.feature_namespace(symbols))
