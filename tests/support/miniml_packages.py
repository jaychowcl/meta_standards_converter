# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Shared builders for miniml packages."""

from __future__ import annotations


def package_v3() -> dict:
    return {
        'miniml_schema_version': '3.0',
        'source': {
            'format': 'MAGE-TAB',
            'version': '1.1',
            'documents': [
                {
                    'kind': 'idf',
                    'name': 'study.idf.txt',
                },
                {
                    'kind': 'sdrf',
                    'name': 'study.sdrf.txt',
                },
            ],
        },
        'database': [],
        'organization': [],
        'contributor': [
            {
                'iid': 'person-1',
                'person': {
                    'first': 'Ada',
                    'last': 'Lovelace',
                },
                'roles': [
                    {
                        'value': 'investigator',
                        'term_source_ref': 'EFO',
                        'term_accession_number': 'EFO:0001739',
                    },
                ],
            },
        ],
        'platform': [],
        'sample': [
            {
                'iid': 'sample-1',
                'channel': [
                    {
                        'source': {
                            'value': 'lung',
                        },
                        'characteristics': [
                            {
                                'name': 'age',
                                'value': '10',
                                'unit': {
                                    'value': 'year',
                                    'term_source_ref': 'UO',
                                    'term_accession_number': 'UO:0000036',
                                    'hz_unit': 'year',
                                    'hz_unit_id': 'UO:0000036',
                                    'hz_unit_onto': 'UO',
                                    'hz_unit_hierarchy_depth': 0,
                                },
                            },
                        ],
                    },
                ],
            },
        ],
        'series': {
            'iid': 'E-MTAB-1',
            'sample_ref': [
                {
                    'ref': 'sample-1',
                },
            ],
            'experiment_date': '2025-01-01',
            'protocols': [
                {
                    'name': 'P-extract',
                    'type': {
                        'value': 'nucleic acid extraction',
                    },
                    'parameters': [
                        'duration',
                    ],
                },
                {
                    'name': 'P-library',
                    'type': {
                        'value': 'library preparation',
                    },
                },
            ],
            'assay_paths': [
                {
                    'document': 'study.sdrf.txt',
                    'steps': [
                        {
                            'kind': 'sample',
                            'name': 'sample-1',
                            'sample_ref': 'sample-1',
                        },
                        {
                            'kind': 'protocol_application',
                            'protocol_ref': 'P-extract',
                            'parameter_values': [
                                {
                                    'name': 'duration',
                                    'value': '30',
                                    'unit': {
                                        'value': 'minute',
                                    },
                                },
                            ],
                        },
                        {
                            'kind': 'protocol_application',
                            'protocol_ref': 'P-library',
                        },
                        {
                            'kind': 'assay',
                            'name': 'assay-1',
                        },
                    ],
                },
            ],
            'quality_controls': [
                {
                    'value': 'biological replicate',
                },
            ],
            'replicate_types': [
                {
                    'value': 'technical replicate',
                },
            ],
            'normalization_types': [
                {
                    'value': 'quantile normalization',
                },
            ],
        },
    }
