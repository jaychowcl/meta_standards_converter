# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Positional IDF accession/source pairs, independent of the semantic overlay."""

VALUE = 'Comment[SecondaryAccession]'
SOURCE = 'Comment[SecondaryAccessionTermSourceRef]'


def _label(row):
    return ''.join(str(row[0]).split()).casefold() if row else ''


def secondary_accession_pairs(rows):
    """Read repeated or multi-column rows without shifting empty source slots.

    A source row belongs to the preceding accession row. A new accession row
    closes that group even when its predecessor omitted the source row entirely.
    """
    pending = None
    pairs = []

    def append(values, sources=()):
        pairs.extend((str(value).strip(), str(sources[i] or '').strip() if i < len(sources) else '')
                     for i, value in enumerate(values) if value is not None and str(value).strip())

    for row in rows:
        label = _label(row)
        if label == VALUE.casefold():
            if pending is not None:
                append(pending)
            pending = row[1:]
        elif label == SOURCE.casefold() and pending is not None:
            append(pending, row[1:])
            pending = None
    if pending is not None:
        append(pending)
    return pairs


def expand_secondary_accession_rows(rows):
    """Final presentation only: overlays still consume a single row per label."""
    labels = {VALUE.casefold(), SOURCE.casefold()}
    pairs = secondary_accession_pairs(rows)
    output = []
    inserted = False
    for row in rows:
        if _label(row) not in labels:
            output.append(row)
        elif not inserted:
            for value, source in pairs:
                output.extend(([VALUE, value], [SOURCE, source]))
            inserted = True
    return output
