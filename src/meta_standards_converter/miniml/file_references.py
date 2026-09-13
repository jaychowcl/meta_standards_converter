"""Provider file absence markers, distinct from biological missing values."""


def is_file_placeholder(value):
    return isinstance(value, str) and value.strip().casefold() == 'none'


def clean_native_file_placeholders(data):
    """Remove legacy phantom file projections on the private native export copy."""
    if data.get('source', {}).get('format') not in {'ENA', 'SRA'}:
        return
    for entity in [data['series'], *data.get('sample', [])]:
        for field in ('raw_data', 'supplementary_data'):
            if field in entity:
                entity[field] = [v for v in entity[field] if not is_file_placeholder(v.get('value'))]
    paths = []
    for path in data['series'].get('assay_paths', []):
        steps = path.get('steps', [])
        cleaned = []
        for step in steps:
            if step.get('kind') in {'array_data_file', 'derived_array_data_file'}:
                value = step.get('link', {}).get('value') or step.get('name')
                if is_file_placeholder(value):
                    continue
            cleaned.append(step)
        # A source plus an absent processed file is not another assay.
        if len(cleaned) != len(steps) and not any(s.get('kind') in
                {'assay', 'scan', 'array_data_file', 'derived_array_data_file'} for s in cleaned):
            continue
        path['steps'] = cleaned
        paths.append(path)
    if 'assay_paths' in data['series']:
        data['series']['assay_paths'] = paths
