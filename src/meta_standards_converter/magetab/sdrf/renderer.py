# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
from __future__ import annotations
from collections import OrderedDict
from .model import SDRFPath, ColumnGroup, SDRFEdge, SDRFAttr

class SDRFRenderer:
    preserve_order = False

    def __init__(self, *, preserve_order=False):
        self.preserve_order = preserve_order

    def plan_columns(self, paths: list[SDRFPath]) -> list[ColumnGroup]:
        columns = OrderedDict()
        constraints = set()
        for path in paths:
            groups = self.path_groups(path=path)
            for group, values in groups:
                self.merge_column_group(columns=columns, group=group)
            constraints.update((a[0].main_key,b[0].main_key) for a,b in zip(groups,groups[1:]))
        if self.preserve_order:
            import heapq
            position = {key:i for i,key in enumerate(columns)}
            edges = {key:set() for key in columns}
            degree = {key:0 for key in columns}
            for before,after in constraints:
                edges[before].add(after); degree[after] += 1
            ready = [(position[key],key) for key in columns if not degree[key]]
            heapq.heapify(ready); ordered = []
            while ready:
                _,key = heapq.heappop(ready); ordered.append(columns[key])
                for after in edges[key]:
                    degree[after] -= 1
                    if not degree[after]: heapq.heappush(ready,(position[after],after))
            if len(ordered) != len(columns):
                blocked = ', '.join(key for key in columns if degree[key])
                raise ValueError('Incompatible native SDRF path ordering: ' + blocked)
            return ordered
        return list(columns.values())


    def merge_column_group(self, columns: OrderedDict, group: ColumnGroup) -> None:
        if group.main_key not in columns:
            columns[group.main_key] = ColumnGroup(
                main_key=group.main_key,
                main_label=group.main_label,
                companions=[],
            )
        target = columns[group.main_key]
        companion_columns = OrderedDict((x.main_key, x) for x in target.companions)
        for companion in group.companions:
            self.merge_column_group(columns=companion_columns, group=companion)
        target.companions = list(companion_columns.values())


    def render_paths(self, columns: list[ColumnGroup], paths: list[SDRFPath]) -> list:
        header = []
        for column in columns:
            header.extend(self.column_labels(column=column))

        rows = [header]
        for path in paths:
            values = {}
            for group, group_values in self.path_groups(path=path):
                values.update(group_values)
            row = []
            for column in columns:
                row.extend(self.column_values(column=column, values=values))
            rows.append(row)
        return rows


    def column_labels(self, column: ColumnGroup) -> list:
        labels = [column.main_label]
        for companion in column.companions:
            labels.extend(self.column_labels(column=companion))
        return labels


    def column_values(self, column: ColumnGroup, values: dict) -> list:
        rendered = [self.render_value(values.get(column.main_key))]
        for companion in column.companions:
            rendered.extend(self.column_values(column=companion, values=values))
        return rendered


    def path_groups(self, path: SDRFPath) -> list[tuple[ColumnGroup, dict]]:
        groups = []
        counts = {}
        node_keys = []
        if self.preserve_order:
            for part in path.parts:
                node_keys.append(None if isinstance(part, SDRFEdge) else self.occurrence_key(counts=counts, label=part.kind))
            counts = {}
        for index, part in enumerate(path.parts):
            if isinstance(part, SDRFEdge):
                anchor = next((k for k in node_keys[index+1:] if k is not None), 'END') if self.preserve_order else None
                key = self.occurrence_key(counts=counts, label="Protocol REF" + ('@' + anchor if anchor else ''))
                groups.append(self.group_with_values(
                    key=key,
                    label="Protocol REF",
                    value=part.protocol_ref,
                    attrs=part.attrs,
                ))
                continue

            key = node_keys[index] if self.preserve_order else self.occurrence_key(counts=counts, label=part.kind)
            groups.append(self.group_with_values(
                key=key,
                label=part.kind,
                value=part.value,
                attrs=part.attrs,
            ))
        return groups


    def group_with_values(self, key: str, label: str, value, attrs: list[SDRFAttr]) -> tuple[ColumnGroup, dict]:
        companions, values = self.attr_columns(parent_key=key, attrs=attrs)
        values[key] = value
        return ColumnGroup(main_key=key, main_label=label, companions=companions), values


    def attr_columns(self, parent_key: str, attrs: list[SDRFAttr]) -> tuple[list[ColumnGroup], dict]:
        groups = []
        values = {}
        counts = {}
        for attr in attrs:
            if attr.value is None and not attr.attrs and not attr.required:
                continue
            key = f"{parent_key}|{self.occurrence_key(counts=counts, label=attr.label)}"
            companions, companion_values = self.attr_columns(parent_key=key, attrs=attr.attrs)
            groups.append(ColumnGroup(
                main_key=key,
                main_label=attr.label,
                companions=companions,
            ))
            values[key] = attr.value
            values.update(companion_values)
        return groups, values


    def occurrence_key(self, counts: dict, label: str) -> str:
        counts[label] = counts.get(label, 0) + 1
        return f"{label}#{counts[label]}"


    def render_value(self, value):
        if value is None:
            return ""
        return value
