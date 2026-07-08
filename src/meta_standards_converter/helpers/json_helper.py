# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
'''
Helper class for JSONs
'''

class JSONHandler:
    def _from_path(self, obj, path_str):
        def parse_path(path):
            return [
                "*" if p == "*" else int(p) if p.isdigit() else p
                for p in path.split(".")
            ]

        def resolve(current, path):
            if not path:
                return [current]

            key, *rest = path

            if current is None:
                return [None]

            if key == "*":
                if not isinstance(current, list):
                    return [None]

                result = []
                for item in current:
                    resolved = resolve(item, rest)

                    if resolved:
                        result.extend(resolved)
                    else:
                        result.append(None)

                return result

            try:
                return resolve(current[key], rest)
            except (KeyError, IndexError, TypeError):
                return [None]

        result = resolve(obj, parse_path(path_str))

        return result

    def _flatten_values(self, value):
        if isinstance(value, list):
            return [
                item
                for element in value
                for item in self._flatten_values(element)
            ]

        if isinstance(value, dict):
            return [
                item
                for element in value.values()
                for item in self._flatten_values(element)
            ]

        return [value]