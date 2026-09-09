#!/usr/bin/env python3
"""Reecrit la ligne binding: d'un YAML BoltzGen avec les hotspots convertis.
Usage: script <yaml> <nouveaux_hotspots_csv>"""
import sys, re

yaml_path, new_hot = sys.argv[1], sys.argv[2]
s = open(yaml_path).read()

# remplacer la valeur apres "binding:" (premiere occurrence)
new_s, n = re.subn(r'(binding:\s*)[0-9,]+', r'\g<1>' + new_hot, s, count=1)
if n != 1:
    print(f"ERREUR: {n} remplacements dans {yaml_path}", file=sys.stderr)
    sys.exit(1)
open(yaml_path, "w").write(new_s)
print(f"  OK {yaml_path} -> binding: {new_hot}")
