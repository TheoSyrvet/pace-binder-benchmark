#!/usr/bin/env python3
"""
Convertit les outputs SolMPNN en FASTA ColabFold
Simple : remplace / par : comme l'ancien prepare_inputs.py
ColabFold gère nativement les X et séquences multi-chaînes
"""
import os, argparse

def parse_mpnn_fa(fa_path):
    sequences = []
    name = None
    with open(fa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif line and name and ("sample=" in name or "T=0.1" in name):
                # Remplacer / par : — ColabFold format
                seq = line.replace("/", ":")
                sequences.append((name, seq))
                name = None
    return sequences

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpnn_dir",    required=True)
    parser.add_argument("--output_dir",  required=True)
    parser.add_argument("--target",      required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    n_generated = n_skipped = n_errors = 0

    design_dirs = sorted([d for d in os.listdir(args.mpnn_dir)
                          if os.path.isdir(os.path.join(args.mpnn_dir, d))])

    for design in design_dirs:
        fa_path = os.path.join(args.mpnn_dir, design, "seqs", f"{design}.fa")
        if not os.path.exists(fa_path):
            n_errors += 1
            continue

        sequences = parse_mpnn_fa(fa_path)

        for i, (header, seq) in enumerate(sequences[:8], 1):
            fasta_name = f"{args.target}_{design}_sample{i}"
            out_path   = os.path.join(args.output_dir, f"{fasta_name}.fasta")

            if os.path.exists(out_path):
                n_skipped += 1
                continue

            with open(out_path, "w") as f:
                f.write(f">{fasta_name}\n")
                f.write(f"{seq}\n")
            n_generated += 1

    print(f"[{args.target}] générés={n_generated} skipped={n_skipped} erreurs={n_errors}")
    print(f"[{args.target}] Total: {len([f for f in os.listdir(args.output_dir) if f.endswith('.fasta')])}")

if __name__ == "__main__":
    main()
