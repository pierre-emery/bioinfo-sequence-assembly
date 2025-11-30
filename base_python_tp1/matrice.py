"""
matrice.py — Génère une matrice 20x20 des scores de chevauchement suffixe(X_i) -> préfixe(X_j)

Ce script lit 20 reads à partir d’un fichier FASTQ et calcule, pour chaque paire ordonnée (i, j),
le score optimal d’un alignement de chevauchement de type « suffixe(X_i) -> préfixe(X_j) ».
Le résultat est écrit dans un CSV (par défaut matrice_20x20.csv) avec en-tête des IDs.

Usage :
  python matrice.py reads.fq
  # options :
  #   --out matrice_20x20.csv
  #   --match 4 --mismatch -4 --gap -8
"""

import argparse
import csv
from utils import read_sequences_from_fastq
from prefixe_suffixe import construire_tables_dp, meilleur_score_derniere_ligne


def overlap_score_only(X_i: str, X_j: str, match: int = 4, mismatch: int = -4, gap: int = -8) -> int:
    """
    Une fonction qui prend en argument deux reads (X_i, X_j) et trois scores (match, mismatch, gap).

    Elle construit les tables de programmation dynamique via construire_tables_dp pour le
    chevauchement suffixe(X_i) -> préfixe(X_j), puis récupère le meilleur score sur la
    dernière ligne à l’aide de meilleur_score_derniere_ligne, sans reconstruire l’alignement.

    Rend un entier qui correspond au score optimal de chevauchement.
    """
    V, _ = construire_tables_dp(X_i, X_j, match, mismatch, gap)
    meilleur, _ = meilleur_score_derniere_ligne(V)
    return int(meilleur)


def main():
    parser = argparse.ArgumentParser(description="Matrice 20x20 des scores de chevauchement (suffixe->préfixe)")
    parser.add_argument("fastq", help="reads.fq (20 reads, FASTQ)")
    parser.add_argument("--out", default="matrice_20x20.csv")
    parser.add_argument("--match", type=int, default=4)
    parser.add_argument("--mismatch", type=int, default=-4)
    parser.add_argument("--gap", type=int, default=-8)
    args = parser.parse_args()

    reads = read_sequences_from_fastq(args.fastq)
    ids = list(reads.keys())
    if len(ids) != 20:
        raise ValueError(f"reads.fq doit contenir 20 reads, trouvé {len(ids)}.")

    n = len(ids)
    M = [[0] * n for _ in range(n)]

    for i, id_i in enumerate(ids):
        X_i = reads[id_i]
        for j, id_j in enumerate(ids):
            if i == j:
                M[i][j] = 0
            else:
                X_j = reads[id_j]
                M[i][j] = overlap_score_only(
                    X_i, X_j,
                    match=args.match, mismatch=args.mismatch, gap=args.gap
                )

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["id"] + ids)
        for i, id_i in enumerate(ids):
            writer.writerow([id_i] + M[i])

    print(f"OK — Matrice de chevauchement écrite dans {args.out}")


if __name__ == "__main__":
    main()
