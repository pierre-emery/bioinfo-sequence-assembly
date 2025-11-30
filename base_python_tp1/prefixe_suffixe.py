"""
prefixe_suffixe.py — Chevauchement suffixe(X_i) -> préfixe(X_j) en O(mn)

Usage:
  python prefixe_suffixe.py two_reads.fastq
  python prefixe_suffixe.py two_reads.fastq --match 4 --mismatch -4 --gap -8

Lit 2 reads en utilisant utils.read_sequences_from_fastq.
Calcule l’alignement de chevauchement optimal suffixe(X_i) -> préfixe(X_j),
et affiche le score, les alignements et la longueur du chevauchement.
"""

import argparse
from utils import read_sequences_from_fastq


def construire_tables_dp(X_i: str, X_j: str, match: int = 4, mismatch: int = -4, gap: int = -8):
    """
    Une fonction qui prend en argument deux reads (X_i, X_j) et trois scores (match, mismatch, gap).

    Elle construit les tables de programmation dynamique pour le chevauchement
    suffixe(X_i) -> préfixe(X_j) :
      - V : table des scores,
      - BT : table des directions pour la remontée ("diag", "up", "left").
    Initialisation :
      V[i][0] = 0                  (préfixe de X_j vide, suffixe de X_i gratuit)
      V[0][j] = V[0][j-1] + gap    (préfixe de X_j pénalisé)
    Récurrence :
      V[i][j] = max(V[i-1][j-1] + s, V[i-1][j] + gap, V[i][j-1] + gap)
      où s = match si X_i[i-1] == X_j[j-1], sinon mismatch.

    Rend un couple (V, BT) où V est la matrice des scores et BT la matrice de backtracking.
    """
    m, n = len(X_i), len(X_j)
    V = [[0] * (n + 1) for _ in range(m + 1)]
    BT = [[''] * (n + 1) for _ in range(m + 1)]

    for j in range(1, n + 1):                  # Initialisation
        V[0][j] = V[0][j - 1] + gap
        BT[0][j] = "left"
    for i in range(1, m + 1):
        V[i][0] = 0
        BT[i][0] = "up"

    for i in range(1, m + 1):                   # Remplissage
        a = X_i[i - 1]
        for j in range(1, n + 1):
            b = X_j[j - 1]
            s = match if a == b else mismatch
            diag_score = V[i - 1][j - 1] + s
            up_score = V[i - 1][j] + gap
            left_score = V[i][j - 1] + gap

            if diag_score >= up_score and diag_score >= left_score:
                V[i][j], BT[i][j] = diag_score, "diag"
            elif up_score >= left_score:
                V[i][j], BT[i][j] = up_score, "up"
            else:
                V[i][j], BT[i][j] = left_score, "left"

    return V, BT


def meilleur_score_derniere_ligne(V):
    """
    Une fonction qui prend en argument la table des scores V.

    Elle parcourt la dernière ligne (i = m) pour trouver le meilleur score et
    la position j_etoile correspondante.

    Rend un couple (meilleur_score, j_etoile).
    """
    derniere = V[-1]
    meilleur = derniere[0]
    j_etoile = 0
    for j in range(1, len(derniere)):
        if derniere[j] > meilleur:
            meilleur = derniere[j]
            j_etoile = j
    return meilleur, j_etoile


def remontee_chevauchement(X_i: str, X_j: str, BT, j_etoile: int):
    """
    Une fonction qui prend en argument les reads X_i et X_j, la table BT et l’indice j_etoile.

    Elle remonte depuis la case (i = len(X_i), j = j_etoile) jusqu’à atteindre j == 0,
    en reconstruisant l’alignement optimal avec des caractères '-' pour les gaps.

    Rend un couple (alignXi, alignXj) correspondant aux deux lignes alignées.
    """
    i, j = len(X_i), j_etoile
    alignXi, alignXj = [], []

    while True:
        if j == 0:
            break  # début du préfixe de X_j atteint
        if i == 0:
            alignXi.append('-')
            alignXj.append(X_j[j - 1])
            j -= 1
            continue

        direction = BT[i][j]
        if direction == "diag":
            alignXi.append(X_i[i - 1])
            alignXj.append(X_j[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            alignXi.append(X_i[i - 1])
            alignXj.append('-')
            i -= 1
        elif direction == "left":
            alignXi.append('-')
            alignXj.append(X_j[j - 1])
            j -= 1

    alignXi.reverse()
    alignXj.reverse()
    return ''.join(alignXi), ''.join(alignXj)


def longueur_chevauchement_depuis_alignement(alignXi: str, alignXj: str) -> int:
    """
    Une fonction qui prend en argument deux alignements (alignXi, alignXj).

    Elle compte le nombre de colonnes « lettre-lettre » (sans '-') communes
    aux deux alignements, ce qui correspond à la longueur du chevauchement.

    Rend un entier représentant la longueur du chevauchement.
    """
    return sum(1 for a, b in zip(alignXi, alignXj) if a != '-' and b != '-')


def calculer_chevauchement(fastq_file: str, match: int = 4, mismatch: int = -4, gap: int = -8):
    """
    Une fonction qui prend en argument le chemin d’un FASTQ (fastq_file) et trois scores (match, mismatch, gap).

    Elle lit exactement deux reads (paire ordonnée) depuis le FASTQ, construit les tables DP,
    récupère le meilleur score et l’indice j_etoile, remonte l’alignement optimal, puis calcule
    la longueur du chevauchement.

    Rend un dictionnaire :
      - 'score' : score optimal (int),
      - 'alignement' : chaîne à deux lignes (Xi puis Xj),
      - 'longueur' : longueur du chevauchement (int),
      - 'ids' : tuple (id_Xi, id_Xj).
    """
    reads = read_sequences_from_fastq(fastq_file)
    ids = list(reads.keys())
    if len(ids) != 2:
        raise ValueError(f"Le FASTQ doit contenir exactement 2 reads (trouvé {len(ids)}).")

    X_i, X_j = reads[ids[0]], reads[ids[1]]

    V, BT = construire_tables_dp(X_i, X_j, match, mismatch, gap)
    score, j_etoile = meilleur_score_derniere_ligne(V)
    alignXi, alignXj = remontee_chevauchement(X_i, X_j, BT, j_etoile)
    L = longueur_chevauchement_depuis_alignement(alignXi, alignXj)

    return {
        "score": int(score),
        "alignement": alignXi + "\n" + alignXj,
        "longueur": int(L),
        "ids": (ids[0], ids[1]),
    }


def main():
    parser = argparse.ArgumentParser(description="Q5 — chevauchement suffixe(X_i) -> préfixe(X_j) (O(mn))")
    parser.add_argument("fastq", help="FASTQ contenant exactement 2 reads (paire ordonnée)")
    parser.add_argument("--match", type=int, default=4)
    parser.add_argument("--mismatch", type=int, default=-4)
    parser.add_argument("--gap", type=int, default=-8)
    args = parser.parse_args()

    res = calculer_chevauchement(args.fastq, args.match, args.mismatch, args.gap)
    id_i, id_j = res["ids"]
    print(f"IDs: {id_i} (X_i) -> {id_j} (X_j)")
    print(f"Score: {res['score']}")
    print(f"Overlap length: {res['longueur']}")
    print("Alignement X_i:")
    print(res["alignement"].splitlines()[0])
    print("Alignement X_j:")
    print(res["alignement"].splitlines()[1])


if __name__ == "__main__":
    main()