"""
sequence_frag.py — Reconstruit la séquence finale à partir de graph_2.dot (graphe réduit),
en réutilisant les fonctions de prefixe_suffixe.py pour calculer le chevauchement optimal.

Imprime uniquement :
- la longueur de la séquence,
- la séquence complète.

Usage:
  python sequence_frag.py
  # options :
  #   --reads reads.fq
  #   --dot   graph_2.dot
  #   --match 4 --mismatch -4 --gap -8
"""

import argparse
import re
from utils import read_sequences_from_fastq
from prefixe_suffixe import (
    construire_tables_dp,
    meilleur_score_derniere_ligne,
    remontee_chevauchement,
)


def extraire_chemin_depuis_dot(chemin_dot: str) -> list[str]:
    """
    Une fonction qui prend en argument le chemin d'un fichier DOT (chemin_dot).

    Elle parcourt le DOT réduit pour reconstruire l'ordre linéaire des nœuds (reads)
    en suivant les arêtes u→v. Le départ est choisi comme le nœud de degré entrant nul
    (ou, à défaut, un nœud avec sortie).

    Rend une liste ordonnée des identifiants de reads (ordre de parcours).
    """
    aretes: dict[str, str] = {}
    indegre: dict[str, int] = {}
    noeuds = set()
    motif_arete = re.compile(r'^\s*"([^"]+)"\s*->\s*"([^"]+)"')
    with open(chemin_dot, "r", encoding="utf-8") as f:
        for ligne in f:
            m = motif_arete.search(ligne)
            if m:
                u, v = m.group(1), m.group(2)
                noeuds.add(u)
                noeuds.add(v)
                aretes[u] = v
                indegre[v] = indegre.get(v, 0) + 1
                indegre.setdefault(u, indegre.get(u, 0))
            else:
                # Enregistrer aussi les nœuds isolés s'ils apparaissent comme  "  "u";"
                m2 = re.match(r'^\s*"([^"]+)"\s*;\s*$', ligne)
                if m2:
                    noeuds.add(m2.group(1))

    departs = [u for u in noeuds if indegre.get(u, 0) == 0 and u in aretes]
    if not departs:
        departs = [u for u in aretes.keys()]
        if not departs:
            raise ValueError("Aucun chemin détecté dans le DOT (aucune arête).")

    u = departs[0]
    ordre = [u]
    vus = {u}
    while u in aretes:
        v = aretes[u]
        if v in vus:  
            break
        ordre.append(v)
        vus.add(v)
        u = v
    return ordre


def indice_prefixe_consomme(X: str, Y: str, match: int, mismatch: int, gap: int) -> int:
    """
    Une fonction qui prend en argument deux reads X (gauche) et Y (droite) ainsi que
    trois scores (match, mismatch, gap).

    Elle construit la DP de chevauchement suffixe(X)→préfixe(Y) via construire_tables_dp,
    récupère j* sur la dernière ligne avec meilleur_score_derniere_ligne, puis remonte
    l’alignement optimal avec remontee_chevauchement. L’indice du préfixe de Y consommé
    est le nombre de caractères non ‘-’ dans l’alignement de Y.

    Rend un entier j (0 ≤ j ≤ len(Y)) correspondant à la longueur du préfixe de Y à ignorer
    lors de la concaténation (on ajoute Y[j:]).
    """
    V, BT = construire_tables_dp(X, Y, match=match, mismatch=mismatch, gap=gap)
    _, j_etoile = meilleur_score_derniere_ligne(V)
    _, alignY = remontee_chevauchement(X, Y, BT, j_etoile)
    j = sum(1 for c in alignY if c != '-')
    return j


def assemble(chemin_reads: str, chemin_dot: str, match: int, mismatch: int, gap: int) -> str:
    """
    Une fonction qui prend en argument le chemin d’un FASTQ (chemin_reads), celui d’un DOT (chemin_dot)
    et trois scores (match, mismatch, gap).

    Elle lit les reads (id → séquence), extrait l’ordre linéaire des reads depuis le DOT réduit,
    puis assemble la séquence finale en calculant, pour chaque paire (X, Y), l’indice j de préfixe
    consommé via les fonctions importées de prefixe_suffixe, et concatène Y[j:].

    Rend la séquence nucléotidique assemblée sous forme de chaîne.
    """
    reads = read_sequences_from_fastq(chemin_reads)
    ordre = extraire_chemin_depuis_dot(chemin_dot)
    if not ordre:
        raise ValueError("Ordre vide extrait du DOT.")

    S = reads[ordre[0]]
    for a, b in zip(ordre, ordre[1:]):
        Xa, Yb = reads[a], reads[b]
        j = indice_prefixe_consomme(Xa, Yb, match, mismatch, gap)
        S += Yb[j:]
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", default="reads.fq")
    ap.add_argument("--dot", default="graph_2.dot")
    ap.add_argument("--match", type=int, default=4)
    ap.add_argument("--mismatch", type=int, default=-4)
    ap.add_argument("--gap", type=int, default=-8)
    args = ap.parse_args()

    S = assemble(args.reads, args.dot, args.match, args.mismatch, args.gap)
    print(len(S))
    print(S)


if __name__ == "__main__":
    main()