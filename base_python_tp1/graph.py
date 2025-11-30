"""
graph.py — Construit les graphes de chevauchement depuis une matrice de scores 20x20.

Le script lit un CSV (matrice_20x20.csv) contenant les scores suffixe(X_i) -> préfixe(X_j),
puis produit :
  - graph_1.dot : graphe filtré au seuil (avant réduction),
  - graph_2.dot : graphe réduit (après suppression des 2-cycles et réduction transitive).

Usage :
  python graph.py matrice_20x20.csv --threshold 80
"""

import argparse
import csv
from collections import defaultdict, deque


def read_scores_matrix(csv_path: str) -> tuple[list[str], dict[str, dict[str, int]]]:
    """
    Une fonction qui prend en argument le chemin d’un fichier CSV (csv_path) représentant une matrice 20x20.

    Elle lit le CSV, vérifie la présence d’un en-tête avec la première colonne 'id', puis
    construit et retourne la liste des identifiants (ids) et un dictionnaire M tel que M[u][v] = score_ij.

    Rend un couple (ids, M) où ids est une liste d’IDs de reads et M une table de scores (dict de dict).
    """
    with open(csv_path, "r", encoding="utf-8") as f:
        rows = [row for row in csv.reader(f) if row]  # ignore lignes vides

    if not rows or len(rows) < 2:
        raise ValueError("CSV invalide: pas assez de lignes.")
    header = rows[0]
    if header[0] != "id" or len(header) < 2:
        raise ValueError("Entête invalide: la première cellule doit être 'id'.")

    ids = header[1:]
    M: dict[str, dict[str, int]] = {u: {} for u in ids}
    for row in rows[1:]:
        rid = row[0]
        vals = row[1:]
        if len(vals) != len(ids):
            raise ValueError(f"Ligne incohérente pour {rid}: {len(vals)} valeurs, attendu {len(ids)}.")
        for j, v in enumerate(ids):
            try:
                M[rid][v] = int(float(vals[j].strip()))
            except ValueError:
                raise ValueError(f"Score non entier en ({rid},{v}): {vals[j]!r}")
    return ids, M


def build_graph_from_matrix(ids: list[str], M: dict[str, dict[str, int]], threshold: int) -> dict[str, dict[str, int]]:
    """
    Une fonction qui prend en argument la liste des IDs (ids), la table des scores M et un seuil (threshold).

    Elle construit un graphe orienté G où l’on garde l’arête u->v si (u != v) et M[u][v] >= threshold.
    En cas de doublon, elle conserve le score maximal pour l’arête.

    Rend un dictionnaire d’adjacence G tel que G[u][v] = score.
    """
    G: dict[str, dict[str, int]] = defaultdict(dict)
    for u in ids:
        for v in ids:
            if u == v:
                continue
            s = int(M[u][v])
            if s >= threshold:
                if v not in G[u] or s > G[u][v]:    # garde le meilleur si doublon
                    G[u][v] = s
    return G


def break_two_cycles(G: dict[str, dict[str, int]]) -> None:
    """
    Une fonction qui prend en argument un graphe orienté G (dictionnaire d’adjacence).

    Elle supprime les 2-cycles : si u->v et v->u existent, elle ne garde que l’arête au score le plus fort
    (en cas d’égalité, elle garde v->u et supprime u->v).

    Rend None (effet de bord : modifie G en place).
    """
    to_del: list[tuple[str, str]] = []
    for u in list(G.keys()):
        for v in list(G[u].keys()):
            if u != v and u in G.get(v, {}):
                su, sv = G[u][v], G[v][u]
                if su >= sv:
                    to_del.append((v, u))  # supprime v->u
                else:
                    to_del.append((u, v))  # supprime u->v
    for a, b in to_del:
        if a in G and b in G[a]:
            del G[a][b]
        if a in G and not G[a]:
            del G[a]


def reachable_excluding_edge(G: dict[str, dict[str, int]], u: str, v: str) -> bool:
    """
    Une fonction qui prend en argument un graphe G et deux nœuds (u, v).

    Elle effectue une BFS depuis u vers v en **ignorant** l’arête directe (u->v),
    afin de savoir si v reste atteignable par un autre chemin.

    Rend True si v est atteignable sans utiliser l’arête (u->v), sinon False.
    """
    visited = {u}
    q = deque([u])
    while q:
        x = q.popleft()
        for y in G.get(x, {}):
            if x == u and y == v:
                continue            # ignorer l'arête directe
            if y not in visited:
                visited.add(y)
                q.append(y)
    return v in visited


def transitive_reduction(G: dict[str, dict[str, int]]) -> None:
    """
    Une fonction qui prend en argument un graphe orienté G.

    Elle applique la réduction transitive : pour toute arête u->v, si v est atteignable
    depuis u par un autre chemin ne passant pas par (u->v), alors l’arête (u->v) est supprimée.

    Rend None (effet de bord : modifie G en place).
    """
    to_del: list[tuple[str, str]] = []
    for u in list(G.keys()):
        for v in list(G[u].keys()):
            if reachable_excluding_edge(G, u, v):
                to_del.append((u, v))
    for u, v in to_del:
        if u in G and v in G[u]:
            del G[u][v]
        if u in G and not G[u]:
            del G[u]


def write_dot(path: str, ids: list[str], G: dict[str, dict[str, int]], title: str) -> None:
    """
    Une fonction qui prend en argument un chemin de sortie (path), la liste d’IDs (ids),
    un graphe G et un titre (title).

    Elle écrit un fichier au format Graphviz DOT avec les nœuds listés dans ids et les arêtes G[u][v]
    annotées par leurs scores (label), en configurant un rendu gauche-droite.

    Rend None (effet de bord : crée/écrit le fichier DOT).
    """
    safe_title = title.replace('"', "'")
    with open(path, "w", encoding="utf-8") as f:
        f.write('digraph "' + safe_title + '" {\n')
        f.write("  rankdir=LR;\n  node [shape=box, fontsize=10];\n  overlap=false;\n  splines=true;\n")
        for u in ids:
            f.write(f'  "{u}";\n')
        for u in ids:
            for v, s in G.get(u, {}).items():
                f.write(f'  "{u}" -> "{v}" [label="{s}"];\n')
        f.write("}\n")


def main():
    ap = argparse.ArgumentParser(description="Génère deux graphes DOT (filtré & réduit) depuis la matrice 20x20")
    ap.add_argument("scores_csv", help="matrice_20x20.csv")
    ap.add_argument("--threshold", type=int, default=80, help="seuil minimum pour garder une arête")
    args = ap.parse_args()

    ids, M = read_scores_matrix(args.scores_csv)

    G_thresh = build_graph_from_matrix(ids, M, threshold=args.threshold)
    write_dot("graph_1.dot", ids, G_thresh, f"Overlap >= {args.threshold}")
    print("Écrit : graph_1.dot")

    G_red = {u: dict(G_thresh.get(u, {})) for u in ids if u in G_thresh}
    break_two_cycles(G_red)
    transitive_reduction(G_red)
    write_dot("graph_2.dot", ids, G_red, "Reduced overlap graph")
    print("Écrit : graph_2.dot")


if __name__ == "__main__":
    main()