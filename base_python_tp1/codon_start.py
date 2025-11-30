"""
codon_start.py — Trouve le cadre de lecture du codon start de la protéine (geneX.fasta)
dans la région génomique (sequence.fasta), en considérant que sequence.fasta
est le brin codant.

Usage:
  python codon_start.py --genome sequence.fasta --protein geneX.fasta
"""
import argparse
import sys
from utils import read_single_fasta_sequence


code_genetique = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def traduire_cadre(adn: str, cadre: int) -> str:
    """
    Une fonction qui prend en argument une séquence d'ADN (adn) et un entier (cadre ∈ {0,1,2}).

    Elle traduit la séquence en acides aminés dans le cadre spécifié selon le code génétique standard
    (les codons inconnus deviennent 'X').

    Rend une chaîne correspondant à la traduction en acides aminés du cadre choisi.
    """
    s = adn[cadre:]
    aa = []
    for i in range(0, len(s) - 2, 3):
        codon = s[i:i+3]
        aa.append(code_genetique.get(codon, 'X'))
    return ''.join(aa)

def lire_sequences(chemin_genome: str, chemin_proteine: str) -> tuple[str, str]:
    """
    Une fonction qui prend en argument deux chemins de fichiers FASTA (chemin_genome, chemin_proteine).

    Elle lit les séquences (ADN et protéine) et normalise : U→T, majuscules, retrait du '*' final de la protéine.

    Rend un couple (adn, proteine) sous forme de chaînes.
    """
    adn = read_single_fasta_sequence(chemin_genome).replace('U', 'T').upper()
    prot = read_single_fasta_sequence(chemin_proteine).replace('U', 'T').upper().rstrip('*')
    return adn, prot

def cadre_par_recherche_complete(adn: str, proteine: str) -> list[tuple[int, int]]:
    """
    Une fonction qui prend en argument l'ADN (adn) et la protéine (proteine).

    Elle cherche la protéine comme sous-chaîne stricte et contiguë dans chaque cadre du brin codant.
    Pour chaque hit, elle collecte (cadre, nt0) où nt0 = position 0-based de l'ATG correspondant.

    Rend une liste de couples (cadre, nt0) triée par nt0 croissant.
    """
    candidats: list[tuple[int, int]] = []
    for f in (0, 1, 2):
        aa = traduire_cadre(adn, f)
        j = aa.find(proteine)
        if j != -1:
            nt0 = f + 3 * j
            candidats.append((f, nt0))
    candidats.sort(key=lambda t: t[1])
    return candidats

def cadre_par_prefixe_met(adn: str, proteine: str, min_prefixe: int = 10) -> int | None:
    """
    Une fonction qui prend en argument l'ADN (adn), la protéine (proteine) et une longueur minimale (min_prefixe).

    Elle cherche, pour chaque cadre, le plus long préfixe de 'proteine' (longueur ≥ min_prefixe) qui apparaît
    comme sous-chaîne dans la traduction du cadre.
    En cas d'égalité de longueur, elle choisit le hit le plus en amont en nucléotides.

    Rend l'entier 'cadre' (0,1,2) si trouvé, sinon None.
    """
    meilleurs: list[tuple[int, int, int]] = []  # (longueur_prefixe, nt0, cadre)
    L = len(proteine)
    for f in (0, 1, 2):
        aa = traduire_cadre(adn, f)
        meilleur_len = 0
        meilleur_nt0 = 10**12
        start = 0                               
        while True:
            j = aa.find('M', start)
            if j == -1:
                break
            max_len = 0
            
            maxL = min(L, len(aa) - j)          # Borne supérieure : ne pas dépasser aa ni proteine
            k = 0
            while k < maxL and aa[j + k] == proteine[k]:
                k += 1
            max_len = k
            if max_len >= min_prefixe:
                nt0 = f + 3 * j
                if (max_len > meilleur_len) or (max_len == meilleur_len and nt0 < meilleur_nt0):
                    meilleur_len = max_len
                    meilleur_nt0 = nt0
            start = j + 1
        if meilleur_len > 0:
            meilleurs.append((meilleur_len, meilleur_nt0, f))
    if not meilleurs:
        return None
    meilleurs.sort(key=lambda t: (-t[0], t[1]))  
    return meilleurs[0][2]

def main():
    p = argparse.ArgumentParser(description="Q1(a) — Cadre du codon start sur le brin codant")
    p.add_argument('--genome', required=True, help='sequence.fasta (ADN)')
    p.add_argument('--protein', required=True, help='geneX.fasta (AA)')
    args = p.parse_args()

    adn, proteine = lire_sequences(args.genome, args.protein)
    hits = cadre_par_recherche_complete(adn, proteine)
    if hits:
        cadre = hits[0][0]
        print(f"Cadre {cadre+1}")
        return

    cadre = cadre_par_prefixe_met(adn, proteine, min_prefixe=10)
    if cadre is None:
        print("Aucun cadre compatible trouvé sur le brin codant.")
        sys.exit(1)

    print(f"Cadre {cadre+1}")


if __name__ == '__main__':
    main()