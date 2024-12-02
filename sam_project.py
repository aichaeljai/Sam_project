import csv # comma seperated values pour pouvoir manipuler des fichiers structurés (lecture & écriture)
import matplotlib.pyplot as plt # Bibliothèque pour les figures
from collections import Counter # Importation de la fonction Counter de la bibl collections pour compter dans un dictionnaire
# import pysam
import argparse # Bibliothèque qui crée un help automatiquement et elle permet de gérer les arguments
import os

def script_call():
    # Créer un parser pour les arguments
    parser = argparse.ArgumentParser(
        description="Ce script réaliste des études sur les fichier SAM et renvoie des informations sur les reads"
    )
    
    # Ajouter des arguments 
    # 1er Arg : sert à donner le chemin du fichier SAM à étudier (parser : récupérer)
    parser.add_argument (
        '-i', '--input', 
        type=str,
        required=True, 
        help="Spécifiez le chemin du fichier SAM."
    )
    #2ème Arg : sert à donner le fichier de sortie contenant les reads filtrés
    parser.add_argument(
        '-o', '--output', 
        type=str, 
        required=True,
        help="Spécifiez le fichier SAM de sortie."
    )
    #3ème Arg : sert à donner le flag size
    parser.add_argument(
        '-fs', '--flag_size', 
        type=int, 
        required=True,
        help="Spécifiez la taille des FLAGS en bits."
    )
    
    # Parser tout les arguments
    args = parser.parse_args()
    
    sam_file_path = args.input
    if os.path.exists (sam_file_path):
        print("fichier SAM trouvé")
    else:
        print("Le chemin est incorrect: Fichier non trouvé")
        print("\n")
    
    output_sam_file = args.output
    flag_size = args.flag_size

    return sam_file_path, output_sam_file, flag_size



# Lecture du fichier SAM
# Définition d'une fonction qui s'appelle sam_reading et qui prend comme argument/param le chemin du fichier SAM
def sam_reading(sam_file_path):
# Ouvrir le fichier SAM en mode read "r"
    with open(sam_file_path, "r") as sam_file:
        # Initialiser le lecteur CSV (camma seperated values) avec le délimiteur de tabulation car le fichier SAM est séparé par des tabulations
        # sam_reader est un tableau (matrice)
        sam_reader = csv.reader(sam_file, delimiter='\t')
        #remplacer avec r split 
        # Parcourir chaque ligne du fichier SAM à partir de la boucle
        # création (déclaration) de listes et dictionnaires vides pour les remplir au fur et à mesure des lignes
        flags=[] #list
        quals=[]
        maps_score=[]
        coverage={} #dictionnaire
        cigars=[]
        references=[]
        start_positions=[]
        seq_lengths=[]
        # Pour chaque ligne de mon fichier SAM je fais une boucle 
        for row in sam_reader:
            # Ignorer les lignes d'en-tête qui commencent par '@' car c'est un header
            if row[0].startswith("@"):
                continue
            
            # Accéder || extraire aux || des informations de chaque colonne
            qname = row[0]    # Nom du read
            flag = int(row[1]) # Flag en entier
            flags.append(flag) # rajouter le flag de chaque ligne dans la liste
            rname = row[2]    # Nom de la séquence de référence
            references.append(rname)
            start_pos = int(row[3])  # Position de début de l'alignement
            start_positions.append(start_pos)
            mapq = int(row[4]) # Qualité de l'alignement
            maps_score.append(mapq) # rajouter la qualité de chaque ligne dans la liste 
            cigar = row[5]    # Chaîne CIGAR (chaîne de caractères qui décrit comment un read s'aligne sur la séquence de référence)
            cigars.append(cigar)
            rnext = row[6]    # Référence du read suivant dans le cas des paires
            pnext = int(row[7]) # Position du read suivant dans le cas des paires
            tlen = int(row[8]) # Longueur du fragment pour les paires
            seq = row[9]      # Séquence de l'ADN
            seq_length = len(seq)
            seq_lengths.append(seq_length)
            qual = row[10]    # Qualité de chaque base dans la séquence
            quals.append(qual) # rajouter la qual de chaque ligne dans la liste

    return flags, quals, maps_score,cigars, references, start_positions, seq_lengths



# Question 1: Combien de reads sont mappés ? # compter le nombre de reads en fonction du flag (colonne #2)
# Je définie une fonction nommée flags_to_binary pour convertir le flag en binaire parceque le flag contient les infos en bit 
# def définit la fonction flags_to_binary qui prend en paramètre (flag_size et flags)
def flags_to_binary(flag_size, flags):
    # Boucle qui va parcourir de 0 à la taille des flags-1
    for i in range(len(flags)):
       # Convertir chaque flag en binaire
        flags[i]=bin(flags[i])
        # Éliminer les 2 premiers caractères du flag (0 et b) & les remplir de 0 sur la gauche -> atteindre la taille souhaitée qui est 12
        flags[i]=flags[i][2:].zfill(flag_size)
        # Faire retourner une liste de flag transformée en binaire
    return flags
        


# Je prend les flags en binaire et on teste la valeur du bit numéro 2 qui indique si le read est mappé ou pas
# Je définie une fonction nommée number_of_mapped_reads qui prend en paramètre la taille des flags et les flags en binaire
def number_of_mapped_reads (flag_size, binary_flags):
    # déclaration d'une liste vide afin d'y mettre les reads mappés
    mapped_reads = []
    flag_size = flag_size - 1 # (0 --> 11)
    nbr = 0
    for i in range(len(binary_flags)):
        # Mettre le flag en binaire de la ligne en question (i) dans la variable flag
        flag=binary_flags[i]
        # le bit 2 code pour l'info "read mappé ou pas", donc on soustrait le chiffre 2 de la taille du flag (12-1=11). Si "0" = mappé 
        if (flag[flag_size-2] == "0"):
        # flag[-2]
        # Rajouter 1 pour compter le nombre de read
            nbr = nbr + 1
            mapped_reads.append(i)
    # Calculer le % des reads mappés par rapport au nbre total de reads
    nbr_prcnt = (nbr / len(binary_flags))*100
    print ("première question:")
    print(f"j'ai {len(binary_flags)} de reads")
    print (f"J'ai {nbr} de read mappés qui représente {nbr_prcnt} % des reads")
    print (f"J'ai {len(binary_flags)-nbr} de read non mappés qui représente {100-nbr_prcnt} % des reads\n")
    return mapped_reads,nbr



# Question 2: Comment les reads et paires de read sont-ils mappés ? Compter le nombre de reads pour chaque flag
# pour les reads mappés, combien sont partiellement/entièrement/incorrectement mappés 

def fully_or_partially_mapped_reads(nbr, cigars,mapped_reads):
    mapped_reads = set(mapped_reads)
    fully_mapped=0
    partially_mapped=0
    incorrectly_mapped=0

    for i in range(len(cigars)):
        if i in mapped_reads:
            if "S" in cigars[i] or "H" in cigars[i] :
                partially_mapped = partially_mapped+1
            elif  "I" in cigars[i] or "D" in cigars[i] or "N" in cigars[i] or "X" in cigars[i] or "P" in cigars[i]:
                incorrectly_mapped = incorrectly_mapped + 1

            else:
                fully_mapped = fully_mapped +1
    fully_mapped_prcnt = (fully_mapped / nbr)*100
    partially_mapped_prcnt = (partially_mapped / nbr)*100
    incorrectly_mapped_prcnt = (incorrectly_mapped / nbr)*100


    print("duxième question:")
    print(f" j'ai {fully_mapped} reads qui sont entièrement mappés ce qui représente {fully_mapped_prcnt}% des reads mappés")
    print(f" j'ai {partially_mapped} reads qui sont partiellement mappés ce qui représente {partially_mapped_prcnt}% des reads mappés")
    print(f" j'ai {incorrectly_mapped} reads qui ne sont pas correctement mappés ce qui représente {incorrectly_mapped_prcnt}% des reads mappés")
    print("\n")


# calculer le nombre de paires de départ et combien de read mappés ; rajouter fonction de la somme dans l'appel de fonction
# Je prend les flags en binaire et on teste la valeur du bit numéro 6 qui indique si c'est le 1er de la paire ainsi que le bit numéro 0 qui montre si c'est une paire ou pas 
def paire_de_reads (flag_size, binary_flags):
    flag_size = flag_size - 1
    nbr_1 = 0 #premier de la paire
    nbr_2= 0 #second de la paire
    for i in range(len(binary_flags)):
        flag=binary_flags[i]
        # Les 2 conditions doivent être VRAI; la 1ère pour dire si le read est le 1er de la paire & la 2ème pour s'assurer que c une paire
        if ((flag[flag_size-6] == "1") and flag[flag_size-0] == "1"):
            nbr_1 = nbr_1 + 1
        elif ((flag[flag_size-7] == "1") and flag[flag_size-0] == "1"):
            nbr_2 = nbr_2 + 1
    
    nbr_de_paires = nbr_1 + nbr_2
    nbr_1_prcnt = (nbr_1 / nbr_de_paires)*100
    nbr_2_prcnt = (nbr_2 / nbr_de_paires)*100        

    print(f"j'ai {nbr_de_paires} paire de reads")
    print (f"J'ai {nbr_1} de reads qui sont les premiers de la paire qui représentent {nbr_1_prcnt}%")
    print (f"J'ai {nbr_2} de reads qui sont les seconds de la paire qui représentent {nbr_2_prcnt}%")
    print("\n")
 

#Question 3: Où les reads sont-ils mappés ? L'alignement est-il homogène le long de la séquence de référence ? compter le nombre de reads par chromosome 
# Faire une boucle qui parcourt les réf et qui compte le nmbre de reads pour chaque réf trouvée 
def repartition_par_chromosome (references):
    chrom_counts = {}
    for i in range(len(references)):
        chrom = references[i]
        if chrom in chrom_counts:
            chrom_counts[chrom] += 1
        else:
            chrom_counts[chrom]=1
    print("3ème question:")
    print("répartition par chromosome:")
    # faire une boucle qui parcourt les réf qui sont les clés de mon dictionnaire et me donner le nbre de reads dans cette ref 
    for chrom in chrom_counts.keys():
        print(f"j'ai {chrom_counts[chrom]} reads sur {chrom}")
    print("\n")

def read_positions (start_positions, seq_lengths, references):
    # Où les reads sont mappés :
    couverture_par_ref={}
    for i in range(len(references)):
        chrom = references[i]
        # si la réf n'existe pas dans le dictionnaire
        if chrom not in couverture_par_ref:
            # accéder à la clé chrom du dict et j'initialise sa valeur avec un 2ème dict vide (j'ai 3 niveaux = réf, pos, nbr de reads)
            couverture_par_ref[chrom]={}
        # parcourir les positions du read en Q     
        for pos in range (start_positions[i], start_positions[i] + seq_lengths[i]):
            # si la position existe dans le 2ème dict imbriqué, j'incrémente la valeur de 1 sinon j'initialise à 1
            if pos in couverture_par_ref[chrom]:
                couverture_par_ref[chrom][pos] +=1
            else:
                couverture_par_ref[chrom][pos] =1 

    # Faire une boucle avec 2 var : ref et couverture (comment mes reads couvrent la réf)
    for ref, coverage in couverture_par_ref.items():
    # Ordonner le dictionnaire coverage selon les positions (keys)
        positions = sorted(coverage.keys())
        # Parcourir toute les positions pour avoir le coverage de chaque position pour les afficher dans le plot
        coverages = list(coverage.values())

        #Afficher un graphique de la couverture
        plt.figure(figsize=(10, 5)) #Initialiser la figure en pouce
        plt.plot(positions, coverages, label=f"couverture sur {ref}") #Déssiner un plot
        plt.xlabel("Position sur la séquence de référence") #Label des axes des X
        plt.ylabel("Nombre de reads (couverture)") #Label des axes des Y
        plt.title("Couverture des reads le long de la séquence de référence") #Titre
        plt.legend() #Afficher la légende
        plt.show() #Afficher la figure

#question 4: Avec quelle qualité les reads sont-ils mappés ? # compter le nombre de reads pour chaque valeur de qualité ou par tranche de valeurs (score de mapping)
# J'avais une liste qui contient les val de qualité et avec counter on a crée un dict qui a la valeur de qualité comme clé et le nbre de read comme valeurs.
def mapping_quality (maps_score): 
    # la liste (maps_score) contient la qualité de mapping pour chaque read. I used Counter pour compter les occurences de chaque score de qualité
    # La fonction  Counter prend en paramètre la liste maps_score et me retourne le dict mapq_counts qui contient les qualités comme clef et les nbre de reads 
    # pour chaque qualité comme val
    mapq_counts = Counter(maps_score)
    print ("quatrième question:")
    for i in mapq_counts.keys():
        print(f"j'ai {mapq_counts[i]} reads avec la qualité {i}")
    

    # défintion du seuil à 30: identification des reads ayant une qualité inf au seuil 30 et les retourne

    seuil_mapq = 30
    read_quality_lower_30 = []
    # ne prendre que les reads avec une qualité inférieure à 30
    for i in range(len(maps_score)):
        if (maps_score[i] <seuil_mapq):
            read_quality_lower_30.append(i)

    return read_quality_lower_30



def filtred_reads (mapped_reads, read_quality_lower_30, sam_file_path, filtred_sam_file): 
    # lignes contenant les reads filtrés <30
    filtered_indices = set(mapped_reads).union(read_quality_lower_30)

    with open(sam_file_path, "r") as sam_file, open(filtred_sam_file, "w", newline='') as output_file:
        # Initialiser le lecteur CSV (camma seperated values) avec le délimiteur de tabulation car le fichier SAM est séparé par des tabulations
        # La fonction csv.reader prend en charge le fichier sam_file et le stocke dans la tableau sam_reader
        sam_reader = csv.reader(sam_file, delimiter='\t')
        output_writer = csv.writer(output_file, delimiter='\t')

        for index, row in enumerate(sam_reader):
        # Conserver les lignes d'en-tête (commencant par '@')
            if row[0].startswith("@"):
                output_writer.writerow(row)
                continue

            # Si l'index du read est dans l'intersection des deux listes (filtré + entièrement mappé), on l'ajoute au fichier de sortie
            if index in filtered_indices:
                output_writer.writerow(row)

# Appel aux fonctions
# gérer la portabilité et généraliser si l'utilisateur décide de réutiliser le flag size ou autre (faire une condition)
sam_file_path,sam_output_file, flag_size= script_call()
# J'appelle la fonction sam_reading qui prend en paramètre le chemin et qui me retourne les flags et les quals ... 
flags, quals, maps_score, cigars, references, start_positions, seq_lengths = sam_reading(sam_file_path)
binary_flags = flags_to_binary(flag_size ,flags)
mapped_reads,nbr = number_of_mapped_reads(flag_size, binary_flags)
fully_or_partially_mapped_reads(nbr, cigars, mapped_reads)
paire_de_reads(flag_size, binary_flags)
repartition_par_chromosome(references)
read_positions(start_positions, seq_lengths, references)
read_quality_greater_30 = mapping_quality(maps_score)
filtred_reads(mapped_reads, read_quality_greater_30, sam_file_path, sam_output_file)

