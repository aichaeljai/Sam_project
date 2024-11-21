import csv # comma seperated values pour pouvoir manipuler des fichiers structurés 
import matplotlib.pyplot as plt # Bibliothèque pour les figures
from collections import Counter # Importation de la fonction Counter pour compter dans un dictionnaire



# Définition d'une fonction qui s'appelle sam_reading et qui prend comme argument le chemin du fichier SAM
def sam_reading(sam_file_path):
# Ouvrir le fichier SAM
    with open(sam_file_path, "r") as sam_file:
        # Initialiser le lecteur CSV (camma seperated values) avec le délimiteur de tabulation car le fichier SAM est séparé par des tabulations
        sam_reader = csv.reader(sam_file, delimiter='\t')
        
        # Parcourir chaque ligne du fichier SAM
        flags=[] #list
        quals=[]
        maps_score=[]
        coverage={} #dictionnaire
        # Pour chaque ligne de mon fichier SAM 
        for row in sam_reader:
            # Ignorer les lignes d'en-tête qui commencent par '@' car c'est un header
            if row[0].startswith("@"):
                continue
            
            # Accéder || extraire aux || des informations de chaque colonne
            qname = row[0]    # Nom du read
            flag = int(row[1]) # Flag en entier
            flags.append(flag) # rajouter le flag de chaque ligne dans la liste
            rname = row[2]    # Nom de la séquence de référence
            start_pos = int(row[3])  # Position de début de l'alignement
            mapq = int(row[4]) # Qualité de l'alignement
            maps_score.append(mapq) # rajouter la qualité de chaque ligne dans la liste
            cigar = row[5]    # Chaîne CIGAR (chaîne de caractères qui décrit comment un read s'aligne sur la séquence de référence)
            rnext = row[6]    # Référence du read suivant dans le cas des paires
            pnext = int(row[7]) # Position du read suivant dans le cas des paires
            tlen = int(row[8]) # Longueur du fragment pour les paires
            seq = row[9]      # Séquence de l'ADN
            seq_length = len(seq)
            qual = row[10]    # Qualité de chaque base dans la séquence
            quals.append(qual) # rajouter la qual de chaque ligne dans la liste

            # pour la question 3 localisation des reads sur la séquence 
            # Une boucle qui va parcourir toute la séquence et qui va remplir le dictionnaire avec l'alignement adécquat sur le début de l'alignement
            for pos in range(start_pos, start_pos + seq_length):
                if pos in coverage:
                    coverage[pos] += 1
                else:
                    coverage[pos] = 1
                    
    return flags, quals, coverage, maps_score



# Je définie une fonction nommée flags_to_binary pour convertir le flag en binaire parceque le flag contient les infos en bit 
# def (mot clé pour définir la fonction) définit la fonction flags_to_binary qui prend en paramètre (flag_size et flags)
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
    mapped_reads = []
    flag_size = flag_size - 1 # (0 --> 11)
    nbr = 0
    for i in range(len(binary_flags)):
        # Mettre le flag en binaire de la ligne en question (i) dans la variable flag
        flag=binary_flags[i]
        # le bit 2 code pour l'info "read mappé ou pas", donc on soustrait le chiffre 2 de la taille du flag (12-1=11). Si "0" = mappé 
        if (flag[flag_size-2] == "0"):
        # Rajouter 1 pour compter le nombre de read
            nbr = nbr + 1
            mapped_reads.append(i)
    print ("J'ai", nbr, "read mappés")
    return mapped_reads



#Je prend les flags en binaire et on teste la valeur du bit numéro 4 qui indique si c'est mappé sur le brin complémentaire
#Je définie une fonction nommée mapped_sur_brin_complementaire qui prend en paramètre la taille des flags et les flags en binaire
def mapped_sur_brin_complementaire (flag_size, binary_flags):
    flag_size = flag_size - 1
    nbr = 0
    for i in range(len(binary_flags)):
        flag=binary_flags[i]
        if (flag[flag_size-4] == "1"):
            nbr = nbr + 1
    print ("J'ai", nbr, "read mappés sur le brin complémentaire")



#Je prend les flags en binaire et on teste la valeur du bit numéro 6 qui indique si c'est le 1er de la paire ainsi que le bit numéro 0 qui montre si c'est une paire ou pas 
def premier_de_la_paire (flag_size, binary_flags):
    flag_size = flag_size - 1
    nbr = 0
    for i in range(len(binary_flags)):
        flag=binary_flags[i]
        # Les 2 conditions doivent être VRAI; la 1ère pour dire si le read est le 1er de la paire & la 2ème pour s'assurer que c une paire
        if ((flag[flag_size-6] == "1") and flag[flag_size-0] == "1"):
            nbr = nbr + 1
    print ("J'ai", nbr, "read qui sont les premiers de la paire")



#Je prend les flags en binaire et on teste la valeur du bit numéro 7 qui indique si c'est le 1er de la paire ainsi que le bit numéro 0 qui montre si c'est une paire ou pas 
def second_de_la_paire (flag_size, binary_flags):
    flag_size = flag_size - 1
    nbr = 0
    for i in range(len(binary_flags)):
        flag=binary_flags[i]
        # Les 2 conditions doivent être VRAI; la 1ère pour dire si le read est le 2ème de la paire & la 2ème pour s'assurer que c une paire
        if ((flag[flag_size-7] == "1") and flag[flag_size-0] == "1"):
            nbr = nbr + 1
    print ("J'ai", nbr, "read qui sont les seconds de la paire")






def read_positions (coverage):
    # Ordonner le dictionnaire coverage selon les positions (keys)
    positions = sorted(coverage.keys())
    # Parcourir toute les positions pour avoir le coverage de chaque position pour les afficher dans le plot
    coverages = [coverage[pos] for pos in positions]

    #Afficher un graphique de la couverture
    plt.figure(figsize=(10, 5)) #Initialiser la figure en pouce
    plt.plot(positions, coverages, label="Couverture") #Déssiner un plot
    plt.xlabel("Position sur la séquence de référence") #Label des axes des X
    plt.ylabel("Nombre de reads (couverture)") #Label des axes des Y
    plt.title("Couverture des reads le long de la séquence de référence") #Titre
    plt.legend() #Afficher la légende
    plt.show() #Afficher la figure


def mapping_quality (maps_score):
    # La fonction Counter crée un dictionnaire avec chaque valeur de qualité qui est associé au nombre de reads
    mapq_counts = Counter(maps_score)
    
    print(mapq_counts)

    # défintion du seuil à 30

    seuil_mapq = 30
    read_quality_greater_30 = []

    for i in range(len(maps_score)):
        if (maps_score[i] > 30):
            read_quality_greater_30.append(i)

    return read_quality_greater_30



def filtred_reads (mapped_reads, read_quality_greater_30, sam_file_path):

    filtred_sam_file = "filtred_reads.sam"
    filtered_indices = set(mapped_reads).intersection(read_quality_greater_30)

    with open(sam_file_path, "r") as sam_file, open(filtred_sam_file, "w", newline='') as output_file:
        # Initialiser le lecteur CSV (camma seperated values) avec le délimiteur de tabulation car le fichier SAM est séparé par des tabulations
        sam_reader = csv.reader(sam_file, delimiter='\t')
        output_writer = csv.writer(output_file, delimiter='\t')

        for index, row in enumerate(sam_reader):
        # Conserver les lignes d'en-tête (commencent par '@')
            if row[0].startswith("@"):
                output_writer.writerow(row)
                continue

            # Si l'index du read est dans l'intersection des deux listes, on l'ajoute au fichier de sortie
            if index in filtered_indices:
                output_writer.writerow(row)



sam_file_path = "/Users/aichaeljai/Desktop/projet_mapping/mapping.sam"
flag_size = 12
# J'appelle la fonction sam_reading qui prend en paramètre le chemin et qui me retourne les flags et les quals
flags, quals, coverage, maps_score = sam_reading(sam_file_path)
binary_flags = flags_to_binary(flag_size ,flags)
mapped_reads = number_of_mapped_reads(flag_size, binary_flags)
mapped_sur_brin_complementaire(flag_size, binary_flags)
premier_de_la_paire(flag_size, binary_flags)
second_de_la_paire(flag_size, binary_flags)
read_positions(coverage)
read_quality_greater_30 = mapping_quality(maps_score)
filtred_reads(mapped_reads, read_quality_greater_30, sam_file_path)


# Enregistrer les qualitès dans un fichier texte avec une en-tête et 'w' pour modifier le fichier
with open("quals_output.txt", "w") as output_file:
    # Écrire l'en-tête en sautant une ligne avec \n
    output_file.write("Les qualitès\n")
    
    # Écrire chaque flag dans une nouvelle ligne
    for qual in quals:
        output_file.write(f"{qual}\n")




