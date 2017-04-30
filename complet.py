#!/ usr/bin/env python
#coding: utf-8
#By: maxime chaput droit:  cc by
#Baser sur le programme de Arnaud FerwQBG69re
#lecture d'un fichier pdb


import string 
import sys # Pour accéder à exit()


def parsePDBMultiConfig(infile) :
    """
        Cette fonction permet de charger un fichier PDB (Protein Data Bank) au format ATOM.
        Puis de parser son contenu (structure 3D d'une molécule) pour le stocker dans une variable Python.
        
        Paramètre(s) :
            - infile : emplacement du fichier à charger et parser
        
        Valeur renvoyée :
            - dddd_PDB : dictionnaire contenant 
        
        Pour plus d'informations sur les fichiers PDB et en particulier le format ATOM, veuillez vous reporter à :
        http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM 
    """

    ### On vérifie que l'ouverture du fichier se passe correctement :
    try:
		f = open(infile, "r")
		lines = f.readlines()
		f.close()
    except:
		print("Le fichier n'a pu être chargé correctement. Vérifiez que le fichier existe bien et relancez votre programme.")
		sys.exit(0) ### Stoppe simplement l'exécution du programme.


    dddddd_PDB = {}
    for line in lines :
		if line[0:5] == "MODEL" :
			config= line[5:14].strip()
			dddddd_PDB[config]={}
			dddddd_PDB[config]["domaine"]=[]
		if line[0:4] == "ATOM" :
			domain=line[72:74].strip()
			if not domain in dddddd_PDB[config]["domaine"]:
				dddddd_PDB[config]["domaine"].append(domain)
				dddddd_PDB[config][domain]={}
				dddddd_PDB[config][domain]["reslist"] = []
			curres = line[22:26].strip()
			if not curres in dddddd_PDB[config][domain]["reslist"] :
				dddddd_PDB[config][domain]["reslist"].append(curres)
				dddddd_PDB[config][domain][curres] = {}
				dddddd_PDB[config][domain][curres]["resname"] = line[17:20].strip()
				dddddd_PDB[config][domain][curres]["atomlist"] = []
			atomtype = line[12:16].strip()
			dddddd_PDB[config][domain][curres]["atomlist"].append(atomtype)
			dddddd_PDB[config][domain][curres][atomtype] = {}
			dddddd_PDB[config][domain][curres][atomtype]["x"] = float(line[30:38])
			dddddd_PDB[config][domain][curres][atomtype]["y"] = float(line[38:46])
			dddddd_PDB[config][domain][curres][atomtype]["z"] = float(line[46:54])
			dddddd_PDB[config][domain][curres][atomtype]["id"] = line[6:11].strip()
<<<<<<< HEAD
    #~ dddddd_PDB[config][domain][curres]["atomlist"]
=======
    #~ print dddddd_PDB[config]["domaine"]
>>>>>>> f616193dd5cfe73ebf3bb0c854e254662fc387a1
    return dddddd_PDB


###
# Ici, vous écrirez vos prochaines fonctions
###

<<<<<<< HEAD
def calculCDM(d_pdb) :
	"""
        Cette fonction permet de calculer le centre de masse de chaque residus de chaque configuration
        
        Paramètre(s) :
            - infile : dictionnaire resultat du parser
        
        Valeur renvoyée :
            - d_pdb : dictionnaire modife
    """
	for confi in d_pdb.keys() :
		for domaine in d_pdb[confi]["domaine"]:
			for res in d_pdb[confi][domaine]["reslist"]:
				mx=0
				my=0
				mz=0
				mtot=0
				for atome in d_pdb[confi][domaine][res]["atomlist"]:
					#donne la ponderation de chaque atome
					if atome[0] == "H" :                                   
						masse=1
					elif atome[0] == "C" :
						masse=12
					elif atome[0] =="O" :
						masse=16
					elif atome[0] == "N" :
						masse=14
					elif atome[0] == "S" :
						masse=32
					elif atome[0] == "P" :
						masse=31
					elif atome[0] == "Z" :
						masse=65
					#faire la somme pondere pour chaque coordonee
					mx+=masse*d_pdb[confi][domaine][res][atome]["x"]     
					my+=masse*d_pdb[confi][domaine][res][atome]["y"]
					mz+=masse*d_pdb[confi][domaine][res][atome]["z"]
					mtot+=masse
				d_pdb[confi][domaine][res]["CDM"]=[mx/mtot,my/mtot,mz/mtot] #calcul des coordonnees du centre de masse
   















=======
>>>>>>> f616193dd5cfe73ebf3bb0c854e254662fc387a1

### Ici, vous pourrez testez vos fonctions :
if __name__ == "__main__":
    
    # Pour afficher une structure de façon un peu plus esthétique :
    #~ import json
    #~ print("Données pur config : \n"+json.dumps(parsePDBMultiConfig("pab21_structure_de_ref.pdb"), indent = 4))
    d_dico=parsePDBMultiConfig("pab21_prod_solute_500frames.pdb")
<<<<<<< HEAD
    calculCDM(d_dico)
=======
>>>>>>> f616193dd5cfe73ebf3bb0c854e254662fc387a1
    #~ print d_dico
    #~ print("\nIdentifiants des chaînes de la protéine 1EJH :")
    #~ print(parsePDBMultiChains("1EJH.pdb")["chains"])
    #~ print("\nIdentifiants des résidus de la chaîne A de la protéine 1EJH :")
    #~ print(parsePDBMultiChains("1EJH.pdb")["A"]["reslist"])
    
	
