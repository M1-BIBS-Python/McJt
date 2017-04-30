#!/ usr/bin/env python
#coding: utf-8
#By: maxime chaput & jonathan torel
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
			
			#verification que les coordonnees sont bien des float:
			dddddd_PDB[config][domain][curres][atomtype]["x"] = float(line[30:38])
			dddddd_PDB[config][domain][curres][atomtype]["y"] = float(line[38:46])
			dddddd_PDB[config][domain][curres][atomtype]["z"] = float(line[46:54])
			dddddd_PDB[config][domain][curres][atomtype]["id"] = line[6:11].strip()
			
    #~ print dddddd_PDB[config]["domaine"]
    return dddddd_PDB


###
# Ici, vous écrirez vos prochaines fonctions
###


### Ici, vous pourrez testez vos fonctions :
if __name__ == "__main__":
    
    # Pour afficher une structure de façon un peu plus esthétique :
    #~ import json
    #~ print("Données pur config : \n"+json.dumps(parsePDBMultiConfig("pab21_structure_de_ref.pdb"), indent = 4))
    d_dico=parsePDBMultiConfig("pab21_prod_solute_500frames.pdb")
    #~ print d_dico
    #~ print("\nIdentifiants des chaînes de la protéine 1EJH :")
    #~ print(parsePDBMultiChains("1EJH.pdb")["chains"])
    #~ print("\nIdentifiants des résidus de la chaîne A de la protéine 1EJH :")
    #~ print(parsePDBMultiChains("1EJH.pdb")["A"]["reslist"])
    
	
