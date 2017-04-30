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

# HEAD
def calculCDM(d_pdb) :
	#auteur : maxime
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
	print d_pdb			
	
#RMSD
#source : https://bioinfo-fr.net/comparaison-de-structures-le-rmsd
   
def RMSDconf(conf1,conf2):
	#fonction permettant de calculer RMSD entre 2 conformations (conf1 étant la reference)
	
	numerateur=float
	n=float
	numerateur=0
	n=0
	
	for resid1 in conf1[""].keys:
		for resid2 in conf2[""].keys:
			if (resid2==resid1):
				for atom in conf2[""][resid2]["atomlist"]: #conf1 valeurs theoriques (ref) conf2 valeurs obs (ref - theo)**2
					numerateur+=((conf1[][resid2][atomtype["x"]-conf2[][resid2][atomtype]["x"])**2)+((conf1[][resid2][atomtype["y"]-conf2[][resid2][atomtype]["y"])**2)+((conf1[][resid2][atomtype["z"]-conf2[][resid2][atomtype]["z"])**2)
					n+=1
	
	rmsd=sqrt(numerateur/n)
	return rmsd
	
	
	
def RMSDdom(conf1,conf2,domaine):
	#fonction permettant de calculer RMSD entre 2 conformations (une étant la reference) pour 1 domaine precis
	#meme fonction que precedement avec juste la condition du domaine
	#source : https://bioinfo-fr.net/comparaison-de-structures-le-rmsd

numerateur=float
	n=float
	numerateur=0
	n=0
	
	for resid1 in conf1[""].keys:
		for resid2 in conf2[""].keys:
			if (resid2==resid1):
				for atom in conf2[""][resid2]["atomlist"]: #conf1 valeurs theoriques (ref) conf2 valeurs obs (ref - theo)**2
					if conf2[][resid2]["domaine"]==domaine:
						numerateur+=((conf1[][resid2][atomtype["x"]-conf2[][resid2][atomtype]["x"])**2)+((conf1[][resid2][atomtype["y"]-conf2[][resid2][atomtype]["y"])**2)+((conf1[][resid2][atomtype["z"]-conf2[][resid2][atomtype]["z"])**2)
						n+=1
	
	rmsd=sqrt(numerateur/n)
	return rmsd
    
	
