#!/ usr/bin/env python
#coding: utf-8
#By: maxime chaput & jonathan torel
#Baser sur le programme de Arnaud FerwQBG69re
#lecture d'un fichier pdb


import string 
import sys # Pour accéder à exit()
import math
import matplotlib.pyplot as plt

def parsePDBMultiConfig(infile) :
    """
        Cette fonction permet de charger un fichier PDB (Protein Data Bank) au format ATOM.
        Puis de parser son contenu (structure 3D de plusieurs molécule pour le stocker dans une variable Python.
        
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

    config = ""
    dddddd_PDB = {}
    for line in lines :
		if line[0:5] == "MODEL" :
			config= line[5:14].strip()
			dddddd_PDB[config]={}
			dddddd_PDB[config]["domaine"]=[]
			
		if line[0:4] == "ATOM" :
			domain=line[72:75].strip()
			if domain[0]=="A" or domain=="B":
				if config=="":
					if not "ref" in dddddd_PDB.keys(): 
						config="ref"
						dddddd_PDB[config]={}
						dddddd_PDB[config]["domaine"]=[]
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
					dddddd_PDB[config][domain][curres]["contact"]=0
					
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
				d_pdb[confi][domaine][res]["CDM"]={}
				d_pdb[confi][domaine][res]["CDM"]["x"]=mx/mtot #calcul des coordonnees du centre de masse		
				d_pdb[confi][domaine][res]["CDM"]["y"]=my/mtot
				d_pdb[confi][domaine][res]["CDM"]["z"]=mz/mtot
	
#RMSD
#source : https://bioinfo-fr.net/comparaison-de-structures-le-rmsd
   
   
def distc(p1,p2):
	distcarre=0
	for i in p1:
		distcarre+=(p1[i]-p2[i])**2
	return distcarre
    
    
def RMSDconf(d_conf1,d_conf2):
	#fonction permettant de calculer RMSD entre 2 conformations (conf1 étant la reference)
	
	numerateur=0
	n=0
	for domain in d_conf1["domaine"]:
		for resid in d_conf1[domain]["reslist"]:
			numerateur+=distc(d_conf1[domain][resid]["CDM"],d_conf2[domain][resid]["CDM"])
			n+=1
	rmsd=math.sqrt(numerateur/n)
	return rmsd

	
	
def RMSDdom(d_conf1,d_conf2,domaine):
	#fonction permettant de calculer RMSD entre 2 conformations (une étant la reference) pour 1 domaine precis
	#meme fonction que precedement avec juste la condition du domaine
	#source : https://bioinfo-fr.net/comparaison-de-structures-le-rmsd

	numerateur=0
	n=0
	for resid in d_conf1[domaine]["reslist"]:
		numerateur+=distc(d_conf1[domaine][resid]["CDM"],d_conf2[domaine][resid]["CDM"])
		n+=1
	rmsd=math.sqrt(numerateur/n)
	return rmsd
	
    
def preSurfaceContact(d_ref):
	selectionner=[]
	for domaine in d_ref["ref"]["domaine"]:
		if domaine!="B":
			for resA in d_ref["ref"][domaine]["reslist"]:
				for resB in d_ref["ref"]["B"]["reslist"]:
					if distc(d_ref["ref"][domaine][resA]["CDM"],d_ref["ref"]["B"][resB]["CDM"])<196:
						selectionner.append(resA)
						break
	return selectionner
					
	
	
	
	
	
	
	
	
	
	
### Ici, vous pourrez testez vos fonctions :
if __name__ == "__main__":
    
    # Pour afficher une structure de façon un peu plus esthétique :
    #~ import json
    #~ print("Données pur config : \n"+json.dumps(parsePDBMultiConfig("pab21_structure_de_ref.pdb"), indent = 4))
    
    d_ref=parsePDBMultiConfig("pab21_structure_de_ref.pdb")
    d_dico=parsePDBMultiConfig("pab21_prod_solute_500frames.pdb")
    #~ d_dico=parsePDBMultiConfig("pab21_structure_de_ref.pdb")
    calculCDM(d_ref)
    calculCDM(d_dico)
    rmsd_conf=[]
    rmsd_dom={}
    for c in sorted(map(int, d_dico.keys())):
		rmsd_conf.append(RMSDconf(d_ref["ref"],d_dico[str(c)]))
		for d in d_dico[str(c)]["domaine"]:
			if not d in rmsd_dom.keys():
				rmsd_dom[d]=[]
			rmsd_dom[d].append(RMSDdom(d_ref["ref"],d_dico[str(c)],d))
    
    ### affiche les graphique de chaque rmsd. Il faut import matplotlib.pylot as plt
    plt.title("RMSD de chaque configuration de la structure")
    plt.plot(rmsd_conf)
    plt.ylabel("RMSD (en Angtrom)")
    plt.xlabel("Configuration")
    plt.show()
    for dom in rmsd_dom.keys():
		plt.title("RMSD de chaque configuration du domaine : "+dom)
		plt.plot(rmsd_dom[dom])
		plt.ylabel("RMSD (en Angtrom)")
		plt.xlabel("Configuration")
		plt.show()
    selected=preSurfaceContact(d_ref)
    #~ for c in d_dico.keys():
		#~ for d in d_dico[c]["domaine"]:
    
    










