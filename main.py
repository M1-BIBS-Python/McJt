#main

#import a definir une fois fonctions finies

#Tout d'abord parsons les données
#le modèle de reference:



ref = parseur("pab21_structure_de_ref.pdb") #nom de fonction a mettre à jour

#les 500 conformations:

conf= ParsePDBMultiConfig("pab21_prod_solute_500frames.pdb")



#calcul CDM : a voir si on calcule RMSD depuis le CDM apres

CDMref=calculCDM("pab21_structure_de_ref.pdb")


CDMconf=("pab21_prod_solute_500frames.pdb")

#calcul RMSD pour chaque conf 

for config in conf.keys():
	RMSD[config]={} 
	RMSD[config]["prot"]=RMSDconf(ref,conf[configu])

#calcul RMSD pour chaque domaine de chaque conf	

for config in conf.keys():	
	
	RMSD[config]["A1"]=RMSDdom(ref,conf[config],"A1")
	RMSD[config]["A2"]=RMSDdom(ref,conf[config],"A2")
	RMSD[config]["A3"]=RMSDdom(ref,conf[config],"A3")
	RMSD[config]["A4"]=RMSDdom(ref,conf[config],"A4")
	RMSD[config]["B"]=RMSDdom(ref,conf[config],"B")
	
#fichier de sortie...
