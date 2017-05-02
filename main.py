#!/ usr/bin/env python
#coding: utf-8
#By: maxime chaput & jonathan torel
#Baser sur le programme de Arnaud FerwQBG69re
#main

#import a definir une fois fonctions finies
import string 
import sys # Pour accéder à exit()
import final_fonction as ps
import math
import matplotlib.pyplot as plt


#Tout d'abord parsons les données
#le modèle de reference:



d_ref=ps.parsePDBMultiConfig("pab21_structure_de_ref.pdb")

#les 500 conformations:

d_conf=ps.parsePDBMultiConfig("pab21_prod_solute_500frames.pdb")



#calcul CDM : a voir si on calcule RMSD depuis le CDM apres

ps.calculCDM(d_ref)
ps.calculCDM(d_conf)
#calcul RMSD pour chaque conf 

rmsd_conf=[]
rmsd_dom={}
for c in sorted(map(int, d_dico.keys())):
	rmsd_conf.append(RMSDconf(d_ref["ref"],d_dico[str(c)]))
#calcul RMSD pour chaque domaine de chaque conf	
	for d in d_conf[str(c)]["domaine"]:
		if not d in rmsd_dom.keys():
			rmsd_dom[d]=[]
		rmsd_dom[d].append(RMSDdom(d_ref["ref"],d_conf[str(c)],d))
	
#graphique resumant les rmsd :
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
	
	
#determine la zone d interface:
##avec une pre-selection avec d_ref pour gagner du temp
selected=ps.preSurfaceContact(d_ref)
ps.surfaceContact(d_dico,selected)

#Modification du dico d_ref pour donner l'occurence de contacte du residus sur l arn
#ainsi que la creation d exemples (1 par domaine)    
d_pair={}
for d in d_ref["ref"]["domaine"]:
	if d!="B":
		d_pair[d]={}
		d_pair[d]["vide"]=1
		for res in d_ref["ref"][d]["reslist"]:
			for c in d_conf.keys():
				if res in d_conf[c]["interface"]:
					d_ref["ref"][d][res]["contact"]+=1
					if d_pair[d]["vide"]==1:
						d_pair[d]["res"]=[d_conf[c][d][res]["resname"],res]
						d_pair[d]["resB"]=d_conf[c][d][res]["avec"]
						d_pair[d]["vide"]=0
