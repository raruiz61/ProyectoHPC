#! /usr/bin/env python
 
import sys, time, os

archivo="/home/proyecto/Documentos/resultados/GardnerellaProcesado/Compilado2.txt"

'''
archivoMHCI="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadomhci.txt"
archivoMHCII="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadomhcii.txt"
archivoImmu="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadoimmunogenicity.txt"
archivoCB="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadocbell.txt"

archivoMHCI="/home/proyecto/Documentos/resultados/Helicocteri/resultadomhci.txt"
archivoMHCII="/home/proyecto/Documentos/resultados/Helicocteri/resultadomhcii.txt"
archivoImmu="/home/proyecto/Documentos/resultados/Helicocteri/resultadoimmunogenicity.txt"
archivoCB="/home/proyecto/Documentos/resultados/Helicocteri/resultadocbell.txt"

'''
'''
GardnerellaProcesado
archivoMHCI2="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadomhci2.txt"
archivoMHCII2="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadomhcii2.txt"
archivoImmu2="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadoimmunogenicity2.txt"
archivoCB2="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadocbell2.txt"
'''
'''
mhci = open(archivoMHCI)
mhcii = open(archivoMHCII)
immu = open(archivoImmu)
cellb = open(archivoCB)
'''
f=open(archivo)

#mhci2 = open(archivoMHCI2,"w")
#mhcii2 = open(archivoMHCII2,"w")
#immu2 = open(archivoImmu2,"w")
#cellb2 = open(archivoCB2,"w")

i=0

for line in f:
	i=i+1
	#line=line.replace('"','')
	epitope=line.split(",")[0]
	#print(line.split(",")[0])
	#epitope=epitope.split("+")[0]
	#gi=line.split(",")[8]
	
	'''
	line=line.rstrip()
	line=line.replace(',',' ')
        line=line.replace('[','')
        line=line.replace("'",'')
        line=line.replace("]",'')
        line=line.replace('-r-n','')
        line=line.replace('(','')
        line=line.replace("'",'')
        line=line.replace(")",'')
	'''
	print(">"+str(i))
	print(epitope)
	#epitope=line.split(" ")[0]
	#valor=line.split(" ")[2]

	if i >6787:
		break

#Cierra el archivo
#file1.close()