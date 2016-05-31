#! /usr/bin/env python
 
import sys, time, os

'''
archivoMHCI="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadomhci.txt"
archivoMHCII="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadomhcii.txt"
archivoImmu="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadoimmunogenicity.txt"
archivoCB="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadocbell.txt"
'''
archivoMHCI="/home/proyecto/Documentos/resultados/Helicocteri/resultadomhci.txt"
archivoMHCII="/home/proyecto/Documentos/resultados/Helicocteri/resultadomhcii.txt"
archivoImmu="/home/proyecto/Documentos/resultados/Helicocteri/resultadoimmunogenicity.txt"
archivoCB="/home/proyecto/Documentos/resultados/Helicocteri/resultadocbell.txt"

'''GardnerellaProcesado
archivoMHCI2="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadomhci2.txt"
archivoMHCII2="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadomhcii2.txt"
archivoImmu2="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadoimmunogenicity2.txt"
archivoCB2="/home/proyecto/Documentos/resultados/mycoplasmaProcesado/resultadocbell2.txt"
'''

mhci = open(archivoMHCI)
mhcii = open(archivoMHCII)
immu = open(archivoImmu)
cellb = open(archivoCB)
'''
#mhci2 = open(archivoMHCI2,"w")
mhcii2 = open(archivoMHCII2,"w")
#immu2 = open(archivoImmu2,"w")
#cellb2 = open(archivoCB2,"w")

file1=mhcii
file2=mhcii2
for line in file1:
	line=line.rstrip()
	line=line.replace(',',' ')
        line=line.replace('[','')
        line=line.replace("'",'')
        line=line.replace("]",'')
        line=line.replace('-r-n','')
        line=line.replace('(','')
        line=line.replace("'",'')
        line=line.replace(")",'')


	epitope=line.split(" ")[0]
	valor=line.split(" ")[2]

	lista=[epitope,valor]
	file2.write('\t'.join(lista)+"\n")



#Cierra el archivo
file1.close()

'''
#Crea un diccionario con el id de la secuancia como llave y la secuencia como valor 
seqs={}
resultados=[]
epitopes=[]
valor=0

for line in mhci:
	line=line.rstrip()
	line=line.replace(',','')
        line=line.replace('(','')
        line=line.replace("'",'')
        line=line.replace(")",'')

	epitope=line.split(" ")[0]
	valor=line.split(" ")[1]

	resultados.append([epitope,valor,'0','0','0','0'])
	epitopes.append(epitope)

#Cierra el archivo
mhci.close()

for line in mhcii:
	line=line.rstrip()
	line=line.replace(',','')
        line=line.replace('(','')
        line=line.replace("'",'')
        line=line.replace(")",'')

	epitope=line.split(" ")[0]
	valor=line.split(" ")[1]	
	
	if epitope in  epitopes:
		num=epitopes.index(epitope)
		resultados[num][2]=valor
	else:
		resultados.append([epitope,'0',valor,'0','0','0'])
		epitopes.append(epitope)

for line in immu:
	line=line.rstrip()
	line=line.replace(',','')
        line=line.replace('(','')
        line=line.replace("'",'')
        line=line.replace(")",'')

	epitope=line.split(" ")[0]
	valor=line.split(" ")[1]	
	if epitope in  epitopes:
		num=epitopes.index(epitope)
		resultados[num][3]=valor
	else:
		resultados.append([epitope,'0','0',valor,'0','0'])
		epitopes.append(epitope)

for line in cellb:
	
	line=line.rstrip()
	line=line.replace(',',' ')
        line=line.replace('[','')
        line=line.replace("'",'')
        line=line.replace("]",'')
        line=line.replace('-r-n','')
	#print(line)

	epitope=line.split('  ')[0]
	valor=line.split('  ')[1]

	#print(epitope)
	#print(valor)	
	
	if epitope in  epitopes:
		num=epitopes.index(epitope)
		resultados[num][4]=valor
	else:
		resultados.append([epitope,'0','0','0',valor,'0'])
		epitopes.append(epitope)

titulo=["Epitope",'',"mhci","mhcii","immun","CellB","Total"]

print('\t'.join(titulo))
t1=0.0
t2=0.0

for r in resultados:
	if not float(r[1]) == 0:	
		t1=500/float(r[1])
	else:
		t1=0.0
	
	if not float(r[2]) == 0:	
		t2=500/float(r[2])
	else:
		t2=0.0

	r[5]=t1+t2+2*float(r[3])+5*float(r[4])


resultados.sort(key=lambda t: t[5], reverse=True)

for r in resultados:

	if not float(r[1]) == 0:	
		t1=500/float(r[1])
	else:
		t1=0.0
	
	if not float(r[2]) == 0:	
		t2=500/float(r[2])
	else:
		t2=0.0

	lista=[str(r[0]),str(round(t1,3)),str(round(t2,3)),str(round(float(r[3]),3)),str(float(r[4])),str(round(r[5],3))]
	print('\t'.join(lista))
#'''
