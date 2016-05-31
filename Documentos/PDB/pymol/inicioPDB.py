#!/usr/bin/python
import sys, time, os

archivo=sys.argv[1]
f = open(archivo)

#Crea un diccionario con el id de la secuancia como llave y la secuencia como valor 
seqs={}
for line in f:
	line=line.rstrip()
	if  line.startswith('>') :#"line[0]=='>': # or line.startswith('>')			
		words=line.split()
		name=words[0][1:]
		seqs[name]=''
	else : # sequence, not header
		seqs[name] = seqs[name] + line	

#Cierra el archivo
f.close()

i=0
for numId,secuencia in seqs.items():
	i=i+1
	t=open("/home/proyecto/Documentos/resultados/tempSeq"+str(i)+".txt","w")

	t.write(">"+numId+"\n")
	t.write(secuencia)
	t.close()

	comando="/home/proyecto/Documentos/PDB/pymol/crearPDB.py /home/proyecto/Documentos/resultados/tempSeq"+str(i)+".txt "+str(i)
	print(comando)	
	resultado = os.system(comando)
	
	