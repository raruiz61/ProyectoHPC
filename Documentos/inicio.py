#!/usr/bin/env python

import sys 
import os


comando=""
#archivo="/home/proyecto/Documentos/mhc_i/examples/input_sequence.fasta"
archivo=sys.argv[1]
resultado="/home/proyecto/Documentos/resultados/resultadomhci.txt"

print(archivo)
#resultado = os.system("python /home/proyecto/Documentos/pri/principal2.py")

comando="/home/proyecto/Documentos/mhc_i/src/predict_binding.py netmhcpan "+"HLA-E*01:01"+" 9 "+archivo+" > "+resultado

#/home/proyecto/Documentos/mhc_i/src/predict_binding.py netmhcpan HLA-E*01:01 9  > /home/proyecto/Documentos/resultados/resultadomhci.txt

resultado = os.system(comando)
