#!/usr/bin/env python

'''
Created on 03.18.2014
@author: Dorjee Tamang
'''

import os
from pylab import * #@UnusedWildImport
from optparse import OptionParser
import re

from src.BCell import BCell

class Prediction():
            
    def predict(self, options, args):
        '''Returns the prediction result.'''

        value = options.method

        import pickle
        
        file_path = os.path.join(os.path.dirname(__file__),'bcell_scales.pickle')
        pickle_file = open(file_path, 'rb')	
        scale_dict = pickle.load(pickle_file)
	#print(scale_dict)
        #print("entrooo")
        # set values for a particular method
        bcell.scale_dict = scale_dict[value]
        
        window_size = bcell.window if options.window_size is False else options.window_size
        center_position = "%d" %round(int(window_size)/2.0)
        
	window_size='9'
  	center_position='5'
	lista=[]
	seqs=''
        if options.swissprot:
            sequence_text = self.get_sequence(args[0])
	    #print(sequence_text)  
        else:
            with open(args[0], "r") as infile:
                #print(infile)
		seqs=infile.read().splitlines( )

	#print(len(seqs))
	for l in seqs:

		sequence_text =l.split(' ')[0]
		
		if len(sequence_text)==9:   
			#print(sequence_text)
			if value == 'Emini':
			    results = bcell.emini_method(value, sequence_text, window_size, center_position)
			elif value == 'Karplus-Schulz':
			    results = bcell.karplusshulz_method(value, sequence_text, window_size, center_position)
			elif value == 'Kolaskar-Tongaonkar':
			    results = bcell.kolaskartongaonkar_method(value, sequence_text, window_size, center_position)
			elif value == 'Bepipred':
			    results = bcell.bepipred_method(value, sequence_text, window_size, center_position)
			else:
			    results = bcell.classical_method(value, sequence_text, window_size, center_position)
			threshold = round(results[1][0], 3)
		
			with open(results[0], "r") as infile:
			    f=open(results[0], "r")
			    l=0
			    for line in f:
				if not l==0:
				    #print infile.read().replace(',','\t')
			    	    #print infile.read()
				    #print(line.split(',')[4]+" "+line.split(',')[5])
			    	    lista.append([line.split(',')[4],line.split(',')[5]])
			        l=l+1
				    
			    lista.sort(key=lambda tup: tup[-1], reverse=True)

			    

			    #comando="/home/proyecto/Documentos/PDB/pymol/crearPDB.py /home/proyecto/Documentos/resultados/resultadocbell.txt > /home/proyecto/Documentos/resultados/resultadocbell.txt"
			    #resultado = os.system(comando)
	
			#print(results)
			'''
			# generate a plot if 'noplot' option is not invoked
			if options.noplot == False:
			    self.create_png(results, value, threshold)
			    print "* A plot has been generated in the 'output' directory."
			'''        

			os.remove(results[0])

	for i in range(0,len(lista)):  
	    	print(lista[i])
	print ""
        
    def get_sequence(self, id):
        from urllib import urlopen
        if re.search( r'^[A-Z][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9]', id):
            url = "http://www.uniprot.org/uniprot/%s.fasta" %id
            fasta = urlopen(url).read()
            input_sequences = fasta.split('\n') 
        else: 
            print "* Please check the input ID."
            exit(0)
        return ''.join(input_sequences[1:])
        
    def create_png(self, results, method_name, threshold):
        
        content = [tuple(content.split(',')) for content in open(results[0]).read().splitlines()]
      
        x=[]
        y=[]
        for ele in content[1:]:
            x.append(int(ele[0]))
            if 'Bepipred' in method_name:
                y.append(float(ele[-2]))
            else:
                y.append(float(ele[-1]))
        
        x = np.array(x)
        y = np.array(y)
        plot(x, y)
        xlabel("Position", fontsize=11)
        ylabel("Score", fontsize=11)
        
        axhline(y=threshold, color='r', linewidth=1.5, label="Threshold")
        legend(loc='upper right', prop={'size':9})
        
        fill_between(x, y, np.float(threshold), where=y>np.float(threshold), color='#FFFF00', interpolate=True)
        fill_between(x, y, np.float(threshold), where=y<np.float(threshold), color='#00CC00', interpolate=True)
        
        grid(True)
        savefig("./output/plot.png")
        
    
    def method_list(self):
        print """
* All available method names:
-----------------------------
1) Chou-Fasman
2) Emini
3) Karplus-Schulz
4) Kolaskar-Tongaonkar
5) Parker
6) Bepipred
"""
        
        
    def commandline_help(self):
        print """
* All available method names:
-----------------------------
Chou-Fasman, Emini, Karplus-Schulz, Kolaskar-Tongaonkar, Parker, Bepipred

1) Make predictions given a file containing a sequence:
-------------------------------------------------------
python predict_bcell.py -m [method-name] [input-file]
Example: python predict_bcell.py -m Chou-Fasman example/test.txt

2) Make predictions given a SwissProt ID:
-----------------------------------------
python predict_bcell.py -m [method-name] --swissprot [swissprot-id]
Example: python predict_bcell.py -m Chou-Fasman --swissprot P02185

You can also use help option (-h or --help) for more information:
*****************************************************************
python predict_bcell.py --help
"""

    def main(self):
        
        parser = OptionParser(usage="usage: %prog -m <method_name> [options] input", version="%prog 1.0")
        
        parser.add_option("-w", "--window",
                      action="store_true",
                      dest="window_size",
                      default=False,
                      help="Sets window size if specified \
                      (default=method specific. Eg: 6,7,...)")
        parser.add_option("-m", "--method",
                          type="choice",
                          dest="method",
                          choices=["Chou-Fasman","Emini","Karplus-Schulz","Kolaskar-Tongaonkar","Parker","Bepipred",],
                          help="Select a method from available method options.",)
        parser.add_option('-s', '--swissprot', 
                          action='store_true', 
                          help="Use when the input is a SwissProt ID (default=file).")
        parser.add_option('-l', '--list', 
                          action='store_true', 
                          help="Show all available method options.")
        parser.add_option('-n', '--noplot', 
                          action='store_true', 
                          default=False, 
                          help="Do not generate a plot.")
        
        (options, args) = parser.parse_args()
        
        if len(args) == 0 and options.list == True: self.method_list()
        elif len(args) == 0: self.commandline_help()
        elif len(args) == 1: self.predict(options, args)
        else: 
            parser.error("Incorrect number of arguments.")
            self.commandline_help()
            
    
if __name__ == '__main__':
    bcell = BCell()
    Prediction().main()

