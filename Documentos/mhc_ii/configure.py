#! /usr/bin/python
import os
import sys
import re
def main():
	curDir = os.getcwd()
	netMHC_bak = curDir + "/netMHCII-1.1/netMHCII_bak"
	smmDir = curDir + "/netMHCII-1.1"
	oldScript = open(netMHC_bak, "r").read()
	smm = re.compile("/home/dorjee/standalone_tools/mhc_ii/netMHCII-1.1")
	newScript = smm.sub(smmDir, oldScript)
	netMHC = curDir + "/netMHCII-1.1/netMHCII"	
	newFile = open(netMHC, "w")
	newFile.write(newScript)
	newFile.close()

	netMHC_bak = curDir + "/netMHCII-2.2/netMHCII-2.2_bak"
	nnDir = curDir + "/netMHCII-2.2"
	oldScript = open(netMHC_bak, "r").read()
	nn = re.compile("/home/dorjee/standalone_tools/mhc_ii/netMHCII-2.2")
	newScript = nn.sub(nnDir, oldScript)
	netMHC = curDir + "/netMHCII-2.2/netMHCII-2.2"	
	newFile = open(netMHC, "w")
	newFile.write(newScript)
	newFile.close()

# 	netMHCIIpan-3.0
	pan_bak = curDir + "/netMHCIIpan-3.0/netMHCIIpan_bak"
	panDir = curDir + "/netMHCIIpan-3.0"
	oldScript = open(pan_bak, "r").read()
	pan = re.compile("/opt/netMHCIIpan-3.0")
	newScript = pan.sub(panDir, oldScript)
	netMHC = curDir + "/netMHCIIpan-3.0/netMHCIIpan"	
	newFile = open(netMHC, "w")
	newFile.write(newScript)
	newFile.close()


	smm_dir = curDir + "/netMHCII-1.1"
	smmFile = curDir + "/smm_dir.txt"
	oFile = open(smmFile, "w")
	oFile.write(smm_dir)
	oFile.close()

	nn_dir = curDir + "/netMHCII-2.2"
	nnFile = curDir + "/nn_dir.txt"
	oFile = open(nnFile, "w")
	oFile.write(nn_dir)
	oFile.close()

	pan_dir = curDir + "/netMHCIIpan-3.0"
	panFile = curDir + "/net2_dir.txt"
	oFile = open(panFile, "w")
	oFile.write(pan_dir)
	oFile.close()

	tool_data = curDir + "/tool_data/MHCII"
	binding_bak = curDir + "/mhc_II_binding.py.temp"	
	oldScript = open(binding_bak, "r").read()	
	newScript = oldScript.replace("tool_data/MHCII", tool_data).replace("smm_dir.txt", smmFile).replace("nn_dir.txt", nnFile).replace("net2_dir.txt", panFile)
	binding_script = curDir + "/mhc_II_binding.py"	
	newFile = open(binding_script, "w")
	newFile.write(newScript)
	newFile.close()


if __name__ == '__main__':   
	main()
