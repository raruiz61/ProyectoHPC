
from Tkinter import *
from pymol import cmd
from pymol import querying
from pymol import computing
from pymol import importing
import sys
import httplib
import Pmw
import random
import string
import os
from treewidgets.constants import *
from treewidgets.texttree import TextTree
from treewidgets.node import TWTreeNode
from treewidgets.node import CustomTreeNode
import json
import traceback
import copy
import base64
from epymol import mae
from epymol.mae import dotmae
from epymol.mae.dotmae import MAEReader, MAEReaderException
import threading

ldir = os.path.dirname(__file__)

local_scripts_data = [ [ "Color By Potential Energy", "color_by_pe", os.path.join(ldir, "color_by_pe.py"), "color_by", [ "ligprep" ] ],
                       [ "Color By Molecular Weight", "color_by_mw", os.path.join(ldir,"color_by_mw.py"), "color_by", [ "Molecule Descriptors"] ],
                       [ "Prime Spectrum Energy", "prime_spectrum_energy", os.path.join(ldir,"prime_spectrum_energy.py"), "prime_spectrum_energy", [ "Prime Analyze Energy" ] ],
                       ]

local_scripts_data_by_name = {}

for k in local_scripts_data:
    local_scripts_data_by_name[k[0]] = k

servers_installed = [ ["localhost", "localhost", 7970, "/resources/", True ],
#                      ["clara", "clara", 7970, "/resources/", True ],
#                      ["june", "june", 7970, "/resources/", True ],
#                      ["VM on clara", "38.96.144.34", 7970, "/resources/", True ],
                      ["nyc-desk-l08", "nyc-desk-l08", 7970, "/resources/", True],
                      ["pymol-test3.dev.bb.schrodinger.com", "pymol-test3.dev.bb.schrodinger.com", 7970, "/resources/", True] ]

def hide_all(nodeobj):
    nodeobj.state = nodeobj.state | NS_PENDING_HIDE
    map(hide_all, nodeobj.children)

def hide_all_nodes(texttree, nodeobj):
    for child in nodeobj.children:
        hide_all_nodes(texttree, child)
    texttree.hideNode(nodeobj)

def file_as_base64(file_name):
    f = open(file_name)
    file_base64 = base64.b64encode("".join(f.readlines()))
    f.close()
    return file_base64

def base64_to_file(file_name, file_base64):
    f = open(file_name, "w")
    data = base64.b64decode(file_base64)
    f.write(data)
    f.close()

def binary_to_file(file_name, file_binary):
    f = open(file_name, "w")
    f.write(file_binary)
    f.close()

def id_generator(size=6, chars=string.ascii_lowercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

server_type_to_type = { "application/vnd.schrodinger.pws.Directory+json" : "Directory",
                        "text/plain" : "Text",
                        "application/vnd.schrodinger.pws.pse" : "PyMOL Session",
                        "application/pdf" : "PDF",
                        "text/pdb" : "PDB",
                        "image/x-png" : "PNG",
                        "image/png" : "PNG",
                        "image/jpeg" : "JPG"}

icons = {
    'folder': """
R0lGODlhDwAPAPIEAAAAADMzM8zMmf//zP///////wAAAAAAACH5BAEBAAUALAAAAAAPAA8A
AAOEWFVVhVVVVVhVVYVVVVVYVVWFVQAACFBVhVVVICgiIoAAAFBYAACAADMzOBAFgzMzMzgz
A4EFMzM4MzODA1EwODMzgzMzEFgwM4MzMzM4EAWDMzMzODMDgQUzMzgzM4MDUQAIAACAAAAQ
WBURgRERERgRVYVVVVVYVVWFVVVVWFVVhVWVADs=""",
    'folderminus': """
R0lGODlhDwAPAPIEAAAAADMzM8zMmf//zP///////wAAAAAAACH5BAEBAAUALAAAAAAPAA8A
AAOEWFVVhVVVVVhVVYVVVVVYVVWFVQAACFBVhVVVICgiIoAAAFBYAACAADMzOBAFgzMzMzgz
A4EFMzM4MzODA1EwOAAAgAAzEFgwM4MzMzM4EAWDMzMzODMDgQUzMzgzM4MDUQAIAACAAAAQ
WBURgRERERgRVYVVVVVYVVWFVVVVWFVVhVWVADs=""",
    'folderplus': """
R0lGODlhDwAPAPIEAAAAADMzM8zMmf//zP///////wAAAAAAACH5BAEBAAUALAAAAAAPAA8A
AAOEWFVVhVVVVVhVVYVVVVVYVVWFVQAACFBVhVVVICgiMoAAAFBYAACAAzMzOBAFgzMzMDgz
A4EFMzM4MDODA1EwOAAAgAAzEFgwM4MDMzM4EAWDMzMwODMDgQUzMzgwM4MDUQAIADKCAAAQ
WBURgRERERgRVYVVVVVYVVWFVVVVWFVVhVWVADs=""",
    'doctype': """
R0lGODlhDwAPAPECABmTAP///wAAAAAAACH5BAEBAAIALAAAAAAPAA8AAAKAlChRokSJEiVK
lChRokSJEiVKlCgRoESAEiVKlAhQokCIEiVKBChRIkCJEiUClChRIESJEgFKlCgRoESJACVK
lCgQokSJACVKlAhQokSJACVKFAhRokSJACVKlChRokSJACUKhChRokSJACUClChRokSJEiVK
lChRokSJEiVKVAEAOw==""",
    'doctypeminus': """
R0lGODlhDwAPAPIEAAAAABmTAI3Jgf///wAAAAAAAAAAAAAAACH5BAEBAAQALAAAAAAPAA8A
AAOESEREhEREREhERIRERERIRESEREREGEEUgUREREgUQYQUQURIRBSBREQRSEREhBFEREgR
RIREIURIRCSCREQRSAAAgABCREgUQoRERCJIRESEEURESBFEhEREEUhERIRERERIFEGEFEFE
SEREhBFEEUhERIRERERIRESERERESEREhESUADs=""",
    'doctypeplus': """
R0lGODlhDwAPAPIEAAAAABmTAI3Jgf///wAAAAAAAAAAAAAAACH5BAEBAAQALAAAAAAPAA8A
AAOESEREhEREREhERIRERERIRESEREREGEEUgUREREgUQoQUQURIRBSBBEQRSEREhBFEQEgR
RIREIURIQCSCREQRSAAAgABCREgUQoQERCJIRESEEURASBFEhEREEUhARIRERERIFEKEFEFE
SEREhBFEEUhERIRERERIRESERERESEREhESUADs=""",
    'page': """
R0lGODlhDwAPAPIDAAAAADMzM5mZ/////////wAAAAAAAAAAACH5BAEBAAQALAAAAAAPAA8A
AAOESEREhEREREhERIAAAAAIAESERDAzODMzgEFEBDgiIoIyEERIBDODMzMDGEREgCMiIigD
QYREMDM4MzOAQUQEOCIigjIQREgEM4MzMwMYRESAIyIiKANBhEQwMzgzM4BBRAQ4IiKCMhBE
SAQzgzMzAxhERIAAAAAIAEGERBQRGBERgUGUADs=""",
    'pageminus': """
R0lGODlhDwAPAPIDAAAAADMzM5mZ/////////wAAAAAAAAAAACH5BAEBAAQALAAAAAAPAA8A
AAOESEREhEREREhERIAAAAAIAECERDAzODMzgxBEBDgiIoIiA0FIBDODMzMzCEFEgCMiIigy
EIREMDM4MzODEEQEOAAAgAADQUgEM4MzMzMIQUSAIyIiKDIQhEQwMzgzM4MQRAQ4IiKCIgNB
SAQzgzMzMwhBRIAAAAAIABCERBQRGBERgRGUADs=""",
    'pageplus': """
R0lGODlhDwAPAPIDAAAAADMzM5mZ/////////wAAAAAAAAAAACH5BAEBAAQALAAAAAAPAA8A
AAOESEREhEREREhERIAAAAAIAECERDAzODMzgxBEBDgiMoIiA0FIBDODAzMzCEFEgCMyMCgy
EIREMDM4MDODEEQEOAAAgAADQUgEM4MDMzMIQUSAIzIwKDIQhEQwMzgwM4MQRAQ4IjKCIgNB
SAQzgzMzMwhBRIAAAAAIABCERBQRGBERgRGUADs=""",
    'pdf' : """
R0lGODlhEAAQAPIAAPI2OOfm57KysrCwsCYmJsfHx4sBA////yH5BAEAAAcALAAAAAAQABAAAANdeBoi4+MFRU+5osSh5xFAKIbGEATmMa7liQbraBADbQJGrutn4QiUywXmGgIVBRexyAgGkjAB82hxmSQuKjRaEJm0J8D2CWgio8qe2QIwKYU+SvqdHAQJ+LyeUEgAADs=
""",
    'pymol' : """
R0lGODlhEAAQAMIHAC5DRDtbiGN6dZei5nav2oSytun4+////yH+GkNyZWF0ZWQgd2l0aCBHSU1QIG9uIGEgTWFjACH5BAEKAAcALAAAAAAQABAAAANZeLbcrSZIQqUNTI5RitmgtAQCAHSFqQpYdBHFZY0qIHSCsLBjAPgeR4xm0jkMvMjKUWz5TEHGLRnYwE4oCmUBGtQA2u1n81SFCYbK10wYUH7HIyFwWNsBhwQAOw=="""
    }

def dummyFunc(a):
    pass

def get_extension(fn):
    filename, ext = os.path.splitext(fn)
    if ext in [ '.gz' ]:
        filename2, ext2 = os.path.splitext(filename)
        if ext2:
            ext = "%s%s" %(ext2, ext)
            filename = filename2
    return filename, ext

def recursivelyRemoveChildren(texttree, nodeobj):
    for child in nodeobj.data['children']:
        tn = texttree.getNodeByID(str(child['id']))
        nodeobj.state = (nodeobj.state | NS_PENDING_HIDE)
        texttree.removeNode(tn)
        recursivelyRemoveChildren(texttree, tn)
        print "child.id=", child['id'], " tn=", tn
    nodeobj.data['children'] = []

def showContentCategory(serverobj, nodeobj):
    if (nodeobj.state & NS_EXPANDED) ^ NS_EXPANDED:
        # not expanded
        nodeobj.expandShallow()
    nodeobj.state = (nodeobj.state | NS_PENDING_SHIFT) 
    serverobj.texttree.scheduleShiftFollowing(nodeobj)
    serverobj.texttree.redisplay()
    serverobj.texttree.focus_set()

def showContentServer(serverobj, nodeobj):
    if (nodeobj.state & NS_EXPANDED) ^ NS_EXPANDED:
        # not expanded
        nodeobj.expandShallow()
        nodeobj.state = (nodeobj.state | NS_PENDING_SHIFT) 
        serverobj.texttree.scheduleShiftFollowing(nodeobj)
        serverobj.texttree.redisplay()
        serverobj.texttree.focus_set()
    else:
        serverobj.addPyMOLServer(edit_server=nodeobj)

def showContentDirectoryImpl(serverobj, nodeobj, services=False):
    sel = False
    fetched = nodeobj.data['fetched']
    if not fetched:
        success = serverobj.connect_and_get(nodeobj, nodeobj.data['server'], nodeobj.data['port'], nodeobj.data['uri'], services)
        if success:
            nodeobj.data['fetched'] = True
        else:
            return
    else:
        if (nodeobj.state & NS_EXPANDED) ^ NS_EXPANDED:
            # not expanded
            nodeobj.expandShallow()
            nodeobj.state = (nodeobj.state | NS_PENDING_SHIFT) 
            serverobj.texttree.scheduleShiftFollowing(nodeobj)
        else:
            nodeobj.collapseAll()
            nodeobj.state = (nodeobj.state | NS_PENDING_SHIFT) 
            serverobj.texttree.scheduleShiftFollowing(nodeobj)

    serverobj.texttree.redisplay()
    serverobj.texttree.focus_set()

def showContentDirectory(serverobj, nodeobj):
    showContentDirectoryImpl(serverobj, nodeobj)

def showContentServices(serverobj, nodeobj):
    showContentDirectoryImpl(serverobj, nodeobj, services=True)

undo_all_molecules_at_once = True
undo_each_molecule = False
already_processed = set([])

def process_molecule_sdf(mol, suspend_undo):
    global already_processed
    global undo_each_molecule
    nameorig = mol[0]
    name = nameorig.strip()
    if name in already_processed:
        return True
    already_processed = set(already_processed | set([name]))
    tmpname = "_temp_%s" % name
    molstr = string.join(mol,'')
    cmd.delete(tmpname)
    try:
        importing.read_sdfstr(molstr, tmpname,object_props='*')
        natoms = cmd.count_atoms(tmpname)
    except:
        return False
    if undo_each_molecule:
        cmd.push_undo( "( %s )" % name, just_coordinates=0)
        cmd.set("suspend_undo", 1, updates=0)
    cmd.remove(name)

    print "before cmd.create(): using copy_properties=True"
    cmd.create(name, tmpname, copy_properties=True)
    cmd.set_title(name, -1, "")

    cmd.delete(tmpname)
    if undo_each_molecule:
        cmd.set("suspend_undo", suspend_undo, updates=0)
        cmd.push_undo("", just_coordinates=0, finish_undo=1)
    return True

def process_output_sdf(outFile, suspend_undo):
    global undo_all_molecules_at_once
    global already_processed
    print "process_output_sdf: outFile=", outFile, " suspend_undo=", suspend_undo
    already_processed = set([])
    if undo_all_molecules_at_once:
        pushstr = ""
        it = iter(querying.get_object_list("enabled"))
        pushstr = it.next()
        for name in it:
            pushstr = "%s or %s" % (pushstr, name)
        cmd.push_undo( "( %s )" % pushstr, just_coordinates=0)
        cmd.set("suspend_undo", 1, updates=0)
    f2 = open(outFile)
    curmol = []
    success = True
    for l in f2.readlines():
        curmol.append(l)
        if l.find("$$$$") >=0:
            success = process_molecule_sdf(curmol, suspend_undo ) and success
            curmol = []
    if len(curmol) > 0:
        success = process_molecule_sdf(curmol, suspend_undo ) and success
    if undo_all_molecules_at_once:
        cmd.set("suspend_undo", suspend_undo, updates=0)
        cmd.push_undo("", just_coordinates=0, finish_undo=1)
    return success


def read_mae_model(maestr):
    mr = MAEReader(1, "*", "*")
    mol_list = mr.listFromStr(maestr)
    len_list = len(mol_list);
    if len_list:
        title = getattr(mol_list[0], 'title', None)
        return title, mol_list[0]

def process_molecule_mae(file_header, mol, suspend_undo):
    global already_processed
    global undo_each_molecule
    tot = []
    tot.extend(file_header)
    tot.extend(mol)
    molstr = string.join(tot,'')

    nameorig, mdl = read_mae_model(molstr)
    name = nameorig.strip()
    if name in already_processed:
        return True
    already_processed = set(already_processed | set([name]))
    tmpname = "_temp_%s" % name
    cmd.delete(tmpname)

    try:
        cmd.load_model(mdl, tmpname)
        natoms = cmd.count_atoms(tmpname)
    except:
        return False
    if undo_each_molecule:
        cmd.push_undo( "( %s )" % name, just_coordinates=0)
        cmd.set("suspend_undo", 1, updates=0)
    cmd.remove(name)
    cmd.create(name, tmpname, copy_properties=True)
    cmd.set_title(name, -1, "")

    cmd.delete(tmpname)
    if undo_each_molecule:
        cmd.set("suspend_undo", suspend_undo, updates=0)
        cmd.push_undo("", just_coordinates=0, finish_undo=1)
    return True

def process_output_mae(outFile, suspend_undo):
    global undo_all_molecules_at_once
    global already_processed
    print "process_output_mae() called : outFile=", outFile
    file_header = None
    already_processed = set([])
    if undo_all_molecules_at_once:
        pushstr = ""
        it = iter(querying.get_object_list("enabled"))
        pushstr = it.next()
        for name in it:
            pushstr = "%s or %s" % (pushstr, name)
        cmd.push_undo( "( %s )" % pushstr, just_coordinates=0)
        cmd.set("suspend_undo", 1, updates=0)
    f2 = open(outFile)
    curmol = []
    success = True
    for l in f2.readlines():
        if l.find("f_m_ct") >=0:
            if not file_header:
                file_header = curmol
            else:
                success = process_molecule_mae(file_header, curmol, suspend_undo ) and success
            curmol = []
        curmol.append(l)
    if len(curmol) > 0:
        success = process_molecule_mae(file_header, curmol, suspend_undo ) and success
    if undo_all_molecules_at_once:
        cmd.set("suspend_undo", suspend_undo, updates=0)
        cmd.push_undo("", just_coordinates=0, finish_undo=1)
    return success

process_output = {}

process_output['.sdf'] = process_output_sdf
process_output['.mae'] = process_output_mae

def showContentService(serverobj, nodeobj):
    server_name = nodeobj.data['server_name']
    service_name = nodeobj.data['name']
    jobinfo = nodeobj.data['jobinfo']
    serverobj._runServiceDialog = Pmw.MessageDialog(serverobj.dialog.component('hull'),
                                               title = 'Run: %s' % service_name,
                                               defaultbutton = 2,
                                               buttons = ('Cancel', 'Run'))
    hull = serverobj._runServiceDialog.component('hull')
    frame = Frame(hull, height=10)
    parent = frame

    source_options = [ "Enabled in PyMOL", "All Entries in PyMOL" ]
    serverobj._sourcevar = StringVar()
    source_combo = Pmw.ComboBox(parent,
                                label_text = 'Source:',
                                labelpos = 'w',
                                entry_textvariable = serverobj._sourcevar,
                                scrolledlist_items = source_options,
                                dropdown = 1)
    source_combo.component('entryfield').component('entry').configure(state = 'readonly')
    source_combo.selectitem(0)
    source_combo.pack(fill='x', expand=1, padx=10, pady=5)


    local_script_options = [ "" ]
    for data in local_scripts_data:
        if service_name in data[4]:
            local_script_options.extend([data[0]])
    serverobj._localscriptvar = StringVar()
    local_script_combo = Pmw.ComboBox(parent,
                                      label_text = 'Script:',
                                      labelpos = 'w',
                                      entry_textvariable = serverobj._localscriptvar,
                                      scrolledlist_items = local_script_options,
                                      dropdown = 1)
    local_script_combo.component('entryfield').component('entry').configure(state = 'readonly')
    local_script_combo.selectitem(0)
    local_script_combo.pack(fill='x', expand=1, padx=10, pady=5)


    parent.pack(expand=YES,fill=BOTH)

    ret = serverobj._runServiceDialog.activate()

    if ret == "Run":
        uid = id_generator()
        tmp_dir = "/tmp/rest_%s" % uid
        ensure_dir(tmp_dir)
        fileName = '%s/_testing.sdf' % tmp_dir

        source = serverobj._sourcevar.get()
        f = open(fileName, "w")
        object_list = None
        if source == "Enabled in PyMOL":
            object_list = querying.get_object_list("enabled")
        elif source == "All Entries in PyMOL":
            object_list = querying.get_object_list()
        for name in object_list:
            input_model = cmd.get_model(name,state=-1)
            (fit_flag, sdf_list) = computing.model_to_sdf_list(cmd,input_model)
            sdf_list[0] = "%s\n" % name
            input_sdf = string.join(sdf_list,'')
            f.write(input_sdf)
        f.close()

        post_data = copy.deepcopy(jobinfo)

        input_file_base64 = file_as_base64(fileName)
        inputfiles = post_data['INPUTFILES']
        for inputf in inputfiles.keys():
            inputfiles[inputf] = [ input_file_base64, inputfiles[inputf] ]

        jobname = "%s_%s" % (service_name, uid)
        post_data = [ [ jobname, post_data ] ]

        print "running %s molecules with service %s on %s" % (len(object_list), service_name, server_name)
        serverobj._statusvar.set("running %s molecules with service %s on %s" % (len(object_list), service_name, server_name))

        class SendAndReceivePostThread(threading.Thread):
            def __init__(self, nodeobj):
                threading.Thread.__init__(self)
                self.nodeobj = nodeobj
            def run(self):
                try:
                    conn = httplib.HTTPSConnection(self.nodeobj.data['server'], self.nodeobj.data['port']) #, timeout=5
                    
                    conn.request("POST", self.nodeobj.data['uri'], str(post_data))
                    response = conn.getresponse()
                    response_data = response.read()
                    conn.close()
                except:
                    print "WARNING: showContentService cannot communicate with server '", self.nodeobj.data['server'], " at port ", self.nodeobj.data['port']
                    return
                try:
                    retobj = json.loads(response_data)
                except:
                    print "WARNING: showContentService: reading bad response_data=", response_data
                    return False

                suspend_undo = cmd.get("suspend_undo")

                success = True
                for jobname, resultsjob in retobj:
                    output_files = resultsjob["OUTPUTFILES"]
                    output_file_base64, output_file_ext = output_files["<output_file_1>"]
                    outFileName = '%s/_outputfile%s' % (tmp_dir, output_file_ext)
                    binary_to_file(outFileName, output_file_base64)
                    try:
                        success = process_output[output_file_ext](outFileName, suspend_undo) and success
                    except:
                        print "Error: process_output: output_file_ext=", output_file_ext
                        success = False
                print "success=", success
                if not success:
                    dialog1 = Pmw.MessageDialog(serverobj.dialog.component('hull'),
                                                title = 'Service: Unreadable Outputs',
                                                message_text = 'Some output molecules were not loaded.' )
                    dialog1.iconname('Service: Unreadable Outputs')
                    dialog1.activate()
                
                if serverobj._localscriptvar.get() in local_scripts_data_by_name:
                    scripts_data = local_scripts_data_by_name[serverobj._localscriptvar.get()]
                    print "should be executing this script scripts_data=", scripts_data, " object_list=", object_list
                    import imp
                    mod = imp.new_module(scripts_data[1])
                    execfile(scripts_data[2], mod.__dict__)
                    at = getattr(mod, scripts_data[3])
                    print "at=", at
                    at(object_list)

                serverobj._statusvar.set("Finished running %s molecules with service %s on %s" % (len(object_list), service_name, server_name))
                print "_runServiceDialog done: ret=", ret
                print "\tserverobj._sourcevar=", serverobj._sourcevar.get()

        t = SendAndReceivePostThread(nodeobj)
        t.start()

def getContentFileAndSave(serverobj, nodeobj, partial=0):
    server = nodeobj.data['server']
    port = nodeobj.data['port']
    uri = nodeobj.data['uri']
    print "server=", server
    conn = httplib.HTTPSConnection(server, port, timeout=5)
    print "uri=", uri
    conn.request("GET", uri)
    response = conn.getresponse()
    response_data = response.read()
    conn.close()
    base, extension = get_extension(uri)
    base = os.path.basename(base)
    tmpfilename = "tmp_%s%s" % (id_generator(5), extension)
    f = open(tmpfilename, "w")
    f.write(response_data)
    f.close()
    return base, tmpfilename

def showContentSession(serverobj, nodeobj, partial=0):
    base, tmpfilename = getContentFileAndSave(serverobj, nodeobj)
    cmd.load(tmpfilename, object=base, partial=partial)

def showContentSessionPartial(serverobj, nodeobj):
    showContentSession(serverobj, nodeobj, partial=1)

def showContentOpenFile(serverobj, nodeobj):
    base, tmpfilename = getContentFileAndSave(serverobj, nodeobj)
    os.system("open %s" % tmpfilename)


def findLastChildOf(nodeobj):
    retnode = nodeobj
    while len(retnode.children)>0:
        retnode = retnode.children[-1]
    return retnode
            
showContentFunctions = { "Category" : showContentCategory, "Server" : showContentServer, "Directory" : showContentDirectory, "Services" : showContentServices, "Service" : showContentService, "PyMOL Session" : showContentSession, "PDB" : showContentSessionPartial, "PDF" : showContentOpenFile, "Text" : showContentOpenFile, "PNG" : showContentOpenFile, "JPG" : showContentOpenFile }

def nfGetType(node):
    return node['type']

def nfGetID(node):
    return str(node['id'])

def nfGetName(node):
    return node['name']

def nfSetName(node, newname):
    node['name'] = newname

def nfGetChildren(node):
    return node['children']

def nfAddChild(node, newchild):
    print "nfAddChild called newchild=", newchild
    node['children'].append(newchild)

def nfDelChild(node, badchild):
    if badchild in node['children']:
        node['children'].remove(badchild)

def showAtts(node):
    print "showAtts node=", node

def glimpse(node):
    pass
#    print "glimpse node=", node

def unGlimpse(node):
    pass
#    print "unGlimpse node=", node

class PrintOne:
    def __init__(self, text):
        self.text = text

    def __call__(self):
        print self.text

def __init__(self):
    # Simply add the menu entry and callback
    print "server.py.__init__ called"
    self.menuBar.addmenuitem('Plugin', 'command',
                             'PyMOL Server',
                             label = 'PyMOL Server',
                             command = lambda s=self : PyMOLServer(s))

pymolServer = None

class PyMOLServer:
    def execute(self, result, refocus=True):
        if result == 'Exit PyMOL Server':
            self.quit()
        elif result == 'Clear':
            try:
                current = self.texttree.getNodeByID(self.texttree.selected)
                if current and (current.state & NS_EXPANDED) & NS_EXPANDED:  # expanded
                    current.data['fetched'] = False
                    current.collapseAll()
                    current.state = (current.state | NS_PENDING_SHIFT) 
                    self.texttree.scheduleShiftFollowing(current)
                    for child in current.data['children']:
                        recursivelyRemoveChildren(self.texttree, child)
                    self.texttree.redisplay()
                    self.showPyMOLContent(current)
            except:
                print "self=", self
                pass

    def quit(self):
        global pymolServer
        pymolServer = None
        self.dialog.destroy()
    def connect_and_get(self, nodeobj, server, port, current_dir, services):
        global server_type_to_type
        nodeobj.data['server']

        try:
            conn = httplib.HTTPSConnection(nodeobj.data['server'], nodeobj.data['port'], timeout=5)
            #        print "\tconnect_and_get: current_dir=", current_dir, " nodeobj=", nodeobj, " nodeobj.id=", nodeobj.id, " nodeobj.data=", nodeobj.data
            print "server=", nodeobj.data['server'], " port=", nodeobj.data['port'] , " current_dir=", current_dir
            conn.request("GET", current_dir)
            response = conn.getresponse()
            response_data = response.read()
            conn.close()
        except:
            print "WARNING: cannot communicate with server '", nodeobj.data['server'], " at port ", nodeobj.data['port']
            return
        try:
            retobj = json.loads(response_data)
        except:
            print "unsuccessful json load: response_data=", response_data
            return False
        if services:
            ret = self.build_children_for_services(nodeobj, retobj)
            return ret
        else:
            ret = self.build_children_for_directory(nodeobj, retobj)
            return ret

    def build_children_for_services(self, nodeobj, retobj):
        lastn = nodeobj
        for servicedata in retobj:
            service_name = servicedata[0]
            self.num_nodes = self.num_nodes + 1
            data = {'type': 'Service',
                    'id' : int(self.num_nodes),
                    'fetched' : False,
                    'name': service_name, 
                    'server_name': nodeobj.data['server_name'],
                    'server' : nodeobj.data['server'],
                    'port'  :  nodeobj.data['port'],
                    'uri' : "/services",
                    'jobinfo' : servicedata[1],
                    'children': [] }
            tn = CustomTreeNode(self.texttree, nodeobj, data, 1,
                                node_funcs={'type': nfGetType,
                                            'id': nfGetID,
                                            'get_name': nfGetName,
                                            'set_name': nfSetName,
                                            'get_children': nfGetChildren,
                                            'add_child': nfAddChild,
                                            'del_child': nfDelChild},
                                state=NS_PENDING_SHOW)
            nodeobj.data["children"].append(data)
            nodeobj.state = (nodeobj.state | NS_HAS_CHILDREN | NS_EXPANDED)
            nodeobj.children.append(tn)
            self.texttree.removeNode(tn)
            self.texttree.insertAfter(lastn, tn)
            lastn = tn
        nodeobj.state = nodeobj.state | NS_PENDING_SHIFT
        self.texttree.scheduleShiftFollowing(nodeobj)
        return True

            
    def build_children_for_directory(self, nodeobj, retobj):
        directory_info = retobj["Directory"]
        dir_entries = directory_info["Entries"]
        dir_url = directory_info["URI"]
        dir_name = directory_info["Name"]
        dir_defaults = directory_info["Defaults"]
        lastn = nodeobj
        for entry in dir_entries:
            entry_type = entry["Type"]
            entry_name = entry["Name"]
            entry_defaults = entry["Defaults"]
            entry_uri = entry["URI"]
            self.num_nodes = self.num_nodes + 1
            datatype = "None"
            if entry_type in server_type_to_type.keys():
                datatype = server_type_to_type[entry_type]
            else:
                print "WARNING: entry_type=", entry_type, " not known, entry_name=", entry_name
            data = {'type': datatype,
                    'id' : int(self.num_nodes),
                    'fetched' : False,
                    'name': entry_name, # .decode('unicode-escape'),
                    'server_name' : nodeobj.data['server_name'],
                    'server' : nodeobj.data['server'],
                    'port'  :  nodeobj.data['port'],
                    'uri' : "/%s" % entry_uri,
                    'children': [] }
            tn = CustomTreeNode(self.texttree, nodeobj, data, 1,
                                node_funcs={'type': nfGetType,
                                            'id': nfGetID,
                                            'get_name': nfGetName,
                                            'set_name': nfSetName,
                                            'get_children': nfGetChildren,
                                            'add_child': nfAddChild,
                                            'del_child': nfDelChild},
                                state=NS_PENDING_SHOW)
            nodeobj.data["children"].append(data)
            nodeobj.state = (nodeobj.state | NS_HAS_CHILDREN | NS_EXPANDED)
            nodeobj.children.append(tn)
            self.texttree.removeNode(tn)
            self.texttree.insertAfter(lastn, tn)
            lastn = tn
        nodeobj.state = nodeobj.state | NS_PENDING_SHIFT
        self.texttree.scheduleShiftFollowing(nodeobj)
        return True

    def focusVScrollOnNode(self, id):
        r = self.texttree.tag_ranges(id)
        start = r[0]
        self.texttree.see(start)

    def _onKeyRight(self, k):
        if self.texttree.selected:
            self._statusvar.set("")
            current = self.texttree.getNodeByID(self.texttree.selected)
            if current and current.data['type'] in [ 'Directory', 'Category', 'Server', 'Services'] and (current.state & NS_EXPANDED) ^ NS_EXPANDED:  # not expanded
                self.texttree._onReturn()
                self.texttree.redisplay()
            return "break"
    def _onKeyLeft(self, k):
        if self.texttree.selected:
            self._statusvar.set("")
            current = self.texttree.getNodeByID(self.texttree.selected)
            if current and (current.state & NS_EXPANDED) & NS_EXPANDED:  # expanded
                current.collapseAll()
                current.state = (current.state | NS_PENDING_SHIFT) 
                self.texttree.scheduleShiftFollowing(current)
                self.texttree.redisplay()
            return "break"
    def _onKeyDown(self, k):
        if self.texttree.selected:
            self._statusvar.set("")
            current = self.texttree.getNodeByID(self.texttree.selected)
            self.texttree._onDown(k)
            self.focusVScrollOnNode(self.texttree.selected)
            return "break"
    def _onKeyUp(self, k):
        if self.texttree.selected:
            self._statusvar.set("")
            current = self.texttree.getNodeByID(self.texttree.selected)
            self.texttree._onUp(k)
            self.focusVScrollOnNode(self.texttree.selected)
            return "break"

    def testPyMOLServer(self):
        print "testPyMOLServer : alias='", self._aliasvar.get(), "' host='", self._hostvar.get(), "' port=", self._portvar.get()
        try:
            conn = httplib.HTTPSConnection(self._hostvar.get(), self._portvar.get(), timeout=5)
            conn.request("GET", '/resources')
            response = conn.getresponse()
            response_data = response.read()
            conn.close()
            self._servervar.set('Reachable')
        except:
            self._servervar.set('Unreachable')
            pass
#            ei = sys.exc_info()
#            print "sys.exc_info=", ei
#            traceback.print_tb(ei[2])
#            print "\tnot connected"
    def addPyMOLServer(self, edit_server=None):
        title_string = 'Add PyMOL Server'
        action_string = 'Add PyMOL Server'
        default_alias = 'Schrodinger Server'
        default_host = 'mobile-pymol.schrodinger.com'
        default_port = 7970
        default_path = '/resources'
        default_encrypt = True
        default_buttons = ('Cancel', action_string)
        if edit_server:
            title_string = 'Edit PyMOL Server'
            action_string = 'Update'
            data = edit_server.data
            default_alias = data['server_name']
            default_host = data['server']
            default_port = data['port']
            default_path = data['uri']
            default_encrypt = data['encrypt']
            default_buttons = ('Cancel', 'Delete', action_string)

        self._addServerDialog = Pmw.MessageDialog(self.dialog.component('hull'),
                                                  title = title_string,
                                                  defaultbutton = 2,
                                                  buttons = default_buttons)
        hull = self._addServerDialog.component('hull')
        frame = Frame(hull, height=10)
        parent = frame
        self._aliasvar = StringVar()
        self._aliasvar.set(default_alias)
        self._alias = Pmw.EntryField(parent,
                                     labelpos = 'w',
                                     label_text = 'Alias:',
                                     validate = None,
                                     entry_textvariable = self._aliasvar)
        self._hostvar = StringVar()
        self._hostvar.set(default_host)
        self._host = Pmw.EntryField(parent,
                                    labelpos = 'w',
                                    label_text = 'Host:',
                                    validate = None,
                                    entry_textvariable = self._hostvar,
                                    entry_width=30)
        self._portvar = IntVar()
        self._portvar.set(default_port)
        self._port = Pmw.EntryField(parent,
                                    labelpos = 'w',
                                    label_text = 'Port:',
                                    validate = Pmw.integervalidator,
                                    entry_textvariable = self._portvar)
        self._pathvar = StringVar()
        self._pathvar.set(default_path)
        self._path = Pmw.EntryField(parent,
                                    labelpos = 'w',
                                    label_text = 'Path:',
                                    validate = None,
                                    entry_textvariable = self._pathvar)
        self._encrypt = BooleanVar()
        self._encrypt.set(default_encrypt)
        self._encryptbutton = Checkbutton(parent, text='Encryption', variable=self._encrypt)
        entries = (self._alias, self._host, self._port, self._path)
        for entry in entries:
            entry.pack(fill='x', expand=1, padx=10, pady=5)
        Pmw.alignlabels(entries)
        self._encryptbutton.pack(fill='x', expand=1, padx=10, pady=5)

        separator = Frame(parent, height=4, bd=1, relief=SUNKEN)
        separator.pack(fill=X, padx=5, pady=5)

        self._servervar = StringVar()
        self._servervar.set('Unknown')
        self._server = Pmw.EntryField(parent,
                                     labelpos = 'w',
                                     label_text = 'Server:',
                                     validate = None,
                                     entry_textvariable = self._servervar)
        self._server.component('entry').configure(state = 'readonly')
        self._networkvar = StringVar()
        self._networkvar.set('WiFi')
        self._network = Pmw.EntryField(parent,
                                    labelpos = 'w',
                                    label_text = 'Network:',
                                    validate = None,
                                    entry_textvariable = self._networkvar,
                                    entry_width=30)
        self._network.component('entry').configure(state = 'readonly')
        entries = (self._server, self._network)
        for entry in entries:
            entry.pack(fill='x', expand=1, padx=10, pady=5)

        self._testbutton = Button(parent, text='Test Server', command=self.testPyMOLServer)
        self._testbutton.pack(fill='x', expand=1, padx=10, pady=5)
        
        parent.pack(expand=YES,fill=BOTH)
        ret = self._addServerDialog.activate()

        if ret == 'Delete':
            dialog1 = Pmw.MessageDialog(self.dialog.component('hull'),
                                        title = 'Delete Server',
                                        defaultbutton = 1,
                                        buttons = ('No', 'Yes'),
                                        message_text = "Are you sure that you want to delete server '%s'" % self._aliasvar.get() )
            dialog1.iconname('Service: Unreadable Outputs')
            ret = dialog1.activate()
            if ret == 'Yes':
                # deleting Server
                hide_all(edit_server)
                self.texttree.redisplay()
                edit_server.state = edit_server.state | NS_PENDING_HIDE
                self.texttree.scheduleShiftFollowing(edit_server)
                self.texttree.redisplay()
                self._statusvar.set("Deleted Server '%s'" % (self._aliasvar.get()))
        elif ret == action_string:
            # adding PyMOL Server
            # TODO : Should check to make sure alias is unique
            if edit_server:
                # updating the server, first remove children
                self.editNodes(edit_server, [ self._aliasvar.get(), self._hostvar.get(), self._portvar.get(), self._pathvar.get(), self._encrypt.get() ])
                edit_server.name = self._aliasvar.get()
                edit_server_resource, edit_server_services = edit_server.children;
                edit_server_resource.collapseAll()
                edit_server_services.collapseAll()
                edit_server_services.data['uri'] = '/services/'
                map(hide_all, edit_server_resource.children)
                map(hide_all, edit_server_services.children)
                self.texttree.redisplay()
                edit_server_services.data['fetched'] = False
                edit_server_resource.data['fetched'] = False
                edit_server.state = edit_server.state | NS_PENDING_SHIFT
                self.texttree.scheduleShiftFollowing(edit_server)
                edit_server_resource.state = (edit_server_resource.state | NS_HAS_CHILDREN) ^ NS_HAS_CHILDREN
                edit_server_services.state = (edit_server_resource.state | NS_HAS_CHILDREN) ^ NS_HAS_CHILDREN

                self.texttree.redisplay()
            else:
                self.addServer([ self._aliasvar.get(), self._hostvar.get(), self._portvar.get(), self._pathvar.get(), self._encrypt.get() ])

        print "_addServerDialog done ret=", ret
        pass
    def editNodes(self, nodeobj, server_data):
        server_name, server_address, server_port, server_uri, server_encrypt = server_data
        nodeobj.data['server_name'] = server_name
        nodeobj.data['name'] = server_name
        nodeobj.data['server'] = server_address
        nodeobj.data['port'] = server_port
        nodeobj.data['encrypt'] = server_encrypt
        nodeobj.data['uri'] = server_uri
        for child in nodeobj.children:
            self.editNodes(child, server_data)

    def addServer(self, server_data, expand=True):
        server_name, server_address, server_port, server_uri, server_encrypt = server_data
        self.num_nodes += 1
        data = {'type': 'Server',
                'id' : self.num_nodes,
                'fetched' : False,
                'name': server_name,
                'server_name': server_name,
                'server' : server_address,
                'port'  :  server_port,
                'encrypt' : server_encrypt,
                'uri' : server_uri,
                'children': [] }
        self.pymol_servers_node.data["children"].append(data)
        tn = CustomTreeNode(self.texttree, self.pymol_servers_node, data, 1,
                            node_funcs={'type': nfGetType,
                                        'id': nfGetID,
                                        'get_name': nfGetName,
                                        'set_name': nfSetName,
                                        'get_children': nfGetChildren,
                                        'add_child': nfAddChild,
                                        'del_child': nfDelChild},
                            state=NS_PENDING_SHOW)
        tn.state = (tn.state | NS_HAS_CHILDREN | NS_EXPANDED)
        self.pymol_servers_node.children.append(tn)
        self.num_nodes += 1
        data = {'type': 'Directory',
                'id' : self.num_nodes,
                'fetched' : False,
                'name': 'Resources',
                'server_name' : server_name,
                'server' : server_address,
                'port'  :  server_port,
                'encrypt' : server_encrypt,
                'uri' : server_uri,
                'children': [] }
        self.pymol_servers_node.data["children"][-1]['children'].append(data)
        tndir = CustomTreeNode(self.texttree, tn, data, 1,
                               node_funcs={'type': nfGetType,
                                           'id': nfGetID,
                                           'get_name': nfGetName,
                                           'set_name': nfSetName,
                                           'get_children': nfGetChildren,
                                           'add_child': nfAddChild,
                                           'del_child': nfDelChild},
                               state=NS_PENDING_SHOW)
        tn.children.append(tndir)
        self.num_nodes += 1
        data = {'type': 'Services',
                'id' : self.num_nodes,
                'fetched' : False,
                'name': 'Services',
                'server_name' : server_name,
                'server' : server_address,
                'port'  :  server_port,
                'encrypt' : True,
                'uri' : '/services/',
                'children': [] }
        self.pymol_servers_node.data["children"][-1]['children'].append(data)
        tnserv = CustomTreeNode(self.texttree, tn, data, 1,
                                node_funcs={'type': nfGetType,
                                            'id': nfGetID,
                                            'get_name': nfGetName,
                                            'set_name': nfSetName,
                                            'get_children': nfGetChildren,
                                            'add_child': nfAddChild,
                                            'del_child': nfDelChild},
                                state=NS_PENDING_SHOW)
        tn.children.append(tnserv)
        self.pymol_servers_node.state = self.pymol_servers_node.state | NS_PENDING_SHIFT
        self.texttree.scheduleShiftFollowing(self.pymol_servers_node)
        if expand:
            self.texttree.redisplay()

    def close(self):
        if pymolServer:
            pymolServer.dialog.withdraw()

    def __init__(self, app):
        global pymolServer
        if pymolServer:
            pymolServer.dialog.show()
            return

        pymolServer = self
        self.parent = app.root

        self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ("Clear", 'Exit PyMOL Server'),
                                 title = 'PyMOL Server',
                                 command = self.execute)
        self.dialog.protocol('WM_DELETE_WINDOW', self.close)
	menuBar = Pmw.MenuBar(self.dialog.component('hull'),
                              hull_relief = 'raised',
                              hull_borderwidth = 1)
        menuBar.pack(fill = 'x')
        self.menuBar = menuBar

        menuBar.addmenu('Settings','remote settings')
        menuBar.addmenuitem('Settings','command','PyMOL Server',
                            command = self.addPyMOLServer,
                            label = 'Add PyMOL Server')

        self.struct = {'type': 'Category',
                       'name': "PyMOL Servers",
                       'fetched' : True,
                       'id' : 1,
                       'properties': NP_ROOT|NP_ALLOW_CHILDREN|NP_AUTOBUILD,
                       'children': []}
        global servers_installed

        self.dialog.title("PyMOL Server Manager - Proof of Concept")
        iconmap = {'type' : 'photo', # bitmap
                   'Category' : ( icons['folder'], icons['folderplus'], icons['folderminus']),
                   'Server' : ( icons['folder'], icons['folderplus'], icons['folderminus']),
                   'Directory' : ( icons['folder'], icons['folderplus'], icons['folderminus']),
                   'Services' : ( icons['folder'], icons['folderplus'], icons['folderminus']),
                   'PDF' : ( icons['pdf'], icons['pdf'], icons['pdf']),
                   'Text' : ( icons['page'], icons['pageplus'], icons['pageminus'] ),
                   'PyMOL Session' : ( icons['pymol'], icons['pymol'], icons['pymol'] ),
                   'PDB' : ('node', 'nodeplus', 'nodeminus'),
                   'PNG' : ('node', 'nodeplus', 'nodeminus'),
                   'JPG' : ('node', 'nodeplus', 'nodeminus'),
                   'default' : ('node', 'nodeplus', 'nodeminus')
                   }

        self.texttree = TextTree(self.dialog.component('hull'),self.dialog,icons=iconmap,
                                funcs={'showContent': self.showPyMOLContent,
                                       'showAtts': showAtts,
                                       'glimpse': glimpse,
                                       'unGlimpse': unGlimpse})
        self.texttree.pack(expand=YES,fill=BOTH)
        self.texttree.showTree(self.struct,DT_CUSTOM,-1,
                               node_funcs={'type': nfGetType,
                                           'id': nfGetID,
                                           'get_name': nfGetName,
                                           'set_name': nfSetName,
                                           'get_children': nfGetChildren,
                                           'add_child': nfAddChild,
                                           'del_child': nfDelChild},
                               props=NP_AUTOBUILD|NP_ALLOW_CHILDREN,
                               state=NS_NONE)
        self.texttree.all_nodes[0].expandShallow()

        self.pymol_servers_node = self.texttree.all_nodes[0]
        self.pymol_servers_node.state = (self.pymol_servers_node.state | NS_HAS_CHILDREN | NS_EXPANDED)

        self.num_nodes = 1
        for server_data in servers_installed:
            # server_data = server_name, server_address, server_port, server_uri, encrypt
            self.addServer(server_data, expand=False)

        for n in self.texttree.all_nodes[0].children:
            n.expandShallow()
            n.state = n.state | NS_PENDING_SHIFT
            self.texttree.scheduleShiftFollowing(n)
        self.texttree.redisplay()        
        self.texttree.bind('<Right>',self._onKeyRight)
        self.texttree.bind('<Left>',self._onKeyLeft)
        self.texttree.bind('<Down>',self._onKeyDown)
        self.texttree.bind('<Up>',self._onKeyUp)


        self._statusvar = StringVar()
        self._statusvar.set('')
        self._status = Pmw.EntryField(self.dialog.component('hull'),
                                      labelpos = 'w',
                                      label_text = 'Status:',
                                      validate = None,
                                      entry_textvariable = self._statusvar)
        self._status.component('entry').configure(state = 'readonly')
        self._status.pack(fill='x', expand=1, padx=10, pady=5)

        self.dialog.show()
    
    def read_dir(self, retobj):
        directory_info = retobj["Directory"]
        dir_entries = directory_info["Entries"]
        dir_url = directory_info["URI"]
        dir_name = directory_info["Name"]
        dir_defaults = directory_info["Defaults"]
        for entry in dir_entries:
            entry_type = entry["Type"]
            entry_name = entry["Name"]
            entry_defaults = entry["Defaults"]
            entry_uri = entry["URI"]

    def showPyMOLContent(self, nodeobj):
        global showContentFunctions
        nodetype = nodeobj.data['type']
        self._statusvar.set("")
        if nodetype in showContentFunctions.keys():
            fn = showContentFunctions[nodetype]
            fn(self, nodeobj)
