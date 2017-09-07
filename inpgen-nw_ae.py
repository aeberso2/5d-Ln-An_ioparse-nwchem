import os, re, string
import sys, argparse
from numpy import *
import pandas as pd
top_direc = os.getcwd()

from inpgen-nw.py import data_sync, init_params, search_basis

def gen_nw_DKH_dft_input(input_set, dft_direc='dft',optimize='n', Nbas_keys):
    '''
    this program generates a set of DFT functional inputs that run in NWChem
    
    this routine is designed for DK3 all-electron calculations, maybe I'll implement zora, if for some reason someone gets stuck doing a benchmark for those. 
    
    Function arguments
    required:
        input_set: the molecule you are making
    optional arguments:
        dft_direc= 'dft', 'sodft', 'both', default='dft'
        optimize= 'n', 'y' ; default='n'
        Nbas_keys[0] = '-aug' ; default=''
        Nbas_keys[1] = 'T','Q','D' ; default='T'
    
    To Do:
        make options for memory specifications
        make further options for optimizations
    '''
    
    # hard coded set of functionals built in for generator
    functionals = ['svwn','bp86','blyp','pw91','pbe','tpss',
                   'm06l','pbe0','b3lyp','bhlyp','b3p86','b97-1',
                   'mpw1k','x3lyp','tpssh','m06', 'm06.2x','b2plyp']
    
    workdir = top_direc
    cmp = input_set
    
    mult_data = pd.read_hdf('input_gendb.h5','multiplicity')
    if cmp not in mult_data.index:
        print('No multiplicity value found for {}'.format(cmp))
        M = str(input('Enter in a multiplicty to use: '))
    else:
        M = str(int(mult_data[cmp]))
    if optimize == None:
        optimize = 'n'
        
    sub_list = re.findall(r'[A-Z][a-z]*|\d+|\(|\)', cmp)
    if len(sub_list) == 1:
        geom_card = 'Pass'
        elem_list = sub_list
        
    else:        
        geom_card = 'Read'
        elem_list = []
        for p in range(1,len(sub_list)):
            if sub_list[p-1].isalpha() and sub_list[p].isalpha():
                if p == len(sub_list)-1:
                    elem_list.append(str(sub_list[p-1]))
                    elem_list.append(str(sub_list[p]))
                else:
                    elem_list.append(str(sub_list[p-1]))               
            elif sub_list[p-1].isalpha() and sub_list[p].isdigit():
                elem_list.append(str(sub_list[p-1]))
        if sub_list[-2].isdigit() and sub_list[-1].isalpha():
            elem_list.append(str(sub_list[-1]))
            
            
    basis_ae_exten = '{}-cc-pVTZ-DK'.format(Nbas_keys)

    basis_ext = basis_ae_exten                           
    
    BASIS = []
    if geom_card == 'Pass':
        atom = elem_list[0]
        GEO = 'zmatrix \n{}\nend\nend'.format(atom)
        BASIS.append(search_basis(atom,basis_ext))
        BASIS.append('\n') 
            
    elif geom_card == 'Read':
        geom_file = pd.read_hdf('input_gendb.h5','geom')
        if cmp not in geom_file.index:
            raise Exception('No geometry entry found for {}')

        GEO = geom_file[cmp]
        for atom in elem_list:
            BASIS.append(search_basis(atom,basis_ext))
            BASIS.append('\n')

    vec_input = ['atomic']
    if optimize == 'n':
        card = 'energy'
        func_span = range(0, len(functionals))
        if dft_direc == 'both':
            DFT_DIRECS = []
            DFT_DIRECS.append('task dft energy')
            DFT_DIRECS.append('task sodft energy')
            exten = []
            exten.append('ae')
            exten.append('so')
            vec_input.append('')
        else:
            DFT_DIRECS = ['task {} energy'.format(dft_direc)]
            exten = ['ae']
            
        new_dir = '{}\\{}'.format(workdir,cmp.lower())
        if os.path.isdir(new_dir) is False:
            os.mkdir(new_dir)
        os.chdir(new_dir) 
        
    elif optimize == 'y':
        card = 'opt+freq'
        if geom_card == 'Pass':
            raise Exception('Atoms do not require optimization')
        if dft_direc == 'both':
            raise Exception('Optimizations are currently limited to one directive block')
        else:
            pass
# default is tpss, no option to change unless you alter the code directly 
        func_span = [5]
        DFT_DIRECS = ['task {} optimize \ntask {} freq\n'.format(dft_direc,dft_direc)]
        exten = ['m' + str(M)]

#------------------------- Main Generator Routine ------------------------#
    for i in func_span:
        Main_file = []
        file_out = open('{}.{}.{}.nw'.format(cmp,functionals[i],exten[0]), 'w',newline='\n')
        Main_file.append('start {}.{}.{}\n'.format(cmp,functionals[i],exten[0]))
        Main_file.append('title "{} {} ae DK-basis {}"\n'.format(cmp,functionals[i],card)) 
        Main_file.append('memory total 6000 stack 2000 mb\necho\n\n')
        
        Main_file.append('geometry\n')
        Main_file.append(GEO)

        Main_file.append('\n\nBASIS spherical nosegment\n')
        for j in range(0, len(BASIS)):
            Main_file.append(BASIS[j])
        Main_file.append('end\n')    
        
        Main_file.append('\nrelativistic\nDouglas-Kroll DK3Full\nend')
        for j in range(0,len(DFT_DIRECS)):
            
            Main_file.append('\n\ndft\n')
            Main_file.append('mult {}\n'.format(M))
            Main_file.append('odft\n{}\n'.format(xckey[functionals[i]]))
            if vec_input[j] == 'atomic':
                Main_file.append('vectors input {} output {}.{}.{}.movecs\n'.format(vec_input[j],
                                                                                    cmp,functionals[i],exten[j])) 
            else:
                inp_arg = '{}.{}.{}'.format(cmp,functionals[i],exten[j-1])
                out_arg = '{}.{}.{}'.format(cmp,functionals[i],exten[j])
                Main_file.append('vectors input {}.movecs output {}.movecs\n'.format(inp_arg,out_arg))
                                 
            Main_file.append('maxiter 200\nend\n\n{}'.format(DFT_DIRECS[j])) 
            
        for k in range(0, len(Main_file)):
            file_out.write(str(Main_file[k]))
            
        file_out.close()
    if os.getcwd() != top_direc:       
        os.chdir(workdir)



parser = argparse.ArgumentParser(description='''
Hi! I make the .nw input files for dft calculations
on the heavy elements. This class is particular for AE calculations
currently set up for DKH, but will add the Zora option later
''',
epilog='\n Usage: inpgen-nw_ae.py -c WCl6 -d sodft -opt y \n')

parser.add_argument('-c','--compound', help="Name of compound", required=True)
parser.add_argument('-opt', '--optimize',help="do an optimization+freq", required=False)
parser.add_argument('-d', '--dft_task',help="task directive for: dft sodft both, default=both", required=False)
parser.add_argument('-aug', '--augmented',help="augmented ligand set?", required=False)
parser.add_argument('-sync', '--db_sync',help="sync database files -sync y", required=False)
args = vars(parser.parse_args())

compnd = args['compound']

info_l = []
if args['optimize']:
    optimize = args['optimize']
else:
    optimize = 'n'
    
if args['dft_task']:
    dft_direc = args['dft_task']
else:
    dft_direc = 'both'    
info_l = []
if args['augmented']:
    Nbas_keys = '-aug'
else:
    Nbas_keys = ''
if args['db_sync']:
	datasync()
else:
    pass 

gen_nw_DKH_dft_input(compnd, dft_direc, optimize, Nbas_keys)





