import os, re, string
import sys, argparse
from numpy import *
import pandas as pd
top_direc = os.getcwd()

def data_sync():
	# commentted out for now until I reimplement the alternative files
#     os.chdir('basis')
#     primitive_basis_db = []
#     indexer = []
#     for i in os.listdir():
#         basis_r = open(f"{i}",'r')
#         indexer.append(re.sub(".txt","",i))
#         entry = basis_r.read()
#         primitive_basis_db.append(entry)
#         basis_r.close()
#     proto_basis_database = pd.Series(primitive_basis_db,index=indexer)
#     os.chdir('..')
    os.chdir('geom')
    primitive_geom_db = []
    indexer = []
    for i in os.listdir():
        geom_r = open(f"{i}",'r')
        indexer.append(re.sub(".geom.txt","",i))
        entry = geom_r.read()
        primitive_geom_db.append(entry)
        geom_r.close()
    proto_database = pd.Series(primitive_geom_db,index=indexer)
    os.chdir('..')
    os.chdir('Bash_scripts')
    mult_info = pd.read_excel('multiplicty.xlsx')
    mult = mult_info.loc[:,'Mult']
    os.chdir('..')

    maindb = pd.HDFStore('input_gendb.h5','a')
#     maindb['basis'] =  proto_basis_database
    maindb['geom']  =  proto_database
    maindb['multiplicity'] = mult
    maindb.close()

def init_params():
    # These elements require a pseudopotential - not inclusive list, just includes the atoms we work with.    
    rel_set = array(['Br','I','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',
                     'Th','U','Pa','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
                     'Md','No','Lr','La','Ce','Pr','Nd','Sm','Eu','Gd','Tb',
                     'Dy','Ho','Er','Tm','Yb','Lu'])

    # functional prefix master list, includes general needed functionals 
    functionals_master = ['svwn',
                          'bp86',
                          'blyp',
                          'pw91',
                          'pbe', 
                          'tpss',
                          'm06l',
                          'pbe0',
                          'b3lyp',
                          'bhlyp',
                          'b3p86',
                          'b97-1',
                          'mpw1k',
                          'x3lyp',
                          'tpssh',
                          'm06',
                          'm06.2x',
                          'm11',
                          'camb3lyp',
                          'b2plyp'
                         ]

    # functional prefix key, order must match above or your file names will be off
    xckey_ms = ['xc slater vwn_5', 
                'xc becke88 perdew86', 
                'xc becke88 lyp', 
                'xc xperdew91 perdew91', 
                'xc xpbe96 cpbe96', 
                'xc xtpss03 ctpss03', 
                'xc m06-l', 
                'xc pbe0', 
                'xc b3lyp', 
                'xc bhlyp',
                'xc vwn_1_rpa perdew86 nonlocal 0.81 HFexch 0.20 slater 0.80 becke88 nonlocal 0.72', 
                'xc becke97-1', 
                'xc mpw1k',  
                'xc vwn_1_rpa 0.129 lyp 0.871 hfexch 0.218 slater 0.782 becke88 nonlocal 0.542 xperdew91 nonlocal 0.167', 
                'xc xctpssh', 
                'xc m06', 
                'xc m06-2x',  
                'xc m11', 
                'xc xcamb88 1.00 lyp 0.81 vwn_5 0.19 hfexch 1.00 \n cam 0.33 cam_alpha 0.19 cam_beta 0.46 \n direct', 
                'xc HFexch 0.53 becke88 0.47 lyp 0.73 mp2 0.27 \n dftmp2 \n direct'
               ]

    #store functional keywords to a dictionary for reference
    temp_frame = pd.Series(xckey_ms,index=functionals_master)
    xckey = temp_frame.to_dict()  
    return xckey, rel_set


def search_basis(atom,basis_ext):
    basis_file = pd.read_hdf('input_gendb.h5','basis')
    basis_entry = '{}{}'.format(atom, basis_ext)
    if basis_entry not in basis_file.index:
        raise Exception('No basis of type {} for {}'.format(basis_ext,atom))
    return basis_file[basis_entry]

def gen_nw_recp_dft_input(input_set, 
                          rel_basis, 
                          dft_direc='dft',
                          optimize='n',
                          level='WB',
                          Nbas_keys=''):
    '''
    this program generates a set of DFT functional inputs that run in NWChem
    First off, before anyone asks, YES, THIS USES THE cc-pV(T+D)Z SET FOR ROW 3 AND BEYOND!
    
    Function arguments
    required:
        input_set: the molecule you are making
        rel_basis: 
                  'ANO': for the atomic natural orbitals Stuggart generated basis sets
                  'SEG': for the Segmented Stuggart generated basis sets
                  '97': for the 1997 RSC Stuggart generated basis set
                  'CC': use the cc-pVTZ-PP basis sets developed from Peterson for 
                        use with the Stuggart pseduopotentials
    optional arguments:
        dft_direc= 'dft', 'sodft', 'both', default='dft'
        optimize= 'n', 'y' ; default='n'
        SG= 'ANO', 'SEG', '97' ; default='ANO'
        level: 'WB', 'DF' ; default='WB'
        Nbas_keys = '-aug' ; default=''
    
    To Do:
        make options for memory specifications
        make further options for optimizations
    '''
    xckey, rel_set = init_params()
    
    # hard coded set of functionals built in for generator
#     functionals = ['svwn','bp86','blyp','pw91','pbe','tpss',
#                    'm06l','pbe0','b3lyp','bhlyp','b3p86','b97-1',
#                    'mpw1k','x3lyp','tpssh','m06', 'm06.2x','b2plyp']
    
    functionals = ['svwn','bp86','blyp','pw91','pbe','tpss',
                   'm06l','pbe0','b3lyp','bhlyp','b3p86','b97-1',
                   'x3lyp','tpssh']   
    
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
                             
    if rel_basis == 'CC':   
        basis_ecp_exten = '-cc-pVTZ-PP'
        SG='cc-pvtz'
    else:
        SG=rel_basis
        basis_ecp_exten = '-{}_ECPscM{}'.format(SG,level)
        
    non_rel_basis  = '{}-cc-pVTZ'.format(Nbas_keys)
    
    ecp_exten = '-ECPscM{}'.format(level) 
    soecp_exten = '-ECPscM{}-SO'.format(level) 
    
    BASIS = []
    ECP = []
    SO = []
    if geom_card == 'Pass':
        atom = elem_list[0]
        if atom in rel_set:
            if atom == 'Br' or atom == 'I':
                basis_ext = '-cc-pVTZ-PP'
            else: 
                basis_ext = basis_ecp_exten
            ECP_req = 'Yes'
        elif atom not in rel_set:
            basis_ext = non_rel_basis
            ECP_req = 'No'
        GEO = 'zmatrix\n{}\nend\nend'.format(atom)

        BASIS.append(search_basis(atom,basis_ext))
        BASIS.append('\n') 
        
        if ECP_req == 'Yes':
            ECP.append(search_basis(atom,ecp_exten))
            ECP.append('\n') 

            SO.append(search_basis(atom,soecp_exten))
            SO.append('\n')
            
    elif geom_card == 'Read':
        geom_file = pd.read_hdf('input_gendb.h5','geom')
        if cmp not in geom_file.index:
            raise Exception('No geometry entry found for {}'.format(cmp))
                             
        GEO = geom_file[cmp]
            
        for atom in elem_list:
            if atom in rel_set:
                ECP_req = 'Yes'
                if atom == 'Br' or atom == 'I':
                    basis_ext = '-cc-pVTZ-PP'
                else: 
                    basis_ext = basis_ecp_exten
                
                BASIS.append(search_basis(atom,basis_ext))
                BASIS.append('\n')

                ECP.append(search_basis(atom,ecp_exten))
                ECP.append('\n') 

                SO.append(search_basis(atom,soecp_exten))
                SO.append('\n')

            elif atom not in rel_set:
                BASIS.append(search_basis(atom,non_rel_basis))
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
            exten.append('rp')
            exten.append('so')
            vec_input.append('')
        else:
            DFT_DIRECS = ['task {} energy'.format(dft_direc)]
            exten = ['rp']
            
        new_dir = '{}\{}'.format(workdir,cmp.lower())
        new_dir = '{}'.format(cmp.lower())
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
        func_span = [5]
        DFT_DIRECS = ['task {} optimize \ntask {} freq\n'.format(dft_direc,dft_direc)]
        exten = ['m' + str(M)]

#------------------------- Main Generator Routine ------------------------#
    for i in func_span:
        Main_file = []
        file_out = open('{}.{}.{}.nw'.format(cmp,functionals[i],exten[0]), 'w', newline='\n')
        Main_file.append('start {}.{}.{}\n'.format(cmp, functionals[i], exten[0]))
        Main_file.append('title "{} {} {} basis {}"\n'.format(cmp, functionals[i], SG, card)) 
        Main_file.append('memory 500 mw\necho\n\n')
        
        Main_file.append('geometry\n')
        Main_file.append(GEO)

        Main_file.append('\n\nBASIS spherical\n')
        for j in range(0, len(BASIS)):
            Main_file.append(BASIS[j])
        Main_file.append('end\n')    
        
        if ECP_req == 'Yes':
            Main_file.append('ECP\n')
            for k in range(0, len(ECP)):
                Main_file.append(ECP[k])
            Main_file.append('end\n')

            Main_file.append('SO\n')
            for s in range(0, len(SO)):
                Main_file.append(SO[s])
            Main_file.append('end\n')

        for j in range(0,len(DFT_DIRECS)):
            
            Main_file.append('\n\ndft\n')
            Main_file.append('mult {}\n'.format(M))
            Main_file.append('odft\n{}\n'.format(xckey[functionals[i]]))
            if vec_input[j] == 'atomic':
                Main_file.append('vectors input {} output {}.{}.{}.movecs\n'.format(vec_input[j],
                                                                                    cmp,
                                                                                    functionals[i],
                                                                                    exten[j])) 
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
on the heavy elements (5d, Ln, and An). This can routines
for dft, sodft, and both. For Recp's there are four basis 
classes: ANO, SEG, correlation consistent, and the old 1997
RSC stuttgart start orbitals. 
''',
epilog='\n Usage: inpgen-nw.py -c UO3 -b ANO -e WB -d sodft -opt y \n')

parser.add_argument('-c','--compound', help="Name of compound", required=True)
parser.add_argument('-b', '--relbasis',help="basis type: CC, SEG, ANO, 97", required=True)
parser.add_argument('-e', '--ecptype',help="relativistic method used to fit ECP: WB (quasi), DF(full)", required=True)
parser.add_argument('-opt', '--optimize',help="y: do an optimization+freq, default=n", required=False)
parser.add_argument('-d', '--dft_task',help="task directive for: dft, sodft, or both, default=both", required=False)
parser.add_argument('-aug', '--augmented',help="augmented ligand basis set?", required=False)
parser.add_argument('-sync', '--db_sync',help="sync database files, default=n", required=False)
args = vars(parser.parse_args())

compnd = args['compound']
rel_basis = args['relbasis']
level = args['ecptype']

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
gen_nw_recp_dft_input(compnd, rel_basis, dft_direc, optimize, level, Nbas_keys)


#Mmove this to separate program 
def gen_nw_DKH_dft_input(input_set, dft_direc='dft',optimize='n', Nbas_keys=['','T']):
    '''
    this program generates a set of DFT functional inputs that run in NWChem
    
    this routine is designed for DK3 all-electron calculations, maybe I'll implement zora, if for some reason
    someone gets stuck doing a benchmark for those. 
    
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
            
            
    basis_ae_exten = '{}-cc-pV{}Z-DK'.format(Nbas_keys[0],Nbas_keys[1])

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


