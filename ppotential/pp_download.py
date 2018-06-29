#!/usr/bin/env python3
import requests



base_path = 'http://www.tcm.phy.cam.ac.uk/~mdt26/pseudo_lib/'

def DiracFock_AREP():
    periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
#    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
#    periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
#    periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
#    periodic += ('Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')
#    for atom in periodic:
#        source_path = base_path + atom.lower() + '/1/pp.data'
#        target_path = 'DiracFock_AREP/' + atom.lower() + '_pp.data'
#        response = requests.get(source_path, stream=True)
#        handle = open(target_path, "wb")
#        for chunk in response.iter_content(chunk_size=512):
#            if chunk:  # filter out keep-alive new chunks
#                handle.write(chunk)

    for atom in periodic:
        source_path = base_path + atom.lower() + '/2/pp_gaussian'
        response = requests.get(source_path)
        flag = 0
        for line in response.text.split('\n'):
            if line.startswith('Element'):
                flag = 1
                print('newecp {}'.format(line.split()[2]))
                print(' N_core 2')
                print(' lmax d')
                continue
            if flag and flag < 3:
                flag += 1
                continue
            if flag == 3:
                print(line)


def HartreeFock():
    periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    for atom in periodic:
        source_path = base_path + atom.lower() + '/4/pp.data'
        target_path = 'HartreeFock/' + atom.lower() + '_pp.data'
        response = requests.get(source_path, stream=True)
        handle = open(target_path, "wb")
        for chunk in response.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)

def SofterDF_AREP():
    periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',       'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
    periodic += ('Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')
    for atom in periodic:
        source_path = base_path + atom.lower() + '/1_2/pp.data'
        target_path = 'SofterDF_AREP/' + atom.lower() + '_pp.data'
        response = requests.get(source_path, stream=True)
        handle = open(target_path, "wb")
        for chunk in response.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)

def SmallCoreDF_AREP():
    periodic += ('Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn')
    periodic += ('Mo', )
    periodic += ('Hf', )
    for atom in periodic:
        source_path = base_path + atom.lower() + '/7/pp.data'
        target_path = 'SmallCoreDF_AREP/' + atom.lower() + '_pp.data'
        response = requests.get(source_path, stream=True)
        handle = open(target_path, "wb")
        for chunk in response.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)

DiracFock_AREP()
