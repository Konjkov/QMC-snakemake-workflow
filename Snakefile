import os
from math import sqrt
from operator import itemgetter
from datetime import timedelta

# 10.1103/PhysRevA.44.7071 (TABLE XI)
exact_atomic_energy = {'h': -0.5, 'be': -14.66736, 'b': -24.65391, 'c': -37.8450, 'n': -54.5892, 'o': -75.0673, 'f': -99.7339,
                       'al': -242.346, 'si': -289.359, 'p':  -341.259, 's': -398.110, 'cl': -460.148}

ATOMS = {'h': 0.0, 'be': 0.0, 'b': 0.0, 'c': 0.0, 'n': 0.0, 'o': 0.0, 'f': 0.0, 'al': 0.0, 'si': 0.0, 'p': 0.0, 's': 0.0, 'cl': 0.0}

# J. Phys. Chem. A, 2008, 112 (50), pp 12868–12886 DOI: 10.1021/jp801805p (TABLE 2 / TAEe^a)
# Zero-point exclusive, nonrelativistic, clamped-nuclei total atomization energies. (kcal/mol)
W4_08_W4 = {'b2h6': 607.02, 'bhf2': 410.97, 'bf3': 470.97, 'c2h6': 713.08, 'h2cn': 343.75, 'nccn': 502.04, 'ch2nh2': 482.28, 'ch3nh': 474.63,
            'ch3nh2': 582.30, 'cf2': 258.78, 'n2h': 224.86, 't-n2h2': 296.53, 'n2h4': 438.28, 'fo2': 134.72, 'foof': 152.37, 'alf3': 430.97,
            'si2h6': 535.89, 'p4': 290.58, 'so2': 260.62, 'so3': 346.94, 'ocs': 335.75, 'cs2': 280.78, 's2o': 208.78, 's3': 168.36,
            's4-c2v': 234.35, 'becl2': 225.27, 'ccl2': 177.36, 'alcl3': 312.65, 'clcn': 285.45, 'oclo': 128.12, 'cloo': 126.39, 'cl2o': 101.46}
W4_08_W44 = {'h2': 109.49, 'oh': 107.21, 'hf': 141.64, 'h2o': 232.97, 'ch': 84.22,  'ch2-trip': 190.74, 'ch3': 307.87, 'ch4': 420.42,
             'cch': 266.16, 'c2h2': 405.52, 'nh3': 298.02, 'c2': 147.02, 'n2': 228.48, 'co': 259.73, 'cn': 181.35, 'no': 152.75,
             'o2': 120.82, 'of': 53.08, 'f2': 39.04, 'ph3': 242.27, 'hs': 87.73, 'h2s': 183.91, 'hcl': 107.50, 'so': 126.47,
             'clo': 65.45, 'clf': 62.80, 'p2': 117.59, 's2': 104.25, 'cl2': 59.75}
W4_08_W43 = {'be2': 2.67, 'b2': 67.46, 'bh': 84.99, 'bh3': 281.29, 'bn': 105.24, 'bf': 182.52, 'nh': 83.10, 'nh2': 182.59, 'hcn': 313.42,
             'hof': 158.65, 'alh': 73.57, 'alh3': 213.17, 'alf': 163.78, 'alcl': 122.62, 'sih': 73.92, 'sih4': 324.95, 'sio': 193.05,
             'sif': 142.71, 'cs': 172.22}
W4_08_W42 = {'bn3pi': 105.82, 'cf': 132.72, 'bef2': 309.10, 'ch2c': 359.93, 'ch2ch': 446.08, 'c2h4': 564.10, 'ch2nh': 439.44, 'hco': 279.42,
             'h2co': 374.66, 'co2': 390.14, 'hno': 205.89, 'no2': 227.88, 'n2o': 270.85, 'o3': 147.43, 'hoo': 175.53, 'hooh': 269.09,
             'f2o': 93.78, 'hocl': 166.23, 'ssh': 165.13}

W4_08 = {}  #  99 total
W4_08.update(W4_08_W4)
W4_08.update(W4_08_W44)
W4_08.update(W4_08_W43)
W4_08.update(W4_08_W42)

# J. Phys. Chem. A, 2009, 113 (29), pp 8434–8447   DOI: 10.1021/jp904369h (TABLE 1 / TAEe - relativ. - spin orbit - DBOC) (kcal/mol)
ALKANES = {'propane': 1007.15 + 0.58 + 0.25 - 0.20 , 'propene': 860.89 + 0.52 + 0.25 - 0.18,
           'propyne': 704.97 + 0.48 + 0.25 - 0.17, 'allene': 703.45 + 0.47 + 0.25 - 0.15}

# Chemical Physics Letters 510 (2011) 165–178      DOI: 10.1016/j.cplett.2011.05.007 (TABLE 4 / TAEe - relativ. - spin orbit - DBOC) (kcal/mol)
W4_11 = {'oxirene': 455.33 + 0.46 + 0.39 - 0.11, 'oxirane': 650.70 + 0.56 + 0.39 - 0.13, 'dioxirane': 409.17 + 0.40 + 0.53 - 0.07,
         'ketene': 532.70 + 0.47 + 0.39 - 0.11, 'acetaldehyde': 677.07 + 0.53 + 0.39 - 0.12,
         'formic': 500.90 + 0.59 + 0.53 - 0.13, 'acetic': 802.82 + 0.79 + 0.62 - 0.21,
         'methanol': 512.86 + 0.46 + 0.31 - 0.13, 'ethanol': 810.39 + 0.65 + 0.39 - 0.19, 'glyoxal': 633.91 + 0.65 + 0.62 - 0.07,
         't-hcoh': 321.87 + 0.33 + 0.31 - 0.04, 'c-hcoh': 317.04 + 0.32 + 0.31 - 0.02,
         'hcnh': 335.90 + 0.32 + 0.08 - 0.05, 'hocn': 409.36 + 0.52 + 0.31 - 0.12,
         'honc': 349.45 + 0.50 + 0.31 - 0.11, 'hnco': 434.00 + 0.54 + 0.31 - 0.11, 'hcno': 364.20 + 0.57 + 0.31 - 0.10,
         'c-n2h2': 290.85 + 0.31 + 0.00 - 0.03, 'hnnn': 331.35 + 0.50 + 0.00 - 0.07, 'hnc': 297.95 + 0.26 + 0.08 - 0.09,
         't-hono': 311.87 + 0.42 + 0.45 - 0.08, 'c-hono': 311.44 + 0.42 + 0.45 - 0.08,
         'nh2cl': 246.92 + 0.39 + 0.84 - 0.09, 't-hooo': 232.38 + 0.31 + 0.67 - 0.06, 'c-hooo': 232.18 + 0.31 + 0.67 - 0.06,
         'ch3f': 422.19 + 0.38 + 0.47 - 0.08, 'ch2f2': 436.34 + 0.54 + 0.85 - 0.07,
         'cf4': 476.36 + 0.85 + 1.63 - 0.07, 'sih3f': 381.00 + 0.95 + 0.81 - 0.01, 'sif4': 573.96 + 1.90 + 1.97 - 0.05,
         'c2h3f': 572.95 + 0.51 + 0.55 - 0.12, 'c2h5f': 720.53 + 0.56 + 0.55 - 0.14, 'fccf': 384.50 + 0.72 + 0.94 - 0.08,
         'hccf': 397.53 + 0.49 + 0.55 - 0.10, 'hcof': 402.62 + 0.49 + 0.69 - 0.06, 'f2co': 418.95 + 0.67 + 1.08 - 0.06,
         'ch2-sing': 181.18 + 0.09 + 0.08 + 0.10}


def atom_charge(symbol):
    periodic = ('X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
    atoms = {v:i for i,v in enumerate(periodic)}
    return atoms[symbol]

def get_max_Z(gwfn_file):
    """Get maximal Z for atoms in molecule."""
    with open(gwfn_file, "r") as gwfn:
        result = ''
        line = gwfn.readline()
        while line and not line.startswith('Atomic numbers for each atom:'):
            line = gwfn.readline()
        if not line:
            raise Exception
        line = gwfn.readline()
        while line and not line.startswith('Valence charges for each atom:'):
            result += line
            line = gwfn.readline()
    return max(map(int, result.split()))

def get_atom_list(molecule):
    with open(os.path.join('..', 'chem_database', molecule+'.in'), 'r') as input_geometry:
        result = dict.fromkeys(ATOMS, 0)
        for line in input_geometry:
            if line.startswith(' '):
                atom_symbol = line.split()[0].lower()
                result[atom_symbol] += 1
        return result

def get_ae_cutoffs(molecule):
    with open(os.path.join('..', 'chem_database', molecule+'.in'), 'r') as input_geometry:
        i = 0
        result = []
        for line in input_geometry:
            if line.startswith(' '):
                i += 1
                result.append('{i}         {i}         0.7                          1'.format(i=i))
        return '\n  '.join(result)

def get_atom_labels(molecule):
    """Returns number of atoms in a set and list of labels for this set.
    Used for generic JASTROW format.
    """
    with open(os.path.join('..', 'chem_database', molecule+'.in'), 'r') as input_geometry:
        i = 0
        result = []
        for line in input_geometry:
            if line.startswith(' '):
                i += 1
                result.append('{i}'.format(i=i))
        return i, ' '.join(result)

def vmc_energy(molecule, method, basis):
    """Get VMC energy without JASTROW optimisation.
     -152.988424660763 +/- 0.003047553900      Correlation time method
    """

    regexp = re.compile(' (?P<energy>[-+]?\d+\.\d+) \+/- (?P<energy_error>[-+]?\d+\.\d+)      Correlation time method')
    with open(os.path.join(molecule, method, basis, 'VMC', '10000000', 'out'), 'r') as vmc_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, vmc_out.read())[-1])
    return value, error

def vmc_opt_energy(molecule, method, basis, opt_method, rank):
    """Get VMC energy with JASTROW optimisation.
     -153.693436512511 +/- 0.003006326588      Correlation time method
    """

    regexp = re.compile(' (?P<energy>[-+]?\d+\.\d+) \+/- (?P<energy_error>[-+]?\d+\.\d+)      Correlation time method')
    with open(os.path.join(molecule, method, basis, opt_method, 'emin', 'casl', rank, '1000000_9', 'out'), 'r') as vmc_opt_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, vmc_opt_out.read())[-1])
    return value, error

def vmc_opt_variance(molecule, method, basis):
    """Get VMC variance with JASTROW optimisation.
      Sample variance of E_L (au^2/sim.cell) : 3.169677109628 +- 0.034986257092
    """

    regexp = re.compile('Sample variance of E_L \(au\^2/sim.cell\) : (?P<variance>[-+]?\d+\.\d+) \+- (?P<variance_error>[-+]?\d+\.\d+)')
    with open(os.path.join(molecule, method, basis, 'VMC_OPT', 'emin', 'casl', '8_8_44', '1000000_9', 'out'), 'r') as vmc_opt_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, vmc_opt_out.read())[-1])
    return value, error

def dmc_energy(molecule, method, basis):
    """Get DMC energy.
          mean:   -153.795024411601 +/-       0.001346260888
    """

    dir = os.path.join(molecule, method, basis, 'VMC_DMC', 'emin', 'casl', '8_8_44', 'tmax_2_1024_1')
    open(os.path.join(dir, '.casino_finished'), 'r').close()
    regexp = re.compile('mean:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
    with open(os.path.join(dir, 'out'), 'r') as dmc_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, dmc_out.read())[-1])
    return value, error

def dmc_stderr(molecule, method, basis):
    """Get DMC standard error.
          stderr:      0.000906128433 +/-       0.000046917552
    """

    dir = os.path.join(molecule, method, basis, 'VMC_DMC', 'emin', 'casl', '8_8_44', 'tmax_2_1024_1')
    open(os.path.join(dir, '.casino_finished'), 'r').close()
    regexp = re.compile('stderr:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
    with open(os.path.join(dir, 'out'), 'r') as dmc_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, dmc_out.read())[-1])
    return value, error

def dmc_ncorr(molecule, method, basis):
    """Get DMC correlation N.
          N_corr:      0.000906128433 +/-       0.000046917552
    """

    dir = os.path.join(molecule, method, basis, 'VMC_DMC', 'emin', 'casl', '8_8_44', 'tmax_2_1024_1')
    open(os.path.join(dir, '.casino_finished'), 'r').close()
    regexp = re.compile('N_corr:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
    with open(os.path.join(dir, 'out'), 'r') as dmc_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, dmc_out.read())[-1])
    return value, error

def dmc_stats_nstep(molecule, method, basis):
    """Get DMC statistic accumulation steps.
          dmc_stats_nstep   : 96000
    """

    dir = os.path.join(molecule, method, basis, 'VMC_DMC', 'emin', 'casl', '8_8_44', 'tmax_2_1024_1')
    regexp = re.compile('dmc_stats_nstep   :\s+(?P<nstep>\d+)')
    with open(os.path.join(dir, 'input'), 'r') as dmc_input:
        # we are only interested in the last occurrence
        value = int(re.findall(regexp, dmc_input.read())[-1])
    return value

def TAE_energy(molecule, method, basis):
    """Total atomization energy (kcal/mol)"""

    atom_list = get_atom_list(molecule)
    energy, energy_error = dmc_energy(molecule, method, basis)

    tae_energy = 630.0 * (energy - sum([atom_list[atom]*dmc_energy(atom, method, basis)[0] for atom in atom_list])) + MOLECULES[molecule]

    if molecule in ATOMS:
        tae_energy_error = 630.0 * energy_error
    else:
        tae_energy_error = 630.0 * sqrt(energy_error**2 + sum([atom_list[atom]*dmc_energy(atom, method, basis)[1]**2 for atom in atom_list]))

    return tae_energy, tae_energy_error

wildcard_constraints:
    molecule='[-\w]+',
    method='[-\w]+',
    basis='[-\w]+',
    jastrow_type='[_\w]+',
    jastrow_rank='[_\w]+',
    backflow_rank='[_\w]+',
    jastrow_opt_method='\w+'

####################################################################################################################

rule VMC_DMC_RUN:
    input:      '{path}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/input',
                '{path}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/gwfn.data',
                '{path}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/parameters.casl',
    output:     '{path}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/out'
    shell:      'cd {wildcards.path}/VMC_DMC/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}/tmax_2_1024_1 && runqmc'


rule VMC_DMC_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/.keep',
                '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/out',
                '{molecule}/{method}/{basis}/VMC/10000000/out'
    output:     '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/input'
    params:
        dt_relative_step = 2.0,
        stderr = 0.001,
        # roughly
        # nstep ~ 1.5 * (E(VMC) - E(DMC))/(nconfig*dtdmc*stderr*stderr)
        # (E(VMC) - E(DMC)) ~ 6 * (E(HF) - E(VMC))
        magic_const = 1.5/6.0
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            gwfn_file = os.path.join(wildcards.molecule, wildcards.method, wildcards.basis, 'gwfn.data')
            hf, _ = vmc_energy(wildcards.molecule, wildcards.method, wildcards.basis)
            vmc, _ = vmc_opt_energy(wildcards.molecule, wildcards.method, wildcards.basis, 'VMC_OPT', wildcards.jastrow_rank)
            dtdmc = 1.0/(get_max_Z(gwfn_file)**2 * 3.0 * params.dt_relative_step)
            nstep = params.magic_const*(hf - vmc)/(int(wildcards.nconfig)*dtdmc*params.stderr*params.stderr)
            nstep=max(50000, int(round(nstep, -3)))
            with open(file_name, 'w') as f:
                f.write(open('../vmc_dmc.tmpl').read().format(
                    neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep // 1000
                ))

rule VMC_DMC_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/out',
                '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/parameters.casl'
    shell:      'ln -s ../../../../../VMC_OPT/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}/10000/parameters.9.casl {output}'

rule VMC_DMC_GWFN:
    input:      '{molecule}/{method}/{basis}/gwfn.data',
                '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_DMC_DIRS:
    input:      '{molecule}/{method}/{basis}/gwfn.data'
    output:     '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_OPT_ENERGY_RUN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/input',
                '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/gwfn.data',
                '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/parameters.casl'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/out'
    shell:      'cd {wildcards.path}/VMC_OPT/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}/1000000_9 && runqmc'

rule VMC_OPT_ENERGY_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc_opt_energy.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule))

rule VMC_OPT_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/out'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/parameters.casl'
    shell:      'ln -s ../10000/parameters.9.casl {output}'

rule VMC_OPT_ENERGY_GWFN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/.keep'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_OPT_ENERGY_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_OPT_RUN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/input',
                '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/gwfn.data',
                '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/parameters.casl',
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/out'
    shell:      'cd {wildcards.path}/VMC_OPT/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}/10000 && runqmc'

rule VMC_OPT_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc_opt.tmpl').read().format(neu=neu, ned=ned, nconfig=10000, method=wildcards.jastrow_opt_method, cycles=9, molecule=wildcards.molecule))

rule VMC_OPT_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/parameters.casl'
    run:
        for file_name in output:
            jastrow = wildcards.jastrow_rank.split('_')
            with open(file_name, 'w') as f:
                f.write(open('../casl.tmpl').read().format(term_2_0=jastrow[0], term_1_1=jastrow[1], term_2_1_1=jastrow[2][0], term_2_1_2=jastrow[2][1]))

rule VMC_OPT_GWFN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/.keep'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_OPT_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_DMC_VARIANCE_PLOT:
    output:     'vmc_dmc_variance.dat'
    params:
        method = 'HF',
        basis = 'cc-pVQZ'
    run:
        with open(output[0], 'w') as output_file:
            print('# molecule      VMC_variance VMC_var_error  DMC_variance  DMC_var_error', file=output_file)
            for molecule in MOLECULES:
                try:
                    vmc_variance, vmc_variance_error = vmc_opt_variance(molecule, params.method, params.basis)
                    dmc_energy_stderr, dmc_energy_stderr_error = dmc_stderr(molecule, params.method, params.basis)
                    n_corr, n_corr_error = dmc_ncorr(molecule, params.method, params.basis)
                    n_step = dmc_stats_nstep(molecule, params.method, params.basis)
                except FileNotFoundError:
                    continue
                print('{:12} {:>13.6f} {:>13.6f} {:>13.6f} {:>13.6f}'.format(
                    molecule,
                    vmc_variance,
                    vmc_variance_error,
                    dmc_energy_stderr**2 * 1024 * n_step/n_corr,
                    (2 * dmc_energy_stderr_error * dmc_energy_stderr/n_corr + dmc_energy_stderr**2 * n_corr_error/n_corr**2) * 1024 * n_step), file=output_file)

rule VMC_DMC_PLOT:
    output:     'dmc_energy.dat'
    params:
        method = 'HF',
        basis = 'cc-pVQZ'
    run:
        with open(output[0], 'w') as output_file:
            print('# molecule\\atoms  H   Be  B   C   N   O   F   Al  Si  P   S   Cl  E(DMC)+TAE(au)  DMC_error(au)  TAE-TAE(DMC)(kcal/mol) TAE(DMC)_error(kcal/mol)', file=output_file)
            dmc = []
            for molecule in MOLECULES:
                atom_list = get_atom_list(molecule)
                try:
                    energy, energy_error = dmc_energy(molecule, params.method, params.basis)
                    tae_energy, tae_energy_error = TAE_energy(molecule, params.method, params.basis)
                except FileNotFoundError:
                    continue
                result = {
                    'molecule': molecule,
                    'energy': energy + MOLECULES[molecule]/630.0,
                    'energy_error': energy_error,
                    'tae_energy': tae_energy,
                    'tae_energy_error': tae_energy_error
                }
                result.update(atom_list)
                dmc.append(result)
            dmc = sorted(dmc, key=itemgetter('tae_energy'))
            for item in dmc:
                print(
                    '{molecule:12}    {h:3} {be:3} {b:3} {c:3} {n:3} {o:3} {f:3} {al:3} {si:3} {p:3} {s:3} {cl:3}  '
                    '{energy:>13.6f} {energy_error:>13.6f}      {tae_energy:>13.6f}          {tae_energy_error:>13.6f}'.format(**item), file=output_file
                )

rule VMC_PLOT:
    output:     'hf_vmc_energy.dat'
    params:
        method = 'HF',
        basis = 'cc-pVQZ'
    run:
        with open(output[0], 'w') as output_file:
            for molecule in MOLECULES:
                print('{:12} {:>13.6f} {:>13.6f} {:>13.6f}'.format(
                    molecule,
                    hf_energy(molecule, params.method, params.basis),
                    vmc_energy(molecule, params.method, params.basis)[0],
                    vmc_energy(molecule, params.method, params.basis)[1]),
                    file=output_file
               )

rule VMC_RUN:
    input:      '{path}/VMC/10000000/input'
    output:     '{path}/VMC/10000000/out'
    shell:      'cd {wildcards.path}/VMC/10000000 && runqmc'

rule VMC_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC/10000000/gwfn.data'
    output:     '{molecule}/{method}/{basis}/VMC/10000000/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule))

rule VMC_GWFN:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC/10000000/gwfn.data'
    shell:      'ln -s ../../gwfn.data {output}'

rule VMC_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC/10000000/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_DMC_BF_RUN:
    input:      '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/input',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/correlation.data',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/parameters.casl',
    output:     '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/out'
    shell:      'cd {wildcards.path}/VMC_DMC_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/tmax_2_1024_1 && runqmc'

rule VMC_DMC_BF_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep',
                '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/out',
                '{molecule}/{method}/{basis}/VMC/10000000/out'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/input'
    params:
        dt_relative_step = 2.0,
        stderr = 0.001,
        # roughly
        # nstep ~ 1.5 * (E(VMC) - E(DMC))/(nconfig*dtdmc*stderr*stderr)
        # (E(VMC) - E(DMC)) ~ 6 * (E(HF) - E(VMC))
        magic_const = 1.5/6.0
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            gwfn_file = os.path.join(wildcards.molecule, wildcards.method, wildcards.basis, 'gwfn.data')
            hf, _ = vmc_energy(wildcards.molecule, wildcards.method, wildcards.basis)
            vmc, _ = vmc_opt_energy(wildcards.molecule, wildcards.method, wildcards.basis, 'VMC_OPT_BF', wildcards.jastrow_rank + '__' + wildcards.backflow_rank)
            dtdmc = 1.0/(get_max_Z(gwfn_file)**2 * 3.0 * params.dt_relative_step)
            nstep = params.magic_const*(hf - vmc)/(int(wildcards.nconfig)*dtdmc*params.stderr*params.stderr)
            nstep=max(50000, int(round(nstep, -3)))
            with open(file_name, 'w') as f:
                f.write(open('../vmc_dmc_bf.tmpl').read().format(
                    neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep // 1000
                ))

rule VMC_DMC_BF_DATA_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out',
                '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/correlation.data'
    shell:      'ln -s ../../../../../VMC_OPT_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/10000/correlation.out.9 {output}'

rule VMC_DMC_BF_CASL_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out',
                '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/parameters.casl'
    shell:      'ln -s ../../../../../VMC_OPT_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/10000/parameters.9.casl {output}'

rule VMC_DMC_BF_GWFN:
    input:      '{molecule}/{method}/{basis}/gwfn.data',
                '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_DMC_BF_DIRS:
    input:      '{molecule}/{method}/{basis}/gwfn.data'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep'
    shell:      'touch {output}'

####################################################################################################################


rule VMC_OPT_BF_ENERGY_RUN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/input',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/gwfn.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/correlation.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/parameters.casl'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/out'
    shell:      'cd {wildcards.path}/VMC_OPT_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/1000000_9 && runqmc'

rule VMC_OPT_BF_ENERGY_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc_opt_energy_bf.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule))

rule VMC_OPT_BF_DATA_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/correlation.data'
    shell:      'ln -s ../10000/correlation.out.9 {output}'

rule VMC_OPT_BF_CASL_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/parameters.casl'
    shell:      'ln -s ../10000/parameters.9.casl {output}'

rule VMC_OPT_BF_ENERGY_GWFN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/.keep'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_OPT_BF_ENERGY_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_OPT_BF_CASL_RUN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/input',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/gwfn.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/correlation.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/parameters.casl',
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out'
    shell:      'cd {wildcards.path}/VMC_OPT_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/10000 && runqmc'

rule VMC_OPT_BF_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc_opt_bf.tmpl').read().format(
                    neu=neu, ned=ned, nconfig=10000,
                    method=wildcards.jastrow_opt_method, cycles=9, molecule=wildcards.molecule))

rule VMC_OPT_BF_DATA_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/correlation.data'
    run:
        for file_name in output:
            jastrow = wildcards.jastrow_rank.split('_')
            with open(file_name, 'w') as f:
                backflow = wildcards.backflow_rank.split('_')
                ae_cutoffs = get_ae_cutoffs(wildcards.molecule)
                number, labels = get_atom_labels(wildcards.molecule)
                f.write(open('../backflow.tmpl').read().format(
                    eta_term=backflow[0],
                    mu_number_of_atoms=number, mu_atom_labels=labels, mu_term=backflow[1],
                    phi_number_of_atoms=number, phi_atom_labels=labels, phi_term_eN=backflow[2][0], phi_term_ee=backflow[2][1],
                    ae_cutoffs=ae_cutoffs))


rule VMC_OPT_BF_CASL_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/parameters.casl'
    run:
        for file_name in output:
            jastrow = wildcards.jastrow_rank.split('_')
            with open(file_name, 'w') as f:
                f.write(open('../casl.tmpl').read().format(term_2_0=jastrow[0], term_1_1=jastrow[1], term_2_1_1=jastrow[2][0], term_2_1_2=jastrow[2][1]))


rule VMC_OPT_BF_GWFN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_OPT_BF_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    shell:      'touch {output}'
