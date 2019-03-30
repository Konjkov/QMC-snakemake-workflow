import os
import csv
from math import sqrt
from collections import defaultdict
from operator import itemgetter
from datetime import timedelta

configfile: "config.yaml"
INPUTS_DIR = '../chem_database'

def atom_charge(symbol):
    periodic = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
    periodic += ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
    periodic += ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']
    periodic += ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']
    periodic += ['Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']
    periodic += ['Fr', 'Ra', 'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']
    periodic[58:58] = ['Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    periodic[90:90] = ['Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
    atoms = {v.lower():i for i,v in enumerate(periodic)}
    return atoms[symbol.lower()]

def get_atomic_symbols(molecule):
    """Get atomic symbols set"""
    with open(os.path.join(INPUTS_DIR, molecule + '.xyz'), 'r') as xyz:
        natoms = int(xyz.readline())
        charge, mult = map(int, xyz.readline().split())
        atomic_symbols = set()
        for atom in range(natoms):
            symbol, x, y, z = xyz.readline().split()
            atomic_symbols.add(symbol)
    return atomic_symbols

def get_XYZ(molecule):
    """Load XYZ-geometry from file."""
    with open(os.path.join(INPUTS_DIR, molecule + '.xyz'), 'r') as xyz:
        natoms = int(xyz.readline())
        charge, mult = map(int, xyz.readline().split())
        geometry = []
        for atom in range(natoms):
            symbol, x, y, z = xyz.readline().split()
            geometry.append((atom_charge(symbol), map(float, (x, y, z))))
    return geometry

def get_XYZ_charge_mul(molecule):
    """Get charge and mutiplicity."""
    with open(os.path.join(INPUTS_DIR, molecule + '.xyz'), 'r') as xyz:
        natoms = int(xyz.readline())
        charge, mult = map(int, xyz.readline().split())
        geometry = []
        for atom in range(natoms):
            symbol, x, y, z = xyz.readline().split()
            geometry.append((atom_charge(symbol), map(float, (x, y, z))))
    return charge, mult

def get_XYZ_nel(molecule):
    """Get number of electrons."""
    with open(os.path.join(INPUTS_DIR, molecule + '.xyz'), 'r') as xyz:
        natoms = int(xyz.readline())
        charge, mult = map(int, xyz.readline().split())
        geometry = []
        nel = 0
        for atom in range(natoms):
            symbol, x, y, z = xyz.readline().split()
            nel += atom_charge(symbol)
    return nel - charge

def get_max_Z(molecule):
    """Get maximal Z for atoms in molecule."""
    return max(Z for Z, _ in get_XYZ(molecule))

def charge_pseudo_atom(Z):
    """Get charge for pseudoatom."""
    if Z <= 2:     # H-He
        return Z
    elif Z <= 10:  # Li-Ne
        return Z - 2
    elif Z <= 18:  # Na-Ar
        return Z - 10
    elif Z <= 30:  # K-Zn
        return Z - 18
    elif Z <= 36:  # Ga-Kr
        return Z - 28
    elif Z <= 48:  # Rb-Cd
        return Z - 36
    elif Z <= 54:  # In-Xe
        return Z - 46
    elif Z <= 71:  # Cs-Lu
        return Z - 54
    elif Z <= 80:  # Hf-Hg
        return Z - 68

def get_lebel_set(molecule):
    """Get set of lebels from xyz-file for every atom type"""
    with open(os.path.join(INPUTS_DIR, molecule + '.xyz'), 'r') as input_geometry:
        natoms = int(input_geometry.readline())
        charge, mult = map(int, input_geometry.readline().split())
        atomic_symbols = []
        for atom in range(natoms):
            symbol, x, y, z = input_geometry.readline().split()
            atomic_symbols.append(symbol)

    result = defaultdict(list)
    for i, item in enumerate(atomic_symbols):
        result[item].append(i+1)
    return result

def get_ae_cutoffs(molecule):
    """Create AE_cutoff initial values.
    Used for Backflow format.
    """
    result = []
    for i, _ in enumerate(get_XYZ(molecule)):
        result.append('{i}         {i}         0.2                          0'.format(i=i+1))
    return '\n  '.join(result)

def casino_time(*path):
    """Get CASINO time.
     Total CASINO CPU time  : : :      378.0500
    """
    regexp = re.compile(' Total CASINO CPU time  : : :\s+(?P<energy_error>\d+\.\d+)')
    try:
        with open(os.path.join(*path, 'out'), 'r') as casino_out:
            # we are only interested in the last occurrence
            return float(re.findall(regexp, casino_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None

def vmc_energy(*path):
    """Get VMC energy without JASTROW optimisation.
     -152.988424660763 +/- 0.003047553900      Correlation time method
    """

    regexp = re.compile(' (?P<energy>[-+]?\d+\.\d+) \+/- (?P<energy_error>[-+]?\d+\.\d+)      Correlation time method')
    try:
        with open(os.path.join(*path, 'out'), 'r') as vmc_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, vmc_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def vmc_variance(*path):
    """Get VMC variance with JASTROW optimisation.
      Sample variance of E_L (au^2/sim.cell) : 3.169677109628 +- 0.034986257092
    """

    regexp = re.compile('Sample variance of E_L \(au\^2/sim.cell\) : (?P<variance>[-+]?\d+\.\d+) \+- (?P<variance_error>[-+]?\d+\.\d+)')
    try:
        with open(os.path.join(*path, 'out'), 'r') as vmc_opt_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, vmc_opt_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def dmc_energy(*path):
    """Get DMC energy.
          mean:   -153.795024411601 +/-       0.001346260888
    """

    dir = os.path.join(*path)
    try:
        open(os.path.join(dir, '.casino_finished'), 'r').close()
        regexp = re.compile('mean:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
        with open(os.path.join(dir, 'out'), 'r') as dmc_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, dmc_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def dmc_stderr(*path):
    """Get DMC standard error.
          stderr:      0.000906128433 +/-       0.000046917552
    """

    dir = os.path.join(*path)
    try:
        open(os.path.join(dir, '.casino_finished'), 'r').close()
        regexp = re.compile('stderr:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
        with open(os.path.join(dir, 'out'), 'r') as dmc_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, dmc_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def dmc_ncorr(*path):
    """Get DMC correlation N.
          N_corr:      0.000906128433 +/-       0.000046917552
    """

    dir = os.path.join(*path)
    try:
        open(os.path.join(dir, '.casino_finished'), 'r').close()
        regexp = re.compile('N_corr:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
        with open(os.path.join(dir, 'out'), 'r') as dmc_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, dmc_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def dmc_stats_nstep(*path):
    """Get DMC statistic accumulation steps.
          dmc_stats_nstep   : 96000
    """

    dir = os.path.join(*path)
    try:
        regexp = re.compile('dmc_stats_nstep   :\s+(?P<nstep>\d+)')
        with open(os.path.join(dir, 'input'), 'r') as dmc_input:
            # we are only interested in the last occurrence
            return int(re.findall(regexp, dmc_input.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None



def get_all_inputs():
    "get file names of all *.xyz input files"
    return sorted((os.path.splitext(filename)[0] for filename in os.listdir(INPUTS_DIR) if os.path.splitext(filename)[1] == '.xyz'))

def pp_basis(basis):
    "check if basis is pseudopotential basis"
    return basis in ('aug-cc-pVDZ-CDF', 'aug-cc-pVTZ-CDF', 'aug-cc-pVQZ-CDF', 'aug-cc-pV5Z-CDF')

wildcard_constraints:
    i = '\d',
    molecule='[-\w+=.]+',
    method='[-\w().]+',
    basis='[-\w]+',
    jastrow_type='[_\w]+',
    jastrow='[_\w]+',
    backflow='[_\w]+',
    opt_plan='\w+',
    nstep='[\d]+',
    nconfig='[\d]+',

####################################################################################################################

rule HF_RESULTS:
    output: 'hf_results.csv'
    run:
        with open(output[0], 'w', newline='') as result_file:
            energy_data = csv.writer(result_file, dialect=csv.unix_dialect, quoting=csv.QUOTE_NONE)
            fieldnames = [
                'method', 'basis', 'molecule', 'hf_energy', 'hf_time',
            ]
            energy_data.writerow(fieldnames)
            for method in METHODS:
                for basis in BASES:
                    for molecule in MOLECULES:
                        path = (method, basis, molecule)
                        try:
                            energy_data.writerow((
                                *path,
                                hf_energy(*path),
                                hf_time(*path),
                            ))
                        except (FileNotFoundError, IndexError) as e:
                            print(e)

rule RESULTS:
    output: 'results.csv'
    run:
        with open(output[0], 'w', newline='') as result_file:
            writer = csv.writer(result_file, dialect=csv.unix_dialect, quoting=csv.QUOTE_NONE)
            fieldnames = [
                'method', 'basis', 'molecule', 'hf_energy', 'hf_time',
                'jastrow',
                'vmc_opt_energy', 'vmc_opt_energy_error', 'vmc_opt_variance', 'vmc_opt_variance_error', 'vmc_opt_time',
                'dmc_energy', 'dmc_energy_error', 'dmc_stderr', 'dmc_stderr_error', 'dmc_ncorr', 'dmc_ncorr_error', 'dmc_time'
            ]
            writer.writerow(fieldnames)
            vmc_path = ('VMC', '10000000')
            for method in METHODS:
                for basis in BASES:
                    for molecule in MOLECULES:
                        path = (method, basis, molecule)
                        for jastrow in JASTROWS:
                            vmc_opt_path = ('VMC_OPT', 'emin', jastrow)
                            dmc_path = ('VMC_DMC', 'emin', jastrow, 'tmax_2_1024_1')
                            try:
                                writer.writerow((
                                    *path,
                                    hf_energy(*path),
                                    hf_time(*path),
                                    jastrow,
                                    *vmc_energy(*path, *vmc_opt_path),
                                    *vmc_variance(*path, *vmc_opt_path),
                                    casino_time(*path, *vmc_opt_path),
                                    *dmc_energy(*path, *dmc_path),
                                    *dmc_stderr(*path, *dmc_path),
                                    *dmc_ncorr(*path, *dmc_path),
                                    casino_time(*path, *dmc_path),
                                ))
                            except (FileNotFoundError, IndexError) as e:
                                print(e)

rule VMC_DMC_RUN:
    input:      '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_{i}/input',
                '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_{i}/%s' % WFN_FILE,
                '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_{i}/parameters.casl',
                '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_{i}/correlation.data',
    output:     '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_{i}/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'


rule DMC_STATS_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_2/%s' % WFN_FILE,
                '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_1/out',
                '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_2/config.in',
    output:     '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_2/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        stderr, _ = dmc_stderr(wildcards.method, wildcards.basis, wildcards.molecule, 'VMC_DMC', wildcards.opt_plan, wildcards.jastrow, 'tmax_{dt}_{nconfig}_1'.format(dt=wildcards.dt, nconfig=wildcards.nconfig))
        max_Z = get_max_Z(wildcards.molecule)
        if pp_basis(wildcards.basis):
            max_Z = charge_pseudo_atom(max_Z)
        dtdmc = 1.0/(max_Z**2 * 3.0 * int(wildcards.dt))
        if STD_ERR > stderr:
            nstep = 10000
        else:
            nstep = ((stderr/STD_ERR)**2 - 1) * 50000
            nstep = int(round(nstep + 5000, -4))
        if pp_basis(wildcards.basis):
            tmove = 'T'
        else:
            tmove = 'F'
        with open(output[0], 'w') as f:
            f.write(open('../dmc_stats.tmpl').read().format(
                neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep//10000,
                tmove=tmove, basis_type=WFN_TYPE, backflow='F'
            ))
        # workaround in pseudopotential
        if pp_basis(wildcards.basis):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')

rule DMC_STATS_CONFIG:
    input:      '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_1/out',
    output:     '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_2/config.in'
    shell:      'ln -rs "$(dirname "{input}")"/config.out "{output}"'

rule VMC_DMC_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_1/%s' % WFN_FILE,
    output:     '{method}/{basis}/{molecule}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_1/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        max_Z = get_max_Z(wildcards.molecule)
        if pp_basis(wildcards.basis):
            max_Z = charge_pseudo_atom(max_Z)
        dtdmc = 1.0/(max_Z**2 * 3.0 * int(wildcards.dt))
        nstep = 50000
        if pp_basis(wildcards.basis):
            tmove = 'T'
        else:
            tmove = 'F'
        with open(output[0], 'w') as f:
            f.write(open('../vmc_dmc.tmpl').read().format(
                neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep//10000,
                tmove=tmove, basis_type=WFN_TYPE, backflow='F'
            ))
        # workaround in pseudopotential
        if pp_basis(wildcards.basis):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')

rule VMC_DMC_DATA_JASTROW:
    input:      '{path}/VMC_OPT/{opt_plan}/{jastrow}/out',
    output:     '{path}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_{i}/correlation.data'
    shell:      'ln -rs "$(dirname "{input}")/correlation.out.final" "{output}"'

rule VMC_DMC_CASL_JASTROW:
    input:      '{path}/VMC_OPT/{opt_plan}/{jastrow}/out',
    output:     '{path}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_{i}/parameters.casl'
    shell:      'ln -rs "$(dirname "{input}")/parameters.final.casl" "{output}"'

rule VMC_DMC_WFN:
    input:      '{path}/%s' % WFN_FILE,
    output:     '{path}/VMC_DMC/{opt_plan}/{jastrow}/tmax_{dt}_{nconfig}_{i}/%s' % WFN_FILE
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_OPT_ENERGY_RUN:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/input',
                '{method}/{basis}/{molecule}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/%s' % WFN_FILE,
                '{method}/{basis}/{molecule}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/parameters.casl',
                '{method}/{basis}/{molecule}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/correlation.data',
    output:     '{method}/{basis}/{molecule}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_OPT_ENERGY_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/%s' % WFN_FILE
    output:     '{method}/{basis}/{molecule}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../vmc_opt_energy.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule, nstep=wildcards.nstep, basis_type=WFN_TYPE, backflow='F'))
        # workaround in pseudopotential
        if pp_basis(wildcards.basis):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')

rule VMC_OPT_ENERGY_DATA_JASTROW:
    input:      '{path}/VMC_OPT/{opt_plan}/{jastrow}/out'
    output:     '{path}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/correlation.data'
    shell:      'ln -rs "$(dirname "{input}")/correlation.out.final" "{output}"'

rule VMC_OPT_ENERGY_CASL_JASTROW:
    input:      '{path}/VMC_OPT/{opt_plan}/{jastrow}/out'
    output:     '{path}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/parameters.casl'
    shell:      'ln -rs "$(dirname "{input}")/parameters.final.casl" "{output}"'

rule VMC_OPT_ENERGY_WFN:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT_ENERGY/{opt_plan}/{jastrow}/{nstep}/%s' % WFN_FILE
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_OPT_RUN:
    input:      '{method}/{basis}/{molecule}/VMC_OPT/{opt_plan}/{jastrow}/input',
                '{method}/{basis}/{molecule}/VMC_OPT/{opt_plan}/{jastrow}/%s' % WFN_FILE,
                '{method}/{basis}/{molecule}/VMC_OPT/{opt_plan}/{jastrow}/parameters.casl',
                '{method}/{basis}/{molecule}/VMC_OPT/{opt_plan}/{jastrow}/correlation.data',
    output:     '{method}/{basis}/{molecule}/VMC_OPT/{opt_plan}/{jastrow}/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'
                '&& ln -s "$(ls parameters.*.casl | sort -t. -k2,2 -g | tail -1)" parameters.final.casl'
                '&& ln -s "$(ls correlation.out.* | sort -t. -k3,3 -g | tail -1)" correlation.out.final'

rule VMC_OPT_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT/{opt_plan}/{jastrow}/%s' % WFN_FILE
    output:     '{method}/{basis}/{molecule}/VMC_OPT/{opt_plan}/{jastrow}/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../opt_plan/{}.tmpl'.format(wildcards.opt_plan)).read().format(
                neu=neu, ned=ned, nstep=VMC_NCONFIG*10, nconfig=VMC_NCONFIG, molecule=wildcards.molecule, basis_type=WFN_TYPE, backflow='F'
            ))
        # workaround in pseudopotential
        if pp_basis(wildcards.basis):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')

rule VMC_OPT_DATA_JASTROW:
    input:      '{path}/correlation.data'
    output:     '{path}/VMC_OPT/{opt_plan}/{jastrow}/correlation.data'
    shell:      'ln -rs "{input}" "{output}"'

rule VMC_OPT_CASL_JASTROW:
    input:      '{path}/VMC_OPT/{opt_plan}/{jastrow}/%s' % WFN_FILE
    output:     '{path}/VMC_OPT/{opt_plan}/{jastrow}/parameters.casl'
    run:
        with open(output[0], 'w') as f:
            f.write(open('../casl/{}.tmpl'.format(wildcards.jastrow)).read())

rule VMC_OPT_WFN:
    input:      '{path}/%s' % WFN_FILE
    output:     '{path}/VMC_OPT/{opt_plan}/{jastrow}/%s' % WFN_FILE
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_RUN:
    input:      '{method}/{basis}/{molecule}/VMC/{nstep}/input',
                '{method}/{basis}/{molecule}/%s' % WFN_FILE,
                '{method}/{basis}/{molecule}/VMC/{nstep}/correlation.data',
    output:     '{method}/{basis}/{molecule}/VMC/{nstep}/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC/{nstep}/%s' % WFN_FILE
    output:     '{method}/{basis}/{molecule}/VMC/{nstep}/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../vmc.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule, nstep=wildcards.nstep, basis_type=WFN_TYPE))
        # workaround in pseudopotential
        if pp_basis(wildcards.basis):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')

rule VMC_DATA_JASTROW:
    input:      '{path}/correlation.data'
    output:     '{path}/VMC/{nstep}/correlation.data'
    shell:      'ln -rs "{input}" "{output}"'

rule VMC_WFN:
    input:      '{path}/%s' % WFN_FILE
    output:     '{path}/VMC/{nstep}/%s' % WFN_FILE
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_DMC_BF_RUN:
    input:      '{path}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_{i}/input',
                '{path}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_{i}/%s' % WFN_FILE,
                '{path}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_{i}/correlation.data',
                '{path}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_{i}/parameters.casl',
    output:     '{path}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_{i}/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule DMC_STATS_BF_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_2/%s' % WFN_FILE,
                '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_1/out',
                '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_2/config.in',
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_2/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        stderr, _ = dmc_stderr(wildcards.method, wildcards.basis, wildcards.molecule, 'VMC_DMC_BF', wildcards.opt_plan, wildcards.jastrow + '__' + wildcards.backflow, 'tmax_{dt}_{nconfig}_1'.format(dt=wildcards.dt, nconfig=wildcards.nconfig))
        dtdmc = 1.0/(get_max_Z(wildcards.molecule)**2 * 3.0 * int(wildcards.dt))
        if STD_ERR > stderr:
            nstep = 10000
        else:
            nstep = ((stderr/STD_ERR)**2 - 1) * 50000
            nstep = int(round(nstep + 5000, -4))
        if pp_basis(wildcards.basis):
            tmove = 'T'
        else:
            tmove = 'F'
        with open(output[0], 'w') as f:
            f.write(open('../dmc_stats.tmpl').read().format(
                neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep//10000,
                tmove=tmove, basis_type=WFN_TYPE, backflow='T'
            ))

rule DMC_STATS_BF_CONFIG:
    input:      '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_1/out',
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_2/config.in'
    shell:      'ln -rs "$(dirname "{input}")"/config.out "{output}"'

rule VMC_DMC_BF_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_1/%s' % WFN_FILE,
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_1/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        dtdmc = 1.0/(get_max_Z(wildcards.molecule)**2 * 3.0 * int(wildcards.dt))
        nstep = 50000
        if pp_basis(wildcards.basis):
            tmove = 'T'
        else:
            tmove = 'F'
        with open(output[0], 'w') as f:
            f.write(open('../vmc_dmc.tmpl').read().format(
                neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep//10000,
                tmove=tmove, basis_type=WFN_TYPE, backflow='T'
            ))

rule VMC_DMC_BF_DATA_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/out',
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_{i}/correlation.data'
    shell:      'ln -rs "$(dirname "{input}")/correlation.out.final" "{output}"'

rule VMC_DMC_BF_CASL_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/out',
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_{i}/parameters.casl'
    shell:      'ln -rs "$(dirname "{input}")/parameters.final.casl" "{output}"'

rule VMC_DMC_BF_WFN:
    input:      '{method}/{basis}/{molecule}/%s' % WFN_FILE,
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{opt_plan}/{jastrow}__{backflow}/tmax_{dt}_{nconfig}_{i}/%s' % WFN_FILE
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_OPT_BF_ENERGY_RUN:
    input:      '{path}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/input',
                '{path}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/%s' % WFN_FILE,
                '{path}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/correlation.data',
                '{path}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/parameters.casl'
    output:     '{path}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_OPT_BF_ENERGY_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/%s' % WFN_FILE
    output:     '{method}/{basis}/{molecule}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../vmc_opt_energy.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule, nstep=wildcards.nstep, backflow='T'))

rule VMC_OPT_BF_DATA_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/out'
    output:     '{path}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/correlation.data'
    shell:      'ln -rs "$(dirname "{input}")/correlation.out.final" "{output}"'

rule VMC_OPT_BF_CASL_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/out'
    output:     '{path}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/parameters.casl'
    shell:      'ln -rs "$(dirname "{input}")/parameters.final.casl" "{output}"'

rule VMC_OPT_BF_ENERGY_WFN:
    input:      '{path}/%s' % WFN_FILE
    output:     '{path}/VMC_OPT_ENERGY_BF/{opt_plan}/{jastrow}__{backflow}/{nstep}/%s' % WFN_FILE
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_OPT_BF_RUN:
    input:      '{path}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/input',
                '{path}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/%s' % WFN_FILE,
                '{path}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/correlation.data',
                '{path}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/parameters.casl',
    output:     '{path}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'
                '&& ln -s "$(ls parameters.*.casl | sort -t. -k2,2 -g | tail -1)" parameters.final.casl'
                '&& ln -s "$(ls correlation.out.* | sort -t. -k3,3 -g | tail -1)" correlation.out.final'

rule VMC_OPT_BF_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/%s' % WFN_FILE
    output:     '{method}/{basis}/{molecule}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../opt_plan/{}.tmpl'.format(wildcards.opt_plan)).read().format(
                neu=neu, ned=ned, nstep=VMC_NCONFIG*10, nconfig=VMC_NCONFIG, molecule=wildcards.molecule, basis_type=WFN_TYPE, backflow='T'
            ))

rule VMC_OPT_BF_DATA_JASTROW:
    input:      '{method}/{basis}/{molecule}/correlation.data'
    output:     '{method}/{basis}/{molecule}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/correlation.data'
    run:
        backflow = wildcards.backflow.split('_')
        mu_set = open('../backflow_mu_set.tmpl').read()
        phi_set = open('../backflow_phi_set.tmpl').read()
        eta_term = open('../backflow_eta_term.tmpl').read()
        mu_term = open('../backflow_mu_term.tmpl').read()
        phi_term = open('../backflow_phi_term.tmpl').read()
        mu_sets = []
        phi_sets = []
        terms = []
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        for nset, (_, labels) in enumerate(get_lebel_set(wildcards.molecule).items()):
            if backflow[1] != '0':
                mu_sets.append(mu_set.format(
                    number_of_atoms=len(labels),
                    spin_dep=0 if neu == ned else 1,
                    atom_labels=' '.join(map(str, labels)),
                    mu_term=backflow[1],
                    nset=nset + 1
                ))
            if backflow[2] != '00':
                phi_sets.append(phi_set.format(
                    number_of_atoms=len(labels),
                    atom_labels=' '.join(map(str, labels)),
                    phi_term_eN=backflow[2][0], phi_term_ee=backflow[2][1],
                    nset=nset + 1
                ))
        ae_cutoffs = get_ae_cutoffs(wildcards.molecule)
        if backflow[0] != '0':
            terms.append(eta_term.format(eta_term=backflow[0]))
        if backflow[1] != '0':
            terms.append(mu_term.format(number_of_mu_sets=nset + 1, mu_sets='\n'.join(mu_sets)))
        if backflow[2] != '00':
            terms.append(phi_term.format(number_of_phi_sets=nset + 1, phi_sets='\n'.join(phi_sets)))
        with open(output[0], 'w') as f:
            f.write(open('../backflow.tmpl').read().format(
                terms='\n'.join(terms),
                ae_cutoffs=ae_cutoffs
            ))
        source_path = os.path.join(wildcards.method, wildcards.basis, wildcards.molecule, 'correlation.data')
        shell('[[ -e "{source_path}" ]] && cat "{source_path}" >> "{output}"; exit 0')


rule VMC_OPT_BF_CASL_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/%s' % WFN_FILE
    output:     '{method}/{basis}/{molecule}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/parameters.casl'
    run:
        with open(output[0], 'w') as f:
            f.write(open('../casl/{}.tmpl'.format(wildcards.jastrow)).read())


rule VMC_OPT_BF_WFN:
    input:      '{path}/%s' % WFN_FILE
    output:     '{path}/VMC_OPT_BF/{opt_plan}/{jastrow}__{backflow}/%s' % WFN_FILE
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'
