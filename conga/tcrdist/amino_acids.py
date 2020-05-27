import string

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

#aana = amino_acids + ['a','c','g','t']

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

short_to_long = {}
for rsd in list(longer_names.keys()):short_to_long[longer_names[rsd]] = rsd


HP = {'I': 0.73, 'F': 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37,
      'M': 0.26, 'A': 0.25, 'G': 0.16, 'C': 0.04, 'Y': 0.02,
      'P': -0.07, 'T': -0.18, 'S': -0.26, 'H': -0.40, 'E': -0.62,
      'N': -0.64, 'Q': -0.69, 'D': -0.72, 'K': -1.10, 'R': -1.76}

HP['X'] = sum(HP.values())/20.

GES = {'F': -3.7, 'M': -3.4, 'I': -3.1, 'L': -2.8, 'V': -2.6,
       'C': -2.0, 'W': -1.9, 'A': -1.6, 'T': -1.2, 'G': -1.0,
       'S': -0.6, 'P': 0.2,  'Y': 0.7,  'H': 3.0,  'Q': 4.1,
       'N': 4.8,  'E': 8.2,  'K': 8.8,  'D': 9.2,  'R': 12.3}

GES['X'] = sum(GES.values())/20.

## KD values (Kyte-Doolittle) taken from http://web.expasy.org/protscale/pscale/Hphob.Doolittle.html

KD = {'A': 1.8, 'C': 2.5, 'E': -3.5, 'D': -3.5, 'G': -0.4, 'F': 2.8, 'I': 4.5, 'H': -3.2, 'K': -3.9, 'M': 1.9, 'L': 3.8, 'N': -3.5, 'Q': -3.5, 'P': -1.6, 'S': -0.8, 'R': -4.5, 'T': -0.7, 'W': -0.9, 'V': 4.2, 'Y': -1.3}
assert len(KD) == 20

aa_charge =  {}
for a in amino_acids: aa_charge[a] = 0.0
aa_charge['K'] = 1.0
aa_charge['R'] = 1.0
aa_charge['D'] = -1.0
aa_charge['E'] = -1.0
aa_charge['X'] = 0.0


if __name__ == '__main__':

    ## covariation between different HP scales
    from scipy import stats

    ges = [ -1*GES[x] for x in amino_acids ]
    kd  = [  -1*KD[x] for x in amino_acids ]
    hp  = [  -1*HP[x] for x in amino_acids ]

    slope, intercept, r_ges_kd, p_value, std_err = stats.linregress(ges,kd)
    slope, intercept, r_ges_hp, p_value, std_err = stats.linregress(ges,hp)
    slope, intercept, r_kd_hp, p_value, std_err = stats.linregress(kd,hp)

    print('r_ges_kd:',r_ges_kd)
    print('r_ges_hp:',r_ges_hp)
    print('r_kd_hp:',r_kd_hp)

