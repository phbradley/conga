## http://python.genedrift.org/2007/04/19/genetic-code-part-i/

gencode = {
    'ATA':'I',    #Isoleucine
    'ATC':'I',    #Isoleucine
    'ATT':'I',    # Isoleucine
    'ATG':'M',    # Methionine
    'ACA':'T',    # Threonine
    'ACC':'T',    # Threonine
    'ACG':'T',    # Threonine
    'ACT':'T',    # Threonine
    'AAC':'N',    # Asparagine
    'AAT':'N',    # Asparagine
    'AAA':'K',    # Lysine
    'AAG':'K',    # Lysine
    'AGC':'S',    # Serine
    'AGT':'S',    # Serine
    'AGA':'R',    # Arginine
    'AGG':'R',    # Arginine
    'CTA':'L',    # Leucine
    'CTC':'L',    # Leucine
    'CTG':'L',    # Leucine
    'CTT':'L',    # Leucine
    'CCA':'P',    # Proline
    'CCC':'P',    # Proline
    'CCG':'P',    # Proline
    'CCT':'P',    # Proline
    'CAC':'H',    # Histidine
    'CAT':'H',    # Histidine
    'CAA':'Q',    # Glutamine
    'CAG':'Q',    # Glutamine
    'CGA':'R',    # Arginine
    'CGC':'R',    # Arginine
    'CGG':'R',    # Arginine
    'CGT':'R',    # Arginine
    'GTA':'V',    # Valine
    'GTC':'V',    # Valine
    'GTG':'V',    # Valine
    'GTT':'V',    # Valine
    'GCA':'A',    # Alanine
    'GCC':'A',    # Alanine
    'GCG':'A',    # Alanine
    'GCT':'A',    # Alanine
    'GAC':'D',    # Aspartic Acid
    'GAT':'D',    # Aspartic Acid
    'GAA':'E',    # Glutamic Acid
    'GAG':'E',    # Glutamic Acid
    'GGA':'G',    # Glycine
    'GGC':'G',    # Glycine
    'GGG':'G',    # Glycine
    'GGT':'G',    # Glycine
    'TCA':'S',    # Serine
    'TCC':'S',    # Serine
    'TCG':'S',    # Serine
    'TCT':'S',    # Serine
    'TTC':'F',    # Phenylalanine
    'TTT':'F',    # Phenylalanine
    'TTA':'L',    # Leucine
    'TTG':'L',    # Leucine
    'TAC':'Y',    # Tyrosine
    'TAT':'Y',    # Tyrosine
    'TAA':'_',    # Stop
    'TAG':'_',    # Stop
    'TGC':'C',    # Cysteine
    'TGT':'C',    # Cysteine
    'TGA':'_',    # Stop
    'TGG':'W',    # Tryptophan
    }


## http://www.cmbi.kun.nl/pythoncourse/spy/index.spy?site=python&action=Grand%20Finale&flag=chap

standard = { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
             'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
             'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
             'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
             'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
             'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
             'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
             'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
             'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
             'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
             'ATA': 'M', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
             'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
             'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
             'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
             'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
             'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

## http://www.pasteur.fr/formation/infobio/python/ch15.html

code = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
        'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
        'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
        'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
        'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
        'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
        'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
        'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
        'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
        'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
        'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
        'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
        'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
        'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
        'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
        'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
        }

genetic_code = {}
reverse_genetic_code = {}

for codon in gencode:
    lowcodon = codon.lower()
    assert code[ lowcodon ] == gencode[ codon ] or ( gencode[ codon ] == '_' and code[ lowcodon ] == '*' )
    aa = code[ lowcodon ]
    genetic_code[ lowcodon ] = aa
    if aa not in reverse_genetic_code: reverse_genetic_code[aa] = []
    reverse_genetic_code[aa].append( lowcodon )

if __name__ == "__main__":
    print('hi')
    for c in code:
        d = c.upper()
        if standard[d] != code[c] or gencode[d] != code[c]:
            print(code[c], standard.get(c.upper(),'?'), gencode.get(c.upper(),'?'))
