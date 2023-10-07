# one-off script to create motif params in each locus

import json
import os

def set_motif_params(locus):
    if locus in ['IGK']:
        J_TRP_MOTIF = 'FGXG'
        J_TRP_OFFSET = 10
        J_SPLICE = 'CGT'
    elif locus in ['IGL']:
        J_TRP_MOTIF = 'FGXG'
        J_TRP_OFFSET = 10
        J_SPLICE = 'GGT'
    elif locus in 'IGH':
        J_TRP_MOTIF = 'WGXG'
        J_TRP_OFFSET = 11
        J_SPLICE = 'GGT'
    elif locus in 'TRB':
        J_TRP_MOTIF = 'FGXG'
        J_TRP_OFFSET = 10
        J_SPLICE = 'GGT'
    elif locus in 'TRA':
        J_TRP_MOTIF = ['FGXG', 'WGXG', 'CGXG', 'FAXG']
        J_TRP_OFFSET = 11
        J_SPLICE = '*GT'
    elif locus in 'TRD':
        J_TRP_MOTIF = 'FGXG'
        J_TRP_OFFSET = 11
        J_SPLICE = '*GT'
    elif locus in 'TRG':
        J_TRP_MOTIF = ['FGXG', 'WGXG', 'CGXG', 'FAXG']
        J_TRP_OFFSET = 10
        J_SPLICE = '*GT'
    else:
        print(f"Error - no J motifs are defined for locus {locus}")
        exit(0)

    if locus in ['IGH', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD']:
        V_RSS_SPACING = 23
    else:
        V_RSS_SPACING = 12

    if locus in ['IGH', 'IGK']:
        J_RSS_SPACING = 23
    else:
        J_RSS_SPACING = 12

    if locus in ['IGH']:
        D_5_RSS_SPACING = 12
        D_3_RSS_SPACING = 12
    elif locus in ['TRB', 'TRD']:
        D_5_RSS_SPACING = 12
        D_3_RSS_SPACING = 23
    else:
        D_5_RSS_SPACING = 0
        D_3_RSS_SPACING = 0

    return {
        'J_TRP_MOTIF': J_TRP_MOTIF,
        'J_TRP_OFFSET': J_TRP_OFFSET,
        'J_SPLICE': J_SPLICE,
        'V_RSS_SPACING': V_RSS_SPACING,
        'J_RSS_SPACING': J_RSS_SPACING,
        'D_5_RSS_SPACING': D_5_RSS_SPACING,
        'D_3_RSS_SPACING': D_3_RSS_SPACING,
    }

for species in ['human', 'rhesus_macaque']:
    for locus in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD']:
        if os.path.exists(f'./src/digger/motifs/{species}/{locus}'):
            params = set_motif_params(locus)
            with open(f'./src/digger/motifs/{species}/{locus}/motif_params.json', 'w') as f:
                json.dump(params, f, indent=4)
