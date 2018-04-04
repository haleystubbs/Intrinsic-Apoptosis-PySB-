

# p53 is sequestered and degraded by mdm2, and p53 transcriptionally upregulates mdm2 forming a negative feedback loop
# Study affects on half life of p53 as well as thresholds for apoptosis activation and cytoplasmic accumulation

from pysb import *
from pysb.util import alias_model_components

from earm import shared
from earm import lopez_modules
from earm import albeck_modules

"""_________________________________________________________________________________________________________________________________________"""

Model()

#p53 is a tetramer present in the nucleus and the cytoplsm, and canbe ubiquitinated, or phosphorylated at 3 different sites
# p53 can be ubiquitinatd on a cluster of six different lys Lys-370, Lys-372, Lys-373, Lys-381, Lys-382 and Lys-386
# p53 can be phosphorylated at Ser-15, Ser-20 and Ser-46 in response to DNA damage
# p53 won't bind mdm2 when phosphorylated
# ser15 phos by ATM, ATR, CHK1, and DNAPK
# ser20 phos by CHK2
# ser46 phos by HIPK2, and PKCdelta
# Cterminal phos at sites 392 (by CKII) or at 371, 376, 378 (by PKC) enhance DNA binding activity
Monomer('p53', ['bf', 'compartment', 'tag', 'phosC', 'phosN'], {'compartment' : ['N', 'C'], 'tag': ['U'], 'phosN': ['None', 'P15', 'P20', 'P46', 'P152046', 'P1520', 'P2046', '1546'], 'phosC': ['P1', 'P2', 'P3', 'P4']})
# mdm2 exists in the nucleus and self ubiquitinates, binds to the transactivation domain of p53
Monomer('mdm2', ['bf', 'state'], {'state' : ['U']})
# ATM phos p53 at ser15
Monomer('ATM', ['bf'])
# ATR phos p53 at ser15
Monomer('ATR', ['bf'])
# CHK1 phos p53 at ser15
Monomer('CHK1', ['bf'])
# DNAPK phos p53 at ser15
Monomer('DNAPK', ['bf'])
# CHK2 phos p53 at ser20
Monomer('CHK2', ['bf'])
# HIPK2 phos p53 at ser46
Monomer('HIPK2', ['bf'])
# PKCdelta phos p53 at 392
Monomer('PKCdelta', ['bf'])


def Define_Parameters()
    Parameter('synthesize_mdm2_ks', ) # mdm2 is synthesized in the cytoplasm, but immediately transported to the nucleus
    Parameter('synthesize_p53_ks', ) # p53 is synthesized in the cytoplasm
    Parameter('equilibrate_p53_compartment_kim', ) #p53 is imported into the nucleus
    Parameter('equilibrate_p53_compartment_kex', )  # p53 is exported to the cytoplasm
    Parameter('ATM_bind_p53_k', )
    Parameter('ATM_phos_p53_kc')                    # ATM catalyzes phosphorylation of p53 at ser15
    Parameter('ATR_bind_p53_k', )
    Parameter('ATR_phos_p53_kc', )                  # ATR catalyzes phosphorylation of p53 at ser15
    Parameter('CHK1_bind_p53_', )
    Parameter('CHK1_phos_p53_kc', )                 # CHK1 catalyzes phosphorylation of p53 at ser15
    Parameter('DNAPK_bind_p53_', )
    Parameter('DNAPK_phos_p53_kc', )                 # DNAPK catalyzes phosphorylation of p53 at ser15
    Parameter('CHK1_bind_p53_', )
    Parameter('CHK1_phos_p53_kc', )  # CHK1 catalyzes phosphorylation of p53 at ser15



def Synthesize_and_Degrade_Monomers()
    Rule('synthesize_mdm2', None >> mdm2(bf=None, state=None), synthesize_mdm2_ks)
    Rule('synthesize_p53', None >> p53(bf=None, tag=None, phos=None, compartment= 'C'), synthesize_p53_ks)
    Rule('equilibrate_p53_compartment', p53(bf=None, tag=None, phos=None, compartment='C') | p53(bf=None, tag=None, phos=None, compartment='N'), _equilibrate_p53_compartment_kim, eqilibrate_p53_compartmen_kex)


def p53_phosphorylated_at_Ser15()

    bind_table([[                                           ATM,                ATR,                CHK1,               DNAPK],
                [p53(bf=None, phosN='None'),                                                                                 ],
                [p53(bf=None, phosN='P15'),            None,   50e-9 * N_A * V, 10e-9 * N_A * V],
                [p53(bf=None, phosN='P20'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P46'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P152046'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P1520'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P2046'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P1546'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],

               kf=1e6 / (N_A * V))


def p53_phosphorylated_at_Ser20()
    bind_table([[                                           CHK2],
                [p53(bf=None, phosN='None'),                                                                                 ],
                [p53(bf=None, phosN='P15'),            None,   50e-9 * N_A * V, 10e-9 * N_A * V],
                [p53(bf=None, phosN='P20'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P46'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P152046'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P1520'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P2046'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P1546'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],



def p53_phosphorylated_at_Ser46()
    bind_table([[                                           HIPK2(bf=None), PKCdelta],
                [p53(bf=None, phosN='None'),                                                                                 ],
                [p53(bf=None, phosN='P15'),            None,   50e-9 * N_A * V, 10e-9 * N_A * V],
                [p53(bf=None, phosN='P20'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P46'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P152046'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P1520'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P2046'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],
                [p53(bf=None, phosN='P1546'), None, 50e-9 * N_A * V, 10e-9 * N_A * V]],



def p53_phosphorylated_at_Ser371()



def p53_phosphorylated_at_Ser376()



def p53_phosphorylated_at_Ser378()



def p53_phosphorylated_at_Ser92()



