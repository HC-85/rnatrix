import os
import pandas as pd
import subprocess
import re
from nupack import Strand, Model, Complex, ComplexSet, SetSpec, complex_analysis
from os import getcwd, remove
from shutil import move

def abe_score(local_min):
    return local_min['deltaG'].astype('float').mean()/len(local_min.iat[0,0])


def accesibility(mRNA):
    a = Strand(mRNA, name='a')
    my_model = Model(material='RNA', celsius=37)
    pairing = Complex([a])
    set1 = ComplexSet(strands=[a], complexes=SetSpec(max_size=1, include=[pairing]))
    complex_results1 = complex_analysis(set1, model=my_model, compute=['pairs'])
    pairing_result = complex_results1[pairing]

    prob_matrix = pairing_result.pairs.to_array()
    matrix_list=[]
    for i in range(len(mRNA)):
        for j in range(i+1,len(mRNA)):
            matrix_list.append(prob_matrix[i][j])

    return sum([1-x for x in matrix_list])


def ale_score(local_min):
    alt = local_min[local_min.deltaG!=0]
    alt['Estructura P/P'] = alt['Estructura P/P'].apply(lambda x: loop_counter(x))
    return (alt['Estructura P/P']/alt['deltaG'].astype('float')).sum()/len(local_min)


def complement(seq):
    return seq.translate(str.maketrans("AUGC", "UACG"))[::-1]


def ind_dG(seq):
    with open('sDeltaGIND.fa', 'w') as fSeqDG:
        fSeqDG.write('>Secuencia-Temporal')
        fSeqDG.write('\n')
        fSeqDG.write(seq)
    os.system("RNAfold -p -d2 --noLP < sDeltaGIND.fa > sDeltaGIND-fold.out")
    fDeltaGIND=open('sDeltaGIND-fold.out', 'r')
    sRenglonDG=fDeltaGIND.readlines()[3]
    aListDG=re.findall(r"[-+]?(?:\d*\.\d+|\d+)", sRenglonDG)
    dDeltaG=float(aListDG[0])
    os.remove('sDeltaGIND.fa')
    os.remove('sDeltaGIND-fold.out')
    return dDeltaG


def local_min(sRNA,dRangoE,workdir):
    """
    Genera un archivo con los mínimos locales
    """
    
    sRangoE="-e "+str(dRangoE)
    if os.path.isfile('RNAsubopt-test.out')== True:
        os.remove('RNAsubopt-test.out')
    if os.path.isfile('barriers-test.out')== True:
        os.remove('barriers-test.out')
    txtsRNA=open('temp.txt', 'w')
    txtsRNA.write(sRNA)
    txtsRNA.close()
    "Creamos lista de estructuras subóptimas"
    pSubopt=subprocess.Popen(['RNAsubopt', sRangoE, '--infile=temp.txt', '--outfile=RNAsubopt-test.out', '-d2', '--noLP', '-s'])
    (output, err)=pSubopt.communicate()
    pSubopt_status=pSubopt.wait()


    "Llamamos a Barriers para encontrar mínimos locales"
    subprocess.run(["barriers -q --rates -G RNA-noLP --max 50 --minh 0.1 < RNAsubopt-test.out > barriers-test.out", '-l'],shell=True, capture_output=True)
    
    "Neceistamos extraer los valores de la energía de los mínmos locales"
    fLocalMinimos=open("barriers-test.out", "r")
    aMatrizBarriers=[]
    for line in fLocalMinimos:
        stripped_line=line.strip()
        line_list=stripped_line.split()
        aMatrizBarriers.append(line_list)
    fLocalMinimos.close()
    move("RNAsubopt-test.out", f"{workdir}/temp/RNAsubopt-test.out")
    move("barriers-test.out", f"{workdir}/temp/barriers-test.out")
    move("temp.txt", f"{workdir}/temp/temp.txt")

    
    "La matriz aMatrizBarriers tiene dentro nuestros datos a partir de la segunda fila"
    "Creamos una nueva con solo los valores de la energía en mínimos locales"
    dfMinLocal=pd.DataFrame(columns=['Estructura P/P', 'deltaG'])

    for i in range(1, len(aMatrizBarriers)):
        dfMinLocal.loc[len(dfMinLocal.index)]=[aMatrizBarriers[i][1], aMatrizBarriers[i][2]]
    return dfMinLocal


def loop_counter(sEPP):
    """
    Obtiene el número de loops dentro de las estructuras con minimos locales
    sEPP es formato punto paréntesis
    """
    count = 0
    while len(re.findall('\(\.+\)', sEPP))>0:
        inner_loops = re.findall('\(\.+\)', sEPP)
        count = count + len(inner_loops)
        sEPP = re.sub('\(\.+\)','', sEPP)
        while len(re.findall('\(\)', sEPP))>0:
            sEPP = re.sub('\(\)','', sEPP)
    return count


def norm(dEntrada, dMin, dMax):
    return (dEntrada-dMin)/(dMax-dMin)


def seq_score(sRNA):
    aDiccionarioPuntosNonGenic={
        "AA":-0.170,
        "UA":-0.302,
        "GA":-0.146,
        "CA":-0.205,
        "AU":-0.191,
        "UU":-0.075,
        "GU":-0.119,
        "CU":-0.106,
        "AG":0.293,
        "UG":0.214,
        "GG":0.193,
        "CG":-0.091,
        "AC":0.187,
        "UC":0.208,
        "GC":0.051,
        "CC":0.211
    }

    aDiccionarioPuntos={
        "AA":-0.021,
        "UA":-0.071,
        "GA":0.109,
        "CA":-0.002,
        "AU":-0.003,
        "UU":0.2,
        "GU":-0.250,
        "CU":0.033,
        "AG":-0.129,
        "UG":0.081,
        "GG":0.143,
        "CG":-0.084,
        "AC":-0.149,
        "UC":0.067,
        "GC":-0.154,
        "CC":0.247
    }

    dSumDinucleótidos=0
    aMatrizTest=[]
    for i in range(len(sRNA)-1):
        dSumDinucleótidos=dSumDinucleótidos+aDiccionarioPuntos[sRNA[i:i+2]]
        aMatrizTest.append(aDiccionarioPuntos[sRNA[i:i+2]])
        
    dSeqScore=dSumDinucleótidos/len(sRNA)
    return dSeqScore


"Cálculo de la energía libre de hibridación mRNA-sRNA -dGAsT"
def srna_mrna_energy(sRNA, mRNA):
    fSeq=open('SeqFxHibridacion.fa', 'w+')
    fSeq.write(sRNA)
    fSeq.write('&')
    fSeq.write(mRNA)
    fSeq.close()
    result = subprocess.run([f'RNAcofold -a -d2 --noLP --output-format=D < SeqFxHibridacion.fa > cofold-hibridacion.csv','-l'], shell=True, capture_output=True)
    dfEnergía = result.stdout
    dfEnergía=pd.read_csv('cofold-hibridacion.csv')
    os.remove('cofold-hibridacion.csv')
    os.remove('SeqFxHibridacion.fa')
    dDeltaGibbs=dfEnergía.iat[0, 6]
    return float(dDeltaGibbs)


def sRNA_to_scores(srna, chaperone, dE, workdir):
    srna_chap = srna + chaperone
    u_score_srna = u_score(srna_chap)
    seq_score_srna = seq_score(srna_chap)
    local_min_srna = local_min(srna_chap, dE,)
    abe_score_srna = abe_score(local_min_srna)
    ale_score_srna = ale_score(local_min_srna)
    return seq_score_srna, u_score_srna, abe_score_srna, ale_score_srna
    

def u_score(sRNA):
    Ucounter10 = sum(map(lambda x: x=='U', list(sRNA[-10:])))
    Ucounter = sum(map(lambda x: x=='U', list(sRNA)))     
    return (Ucounter10/10)/(Ucounter/len(sRNA))

def a_density(sRNA):
    sRNA = sRNA.replace('a','A')
    return (sRNA.count('A'))/(len(sRNA))

def u_density(sRNA):
    sRNA = sRNA.replace('u','U')
    return (sRNA.count('U'))/(len(sRNA))

def g_density(sRNA):
    sRNA = sRNA.replace('g','U')
    return (sRNA.count('G'))/(len(sRNA))

def c_density(sRNA):
    sRNA = sRNA.replace('c','U')
    return (sRNA.count('C'))/(len(sRNA))

def BSP(sequence):
    s = subprocess.check_output(["RNAup", "-b" , '-d2'], text = True ,input= sequence)
    s = re.findall('([0-9]+),([0-9]+)', s)
    os.remove('RNA_w25_u1.out')
    return(s)

def harmonic_mean(x):
        return 2/(1/x[0] + 1/x[1])