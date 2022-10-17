import os
import numpy as np
import pandas as pd
from nupack import *
from functions import *
from tqdm import tqdm
from os import listdir, remove
from re import match
from keras.models import load_model
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
import argparse



parser = argparse.ArgumentParser(description="BELINDA", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("win_length", help="Desired sequence length")
parser.add_argument("dE", help="Free energy threshold")
parser.add_argument("--nbest", '-n', type=int, default=5, help="Number of predicted sequences for each desired behavior")
parser.add_argument("--train",  action='store_true')
parser.add_argument("--full",'-f',action='store_true' )
workdir = os.getcwd()
args = parser.parse_args()
argvars = vars(args)



try:
    genefile = listdir(f'{workdir}/input/mRNA')[0]
    with open(f'{workdir}/input/mRNA/{genefile}') as file:
        gene = file.read().upper().replace("T","U")
        gene = gene.replace('\n','')
        gene = gene.replace('\r','')
    Design = 0
except:
    raise Exception('No file was found in the mRNA folder. Move your your desire sRNA to the folder and try again')


try:
    input_sRNA = listdir(f'{workdir}/input/sRNA')[0]
    with open(fr'{workdir}/input/sRNA/{input_sRNA}') as file:
        sRNA = file.read().upper().replace("T","U")
        sRNA = sRNA.replace('\n','')
        sRNA = sRNA.replace('\r','')

    Design = 1
except:
    print('There is nothing in the input/sRNA file, the program is desingning and predicting sRNA´s')
    Desing = 0



shfq='UUUCUGUUGGGCCAUUGCAUUGCCACUGAUUUUCCAACAUAUAAAAAGACAAGCCCGAACAGUCGUCCGGGCUUUUUUUCUCGAG'
sHFQ_DNA='tttctgttgggccattgcattgccactgattttccaacatataaaaagacaagcccgaacagtcgtccgggctttttttctcgag'

print("BELINDAsRNA v1.0")


win_length = int(argvars['win_length'])
dE = int(argvars['dE']) 
Display = argvars['nbest']
train = argvars['train']
full = argvars['full']

if train:
    subprocess.run([f'python' ,f'{workdir}/regulation_predictor/NeuralNetworktraining.py'])

Loadingmodel = tqdm(total = 1 )
Loadingmodel.set_description("Loading the trained model for sRNA prediction")
sRNApredictor = load_model(f'{workdir}/regulation_predictor/sRNAprediction.h5')
Loadingmodel.update(1)

with open(fr'{workdir}/input/mRNA/{genefile}') as file:
    gene = file.read().upper().replace("T","U")

if Design == 1:
    print(f'Processing {genefile}')
    sRNA = sRNA.replace('/n','')
    sRNA = sRNA.replace('/r','')
    sRNA = sRNA.replace('T','U')
    mRNA = gene
    sequence = [sRNA,'&',mRNA]
    sequence = "".join(sequence)

    
    positionBS = BSP(sequence)
    local_min_sRNA = local_min(sRNA,dE,workdir)

    preddata = [None]*(1)
    sRNAalescore = ale_score(local_min_sRNA)
    sRNAabescore = abe_score(local_min_sRNA)
    sRNAuscore = u_score(sRNA)
    sRNAseqscore = seq_score(sRNA)
    sRNAudensity = u_density(sRNA)
    sRNAadensity = a_density(sRNA)
    sRNAgdensity = g_density(sRNA)
    sRNAcdensity = c_density(sRNA)
    sRNABSP = int(positionBS[0][0])
    sRNABSlen = int(positionBS[0][1])-int(positionBS[0][0])
    sRNAMFE = ind_dG(sRNA[int(positionBS[1][0])-1:int(positionBS[1][1])])
    mRNAMFE = ind_dG(mRNA[int(positionBS[0][0])-1:int(positionBS[0][1])])
    mRNAFA = accesibility(mRNA[int(positionBS[0][0])-1:int(positionBS[0][1])])

    bfsRNABSP = sRNA[0:int(positionBS[1][0])]
    aftersRNABSP = sRNA[int(positionBS[1][1]):len(sRNA)]


    if bfsRNABSP == '':
        Before_BSP_sRNA_ABE = 1
        Before_BSP_sRNA_ALE = 1
        Before_BSP_sRNA_len = 0
    else:
        localminsRNABf = local_min(bfsRNABSP,5,workdir)
        Before_BSP_sRNA_ABE = abe_score(localminsRNABf)
        Before_BSP_sRNA_ALE = ale_score(localminsRNABf)
        Before_BSP_sRNA_len = len(bfsRNABSP)

    if aftersRNABSP == '':
        After_BSP_sRNA_ABE = 1
        After_BSP_sRNA_ALE = 1
        After_BSP_sRNA_len = 0
    else:
        localminsRNAafter = local_min(aftersRNABSP,5,workdir)
        After_BSP_sRNA_ABE = abe_score(localminsRNAafter)
        After_BSP_sRNA_ALE = ale_score(localminsRNAafter)
        After_BSP_sRNA_len = len(aftersRNABSP)
    
    HE = srna_mrna_energy(sRNA[int(positionBS[1][0])-1:int(positionBS[1][1])],mRNA[int(positionBS[0][0])-1:int(positionBS[0][1])])



    preddata[0] = [sRNAalescore, sRNAabescore, sRNAuscore, sRNAseqscore, sRNAudensity ,sRNAadensity,sRNAgdensity, sRNAcdensity,
                sRNABSP, sRNABSlen,sRNAMFE, mRNAMFE, mRNAFA,After_BSP_sRNA_ALE,After_BSP_sRNA_ABE,Before_BSP_sRNA_ALE,Before_BSP_sRNA_ABE,After_BSP_sRNA_len,Before_BSP_sRNA_len, HE]

    onlysRNApred=pd.DataFrame(data = preddata ,columns=['ale_score','abe_score','u_score','seq_score','U_density','A_density','G_density','C_density',
        'BSP','BSlen','sRNA_MFE','mRNA_MFE','FA','After_BSP_sRNA_Ale','After_BSP_sRNA_Abe','Before_BSP_sRNA_Ale','Before_BSP_sRNA_Abe',
        'After_BSP_sRNA_len','Before_BSP_sRNA_len','HE'])
    
    onlysRNApred.to_csv(f'{workdir}/output/sRNAtrainingdata.csv')

    sRNAprediction = tqdm(total = 1 )
    sRNAprediction.set_description("Doing some predictions")

    x = onlysRNApred
    for k in [2]:
        for col in x.columns:
            x[col+'_'+str(k)] = x[col].apply(lambda x: x**k)
    
    ct = ColumnTransformer([('numeric',StandardScaler(),x.columns)])
    x = ct.fit_transform(x)
    y = sRNApredictor.predict(x)
    probabilidad = (np.max(y, axis = 1))
    predecido = (np.argmax(y,axis=1))
    if predecido == 1:
        predecido = 'Up regulation'
    else:
        predecido = 'Down regulation'


    sRNAprediction.update(1)

    print(f'The {genefile} and {input_sRNA} interaction is a {predecido} with a probability of {probabilidad}')


if Design == 0:
    print(f'Processing {genefile}')
    datalist = [None]*(len(gene) - win_length + 1)
    predictionlist = [None]*(len(gene) - win_length + 1)

    after_local_min = local_min(shfq,5,workdir)
    after_abe = abe_score(after_local_min)
    after_ale = ale_score(after_local_min)

    DesigningsRNAs = tqdm(total = len(gene) - win_length + 1)

    bindingsitep = 1
    for i in range(len(gene) - win_length + 1):
        DesigningsRNAs.set_description("We are designing some amazing sRNA's %s" % i)
        mrna = gene[i: i+win_length]
        if i==0:
                mRNAtF = gene[i: i+win_length+1]
        
        elif i==(len(gene)-win_length):
            mRNAtF = gene[i-1: i+win_length]
            
        else:
            mRNAtF = gene[i-1: i+win_length+1]
        srna = complement(mrna)
        srna_hfq = srna + shfq
        u_score_srna = u_score(srna_hfq)
        seq_score_srna = seq_score(srna_hfq)
        local_min_srna = local_min(srna_hfq, dE, workdir)
        abe_score_srna = abe_score(local_min_srna)
        ale_score_srna = ale_score(local_min_srna)
        mod1_score = seq_score_srna + u_score_srna - abe_score_srna - ale_score_srna
        hybrid_free_energy = srna_mrna_energy(srna, mrna)
        mrna_free_energy = ind_dG(mrna)
        srna_free_energy = ind_dG(srna)
        acc_factor = accesibility(mrna)
        dG = hybrid_free_energy - mrna_free_energy - srna_free_energy
        mod2_score = acc_factor * mrna_free_energy + mrna_free_energy * hybrid_free_energy + \
                    hybrid_free_energy + mrna_free_energy + acc_factor
        Adensity = a_density(srna_hfq)
        Udensity = u_density(srna_hfq)
        Gdensity = g_density(srna_hfq)
        Cdensity = c_density(srna_hfq)

        
        datalist[i] = [mrna, srna_hfq, seq_score_srna, abe_score_srna, ale_score_srna, u_score_srna, mod1_score,
                    hybrid_free_energy, mrna_free_energy, srna_free_energy, acc_factor, dG, mod2_score, bindingsitep]
        predictionlist[i] = [ale_score_srna, abe_score_srna, u_score_srna, seq_score_srna, Udensity,Adensity,Gdensity,Cdensity,
                    i, win_length,srna_free_energy, mrna_free_energy, acc_factor, after_ale, after_abe, 1,1, 0 ,len(shfq), hybrid_free_energy ]
        
        DesigningsRNAs.update(1)
        bindingsitep = bindingsitep +  1


    for file in listdir():
        if match('.*\.ps$', file) is not None:
            remove(file)
        elif match('^rates', file) is not None:
            remove(file)

    dfsRNA=pd.DataFrame(data = datalist, columns=['mRNA', 'sRNA', 'seq_score', 'abe_score', 'ale_score','u_score', 'mod1_score',
            'hybrid_free_energy', 'mrna_free_energy', 'srna_free_energy', 'acc_factor', 'dG', 'mod2_score', 'BSP'])

    dfPrediction=pd.DataFrame(data = predictionlist, columns=['ale_score','abe_score','u_score','seq_score','U_density','A_density','G_density','C_density',
            'BSP','BSlen','sRNA_MFE','mRNA_MFE','FA','After_BSP_sRNA_Ale','After_BSP_sRNA_Abe','Before_BSP_sRNA_Ale','Before_BSP_sRNA_Abe',
            'After_BSP_sRNA_len','Before_BSP_sRNA_len','HE'])

    sRNAprediction = tqdm(total = 1 )
    sRNAprediction.set_description("Doing some predictions")
    x = dfPrediction
    for k in [2]:
        for col in x.columns:
            x[col+'_'+str(k)] = x[col].apply(lambda x: x**k)

    ct = ColumnTransformer([('numeric',StandardScaler(),x.columns)])
    x = ct.fit_transform(x)
    y = sRNApredictor.predict(x)
    sRNAprediction.update(1)

    dfsRNA['mod1_score_norm'] = dfsRNA['mod1_score'].apply(lambda x: norm(x, dfsRNA['mod1_score'].min(), dfsRNA['mod1_score'].max()))
    dfsRNA['mod2_score_norm'] = dfsRNA['mod2_score'].apply(lambda x: norm(x, dfsRNA['mod2_score'].min(), dfsRNA['mod2_score'].max()))
    dfsRNA['final_original'] = (dfsRNA['mod1_score_norm'] + dfsRNA['mod2_score_norm'])

    dfsRNA = dfsRNA.reset_index()

    dfsRNA = dfsRNA.sort_values(by='mod1_score_norm')
    dfsRNA['mod1_score_reg'] = (dfsRNA.mod1_score_norm/dfsRNA.mod1_score_norm.sum()).cumsum()

    dfsRNA = dfsRNA.sort_values(by='mod2_score_norm')
    dfsRNA['mod2_score_reg'] = (dfsRNA.mod2_score_norm/dfsRNA.mod2_score_norm.sum()).cumsum()

    dfsRNA = dfsRNA.sort_values(by='index').drop(columns='index')

    dfsRNA['sRNA Final score'] = (dfsRNA['mod1_score_reg'] + dfsRNA['mod2_score_reg'])/2
    dfsRNA['Prediction'] = np.argmax(y, axis = 1)
    dfsRNA['Prediction Confidence'] = np.max(y, axis = 1)


    dfsRNA['Harmonic Mean'] = dfsRNA[['sRNA Final score', 'Prediction Confidence']].apply(harmonic_mean, axis=1)
    dfsRNA = dfsRNA.drop(columns = ['final_original','mod2_score_norm','mod1_score_norm','mod2_score','mod1_score'])

    

    best_down = dfsRNA.query('Prediction == 0')[['mRNA', 'sRNA','BSP','sRNA Final score', 'Prediction','Prediction Confidence' ,'Harmonic Mean']].sort_values('Harmonic Mean', ascending=False).head(Display)
    best_up = dfsRNA.query('Prediction == 1')[['mRNA', 'sRNA','BSP','sRNA Final score', 'Prediction','Prediction Confidence' ,'Harmonic Mean']].sort_values('Harmonic Mean', ascending=False).head(Display)
    belinda = pd.concat([best_up, best_down])

    if full:
        dfsRNA.to_csv(f'{workdir}/output/Results.csv')
        
    else:
         belinda.to_csv(f'{workdir}/output/Results.csv', index = False)

    print(f'The results file was created succesfully in {workdir}/output/Results.csv')