import os
import shutil
import numpy as np
import pandas as pd
import time
import matplotlib
import matplotlib.pyplot as plt
from easygui.boxes.fileopen_box import tk
from matplotlib import patches
import easygui


#### Verify the version of the pipeline used and the name of the folder containing your probes design:
version_pipeline = "V4"
Probes_Design = 'Probes_Design'



##### Functions ######


def timer(start_time=None):
    if not start_time:
        start_time = time.time()
        return start_time
    elif start_time:
        #print("---- %s seconds ----" % (time.time() - start_time))
        return time.time() - start_time


start_time = timer()


def select_repertoire():
    repertoire_run = easygui.diropenbox(msg="Select the folder of your run inside the out folder", title= "Folder of the run to analyse")
    repertoire = os.path.dirname(os.path.dirname(repertoire_run))
    return repertoire, repertoire_run


def get_run_name(repertoire_run):
    run_name = os.path.basename(repertoire_run)
    return run_name

repertoire, repertoire_run = select_repertoire()
run_name = get_run_name(repertoire_run)
print('Your working folder is: ' + repertoire)
print('The run name is: ' + run_name)
def Creation_Matrices_Splice(Data):
    if Data['Fusion'].empty:
        print(f"⚠️ Warning : No read for the current gene. Moving to the next.")
    else:
        sondes_G, sondes_D = zip(*[chaine.split('G_') for chaine in Data['Fusion']])
        sondes_G = [chaine + 'G' for chaine in list(sondes_G)]
        Data['Sd_G'] = sondes_G
        Data['Sd_D'] = sondes_D
        splices = pd.DataFrame(0, index = list(set(sondes_G)), columns = list(set(sondes_D)))
        splices = splices.astype(float)
        for i in Data.index:
            splices.loc[Data['Sd_G'][i], Data['Sd_D'][i]] = Data['Count_Full'][i]
        return splices
comment = ''

def Figure_Splice(gene):
    design = pd.read_csv(repertoire + os.path.sep + Probes_Design + os.path.sep +  gene + '_Design.csv', sep = ';', index_col = 'Probe')
    splices = pd.read_csv(repertoire + os.path.sep + 'out' + os.path.sep + run_name + os.path.sep + patient + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices' + os.path.sep + '' + patient + '_Counts_' + gene + '.csv', sep = ';', index_col = 0)
    figure = plt.figure()
    axes = figure.add_subplot(111)
    axes.set_xlim(0, max([int(x) for x in design.loc['Pos']]) * 1.1)
    axes.set_ylim(0, max([int(x) for x in design.loc['Ypos']]) + 50)
    axes.set_frame_on(False)
    axes.xaxis.set_visible(False)
    axes.yaxis.set_visible(False)
    for i in design.columns:
        if design[i]['Type'] == 'I':
            axes.add_artist(matplotlib.lines.Line2D((int(design[i]['Pos']), int(design[i]['Pos']) + int(design[i]['Size'])), (60, 60),linewidth = 0.5, color = 'darkorange'))
        if design[i]['Num'] != '0'and len(str(design[i]['Num'])) < 3:
            axes.text(int(design[i]['Pos']) + int(design[i]['Size'])/2,55, design[i]['Num'], fontsize = 4)
        if design[i]['Type'] == 'G':
            axes.add_artist(patches.Rectangle((int(design[i]['Pos']),58), int(design[i]['Size']), 4, edgecolor = 'black', facecolor = 'lightgrey', fill = True, linestyle = 'solid', linewidth = 0.5, zorder = 1))
        if design[i]['Type'] == 'D' and design[i]['Size'] != '0' :
            axes.add_artist(patches.Rectangle((int(design[i]['Pos']),58), int(design[i]['Size']), 4, edgecolor = 'black', facecolor = 'lightgrey', fill = True, linestyle = 'solid', linewidth = 0.5, zorder = 1))
            lower = 0
    splice_2 = splices.copy()
    left = []
    right = []
    for k in design.columns.tolist():
        if design[k]['Type'] in ['S', 'N', 'P', 'Geno'] and k[-1] == 'G' and k in splice_2.index.tolist():
            splice_2 = splice_2.drop(k, axis = 0)
        if design[k]['Type'] in ['S', 'N', 'P', 'Geno'] and k[-1] == 'D' and k in splice_2.columns.tolist():
            splice_2 = splice_2.drop(k, axis = 1)
    higher = splice_2.max().max() * 1.1
    for droite in splices.columns:
            for gauche in splices.index:
                    YSplice = splices[droite][gauche]
                    start = int(design[gauche]['Pos']) + int(design[gauche]['Size'])
                    stop = int(design[droite]['Pos']) + int(design[droite]['Size'])
                    YposStart = int(design[gauche]['Ypos'])
                    YposStop = int(design[droite]['Ypos'])
                    Nat = design[gauche]['Type']
                    NatD = design[droite]['Type']
                    if stop >= start and YposStart != 0 and Nat != 'S':
                        coul = 'green'
                        ep = 0.4
                        if Nat == 'Geno':
                            YSplice = 0
                            coul = 'blue'
                            ep = 0.5
                        if YSplice >= (higher * 0.2 / 100):
                            high = YposStart + 2 + (YSplice / higher * 40)
                            axes.add_artist(matplotlib.lines.Line2D((start, start+(stop-start)/2), (YposStart+2, high), color = coul,linewidth = ep))
                            axes.add_artist(matplotlib.lines.Line2D((start+(stop-start)/2, stop), (high, YposStart+2), color = coul,linewidth = ep))
                    if stop >= start and YposStart != 0 and Nat == 'S':
                        coul = 'red'
                        ep = 1.3
                        if NatD == 'N':
                            coul = 'green'
                        if YSplice >= (higher * 0.5 / 100):
                            high = YposStart+2 + (YSplice / higher * 45)
                    if start >= stop :
                        ep = 0.3
                        if YSplice >= (higher * 0.5 / 100):
                            high = YposStart-2 - (YSplice / higher * 40)
                            axes.add_artist(matplotlib.lines.Line2D((start, start+(stop-start)/2), (YposStart-2, high), color = 'red',linewidth = ep))
                            axes.add_artist(matplotlib.lines.Line2D((start+(stop-start)/2, stop), (high, YposStop-2), color = 'red',linewidth = ep))
    col = 'black'
    if patient != 'Mean':
        if splice_2.sum().sum() / len(splice_2) < 100 :
            col = 'red'
            TotReads = patient + ' ' + gene + ' n Reads = ' + str(int(splice_2.sum().sum())) + ' Warning : Low-Reads/junctions Counts (' + str(int(splice_2.sum().sum() / len(splice_2))) + ')'
        else:
            TotReads = patient + ' ' + gene + ' n Reads = ' + str(int(splice_2.sum().sum()))
    else:
        TotReads = patient + ' ' + gene
    axes.text(50, 1,TotReads, fontsize = 8, color = col)
    axes.text(10, 8,comment, fontsize = 5, color = col)
    if patient != 'Mean':
        orange_low = False
        blue_high = False
        anos = pd.read_csv(repertoire_run + 'Results' + os.path.sep + 'Ano_Splice_' + patient +  '.csv', index_col = 0, sep = ';')
        anos = anos[anos['Gene'] == gene]
        anos = anos[anos['Norm_Back'] == 'Normal_Splice']
        anos = anos[anos['Jonction_Forte'] == True]
        anos = anos[['Tested_Junction','Start', 'Stop', 'Too_Low_Global', 'Too_High_Global']]
        for i in anos.index:
            if anos['Too_Low_Global'][i] == 'True':
                orange_low = True
            if anos['Too_High_Global'][i] == 'True':
                blue_high = True
        for i in anos.index:
            if anos['Too_Low_Global'][i] == 'True':
                axes.add_artist(matplotlib.lines.Line2D(((anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2), (anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2)), (45, 50), color = 'red',linewidth = 0.6))
                axes.text((anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2-2), 25, anos['Tested_Junction'][i], rotation = 90, fontsize = 3)            
            if anos['Too_Low_Global'][i] == 'q25' and (orange_low == True or Exp_Loss == True):
                axes.add_artist(matplotlib.lines.Line2D(((anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2), (anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2)), (45, 50), color = 'orange',linewidth = 0.3))
                axes.text((anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2-2), 25, anos['Tested_Junction'][i], rotation = 90, fontsize = 3)            
            if anos['Too_High_Global'][i] == 'True':
                axes.add_artist(matplotlib.lines.Line2D(((anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2), anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2), (45, 50), color = 'green',linewidth = 0.6))
                axes.text((anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2-2), 25, anos['Tested_Junction'][i], rotation = 90, fontsize = 3)
            if anos['Too_High_Global'][i] == 'q75' and blue_high == True:
                axes.add_artist(matplotlib.lines.Line2D(((anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2), (anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2)), (45, 50), color = 'green',linewidth = 0.2))
                axes.text((anos['Start'][i] + (anos['Stop'][i] - anos['Start'][i])/2-2), 25, anos['Tested_Junction'][i], rotation = 90, fontsize = 3)        
    return(figure)

## Creation of the depositories
repertoire_run = repertoire + os.path.sep + 'out' + os.path.sep + run_name + os.path.sep
if os.path.exists(repertoire_run + 'Results'):
    shutil.rmtree(repertoire_run + 'Results')

if os.path.exists(repertoire_run + 'Matrices'):
    shutil.rmtree(repertoire_run + 'Matrices')

if os.path.exists(repertoire_run + 'Mean'):
    shutil.rmtree(repertoire_run + 'Mean')

## Patient list
list_patients = os.listdir(repertoire + os.path.sep + 'out' + os.path.sep + run_name)
#remove from list_patients the folders begining with QC, Res, Fail
list_patients = [element for element in list_patients if not element.startswith('QC')]
list_patients = [element for element in list_patients if not element.startswith('Res')]
list_patients = [element for element in list_patients if not element.startswith('Fail')]
list_patients = [element for element in list_patients if not element.startswith('~$')]
list_genes = os.listdir(repertoire + os.path.sep + 'out' + os.path.sep + run_name + os.path.sep + list_patients[0])
list_genes = [element for element in list_genes if not element.startswith('Count_')]
os.mkdir(repertoire_run + 'Matrices')
os.mkdir(repertoire_run + 'Mean')
for gene in list_genes:
    os.mkdir(repertoire_run + 'Mean' + os.path.sep + gene)
    os.mkdir(repertoire_run + 'Mean' + os.path.sep + gene + os.path.sep + 'Counts')
    os.mkdir(repertoire_run + 'Mean' + os.path.sep + gene + os.path.sep + 'Figures')
    os.mkdir(repertoire_run + 'Mean' + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices')
## Matrices
print('Begining normalisation of matrices and generation of reference profiles')   
for gene in list_genes:
    print('Gene :', gene)
    list_columns = []
    for patient in list_patients:
        Data = pd.read_csv(repertoire_run + patient + os.path.sep + gene + os.path.sep + 'Counts' + os.path.sep + patient + '_Counts_' + gene + '.csv', sep = ';')
        for fusion in Data['Fusion']:
            if fusion not in list_columns:
                list_columns.append(fusion)
    Matrice_Full = pd.DataFrame(0, index = list_patients, columns = list_columns)
    Matrice_UMI = pd.DataFrame(0, index = list_patients, columns = list_columns)
    for patient in list_patients:
        Data = pd.read_csv(repertoire_run + patient + os.path.sep + gene + os.path.sep + 'Counts' + os.path.sep+  patient + '_Counts_' + gene + '.csv', sep = ';')
        for i in Data.index:
            Matrice_Full.loc[patient, Data['Fusion'][i]] = Data['Count_Full'][i]

            Matrice_UMI.loc[patient, Data['Fusion'][i]] = Data['Count_UMI'][i]
    # REmove probes for SNP and genomic contamination from normalisation 
    design = pd.read_csv(repertoire + os.path.sep + Probes_Design + os.path.sep + gene + '_Design.csv', sep = ';', index_col = 'Probe')
    SNPGEno = design.loc['Type'][design.loc['Type'].isin(['S', 'N', 'P', 'Geno'])].index.tolist()
    for chaine in SNPGEno:
        colonnes_a_supprimer = [col for col in Matrice_Full.columns if chaine in col]
        Matrice_Full.drop(colonnes_a_supprimer, axis=1, inplace=True)
        Matrice_UMI.drop(colonnes_a_supprimer, axis=1, inplace=True)
    Matrice_Full.insert(0, 'Tot_Reads', Matrice_Full.sum(axis = 1))
    Matrice_UMI.insert(0, 'Tot_UMI', Matrice_UMI.sum(axis = 1))
    Matrice_Full.to_csv(repertoire_run + 'Matrices' + os.path.sep + 'Matrice_Full_' + gene + '_' + run_name + '.csv', sep = ';')
    Matrice_UMI.to_csv(repertoire_run + 'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_' + run_name + '.csv', sep = ';')
    #Matrice_Full.iloc[:, 1:] = Matrice_Full.iloc[:,1:].div(Matrice_Full['Tot_Reads'] / 100, axis = 0)
    #Matrice_UMI.iloc[:, 1:] = Matrice_UMI.iloc[:,1:].div(Matrice_UMI['Tot_UMI'] / 100, axis = 0)
    Matrice_UMI = Matrice_UMI.astype(float)
    Matrice_UMI.loc[:, Matrice_UMI.columns[1:]] = Matrice_UMI.loc[:, Matrice_UMI.columns[1:]].div(Matrice_UMI['Tot_UMI'] / 100, axis=0)
    Matrice_Full= Matrice_Full.astype(float)
    Matrice_Full.loc[:, Matrice_Full.columns[1:]] = Matrice_Full.loc[:, Matrice_Full.columns[1:]].div(Matrice_Full['Tot_Reads'] / 100, axis=0)
    Matrice_Full.to_csv(repertoire_run + 'Matrices' + os.path.sep + 'Matrice_Full_' + gene + '_Norm_' + run_name + '.csv', sep = ';')
    Matrice_UMI.to_csv(repertoire_run + 'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_Norm_' + run_name + '.csv', sep = ';')
    # Profil de reference
    patient = 'Mean'
    Fusion_Mean = []
    Count_Full_Mean = []
    Count_UMI_Mean = []
    for junction in Matrice_Full.columns[1:]:
        Fusion_Mean.append(junction)
        Count_Full_Mean.append(np.mean(Matrice_Full[junction]))
        Count_UMI_Mean.append(np.mean(Matrice_UMI[junction]))
    splices = pd.DataFrame()
    splices['Fusion'] = Fusion_Mean
    splices['Count_Full'] = Count_Full_Mean
    splices['Count_UMI'] = Count_UMI_Mean
    splices = Creation_Matrices_Splice(splices)
    if splices is not None:
        splices.to_csv(repertoire_run + patient + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices' + os.path.sep + patient + '_Counts_' + gene + '.csv',
        sep=';', index=True)
        figure = Figure_Splice(gene)
        figure.savefig(repertoire_run + os.path.sep + patient + os.path.sep + gene + os.path.sep +'Figures'+ os.path.sep + patient + '_Splice_' + gene + '.png', bbox_inches = 'tight', dpi = 400)
        plt.close()
    else:
        print(f"⚠️ Warning : 'splices' is empty for the gene {gene}. No matrice saved.")
    

## Making a list of highly expressed junctions ("Jun_fortes") contaning all junction with a median superior to 50 UMI reads
J_Fortes_UMI = pd.DataFrame(columns = ['Sd_G', 'Sd_D', 'Pos_G', 'Pos_D'])
for gene in list_genes:
    design = pd.read_csv(repertoire + os.path.sep + Probes_Design + os.path.sep +  gene + '_Design.csv', sep = ';', index_col = 'Probe')
    Matrice_UMI = pd.read_csv(repertoire_run + 'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_' + run_name + '.csv', sep = ';', index_col = 0)
    Matrice_UMI = Matrice_UMI.drop(['Tot_UMI'], axis = 1)
    mediane_colonnes = Matrice_UMI.median()
    list_junc_Fortes_UMI = []
    for junction in Matrice_UMI.columns:
        pos_G_ = junction.rfind('G_')
        pos_Eg_ = junction[:pos_G_].rfind('E')
        pos_Ed_ = junction.rfind('E')
        pos_D = junction.rfind('D')
        exon_g = junction[pos_Eg_ + 1:pos_G_]
        exon_d = junction[pos_Ed_ + 1:pos_D]
        median_junction = np.quantile(Matrice_UMI[junction], 0.50)
        if (exon_g.isdigit() and exon_d.isdigit() and (int(exon_d) == int (exon_g) + 1) and median_junction > 50 ) :
            list_junc_Fortes_UMI.append(junction)
    sd_G = []
    sd_D = []
    pos_G = []
    pos_D = []
    val_exp = []
    for i in list_junc_Fortes_UMI:
        sd_gauche = i.split('_')[0]
        sd_G.append(sd_gauche)
        sd_droite = i.split('_')[1]
        sd_D.append(sd_droite)
        pos_G.append(int(design[sd_gauche]['Pos']))
        pos_D.append(int(design[sd_droite]['Pos']))
    J_Fortes = pd.DataFrame(index = list_junc_Fortes_UMI)
    J_Fortes['Sd_G'] = sd_G
    J_Fortes['Sd_D'] = sd_D
    J_Fortes['Pos_G'] = pos_G
    J_Fortes['Pos_D'] = pos_D
    J_Fortes = J_Fortes.sort_values(by = ['Pos_G', 'Pos_D'])
    J_Fortes_UMI = pd.concat([J_Fortes_UMI, J_Fortes])

## Global Matrices
if os.path.exists(repertoire_run + 'Results'):
    shutil.rmtree(repertoire_run + 'Results')
os.mkdir(repertoire_run + 'Results')
os.mkdir(repertoire_run + 'Results' + os.path.sep + 'Exon_Skipping' )
os.mkdir(repertoire_run + 'Results' + os.path.sep + 'Drop' )
os.mkdir(repertoire_run + 'Results' + os.path.sep + 'Dup_or_circRNA' )
os.mkdir(repertoire_run + 'Results' + os.path.sep + 'Dup+++_or_circRNA' )
os.mkdir(repertoire_run + 'Results' + os.path.sep + 'Low_Gene_Expression' )
Matrice_Norm_UMI_Genes = pd.DataFrame(index = list_patients)
for gene in list_genes:
    Matrice_Norm = pd.read_csv(repertoire_run + 'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_' + run_name + '.csv', sep = ';', index_col = 0)
    Matrice_Norm = Matrice_Norm.astype(float)
    for i in Matrice_Norm.columns:
        Matrice_Norm_UMI_Genes[i] = Matrice_Norm[i].tolist()
        Matrice_Norm_UMI_Genes = Matrice_Norm_UMI_Genes.copy()
Matrice_Norm_UMI_Genes = Matrice_Norm_UMI_Genes.astype(float)
Matrice_Norm_UMI_Genes['Tot_Reads'] = 0
Matrice_Norm_UMI_Genes['Tot_Reads'] = Matrice_Norm_UMI_Genes.sum(axis=1)
Matrice_Norm_UMI_Genes.iloc[:, 1:] = Matrice_Norm_UMI_Genes.iloc[:,1:].div(Matrice_Norm_UMI_Genes['Tot_Reads'] / 100, axis = 0)
Matrice_Norm_UMI_Genes.to_csv(repertoire_run + 'Matrices' + os.path.sep + 'Matrice_Norm_All_Genes.csv', sep = ';')
for patient in list_patients:
    print('Patient :', patient)
    G_e = []
    T_j = []
    V_p = []
    mean = []
    q25 = []
    q50 = []
    q75 = []
    outlier_H = []
    outlier_L = []
    start = []
    stop = []
    V_p_Global = []
    q25_Global = []
    q50_Global = []
    q75_Global = []
    outlier_H_Global = []
    outlier_L_Global = []
    # Design et normalised data : Table Expression
    for gene in list_genes:
        design = pd.read_csv(repertoire + os.path.sep + Probes_Design + os.path.sep +  gene + '_Design.csv', sep = ';', index_col = 'Probe')
        Matrice_Norm = pd.read_csv(repertoire_run + 'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_Norm_' + run_name + '.csv', sep = ';', index_col = 0)
        new_index = {old_index: str(old_index) for old_index in Matrice_Norm.index}
        Matrice_Norm.rename(index = new_index, inplace=True)
        Sondes_Gauche = design.columns[design.loc['Type'].str.contains('G')].tolist()
        Sondes_Droite = design.columns[design.loc['Type'].str.contains('D')].tolist()
        for Sd_G in Sondes_Gauche:
            for Sd_D in Sondes_Droite:
                Tested_Jonction = Sd_G + '_' + Sd_D
                if Tested_Jonction in Matrice_Norm.columns:
                    G_e.append(gene)
                    start.append(int(design[Sd_G]['Pos']))
                    stop.append(int(design[Sd_D]['Pos']))
                    T_j.append(Tested_Jonction)
                    V_p.append(Matrice_Norm[Tested_Jonction][patient])
                    mean.append(np.mean(Matrice_Norm[Tested_Jonction]))
                    q25.append(np.quantile(Matrice_Norm[Tested_Jonction], 0.25))
                    q50.append(np.quantile(Matrice_Norm[Tested_Jonction], 0.50))
                    q75.append(np.quantile(Matrice_Norm[Tested_Jonction], 0.75))
                    outlier_H.append(np.quantile(Matrice_Norm[Tested_Jonction], 0.75) + 1.5*(np.quantile(Matrice_Norm[Tested_Jonction], 0.75)-np.quantile(Matrice_Norm[Tested_Jonction], 0.25)))
                    outlier_L.append(np.quantile(Matrice_Norm[Tested_Jonction], 0.25) - 1.5*(np.quantile(Matrice_Norm[Tested_Jonction], 0.75)-np.quantile(Matrice_Norm[Tested_Jonction], 0.25)))
                    V_p_Global.append(Matrice_Norm_UMI_Genes[Tested_Jonction][patient])
                    q25_Global.append(np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction],0.25))
                    q50_Global.append(np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction],0.50))
                    q75_Global.append(np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction],0.75))
                    outlier_H_Global.append(np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction], 0.75) + 1.5*(np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction], 0.75)-np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction], 0.25)))
                    outlier_L_Global.append(np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction], 0.25) - 1.5*(np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction], 0.75)-np.quantile(Matrice_Norm_UMI_Genes[Tested_Jonction], 0.25)))

    Table_Exp = pd.DataFrame()
    Table_Exp['Gene'] = G_e
    Table_Exp['Tested_Junction'] = T_j
    Table_Exp['Start'] = start
    Table_Exp['Stop'] = stop
    Table_Exp['Q25-1.5IQR'] = outlier_L    
    Table_Exp['Exp_Val'] = V_p
    Table_Exp['Q75+1.5IQR'] = outlier_H
    Table_Exp['Mean'] = mean
    Table_Exp['q25'] = q25
    Table_Exp['q50'] = q50
    Table_Exp['q75'] = q75
    Table_Exp['Q25-1.5IQR_Global'] = outlier_L_Global
    Table_Exp['Exp_Val_Global'] = V_p_Global
    Table_Exp['Q75+1.5IQR_Global'] = outlier_H_Global  
    Table_Exp['q25_global'] = q25_Global
    Table_Exp['q50_global'] = q50_Global
    Table_Exp['q75_global'] = q75_Global
    Table_Exp = Table_Exp.drop(Table_Exp[(Table_Exp['Exp_Val'] == 0) & (Table_Exp['q25'] == 0) & (Table_Exp['q50'] == 0) & (Table_Exp['q75'] == 0)].index)
    J_Fortes = []
    Norm_Back = []
    for i in Table_Exp.index:
        if Table_Exp['Tested_Junction'][i] in J_Fortes_UMI.index:
            J_Fortes.append(True)
        else:
            J_Fortes.append(False)
        if Table_Exp['Stop'][i] > Table_Exp['Start'][i]:
            Norm_Back.append('Normal_Splice')
        else:
            Norm_Back.append('Back_Splice')
    Table_Exp['Jonction_Forte'] = J_Fortes
    Table_Exp['Norm_Back'] = Norm_Back
    Table_Exp['Too_High'] = Table_Exp.apply(lambda row: 'q75' if row['Exp_Val'] > row['q75'] else '', axis=1)
    Table_Exp['Too_Low'] = Table_Exp.apply(lambda row: 'q25' if row['Exp_Val'] < row['q25'] else '', axis=1)
    Table_Exp['Too_High'] = Table_Exp.apply(lambda row: True if row['Exp_Val'] > row['Q75+1.5IQR'] else row['Too_High'],axis=1)
    Table_Exp['Too_Low'] = Table_Exp.apply(lambda row: True if row['Exp_Val'] < row['Q25-1.5IQR'] else row['Too_Low'],axis=1)
    Table_Exp['Too_High_Global'] = Table_Exp.apply(lambda row: 'q75' if row['Exp_Val_Global'] > row['q75_global'] else '', axis=1)
    Table_Exp['Too_Low_Global'] = Table_Exp.apply(lambda row: 'q25' if row['Exp_Val_Global'] < row['q25_global'] else '', axis=1)
    Table_Exp['Too_High_Global'] = Table_Exp.apply(lambda row: True if row['Exp_Val_Global'] > row['Q75+1.5IQR_Global'] else row['Too_High_Global'], axis=1)
    Table_Exp['Too_Low_Global'] = Table_Exp.apply(lambda row: True if row['Exp_Val_Global'] < row['Q25-1.5IQR_Global'] else row['Too_Low_Global'], axis=1)    
    Table_Exp.to_csv(repertoire_run + 'Results' + os.path.sep + 'Ano_Splice_' + patient +  '.csv', index = True, sep = ';')



###############################################################################################
## Search for splicing anomalies
print('Searching for splicing anomalies')
episs_patient = []
episs_gene = []
episs_junction = []
episs_junction_L = []
episs_junction_R = []
episs_type = []
episs_score = []
episs_FC = []
episs_percentage_skipping_or_drop_code = []
episs_percentage_of_internal_junction_inf_q25 = []
episs_patient_UMI = []
episs_mean_UMI = []
episs_tot_UMI = []
episs_mean_tot_UMI = []

for patient in list_patients:
    Exp_Loss = False
    print()
    print('Patient : ', patient)
    Table_Exp=pd.read_csv(repertoire_run + 'Results' + os.path.sep + 'Ano_Splice_' + patient +  '.csv', index_col = 0, sep = ';')

# Exon skipping and duplication/circRNA
    Ex_Skip_1 = []
    Back_Splice = []
    add_exS = 2
    for ano in ['Normal_Splice', 'Back_Splice']:
        if ano == 'Normal_Splice':
            event = 'Exon_Skipping'
        else:
            event = 'Dup_or_circRNA'
        for i in Table_Exp.index:
            gene = Table_Exp['Gene'][i]
            sub = Table_Exp[(Table_Exp['Gene'] == Table_Exp['Gene'][i]) & (Table_Exp['Jonction_Forte'] == True) & (Table_Exp['Norm_Back'] == 'Normal_Splice')]
            min_exp = np.mean(sub['q50'])
            if (((gene == 'PMS2' and Table_Exp['Exp_Val'][i] > min_exp / 50) or (gene != 'PMS2' and Table_Exp['Exp_Val'][i] > min_exp / 25)) and (Table_Exp['Too_High'][i] == 'True' and Table_Exp['Norm_Back'][i] == ano and ano == 'Normal_Splice') or Table_Exp['Exp_Val'][i] > min_exp / 3 and Table_Exp['Too_High'][i] == 'True' and Table_Exp['Norm_Back'][i] == ano and ano == 'Back_Splice'):
                Pass = True
                if ano == 'Normal_Splice':
                    sub = Table_Exp[(Table_Exp['Gene'] == Table_Exp['Gene'][i]) & (Table_Exp['Norm_Back'] == ano) & (Table_Exp['Start'] >= Table_Exp['Start'][i]) & (Table_Exp['Stop'] <= Table_Exp['Stop'][i])]
                    sub = sub.drop(index = i)
                    sub_for_mean = Table_Exp[(Table_Exp['Gene'] == Table_Exp['Gene'][i]) & (Table_Exp['Norm_Back'] == ano) & (Table_Exp['Start'] >= Table_Exp['Start'][i]) & (Table_Exp['Stop'] <= Table_Exp['Stop'][i]) & Table_Exp['Jonction_Forte'] == True]
                    inf_q25 = 0
                    inf_IQR = 0
                    mean_sub = np.mean(sub_for_mean['Exp_Val'])
                    for j in sub.index:
                        if sub['Exp_Val'][j] < sub['Q25-1.5IQR'][j]:
                            inf_IQR += 1
                        if sub['Exp_Val'][j] < sub['q25'][j]:
                            inf_q25 += 1
                if ano == 'Back_Splice':
                    sub = Table_Exp[(Table_Exp['Gene'] == Table_Exp['Gene'][i]) & (Table_Exp['Norm_Back'] == 'Normal_Splice') & (Table_Exp['Start'] >= Table_Exp['Stop'][i]) & (Table_Exp['Stop'] <= Table_Exp['Start'][i]) & (Table_Exp['Jonction_Forte'] == True)]
                    inf_q25 = 0
                    inf_IQR = 0
                    sub_for_mean = Table_Exp[
                        (Table_Exp['Gene'] == Table_Exp['Gene'][i]) & (Table_Exp['Norm_Back'] == 'Back_Splice') & (
                                    Table_Exp['Start'] >= Table_Exp['Stop'][i]) & (
                                    Table_Exp['Stop'] <= Table_Exp['Start'][i])]
                    mean_sub = np.mean(sub_for_mean['Exp_Val'])
                    for j in sub.index:
                        if sub['Exp_Val'][j] > sub['Q75+1.5IQR'][j]:
                            inf_IQR += 1
                        if sub['Exp_Val'][j] > sub['q75'][j]:
                            inf_q25 += 1
                junction = Table_Exp['Tested_Junction'][i]
                if ano == 'Normal_Splice':
                    pos_G_ = junction.rfind('G_')
                    pos_Eg_ = junction[:pos_G_].rfind('E')
                    exon_g = ''.join([c for c in junction[pos_Eg_+1:pos_G_] if c.isdigit()])
                    pos_Ed_ = junction.rfind('E')
                    pos_D = junction.rfind('D')
                    exon_d = ''.join([c for c in junction[pos_Ed_+1:pos_D] if c.isdigit()])
                    if exon_g.isdigit() and exon_d.isdigit() and int(exon_d) == int(exon_g) + 1:
                        Pass = False
                Matrice_UMI = pd.read_csv(
                    repertoire_run + os.path.sep + 'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_' + run_name + '.csv',
                    index_col=0, sep=';')
                patient_UMI = Matrice_UMI.loc[patient, Table_Exp['Tested_Junction'][i]]
                tot_UMI = Matrice_UMI['Tot_UMI'][patient]
                mean_UMI = np.mean(Matrice_UMI[Table_Exp['Tested_Junction'][i]])
                mean_tot_UMI = np.mean(Matrice_UMI['Tot_UMI'])
                if Table_Exp['q50'][i] > 0:
                    seuil = (Table_Exp['Exp_Val'][i] / Table_Exp['q50'][i])
                    ratio_patient_q50 = ' (' + event + ') Fold_Change:_'  + str(round(Table_Exp['Exp_Val'][i] / Table_Exp['q50'][i], 2)) + '_or_' + str(round(Table_Exp['Exp_Val'][i] / mean_sub * 100, 2)) + '%_(' + str(round(Table_Exp['Exp_Val'][i], 2)) + '_vs_meanGeneExp_' + str(round(min_exp,2)) + '_vs_mean_sub_'+ str(round(mean_sub, 2)) + ')_('
                    Fold_change = str(round(Table_Exp['Exp_Val'][i] / Table_Exp['q50'][i], 2))
                    percentage_skipping = str(round(Table_Exp['Exp_Val'][i] / mean_sub * 100, 2)) 
                else:
                    ratio_patient_q50 = ' (' + event +') Fold_Change:_inf_or_' + str(round(Table_Exp['Exp_Val'][i] / mean_sub * 100, 2)) + '%_(' +  str(round(Table_Exp['Exp_Val'][i], 2)) + '_vs_meanGeneExp_' + str(round(min_exp,2)) + '_vs_mean_sub _' + str(round(mean_sub,2)) +')_('
                    seuil = 3
                    Fold_change = 'inf'
                    percentage_skipping = str(round(Table_Exp['Exp_Val'][i] / mean_sub * 100, 2)) 
                    if ano == 'Back_Splice':
                        event = 'Dup+++_or_circRNA'
                if patient_UMI < 5 or (5 <= patient_UMI <= 20 and Table_Exp['Exp_Val'][i] < Table_Exp['Mean'][i] * 10) :
                    Pass = False
                if (len(sub) > 0
                        and (inf_q25 / len(sub) >= 0.10 or gene == 'PMS2')
                        and Pass == True
                        and seuil > 2
                        or (len(sub) == 0 and ano == 'Back_Splice' and Pass == True and seuil > 2)):
                    if ano == 'Normal_Splice':
                        Ex_Skip_1.append(ratio_patient_q50 + str(len(sub)) + '-IntJunc : ' + str(inf_IQR)+'_inf-IQR, ' + str(inf_q25)+'_inf-q25)')
                        percentage_of_internal_junctions_inf_q25 = str(round(inf_q25 / len(sub) * 100, 2))

                    else:
                        Ex_Skip_1.append(ratio_patient_q50 + str(len(sub)) + '-IntJunc : ' + str(inf_IQR)+'_sup-IQR, ' + str(inf_q25)+'_sup-q75)')
                        
                        percentage_of_internal_junctions_inf_q25 = "NA"
                    shutil.copy2(repertoire_run + os.path.sep + 'Mean' + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Mean_Splice_' + gene + '.png', repertoire_run + os.path.sep +'Results' + os.path.sep + event + os.path.sep + gene + '__Mean.png')
                    comment = str(Table_Exp['Tested_Junction'][i]) + '  ' + Ex_Skip_1[-1] + '  ' + str(patient_UMI) + '  ' + str(mean_UMI)
                    print(patient, comment)
                    figure = Figure_Splice(gene)
                    if gene + '_' + patient + '.png' in os.listdir(repertoire_run + os.path.sep +'Results' + os.path.sep + event):
                        figure.savefig(repertoire_run + os.path.sep +'Results' + os.path.sep + event + os.path.sep + gene + '_' + patient + '_' + str(add_exS) + '.png', bbox_inches = 'tight', dpi = 400)
                        add_exS += 1
                        print('n :', add_exS)
                        plt.close()
                    else:
                        figure.savefig(repertoire_run + os.path.sep +'Results' + os.path.sep + event + os.path.sep + gene + '_' + patient + '.png', bbox_inches = 'tight', dpi = 400)
                        plt.close()
                    episs_patient.append(patient)
                    episs_gene.append(gene)
                    episs_junction.append(Table_Exp['Tested_Junction'][i])
                    episs_type.append(event)
                    episs_score.append(Ex_Skip_1[-1])
                    episs_FC.append (Fold_change)
                    episs_percentage_skipping_or_drop_code.append (percentage_skipping)
                    episs_percentage_of_internal_junction_inf_q25.append (percentage_of_internal_junctions_inf_q25)
                    episs_patient_UMI.append(patient_UMI)
                    episs_mean_UMI.append(mean_UMI)
                    episs_tot_UMI.append(tot_UMI)
                    episs_mean_tot_UMI.append(mean_tot_UMI)
                else:
                    Ex_Skip_1.append('')
            else:
                Ex_Skip_1.append('')        
    Table_Exp['Ex_Skip_1'] = Ex_Skip_1[0:int((len(Ex_Skip_1)/2))]
    Table_Exp['Back_Splice'] = Ex_Skip_1[int((len(Ex_Skip_1)/2)):]
    print()
    
    # Drops (unique junctions)
    Drops = []
    add_other = 2
    Table_Exp['Score_Drop'] = ''
    low_junc_row = Table_Exp[(Table_Exp['Too_Low'] == 'True') & (Table_Exp['Jonction_Forte'] == True) & (Table_Exp['Norm_Back'] == 'Normal_Splice')].index
    sub_ini = Table_Exp[(Table_Exp['Jonction_Forte'] == True) & (Table_Exp['Norm_Back'] == 'Normal_Splice')]
    for i in low_junc_row:
        junc = sub_ini['Tested_Junction'][i]
        gene = Table_Exp['Gene'][i]
        Matrice_Norm = pd.read_csv(repertoire_run + os.path.sep +'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_Norm_' + run_name + '.csv', sep = ';', index_col = 0)
        Matrice_Norm.index = Matrice_Norm.index.astype(str)
        sub = sub_ini[sub_ini['Gene'] == gene]
        list_junc_ok = []
        for junction in sub['Tested_Junction']:
            pos_G_ = junction.rfind('G_')
            pos_Eg_ = junction[:pos_G_].rfind('E')
            pos_Ed_ = junction.rfind('E')
            pos_D = junction.rfind('D')
            exon_g = junction[pos_Eg_+1:pos_G_]
            exon_d = junction[pos_Ed_+1:pos_D]
            if (exon_g.isdigit() and exon_d.isdigit()) or (junction == junc):
                list_junc_ok.append(junction)
        sub = sub[sub['Tested_Junction'].isin(list_junc_ok)]
        j = sub.index.tolist().index(i)        
        start_j = sub['Start'][sub.index[j]]
        sub = sub[(sub['Start'] != start_j) | ((sub['Start'] == start_j) & (sub['Tested_Junction'] == junc))]        
        j = sub.index.tolist().index(i)
        if sub.index.tolist()[j] != sub.index.tolist()[0]:
            pre_junc = sub['Tested_Junction'][sub.index.tolist()[j-1]]
        else:
            pre_junc = []
        if sub.index.tolist()[j] != sub.index.tolist()[-1]:
            post_junc = sub['Tested_Junction'][sub.index.tolist()[j+1]]
        else:
            post_junc = []
        Vals_junc = np.array(Matrice_Norm[junc])
        Val_patient = Matrice_Norm[junc][patient]
        Drop = (Val_patient / np.quantile(Vals_junc, 0.50) * 100)
        if Drop <= 75:
            score = 200
        else:
            score = 100
        if len(pre_junc) > 0:
            Vals_pre = np.array(Matrice_Norm[pre_junc])
            ratio_pre = Vals_junc / Vals_pre
            min_ratio_pre = np.quantile(ratio_pre, 0.25) - 1.5 * (np.quantile(ratio_pre, 0.75) - np.quantile(ratio_pre, 0.25))
            ratio_pre_junc = Matrice_Norm[junc][patient] / Matrice_Norm[pre_junc][patient]
            if ratio_pre_junc < min_ratio_pre:
                score += 10
        if len(post_junc) > 0:
            Vals_post = np.array(Matrice_Norm[post_junc])
            ratio_post = Vals_junc / Vals_post
            min_ratio_post = np.quantile(ratio_post, 0.25) - 1.5 * (np.quantile(ratio_post, 0.75) - np.quantile(ratio_post, 0.25))
            ratio_post_junc = Matrice_Norm[junc][patient] / Matrice_Norm[post_junc][patient]
            if ratio_post_junc < min_ratio_post:
                score += 1
        
        sd_G = junc.split('_')[0]
        sd_D = junc.split('_')[1]
        Sub_G = Table_Exp[(Table_Exp['Gene'] == gene) & (Table_Exp['Start'] == Table_Exp['Start'][i]) & (Table_Exp['Norm_Back'] == 'Normal_Splice') & (Table_Exp['q25'] > 0)]
        Sub_D = Table_Exp[(Table_Exp['Gene'] == gene) & (Table_Exp['Stop'] == Table_Exp['Stop'][i]) & (Table_Exp['Norm_Back'] == 'Normal_Splice') & (Table_Exp['q25'] > 0)]
        cpt_G_low = 0
        for l in Sub_G.index:
            if (Table_Exp['Exp_Val'][l] < Table_Exp['q25'][l]):
                cpt_G_low += 1
        cpt_D_low = 0
        for l in Sub_D.index:
            if (Table_Exp['Exp_Val'][l] < Table_Exp['q25'][l]):
                cpt_D_low += 1
        to_add_G = '_(' + sd_G + '-x <q25 : ' + str(cpt_G_low) + '/' + str(len(Sub_G)) + '_'  
        to_add_D = 'x-' + sd_D + ' <q25 : ' + str(cpt_D_low) + '/' + str(len(Sub_D)) + ')'
        code = str(round(Drop,1)) + '%' + '_(code_' + str(score)[1:] + ')' + to_add_G + to_add_D
        Matrice_UMI = pd.read_csv(
            repertoire_run + os.path.sep +'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_' + run_name + '.csv', index_col=0,
            sep=';')
        
        patient_UMI = Matrice_UMI.loc[patient, Table_Exp['Tested_Junction'][i]]
      
        tot_UMI = Matrice_UMI['Tot_UMI'][patient]
        mean_UMI = np.mean(Matrice_UMI[Table_Exp['Tested_Junction'][i]])
        mean_tot_UMI = np.mean(Matrice_UMI['Tot_UMI'])
        if score > 200:
            Table_Exp.loc[i, 'Score_Drop'] = code
            episs_patient.append(patient)
            episs_gene.append(gene)
            episs_junction.append(junc)
            episs_type.append('Low_Junction')
            episs_score.append('drop to : ' + code)
            episs_FC.append (Val_patient / np.quantile(Vals_junc, 0.50))
            episs_percentage_skipping_or_drop_code.append (str(score)[1:])
            episs_percentage_of_internal_junction_inf_q25.append ("NA")
            episs_patient_UMI.append(patient_UMI)
            episs_mean_UMI.append(mean_UMI)
            episs_tot_UMI.append(tot_UMI)
            episs_mean_tot_UMI.append(mean_tot_UMI)
            comment = junc + ' drop to : ' + code
            print(patient, comment)
            shutil.copy2(repertoire_run + os.path.sep +'Mean' + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Mean_Splice_' + gene + '.png', repertoire_run +os.path.sep + 'Results' + os.path.sep + 'Drop' + os.path.sep + gene + '__Mean.png')
            figure =Figure_Splice(gene)
            if gene + '_' + patient + '.png' in os.listdir(repertoire_run + os.path.sep +'Results' + os.path.sep + 'Drop'):
                figure.savefig(repertoire_run + os.path.sep +'Results' + os.path.sep + 'Drop' + os.path.sep + gene + '_' + patient + '_' + str(add_other) + '.png', bbox_inches = 'tight', dpi = 400)
                add_other += 1
                print('n :', add_other)
                plt.close()
            else:
                figure.savefig(repertoire_run + os.path.sep +'Results' + os.path.sep + 'Drop' + os.path.sep + gene + '_' + patient + '.png', bbox_inches = 'tight', dpi = 400)
                plt.close()
    Table_Exp.to_csv(repertoire_run + os.path.sep +'Results' + os.path.sep + 'Ano_Splice_' + patient +  '.csv', index = True, sep = ';')

    # Low gene expression
    for gene in list_genes:
        comment = ''
        sub = Table_Exp[Table_Exp['Gene'] == gene]
        sub = sub[sub['Norm_Back'] == 'Normal_Splice']
        sub = sub[sub['Jonction_Forte'] == True]
        
        ind_dG = []
        testj_dG = []
        valpt_dG = []
        q50_dG = []
        for i in sub.index:
            ind_dG.append(i)
            testj_dG.append(sub['Tested_Junction'][i])
            valpt_dG.append(Matrice_Norm_UMI_Genes[sub['Tested_Junction'][i]][patient])
            q50_dG.append(np.quantile(Matrice_Norm_UMI_Genes[sub['Tested_Junction'][i]], .50))
            ratios = [denom / num for num, denom in zip(q50_dG, valpt_dG)]
        
        #sub_low = sub[sub['Too_Low_Global'] == 'True']
        sub_low = sub[sub['Too_Low_Global'] == sub['Too_Low_Global']]
        if len(sub) > 0:
            if 100 * len(sub_low) / len(sub) > 80 and np.median(ratios)*100 < 200/3:
                Exp_Loss = True
                episs_patient.append(patient)
                episs_gene.append(gene)
                episs_junction.append('Global')
                episs_type.append('Low_Gene_Expression')
                Matrice_UMI = pd.read_csv(
                    repertoire_run + 'Matrices' + os.path.sep + 'Matrice_UMI_' + gene + '_' + run_name + '.csv',
                    index_col=0, sep=';')
                patient_UMI = Matrice_UMI.loc[patient, Table_Exp['Tested_Junction'][i]]
                tot_UMI = Matrice_UMI['Tot_UMI'][patient]
                mean_UMI = np.mean(Matrice_UMI[Table_Exp['Tested_Junction'][i]])
                mean_tot_UMI = np.mean(Matrice_UMI['Tot_UMI'])
                episs_patient_UMI.append(patient_UMI)
                episs_mean_UMI.append(mean_UMI)
                episs_tot_UMI.append(tot_UMI)
                episs_mean_tot_UMI.append(mean_tot_UMI)
                comment = 'Gene Expression Ratio / q50 : ' + str(round(np.median(ratios)*100,1)) + '%' + ' & ' + str(round(100 * len(sub_low) / len(sub),1)) + '% Junc. of ' + gene + ' inf q25'
                episs_score.append(comment)
                print(patient, gene, 'Expression_Loss', comment)
                figure =Figure_Splice(gene)
                figure.savefig(repertoire_run + 'Results' + os.path.sep + 'Low_Gene_Expression' + os.path.sep + gene + '_' + patient + '.png', bbox_inches = 'tight', dpi = 400)
                shutil.copy2(repertoire_run + 'Mean' + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Mean_Splice_' + gene + '.png', repertoire_run + 'Results' + os.path.sep + 'Low_Gene_Expression' + os.path.sep + gene + '__Mean.png')
                plt.close()
                
Res_Final = pd.DataFrame()
Res_Final['Patient'] = episs_patient
Res_Final['Gene'] = episs_gene
Res_Final['Jonction'] = episs_junction
Res_Final['Type'] = episs_type
Res_Final['Score'] = episs_score
Res_Final['FC'] = episs_FC
Res_Final['Percentage_Skipping_or_Drop_Code'] = episs_percentage_skipping_or_drop_code
Res_Final['Percentage_of_Internal_Junction_inf_q25'] = episs_percentage_of_internal_junction_inf_q25
Res_Final['UMI_Count'] = episs_patient_UMI
Res_Final['Mean_UMI'] = episs_mean_UMI
Res_Final['Tot_UMI'] = episs_tot_UMI
Res_Final['Mean_Tot_UMI'] = episs_mean_tot_UMI

# For patients with no anomalies detected, adding of a No_Alarm result in the "Type" field:
no_alarm_patients = []
for patient in list_patients:
    if patient not in Res_Final['Patient'].tolist():
        no_alarm_patients.append({'Patient': patient, 'Gene': 'None', 'Jonction': 'None', 'Type': 'No_Alarm', 'Score': '.', 'UMI_Count': '.', 'Mean_UMI': '.'})
# Add no_alarm_patients to the results:
if no_alarm_patients:
    no_alarm_df = pd.DataFrame(no_alarm_patients)
    Res_Final = pd.concat([Res_Final, no_alarm_df], ignore_index=True)

Res_Final.to_csv(repertoire_run + 'Results' + os.path.sep + 'Results_' + run_name +'_UMI_' + version_pipeline +'.csv', sep = ';', index = False)


run_time = time.time() - start_time

#message box pour dire que c'est fini et donner le temps en minute pour réaliser l'analyse
easygui.msgbox('Run analysis is complete! \n'
               'It lasted: ' + str(round(run_time/60, 1)) + ' minutes to analyse the ' + str(len(list_patients)) + ' patients of the run: ' + run_name + ' !')
