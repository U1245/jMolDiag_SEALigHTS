import os
import shutil
import gzip
import time
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import transforms
from matplotlib.patches import Rectangle
import easygui

######################### Instructions #############################

######## - Store your FASTQ in a directory named after the run_name

######## - Copy the directory provided containing the probes design next to the FASTQ directory 

####### - Indicate your list of genes:
list_genes = ['MRE11', 'SMAD4', 'PMS2', 'MSH2', 'TP53', 'APC', 'NBN', 'MBD4', 'MSH3', 'ATM', 'STK11', 'MLH1', 'RAD50', 'RPS20', 'POLE', 'CDH1', 'POT1',  'RINT1', 'BARD1', 'RAD51B', 'BAP1', 'BRCA1', 'MSH6', 'FANCM', 'BUB1', 'CHEK2', 'FAM175A', 'RAD51C', 'BRIP1', 'FAN1', 'PALB2', 'PTEN', 'BMPR1A', 'BRCA2', 'GALNT12', 'RAD51D', 'POLD1', 'MUTYH',  'AXIN2', 'XRCC2', 'RNF43', 'EPCAM', 'NTHL1', 'GREM1']

####### - Indicate your probe design:
Probes_Design = 'Probes_Design'

"""
Directories architecture:
- FASTQ\Run_name
- Probes_Design
- out\Run_name (created by the script if it doesn't already exist)

"""
#----------------------------------------------------------------------------------------------------------------------

## Fonctions
#Fonction to determine how long did the script run for
def timer(start_time=None):
    if not start_time:
        start_time = time.time()
        return start_time
    elif start_time:
        #print("---- %s seconds ----" % (time.time() - start_time))
        return time.time() - start_time


def select_repertoire():
    repertoire_fastq = easygui.diropenbox(msg="Sélectionner le dossier contenant les fichiers \"FastQ\" que vous souhaitez analyser", title="Dossier du run à analyser")
    repertoire = os.path.dirname(os.path.dirname(repertoire_fastq))
    return repertoire, repertoire_fastq


def get_run_name(repertoire_fastq):
    run_name = os.path.basename(repertoire_fastq)
    return run_name


# Extraction of the reads from the FastQ
def From_FastQ_To_Reads(FastQ):
    list_reads_fastq = []
    with gzip.open(FastQ) as f:
        for lineno, line in enumerate(f):
            if (lineno-1) % 4 == 0:
                list_reads_fastq.append(str(line)[2:-3])
    return list_reads_fastq

# Fonction Demultiplexing and counting the junctions from the reads
def Demultiplex_Reads(list_reads):
    Data_Final = pd.DataFrame()
    Data_Final['Sequences'] = list_reads
    Data_Final['Left_Probes'] = 'na'
    Data_Final['len_Left_probes'] = 'na'
    Sequences = Data_Final['Sequences']
    # Sondes G : full genes
    for probe in probes_G.index:
        print(patient, probes_G['probe'][probe], len(Sequences))
        Left_Probe_13bp = probes_G['seq'][probe][0:13]
        indexs_probe = Sequences[(Sequences.str.find(Left_Probe_13bp) == 7) == True].index
        Data_Final['Left_Probes'][indexs_probe] = probes_G['probe'][probe]
        Data_Final['len_Left_probes'][indexs_probe] = len(probes_G['seq'][probe])
        Sequences = Sequences.drop(indexs_probe)
    # Exceptions : sondes G identiques en 5'
    Data = Data_Final[Data_Final['Left_Probes'] == 'ATME63aG']
    Data_1 = Data[Data['Sequences'].str.slice(26).str.startswith('A')]
    Data_2 = Data[Data['Sequences'].str.slice(26).str.startswith('T')]
    Data_Final.loc[Data.index.tolist(), 'Left_Probes'] = 'na'
    Data_Final.loc[Data_1.index.tolist(), 'Left_Probes'] = 'CHEK2Eins2aG'
    Data_Final.loc[Data_2.index.tolist(), 'Left_Probes'] = 'ATME63aG'
    # Sondes D : Gene spe
    Data_Final['Right_Probes'] = 'na'
    for gene in list_genes:
        print(patient, 'right probes', gene)
        probes_Gene = probes_D[probes_D['probe'].str.startswith(gene)]
        Data_Gene = Data_Final[Data_Final['Left_Probes'].str.startswith(gene)].copy()
        if len(Data_Gene) > 0:
            Data_Gene['Droite'] = Data_Gene.apply(lambda row: row['Sequences'][row['len_Left_probes']+7:], axis=1)
            Sequences = Data_Gene['Droite']
            for probe in probes_Gene.index:
                print(patient, probes_Gene['probe'][probe], len(Sequences))
                Right_Probe_13bp = probes_Gene['seq'][probe][0:13]
                indexs_probe = Sequences[Sequences.str.startswith(Right_Probe_13bp) == True].index
                Data_Final['Right_Probes'][indexs_probe] = probes_Gene['probe'][probe]
                Sequences = Sequences.drop(indexs_probe)
    Data_Final = Data_Final.drop('len_Left_probes', axis = 1)
    Data_Final['Right_Probes'] = Data_Final['Right_Probes'].str.replace('_2D', 'D', regex=True)
    Data_Final['Left_Probes'] = Data_Final['Left_Probes'].str.replace('_2G', 'G', regex=True)
    Data_Final['UMI'] = Data_Final['Sequences'].str.slice(0, 7)
    return Data_Final

# Fonction Creation of counting files
def Creation_Matrices(global_count):
    global_gene = global_count[global_count['Left_Probes'].str.contains(gene) & global_count['Right_Probes'].str.contains(gene)].copy()
    global_gene['UMI_Fusion'] = '_'
    global_gene['UMI_Fusion'] = global_gene['UMI'] + global_gene['UMI_Fusion'] + global_gene['Left_Probes'] + global_gene['UMI_Fusion'] + global_gene['Right_Probes']
    global_gene['Fusion'] = '_'
    global_gene['Fusion'] = global_gene['Left_Probes'] + global_gene['Fusion'] + global_gene['Right_Probes']
    UMI_Count = pd.Series(global_gene['UMI_Fusion']).value_counts()
    probes_Gene_D = probes_D[probes_D['probe'].str.startswith(gene)]
    probes_Gene_G = probes_G[probes_G['probe'].str.startswith(gene)]
    list_junctions = []
    for G in probes_Gene_G['probe']:
        for D in probes_Gene_D['probe']:
            list_junctions.append(G + '_' + D)
    list_counts_UMI = [0] * len(list_junctions)
    list_counts_Reads = [0] * len(list_junctions)
    for hit in UMI_Count.index.tolist():
        junction = hit[hit.find('_') + 1:]
        position = list_junctions.index(junction)
        list_counts_UMI[position] += 1
    for hit in global_gene['Fusion'].tolist():
        position = list_junctions.index(hit)
        list_counts_Reads[position] += 1
    Data = pd.DataFrame()
    Data['Fusion'] = list_junctions
    Data['Count_Full'] = list_counts_Reads
    Data['Count_UMI'] = list_counts_UMI
    Data = Data.sort_values(by = 'Count_Full', ascending = False)
    Data = Data.loc[(Data['Count_Full'] != 0) & (Data['Count_UMI'] != 0)]
    return Data

# Fonction Creation of splice matrices 
def Creation_Matrices_Splice(Data):
    sondes_G, sondes_D = zip(*[chaine.split('G_') for chaine in Data['Fusion']])
    sondes_G = [chaine + 'G' for chaine in list(sondes_G)]
    Data['Sd_G'] = sondes_G
    Data['Sd_D'] = sondes_D
    splices = pd.DataFrame(0, index = list(set(sondes_G)), columns = list(set(sondes_D)))
    for i in Data.index:
        splices.loc[Data['Sd_G'][i], Data['Sd_D'][i]] = Data['Count_Full'][i]
    return splices

# Fonction Quality check figures
def Figure_QC(global_count):

    pastel_green = (0.5, 0.7, 0.5)
    fig, axs = plt.subplots(1, 2, figsize=(21, 5), gridspec_kw = {'width_ratios': [1, 4]})
    total_reads = len(global_count)
    total_reads_left_ok = len(global_count[global_count['Left_Probes'] != 'na'])
    total_reads_right_ok = len(global_count[global_count['Right_Probes'] != 'na'])
    pie_data = [total_reads - total_reads_left_ok, total_reads_left_ok - total_reads_right_ok, total_reads_right_ok]
    axs[0].set_title('Probes Alignment :', fontsize = 9)
    axs[0].pie(pie_data, labels=['Fail', 'L Only', 'L+R Ok'], startangle = 90, colors = ['red', 'lightcoral', pastel_green])
    axs[0].text(0.5, 1.2, 'Total : ' + str(total_reads) + ' Reads', transform=axs[0].transAxes, ha='center', fontsize = 12)
    axs[0].text(0.5, 1.31, 'Sample : ' + patient, transform=axs[0].transAxes, ha='center', fontsize = 12)
    axs[0].text(0.5, -0.05, 'Fail : ' + str(round(100 * (total_reads - total_reads_left_ok) / total_reads,2)) + '%', transform=axs[0].transAxes, ha='center', fontsize = 8)
    axs[0].text(0.5, -0.15, 'L Only : ' + str(round(100 * (total_reads_left_ok - total_reads_right_ok) / total_reads,2)) + '%', transform=axs[0].transAxes, ha='center', fontsize = 8)
    if round(100 * (total_reads_right_ok) / total_reads,2) < 70:
        col = 'red'
    else:
        col = pastel_green
    axs[0].text(0.5, -0.25, 'L+R Ok : ' + str(round(100 * (total_reads_right_ok) / total_reads,2)) + '%', transform=axs[0].transAxes, ha='center', fontsize = 10, color = col)
    if total_reads_right_ok < 3000000:
        col = 'red'
    else:
        col = pastel_green
    axs[0].text(0.5, -0.35, 'Reads L+R Ok : ' + str(total_reads_right_ok), transform=axs[0].transAxes, ha='center', fontsize = 10, color = col)
    left_Ok = []
    right_Ok = []
    for gene in list_genes:
        nb_probes = len(pd.read_csv(repertoire + os.path.sep + Probes_Design + os.path.sep +  gene + '_probes.csv', sep = ';')) / 2
        left_Ok.append(len(global_count[global_count['Left_Probes'].str.startswith(gene)]) / nb_probes)
        right_Ok.append(len(global_count[global_count['Left_Probes'].str.startswith(gene) & global_count['Right_Probes'].str.startswith(gene)]) / nb_probes)
    axs[1].bar(np.arange(len(list_genes)), left_Ok, 0.3, color = 'lightcoral', edgecolor='black')
    axs[1].bar(np.arange(len(list_genes)) + 0.3, right_Ok, 0.3, color = pastel_green, edgecolor='black')
    axs[1].set_ylim(0, max(left_Ok) * 1.1)
    positions = np.arange(len(list_genes))
    axs[1].set_xticks(positions + 0.15)
    axs[1].set_xticklabels(list_genes, rotation = 70, fontsize = 8)
    axs[1].set_title('Number of reads per gene normalised by junctions number (Probes Left / Left + Right)')
    return plt

# Fonction Figure Splice
def Figure_Splice(gene):
    design = pd.read_csv(repertoire + os.path.sep + Probes_Design + os.path.sep + gene + '_Design.csv', sep = ';', index_col = 'Probe')
    splices = pd.read_csv(repertoire + os.path.sep + 'out' + os.path.sep + run_name + os.path.sep + patient + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices'+ os.path.sep + patient + '_Counts_' + gene + '.csv', sep = ';', index_col = 0)
    figure = plt.figure()
    axes = figure.add_subplot(111)
    axes.set_xlim(0, max([int(x) for x in design.loc['Pos']]) * 1.1)
    axes.set_ylim(0, max([int(x) for x in design.loc['Ypos']]) + 50)
    axes.set_frame_on(False)
    axes.xaxis.set_visible(False)
    axes.yaxis.set_visible(False)
    for i in design.columns:
        if design[i]['Type'] == 'I':
            axes.add_artist(matplotlib.lines.Line2D((int(design[i]['Pos']), int(design[i]['Pos']) + int(design[i]['Size'])), (60, 60),linewidth = 0.75, color = 'darkorange'))
        if design[i]['Num'] != '0'and len(str(design[i]['Num'])) < 3:
            axes.text(int(design[i]['Pos']) + int(design[i]['Size'])/2,55, design[i]['Num'], fontsize = 4)
        if design[i]['Num'] != '0'and len(str(design[i]['Num'])) > 3:
            axes.text(int(design[i]['Pos']) + int(design[i]['Size'])/2,60 - (len(design[i]['Num']) + 1), design[i]['Num'], rotation = 90, fontsize = 4)
        if design[i]['Type'] == 'G':
            axes.add_artist(patches.Rectangle((int(design[i]['Pos']),58), int(design[i]['Size']), 4, edgecolor = 'black', facecolor = 'lightgrey', fill = True, linestyle = 'solid', linewidth = 0.5, zorder = 1))
        if design[i]['Type'] == 'D' and design[i]['Size'] != '0' :
            axes.add_artist(patches.Rectangle((int(design[i]['Pos']),58), int(design[i]['Size']), 4, edgecolor = 'black', facecolor = 'lightgrey', fill = True, linestyle = 'solid', linewidth = 0.5, zorder = 1))
            lower = 0
    higher = splices.max().max() * 1.1
    for droite in splices.columns:
            for gauche in splices.index:
                    YSplice = splices[droite][gauche]
                    start = int(design[gauche]['Pos']) + int(design[gauche]['Size'])
                    stop = int(design[droite]['Pos']) + int(design[droite]['Size'])
                    YposStart = int(design[gauche]['Ypos'])
                    YposStop = int(design[droite]['Ypos'])
                    Nat = design[gauche]['Type']
                    NatD = design[droite]['Type']
                    # Normal Splice
                    if stop >= start and YposStart != 0 and Nat != 'S':
                        coul = 'green'
                        ep = 0.4
                        if Nat == 'Geno':
                            coul = 'blue'
                            ep = 0.5
                        if YSplice >= (higher * 0.2 / 100):
                            high = YposStart + 2 + (YSplice / higher * 40)
                            axes.add_artist(matplotlib.lines.Line2D((start, start+(stop-start)/2), (YposStart+2, high), color = coul,linewidth = ep))
                            axes.add_artist(matplotlib.lines.Line2D((start+(stop-start)/2, stop), (high, YposStart+2), color = coul,linewidth = ep))
                    # SNPs
                    if stop >= start and YposStart != 0 and Nat == 'S':
                        coul = 'red'
                        ep = 1.3
                        if NatD == 'N':
                            coul = 'green'
                        if YSplice >= (higher * 0.5 / 100):
                            high = YposStart+2 + (YSplice / higher * 45)
                            axes.add_artist(matplotlib.lines.Line2D((stop, stop), (YposStart+2, high), color = coul,linewidth = ep))
                    # Back Splices
                    if start >= stop :
                        ep = 0.3
                        if YSplice >= (higher * 0.5 / 100):
                            high = YposStart-2 - (YSplice / higher * 40)
                            axes.add_artist(matplotlib.lines.Line2D((start, start+(stop-start)/2), (YposStart-2, high), color = 'red',linewidth = ep))
                            axes.add_artist(matplotlib.lines.Line2D((start+(stop-start)/2, stop), (high, YposStop-2), color = 'red',linewidth = ep))
    TotReads = patient + ' ' + gene + ' n Reads = ' + str(splices.sum().sum())
    axes.text(50, 1,TotReads, fontsize = 8)
    return(figure)

##################### Analyse ####################

start_time = timer()
repertoire, repertoire_fastq = select_repertoire()
run_name = get_run_name(repertoire_fastq)
print('Your working directory is: ' + repertoire)
print('The run name is: ' + run_name)

# Patient list:
list_fastq = os.listdir(repertoire + os.path.sep + 'FastQ'+ os.path.sep  + run_name)
list_fastq = [fastq for fastq in list_fastq if '_R1' in fastq]
list_patients = list(set([fichier.split("_")[0] for fichier in list_fastq]))
list_patients.sort()
print('Patients :', len(list_patients))

# Creation of the "out" depositories
if 'out' not in os.listdir(repertoire):
    os.mkdir(repertoire + os.path.sep + 'out')

if os.path.exists(repertoire + os.path.sep + 'out' + os.path.sep + run_name):
    shutil.rmtree(repertoire + os.path.sep + 'out' + os.path.sep  + run_name)
os.mkdir(repertoire + os.path.sep + 'out' + os.path.sep  + run_name)
os.mkdir(repertoire + os.path.sep + 'out' + os.path.sep  + run_name + os.path.sep + 'QC_patients')
for patient in list_patients:
    os.mkdir(repertoire + os.path.sep + 'out' + os.path.sep  + run_name + os.path.sep + patient)
    for gene in list_genes:
        os.mkdir(repertoire + os.path.sep + 'out' + os.path.sep  + run_name + os.path.sep + patient + os.path.sep + gene)
        os.mkdir(repertoire + os.path.sep + 'out' + os.path.sep  + run_name + os.path.sep + patient + os.path.sep + gene + os.path.sep + 'Counts')
        os.mkdir(repertoire + os.path.sep + 'out' + os.path.sep  + run_name + os.path.sep + patient + os.path.sep + gene + os.path.sep + 'Figures')
        os.mkdir(repertoire + os.path.sep + 'out' + os.path.sep  + run_name + os.path.sep + patient + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + 'Matrices')


list_probes_names = []
list_probes_seq = []
for gene in list_genes:
    probes = pd.read_csv(repertoire + os.path.sep + Probes_Design + os.path.sep +  gene + '_Probes.csv', sep = ';')
    list_probes_names.extend(probes['Sonde'])
    list_probes_seq.extend(probes['Seq'])
probes = pd.DataFrame()
probes['probe'] = list_probes_names
probes['seq'] = list_probes_seq
probes_G = probes[probes['probe'].str.endswith('G')]
probes_D = probes[probes['probe'].str.endswith('D')]




######### Main loop  ########

for patient in list_patients:
    print()
    print(patient)
    print('Extraction of the reads')
    os.chdir(repertoire + os.path.sep + 'FastQ'+ os.path.sep + run_name)
    fastqs_patients = [fastq for fastq in list_fastq if patient in fastq]
    list_reads_patient = []
    for fastq in fastqs_patients:
        print(fastq)
        list_reads_patient.extend(From_FastQ_To_Reads(fastq))
    print('Reads :',len(list_reads_patient))
    
    print('Searching for probes')
    list_reads_patient = list_reads_patient[0:]
    global_count = Demultiplex_Reads(list_reads_patient)
    os.chdir(repertoire + os.path.sep + 'out' + os.path.sep+ run_name + os.path.sep + patient)
    
    print('Generation of the counting matrices')
    for gene in list_genes:
        Data = Creation_Matrices(global_count)
        Data.to_csv(repertoire + os.path.sep + 'out' + os.path.sep + run_name + os.path.sep + patient + os.path.sep + gene + os.path.sep + 'Counts' + os.path.sep + patient + '_Counts_' + gene + '.csv', sep = ';', index = False)
        if len(Data) > 0:
            splices = Creation_Matrices_Splice(Data)
        else:
            splices = pd.DataFrame()
        splices.to_csv(repertoire + os.path.sep + 'out' + os.path.sep + run_name + os.path.sep + patient + os.path.sep + gene + os.path.sep + 'Figures'+ os.path.sep + 'Matrices' + os.path.sep + patient + '_Counts_' + gene + '.csv', sep = ';', index = True)

    # Figures Quality Control
    print('QC Figures')
    plt = Figure_QC(global_count)
    plt.savefig(repertoire + os.path.sep + 'out' + os.path.sep + run_name + os.path.sep + 'QC_patients' + os.path.sep +'QC_' + patient + '.png', bbox_inches='tight')
    plt.close()
    # Figures Splice
    print('Splice Figures')
    for gene in list_genes:
        figure = Figure_Splice(gene)
        figure.savefig(repertoire + os.path.sep + 'out' + os.path.sep + run_name + os.path.sep + patient + os.path.sep + gene + os.path.sep + 'Figures' + os.path.sep + patient + '_Splice_' + gene + '.png', bbox_inches = 'tight', dpi = 400)
        plt.close()

run_time = time.time() - start_time

easygui.msgbox('Run demultiplexing is over! You can now run the 2nd script to search for anomalies \n'
               'Demultiplexing time: ' + str(round(run_time/60, 1)) + ' minutes \n'
               'Count files are in the following folder: ' + repertoire + os.path.sep + 'out' + os.path.sep + run_name)
