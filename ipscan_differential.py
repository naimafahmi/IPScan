import csv
import glob, os
from operator import attrgetter
import numpy as np
import pandas as pd
import time
import bisect
import sys
from bisect import bisect_left
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def bi_contains(lst, item):
    return bisect_left(lst, item)

class Stack:
	def __init__(self):
		self.items = []

	def size(self):
		return len(self.items)

	def isEmpty(self):
		return self.items == []

	def push(self, val):
		self.items.append(val)

	def top(self):
		if self.isEmpty():
			return None
		else:
			return self.items[self.size()-1]

	def pop(self):
		if self.isEmpty():
			return None
		else:
			return self.items.pop()

def CountTotalReadCount(chrom, start, end, bam_list, position_row):
	totalCount = 0
	length = end - start + 1
	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)

	if(pos1 < len(bam_list) and pos2 < len(bam_list)):
		if(int(bam_list[pos2][0]) != end):
			pos2 = pos2 - 1

		for t in range(pos1, pos2+1):
			read = int(bam_list[t][1])
			totalCount += read
		
	return totalCount, length

def InsertIntoOldChromDict(ChromDict, chrom, exListNew, geneName, strand):
	GeneDict = ChromDict[chrom]
	if geneName not in GeneDict.keys():
		GeneDict[geneName] = (exListNew, strand)
	else:
		(exList, strand) = GeneDict[geneName]
		exList.extend(exListNew)
		GeneDict[geneName] = (exList, strand)
	return GeneDict

def MakeFullDictionary(ann_df):
	ChromDict = {}
	
	for a_row in ann_df.itertuples():
		chrom = a_row[3]
		geneName = a_row[13].strip()
		exonCount = int(a_row[9])
		exonStartList = a_row[10].split(',')
		exonEndList = a_row[11].split(',')
		strand = a_row[4]

		exList = []
		for i in range(exonCount):
			exonStart = int(exonStartList[i])
			exonEnd = int(exonEndList[i])
			exList.append((exonStart, exonEnd))
		
		if chrom not in ChromDict.keys():
			GeneDict = {}
			GeneDict[geneName] = (exList, strand)
		else:
			GeneDict = InsertIntoOldChromDict(ChromDict, chrom, exList, geneName, strand)

		ChromDict[chrom] = GeneDict

	return ChromDict

def MergeIntervals(inputlist):
	n = len(inputlist)
	inputlist.sort(key = lambda x: x[0], reverse = False)

	st = Stack()
	st.push(inputlist[0])

	for i in range(1,n):
		stacktop_start, stacktop_end = st.top()
		if inputlist[i][0]<= stacktop_end:
			st.pop()
			stacktop_end = max(stacktop_end, inputlist[i][1])
			st.push((stacktop_start, stacktop_end))
		else:
			st.push(inputlist[i])

	mergedList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		mergedList.append(stacktop)
		st.pop()

	mergedList.sort(key = lambda x: x[0], reverse = False)
	return mergedList

def GetIntronList(exonList):
	intronList = []
	start = exonList[0][1] + 1
	for item in exonList[1:]:
		end = item[0] - 1
		intronList.append((int(start), int(end)))
		start = item[1] + 1

	return intronList


def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]


def Generate(ChromDict, input_dir, sample, outdir):
	#peak_df = pd.read_csv("/home/naima/codes/IPScan/generate_peaks_and_motifs/M_Minus_Plus/Peak_positions.csv", delimiter='\t')
	peak_df = pd.read_csv("/home/naima/codes/IPScan/generate_peaks_and_motifs/M_Minus_Plus/Signal_positions.csv", delimiter='\t')

	tt = time.time()
	with open(outdir+"Read_coverages.csv",'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Chrom', 'Gene Name', 'Strand', 'Intron Start', 'Intron End', 'Peak position', 'Ratio_1', 'Ratio_2', 'Ratio_1/Ratio_2', 'R1', 'R2', 'R3', 'L1', 'L2', 'L3', 'C1', 'C2', 'C3'])
		intron_motif_count, found = 0, 0
		intron_flag = []
		count = 0
		for chrom in chromosomes:
			print("Starting:",chrom)
			
			bam_file_reader1 = open(input_dir+'/'+sample+'/'+chrom+".txt", "rt")
			bam_read1 = csv.reader(bam_file_reader1, delimiter="\t")
			bam_list1 = list(bam_read1)
			position_row1 = [int(bam_list1[i][0]) for i in range(len(bam_list1))]

			GeneDict = ChromDict[chrom]
			#print(geneDict)
			for geneID in GeneDict.keys():
				#if geneID == 'Zfp865':
				#print("geneID GeneDict[geneID]", geneID, GeneDict[geneID])
				selected_rows = peak_df.loc[(peak_df['Chrom']==chrom) & (peak_df['Gene Name']==geneID)]
				if len(selected_rows) > 0:
					(exonList, strand) = GeneDict[geneID]
					#print("geneID, exonList", geneID, exonList)
					mergedExonList = MergeIntervals(exonList)
					intronList = GetIntronList(mergedExonList)
					#print("intronList ", mergedExonList, intronList, strand)
					
					for row in selected_rows.itertuples():
						gene_peaks = (str(row[6]).strip('[ ]')).split(',')
						#print(gene_peaks)
						for signal_pos in gene_peaks:
							signal_pos = int(signal_pos)
							for ind in range(len(intronList)):
								(intron_start, intron_end) = intronList[ind]
								#print("intron_start, intron_end", intron_start, intron_end, ind)
								if signal_pos in range(intron_start, intron_end+1):
									intron_motif_count += 1
									R1, L1 = 0, 0
									for (exonStart, exonEnd) in mergedExonList:
										Re, Le = CountTotalReadCount(chrom, exonStart, exonEnd, bam_list1, position_row1)
										R1 += Re
										L1 += Le
									if strand == '+':
										R2, L2 = CountTotalReadCount(chrom, intron_start, signal_pos, bam_list1, position_row1)
										R3, L3 = CountTotalReadCount(chrom, signal_pos+1, intron_end, bam_list1, position_row1)
									else:
										R2, L2 = CountTotalReadCount(chrom, signal_pos, intron_end, bam_list1, position_row1)
										R3, L3 = CountTotalReadCount(chrom, intron_start, signal_pos-1, bam_list1, position_row1)
									
									#print("R1, R2, R3, L1, L2, L3", R1, R2, R3, L1, L2, L3)
									if 0 not in [R1,R2,L1,L2,L3]:
										C1 = float(R1)/L1
										C2 = float(R2)/L2
										C3 = float(R3)/L3
										ratio_1 = float(C2)/C1
										ratio_2 = float(C3)/C1
										#print("dhuksiii#################", C1, C2, C3, ratio_1, ratio_2)
										if ratio_1 > 0 and ratio_2 > 0 and C1 > C2 and C1 > C3 and C2 > C3:
											found += 1
											writer.writerow([chrom, geneID, strand, intron_start, intron_end, signal_pos, ratio_1, ratio_2, round(ratio_1/ratio_2, 3), R1, R2, R3, L1, L2, L3, C1, C2, C3])

										
		f.close()

	print("Time:", time.time() - tt)

	writer_out = pd.ExcelWriter(outdir+"Read_coverages.xlsx", engine='xlsxwriter')
	df = pd.read_csv(outdir+"Read_coverages.csv", delimiter = '\t')
	df.to_excel(writer_out, sheet_name='Sheet1', index = None, header=True)
	os.remove(outdir+"Read_coverages.csv")
	writer_out.save()

	print("intron_motif_count, found", intron_motif_count, found)


def Generate_overlap_with_3seq(ChromDict, input_dir, sample, filename, outdir, ann):
	seq_df = pd.read_csv(outdir+"Peak_positions.csv", delimiter='\t')
	signal_df = pd.read_csv(ann+"_Signal_positions.csv", delimiter='\t')

	with open(outdir+"Read_coverages.csv",'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Chrom', 'Gene Name', 'Strand', 'Intron Start', 'Intron End', 'Peak position', 'Ratio_1', 'Ratio_2', 'Ratio_1/Ratio_2', 'R1', 'R2', 'R3', 'L1', 'L2', 'L3', 'C1', 'C2', 'C3'])
		intron_motif_count, overlap, found = 0, 0, 0
		count = 0
		for chrom in chromosomes:
			print("Starting:",chrom)
			
			bam_file_reader1 = open(input_dir+sample+'/'+chrom+".txt", "rt")
			bam_read1 = csv.reader(bam_file_reader1, delimiter="\t")
			bam_list1 = list(bam_read1)
			position_row1 = [int(bam_list1[i][0]) for i in range(len(bam_list1))]
			
			GeneDict = ChromDict[chrom]
			#print(geneDict)
			for geneID in GeneDict.keys():
				#print("chrom, geneID", chrom, geneID)
				#if geneID == 'PWP2':
					#print(geneID)
				pos_flag, pos_flag2 = [], []
				signal_rows = signal_df.loc[(signal_df['Chrom']==chrom) & (signal_df['Gene Name']==geneID)]
				seq_rows = seq_df.loc[(seq_df['Chrom']==chrom) & (seq_df['Gene Name']==geneID)]
				#print(seq_rows)
				if len(signal_rows) > 0 and len(seq_rows) > 0:
					(exonList, strand) = GeneDict[geneID]
					mergedExonList = MergeIntervals(exonList)
					intronList = GetIntronList(mergedExonList)
					#print("intronList", intronList)

					for row in signal_rows.itertuples():
						gene_pos = (str(row[4]).strip('[ ]')).split(',')
						gene_pos.sort()
						#print("gene_pos", gene_pos)
						for signal_pos in gene_pos:
							signal_pos = int(signal_pos)
							#print("signal_pos", signal_pos)
							if signal_pos not in pos_flag:
								pos_flag.append(signal_pos)
								overlap_flag = 0
								for ind in range(len(intronList)):
									(intron_start, intron_end) = intronList[ind]
									#print("(intron_start, intron_end)", intron_start, intron_end)
									if signal_pos in range(intron_start, intron_end+1):
										intron_motif_count += 1

										for row2 in seq_rows.itertuples():
											positionsList = (str(row2[4]).strip('[ ]')).split(',')
											positionsList.sort()
											#print("positionsList", positionsList)

											for peak_pos in positionsList:
												diff = abs(signal_pos-int(peak_pos))
												if diff<=50:
													overlap_flag = 1

										if overlap_flag == 1:
											#print("overlap_flag ###############")
											overlap += 1
											R1, L1 = 0, 0
											
											for (exonStart, exonEnd) in mergedExonList:
												Re, Le = CountTotalReadCount(chrom, exonStart, exonEnd, bam_list1, position_row1)
												R1 += Re
												L1 += Le
											if strand == '+':
												R2, L2 = CountTotalReadCount(chrom, intron_start, signal_pos, bam_list1, position_row1)
												R3, L3 = CountTotalReadCount(chrom, signal_pos+1, intron_end, bam_list1, position_row1)
											else:
												R2, L2 = CountTotalReadCount(chrom, signal_pos, intron_end, bam_list1, position_row1)
												R3, L3 = CountTotalReadCount(chrom, intron_start, signal_pos-1, bam_list1, position_row1)
											
											#print("R1, R2, R3, L1, L2, L3", R1, R2, R3, L1, L2, L3)
											
											if 0 not in [R1,R2,L1,L2,L3]:
												#print([R1,R2,L1,L2,L3])
												C1 = float(R1)/L1
												C2 = float(R2)/L2
												C3 = float(R3)/L3
												ratio_1 = float(C2)/C1
												ratio_2 = float(C3)/C1
												#print("dhuksiii################# C1, C2, C3, c2/c1, c3/c1", C1, C2, C3, ratio_1, ratio_2)
												if ratio_1 > 0 and ratio_2 > 0 and C1 > C2 and C1 > C3 and C2 > C3:
													found += 1
													writer.writerow([chrom, geneID, strand, intron_start, intron_end, signal_pos, ratio_1, ratio_2, round(ratio_1/ratio_2, 3), R1, R2, R3, L1, L2, L3, C1, C2, C3])
													#print("Found", chrom, geneID, strand, intron_start, intron_end, signal_pos, ratio_1, ratio_2, round(ratio_1/ratio_2, 3), R1, R2, R3, L1, L2, L3, C1, C2, C3)
													#sys.exit()
											
		f.close()

	writer_out = pd.ExcelWriter(outdir+filename, engine='xlsxwriter')
	df = pd.read_csv(outdir+"Read_coverages.csv", delimiter = '\t')
	df.to_excel(writer_out, sheet_name='Sheet1', index = None, header=True)
	os.remove(outdir+"Read_coverages.csv")
	writer_out.save()

	print("intron_motif_count, overlap, found", intron_motif_count, overlap, found)


def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, start, end, peak_pos):
	bam_file_reader= open(pathin+'/'+sample+'/'+chrom+".txt", "rt")
	bam_read = csv.reader(bam_file_reader, delimiter="\t")
	bam_list = list(bam_read)
	position_row = [int(bam_list[i][0]) for i in range(len(bam_list))]

	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)
	#print(chrom, start, end, pos1, pos2)
	if(int(bam_list[pos2][0]) != end):
		pos2 = pos2 - 1

	p = []
	c = []
	read = 0
	length = end - start + 1
	
	for t in range(length):
		p.append(t+start)
		c.append(0)
		
	for t in range(pos1, pos2+1):
		position = int(bam_list[t][0])
		read = int(bam_list[t][1])
		index = p.index(position)
		c[index] = read

	p = np.array(p)
	c = np.array(c)

	caption = ax.fill_between(p,c, color="blue", alpha=0.9)
	#ax.legend(handles = [caption])
	ax.vlines(x=int(peak_pos), ymin=0, ymax=max(c), colors='crimson', linestyles='solid', linewidth=1)
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.set_xlim(start, end)
	ax.autoscale(enable = True)
	ax.set_xticklabels([])
	ax.tick_params(axis='both', bottom=False, which='major', labelsize=8)

def Generate_annotation_plot(ax, isoforms, exonCountList, exonStartList, exonEndList, region_start, region_end):
	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)
	print("isoforms are: ",isoforms)

	ystart = 0
	height = 3
	for i in range(isoforms):
		if i>=15:
			print("15 isoforms of this Gene is plotted.")
			break;
		else:
			ax.hlines(y=(ystart+ystart+height)/2, xmin=region_start, xmax=region_end, linewidth=1, color='skyblue', linestyle = '--')

			ecount = int(exonCountList[i])
			stList = exonStartList[i]
			enList = exonEndList[i]
			for p in range(ecount):
				ex_s = int(stList[p])
				ex_e = int(enList[p])
				if (ex_s >= region_start) and (ex_e <= region_end):
					width = int(enList[p]) - int(stList[p]) + 1
					
					rect = patches.Rectangle((ex_s,ystart), width, height, color = 'skyblue', fill = True)
					ax.add_patch(rect)
		
			ystart +=5

	ax.set_xlim(region_start, region_end)
	ax.autoscale(enable = True)
	return

def Generate_plots_whole_gene(chrom, geneID, strand, peak_pos, input_dir, sample, ann_df, ChromDict):
	GeneDict = ChromDict[chrom]
	exList, strand = GeneDict[geneID]
	mergedExonList = MergeIntervals(exList)
	region_start = int(exList[0][0])
	region_end = int(exList[-1][1])
	#print(region_start, region_end)

	ann_tt = ann_df.loc[(ann_df['name2']==geneID) & (ann_df['chrom']==chrom)]

	exonStartList = {}
	exonEndList = {}
	exonCountList = {}

	isoforms = 0
	for a_row in ann_tt.itertuples():
		exonCount = int(a_row[9])
		exonCountList[isoforms] = exonCount
		exonStartList[isoforms] = a_row[10].split(',')
		exonEndList[isoforms] = a_row[11].split(',')
		isoforms+=1

	#figsize=(12,4)
	fig = plt.figure()
	ax1 = fig.add_subplot(2, 1, 1)
	
	Generate_read_coverate_plot(ax1, input_dir, sample, chrom, geneID, int(region_start), int(region_end), peak_pos)
	ax1.set_ylabel('Read coverage')
	ax2 = fig.add_subplot(2, 1, 2)
	Generate_annotation_plot(ax2, isoforms, exonCountList, exonStartList, exonEndList, region_start, region_end)
	ax2.set_ylabel('Annotation')
	ax2.set_xlabel('Position')
	ax2.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	ax2.spines["left"].set_visible(False)
	ax2.set_yticklabels([])
	plt.setp(ax2.get_xticklabels(), rotation=20, horizontalalignment='right')
	ax2.tick_params(left=False, axis='both', which='major')
	plt.savefig("Figures_wholegene_ratio_5.0/"+geneID+"_"+str(peak_pos)+'.png')


def Generate_plots(chrom, geneID, intron_start, intron_end, peak_pos, input_dir, sample, ann_df, ChromDict, plotout):
	GeneDict = ChromDict[chrom]
	exList, strand = GeneDict[geneID]
	mergedExonList = MergeIntervals(exList)
	
	region_start, region_end = 0, 0
	for ind in range(1, len(mergedExonList)):
		exon_end_prev, exon_start_next = int(mergedExonList[ind-1][1]), int(mergedExonList[ind][0])
		if (exon_end_prev <= intron_start) and (exon_start_next >= intron_end):
			region_start = int(mergedExonList[ind-1][0])
			region_end = int(mergedExonList[ind][1])

	ann_tt = ann_df.loc[(ann_df['name2']==geneID) & (ann_df['chrom']==chrom)]

	exonStartList = {}
	exonEndList = {}
	exonCountList = {}

	isoforms = 0
	for a_row in ann_tt.itertuples():
		exonCount = int(a_row[9])
		exonCountList[isoforms] = exonCount
		exonStartList[isoforms] = a_row[10].split(',')
		exonEndList[isoforms] = a_row[11].split(',')
		isoforms+=1

	#figsize=(12,4)
	fig = plt.figure()
	ax1 = fig.add_subplot(2, 1, 1)
	
	Generate_read_coverate_plot(ax1, input_dir, sample, chrom, geneID, int(region_start), int(region_end), peak_pos)
	ax1.set_ylabel('Read coverage')
	ax2 = fig.add_subplot(2, 1, 2)
	Generate_annotation_plot(ax2, isoforms, exonCountList, exonStartList, exonEndList, region_start, region_end)
	ax2.set_ylabel('Annotation')
	ax2.set_xlabel('Position')
	ax2.spines["top"].set_visible(False)
	ax2.spines["right"].set_visible(False)
	ax2.spines["left"].set_visible(False)
	ax2.set_yticklabels([])
	plt.setp(ax2.get_xticklabels(), rotation=20, horizontalalignment='right')
	ax2.tick_params(left=False, axis='both', which='major')
	plt.title(sample)
	name = geneID+"_"+str(peak_pos)
	plt.savefig(plotout+name+'.png')
	plt.savefig(plotout+name+'.eps', format = 'eps', dpi = 100)
	print(name+" printed successfully")






########################## Main program starts here ########################

ann = 'mm10'
chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
annotation = '/home/naima/codes/IPScan/generate_peaks_and_motifs/mm10_refseq_2019_June20.txt'
input_dir = "/home/naima/input/mouse_M-_M+/RNA-seq_bam/"
sample_namelist = ['Minus_M', 'Plus_M']


"""
ann = 'hg38'
chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19', 'chr20', 'chr21', 'chr22', 'chrX','chrY']
annotation = '/home/naima/codes/IPScan/generate_peaks_and_motifs/hg38_refseq_2018May1.txt'
"""
"""
input_dir = "/home/naima/input/Breast_Cancer_Cell_line/MCF7_mock_torin/RNAseq"
#input_dir = "/home/naima/input/Breast_Cancer_Cell_line/MCF7_mock_torin/RNAseq/STAR_hg38/"
sample_namelist = ['MCF-7_S3', 'MCF-7_Torin_S4']
"""
"""
input_dir = "/home/naima/input/Breast_Cancer_Cell_line/BT549_mock_torin/RNA_seq/hisat2_hg38/"
#input_dir = "/home/naima/input/Breast_Cancer_Cell_line/BT549_mock_torin/RNA_seq/STAR_hg38/"
sample_namelist = ['BT549_Mock', 'BT549_Torin']
"""

ann_df = pd.read_csv(annotation, delimiter = '\t')
ChromDict = MakeFullDictionary(ann_df)

########## These are the main functions which generate the final events in Result_analysis.xlsx file
#for sample in sample_namelist:
sample = sample_namelist[0]    #### Only need to change here
print("Running sample: ", sample)
outdir = sample+"/"
os.makedirs(outdir, exist_ok=True)

######Generate(ChromDict, input_dir, sample, outdir) ####Old method, not using anymore
filename = "Read_coverages_signal_peak_both.xlsx"
final_filename = sample+"_Result_analysis.xlsx"

"""
Generate_overlap_with_3seq(ChromDict, input_dir+sample+"/", sample, filename, outdir, ann)
C2_cutoff = 10
writer_out = pd.ExcelWriter(outdir+final_filename, engine='xlsxwriter')
df = pd.read_excel(outdir+filename)
newdf = df.loc[(df['C2']>C2_cutoff)&(df['Ratio_1/Ratio_2']>=5.0)]
newdf.to_excel(writer_out)
writer_out.save()
"""

########## Now we will generate some plots

"""
plotout = outdir+"Plots/" 
#os.makedirs(plotout, exist_ok=True)
count = 0
plt_pd = pd.read_excel(outdir+final_filename, sheet_name="Sheet1")
plt_list = plt_pd.values.tolist()
#chrom, geneID, start, end, peak_pos, strand = 'chr16', '1810013L24Rik', 8843244, 8855792, 8851681, '+'
chrom, geneID, start, end, peak_pos, strand = 'chr11', 'Vmp1', 86611423, 86635173, 86615305, '-'
Generate_plots(chrom, geneID, start, end, peak_pos, input_dir+sample+"/", sample, ann_df, ChromDict, plotout)
"""

"""
for chrom in chromosomes:
	chr_rows = plt_pd.loc[plt_pd['Chrom']==chrom]
	#print(chrom, chr_rows)
	for eachrow in chr_rows.itertuples():
		geneID, start, end, peak_pos, ratio, strand = eachrow[3], eachrow[5], eachrow[6], eachrow[7], eachrow[10], eachrow[4]
		if (end-start) <= 1000000:
			#print("----------------", chrom, geneID, start, end, peak_pos, ratio, strand)
			#if float(ratio) > 5.0:
				#count += 1
			print(chrom, geneID, start, end)
			Generate_plots(chrom, geneID, start, end, peak_pos, input_dir+sample+"/", sample, ann_df, ChromDict, plotout)
			#Generate_plots_whole_gene(chrom, geneID, strand, peak_pos, input1_dir, s1_namelist, ann_df, ChromDict)
"""