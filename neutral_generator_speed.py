
import pandas as pd
from pathlib import Path
from itertools import chain
import numpy as np
import json
#think should redo seq etraction, combining eon and exon end first-
# as causing contiuous extraction in extracted exon seqs (ie exon plus codon), 
# but in seqs not combining + adding headache to merging exons
# so preprocess gene coords, combining stop codons with exon 
# sort df and groupby gene, combine lastrow end = row start , apply to group, like merge func here
'''
NOTE - need to edit intron file input so i first intersect with genes and strand info using intron coords 
and bedtools, then remove duplicates

DATA
exon_file - species data with reconstructed windows, all seq on + strand, 
bedtools intersected with full exon coords (without dup or cpg removed, but with overlapping gene segments removed)
all seq and coords contain flanks (added +- 1 to coords )
gene_annotations - contains full gene coords (without duplications removed) 
to add any missing exons in genes in the exon_file-without flanks in coords
OUTPUT
binned exome into equal l_s lengths, and caclulated s_obs per bin

Note - main difference is adding tracking of positions of all possible mutaions
#all added positions are with respect to the POSITIVE STRAND
# so when on - strand start counting fromm last positions  
note never transform coordinates when reverse compliemnt within the intron 


PIPELINE-EXONS
1) load data for seq (ancestral) and target seq
2) remove flanks from exon coords (as not using trinucs in this pipeline)
2) merge segments of an exon if multiple exist in exon file (for seq &target_seq)
3) fill ends of exons, so correct lengths (with flanks still, for seq&target_seq)
4)add list of positions of all bases (so can order with introns)
5)  reverse complement whole exons on negative strand (seq, target_seq)
 + reverse positions
Note - main difference is adding tracking of positions of all possible mutaions
#all added positions are with respect to the POSITIVE STRAND
# so when on - strand start counting fromm last positions  
6) for each exon record list of positions of each subsitions between seq & target seq 
+ alt (to determine if synonymous later)
7) drop target seq-from now only work with seq
7) add any extra exons
8) reconstruct gene seqs ( and positions) by merging exons in genes 
per gene: 
9) iterate through codons and identify if all possible muts are syn or non-syn
10) for syn target extract correpsonding position and alt


#PIPELINE-INTRONS
add positions
rc seq, target seq and reverse order of positions for - strand
for each intron record list of positions of each subsitions between seq & target seq
discard target seq
per gene:
1) iterate through codons and identify all possible muts
2) save position and alt of all possible muts (all synonymous) 
#note as all synonymous do not really need to save alt to match to subs for introns

combined per gene:
identify synoyn subs per gene as overlap between subs and possible muts
calculate l_s = intron and exon possible syn muts per gene, 
if l_s < bin size: s_obs += subs in gene and process next gene , 
else find pos where l_s = bin_size, and add subs before this to s_obs 
then reset counters

need to be careful about saving location of bins, when have a mix of positive and negative strand genes
possibly save list of all ranges in bin?
or as care about genomic proximity, bin negative genes also with same order of positive positions
ie do not reverse order of subs_pos, then get simpler regions
'''

'''
idea2
very slow when running with whle itnron table
instead for exons run normal pipeline + positions and save with trinuc and alt as normal

for intron
save only subs positions (do not need to differentiate alt)
do not need to iterate through possible muts, just if valid base, thn add pos x3


'''
#keep secon sequence up to point of adding positions and reverse compleimenting - then if sub add col of sub
#then when check if for mut types only cosider ancestral sequence
#then in gene mut df still have col of alt
#then add muts that actually happen to the mut matrix!


#for reverse compliment must apply position BEFORE reverse complimenent 
#and sort in descending order
#or do proper transformation fo coordinates

from MutationMatrixGenerator import MutationMatrixGenerator

class NeutralGenerator(MutationMatrixGenerator):
    def __init__(self, species,target_species,  possible_species, exon_file, gene_annotations, intron_file, output_dir):
        super().__init__(species,  possible_species, exon_file, gene_annotations, output_dir)
        self.intron_file = intron_file
        self.target_species = target_species
        
    #edititng to add position data before reverse complimenting
    def load_data(self):
        #drop all species but choosen input
        header = ['chr', 'start', 'end'] + self.possible_species + ['chr_copy', 'exon_start', 'exon_end', 'gene', 'strand']
        
        self.reconstructed_df = pd.read_csv(self.exon_file, sep='\t', names =header)
        self.reconstructed_df.drop(columns='chr_copy', inplace=True)
        #test 
        #self.reconstructed_df = self.reconstructed_df[:1000]

        for speci in self.possible_species:
            if speci == self.species:
                self.reconstructed_df.rename(columns={speci : 'seq'}, inplace=True)
            elif speci == self.target_species:
                self.reconstructed_df.rename(columns={speci : 'target_seq'}, inplace=True)
            else:
                self.reconstructed_df.drop(columns=[speci], inplace=True)
                
        print(f'reconstructed genes = {len(self.reconstructed_df['gene'].unique())}')

        coord_df = pd.read_csv(self.gene_annotations, sep='\t', names=['chr', 'exon_start', 'exon_end', 'gene', 'gene_name', 'strand'])
        coord_df.drop(columns=['gene_name'], inplace=True)
       
        gene_list = self.reconstructed_df['gene'].unique()
        cut_coord_df = coord_df[coord_df['gene'].isin(gene_list)]
        del coord_df

        self.reconstructed_df = self.remove_flanks(self.reconstructed_df)
        #input segments of exons, output full exons and corresponding trinucs 
        exon_merged_df = self.apply_segment_merges()
        full_exons_df = self.fill_gaps_in_exons(exon_merged_df)
        full_exons_df['positions'] = full_exons_df.apply(lambda row : range(row['exon_start'], row['exon_end']), axis=1)
        #reverse comp seqs
        rc_df = self.apply_reverse_complement_seq_target(full_exons_df)
        #reverse positions
        rc_df['positions'] = rc_df.apply(lambda row: row['positions'][::-1] if row['strand'] == '-' else row['positions'],axis=1)
        #subs columns contains lists of positions of subs, subs_alts contains alts
        #to later check if subs are synonymous or not 
        rc_df[['subs', 'subs_alts']] = rc_df.apply(lambda row :pd.Series(self.find_subs(row['seq'], row['target_seq'], row['positions'])),axis=1)
        rc_df.drop(columns='target_seq', inplace=True)
        self.all_exons_df = self.add_extra_exons(rc_df, cut_coord_df) #also fill mutations column
       
    def load_intron_data(self):
        header = ['chr', 'start', 'end'] + self.possible_species + ['gene', 'strand']
        intron_df = pd.read_csv(self.intron_file, sep="\t", names=header)
        #test 
        #intron_df = intron_df[intron_df['gene'].isin(self.all_exons_df['gene'].values)]
    
        for speci in self.possible_species:
            if speci == self.species:
                intron_df.rename(columns={speci : 'seq'}, inplace=True)
            elif speci == self.target_species:
                intron_df.rename(columns={speci : 'target_seq'}, inplace=True)
            else:
                intron_df.drop(columns=[speci], inplace=True)
        #intron_df['positions'] = intron_df.apply(lambda row : range(row['start'], row['end']),axis=1)
        
        rc_intron_df = self.apply_reverse_complement_seq_target(intron_df)
        #rc_intron_df['positions'] = rc_intron_df.apply(lambda row: row['positions'][::-1] if row['strand'] == '-' else row['positions'],axis=1)
        rc_intron_df['subs'] = rc_intron_df.apply(lambda row : self.find_subs_intron(row['seq'], row['target_seq'], row['start'], row['end'], row['strand']),axis=1)
        rc_intron_df.drop(columns='target_seq',inplace=True)
        self.intron_df_group = rc_intron_df.groupby('gene')


    def find_subs(self, seq, target_seq, positions):
        df = pd.DataFrame({'seq':list(seq), 'target_seq':list(target_seq), 'positions':list(positions)})
        ok = df[['seq', 'target_seq']].isin(['A', 'C', 'G', 'T'])
        # rows where ALL species are valid
        valid_rows = ok.all(axis=1)
        df =df.loc[valid_rows, :]
       
        df = df[df['seq']!=df['target_seq']]
        vals = df[['positions','target_seq']].to_numpy()
        return vals[:,0].tolist(), vals[:,1].tolist()
    
    
    def find_subs_intron(self, seq, target_seq, start, end, strand):
        df = pd.DataFrame({'seq':list(seq), 'target_seq':list(target_seq)})
        ok = df[['seq', 'target_seq']].isin(['A', 'C', 'G', 'T'])
        # rows where ALL species are valid
        valid_rows = ok.all(axis=1)
        df =df.loc[valid_rows, :]
        df = df[df['seq']!=df['target_seq']]
        if strand == '+':
            pos = [start + x for x in list(df.index)]
        if strand == '-':
            pos = [end - x for x in list(df.index)]
        return pos
    
    
    #chatgpt sugestion, not convinced
    def find_subs_intron_array(self, seq, target_seq, start, end, strand):
        # Convert sequences to numpy arrays of single characters
        seq_arr = np.frombuffer(seq.encode("ascii"), dtype="S1")
        tgt_arr = np.frombuffer(target_seq.encode("ascii"), dtype="S1")

        # Mask: only valid bases AND different
        valid = np.isin(seq_arr, [b"A", b"C", b"G", b"T"]) & np.isin(tgt_arr, [b"A", b"C", b"G", b"T"])
        mask = valid & (seq_arr != tgt_arr)

        idx = np.nonzero(mask)[0]

        if strand == '+':
            return (start + idx).tolist()
        else:  # strand == '-'
            return (end - idx).tolist()

        '''
        if len(df)>1:
            return df['positions'].values, df['target_seq'].values
        elif len(df)==1:
            return [df['positions']], df['target_seq']
        else:
            return [],[]
        '''

    #editting so apply merge to seq and target seq
    def merge_segments_of_exons(self,exon_df):
        if len(exon_df) == 1:
            return exon_df.iloc[[0]]
        else:
            last_row = None
            for idx, row in exon_df.iterrows():
                if last_row is None:
                    last_row = row.copy()
                else:
                    if (row == last_row).all():  #if have overlapping exons
                        pass
                    else:
                        diff = row['start'] - last_row['end']
                        if diff >= 0:
                            new_row = row.copy()
                            new_row['seq'] = last_row['seq'] + 'K' * diff + row['seq']
                            new_row['target_seq'] = last_row['target_seq'] + 'K' * diff + row['target_seq']
                            new_row['start'] = last_row['start']
                            new_row['end'] = row['end']
                            last_row = new_row
                        
                        #think will never happen, as have sorted coords
                        elif row['start'] <= last_row['start'] and row['end'] >= last_row['end']:
                            last_row = row.copy() #new row fully contains old row so replace, not merge
                            
                        
                        elif row['start'] >= last_row['start'] and row['end'] <= last_row['end']:
                            pass  #keep with last row, as fully contains new row
                        #overlapping, eg due to exon and stop codon adding flanks
                        elif diff < 0: # then trim first row, and merge
                            last_row_len = row['start'] - last_row['start']
                            new_row=row.copy()
                            new_row['seq'] = last_row['seq'][:last_row_len] + row['seq']
                            new_row['target_seq'] = last_row['target_seq'][:last_row_len] + row['target_seq']
                            new_row['start'] = last_row['start']
                            last_row = new_row
                         
                        else:
                            raise ValueError(f"Overlapping segments in group {exon_df.name}")
                    
            return pd.DataFrame([last_row])

    #added to pad second sequence also
    def fill_gaps_in_exons(self, df):
        pad_start = (df['start'] - df['exon_start']).astype(int)
        pad_end = (df['exon_end'] - df['end']).astype(int)
        if (pad_start < 0).any() or (pad_end < 0).any():
            print(df[(df['start'] - df['exon_start'])<0])
            print(df[(df['exon_end'] - df['end'])<0])
            raise ValueError("padding problem: negative lengths found")
        
        df['seq'] = pad_start.apply(lambda x: 'K' * x) + df['seq'] + pad_end.apply(lambda x: 'K' * x)
        df['target_seq'] = pad_start.apply(lambda x: 'K' * x) + df['target_seq'] + pad_end.apply(lambda x: 'K' * x)
        #df['trinucs'] = pad_start.apply(lambda x: [-1] * x) + df['trinucs'] + pad_end.apply(lambda x: [-1] * x)
        df.drop(columns=['start', 'end'], inplace=True)
        return df
    
    def apply_reverse_complement_seq_target(self, df):
        df['seq'] = df.apply(
        lambda row: 
            self.reverse_complement(row['seq']) if row['strand'] == '-' else row['seq'],
        axis=1)
        df['target_seq'] = df.apply(
        lambda row: 
            self.reverse_complement(row['target_seq']) if row['strand'] == '-' else row['target_seq'],
        axis=1)
        return df

    #editting to also add positions for extra exons
    def add_extra_exons(self, df, coord_df):
        merged = coord_df.merge(df[['chr', 'exon_start', 'exon_end', 'gene', 'strand']],
                                on=['chr', 'exon_start', 'exon_end', 'gene', 'strand'],
                                how='left', indicator=True)
        missing = merged[merged['_merge'] == 'left_only'].copy()
        if missing.empty:
            return df
        else:
            missing['seq'] = missing.apply(lambda row: 'K' * (row['exon_end'] - row['exon_start']), axis=1)
            missing['positions'] = missing.apply(lambda row : range(row['exon_start'], row['exon_end']),axis=1)
            missing['subs'] = missing.apply(lambda row: [], axis=1)
            missing['subs_alts'] = missing.apply(lambda row: [], axis=1)
            missing.drop(columns=['_merge'], inplace=True)
            return pd.concat([df, missing], ignore_index=True)

    #if using same input as target size
    #do not need flanks here is not considering trinuc contexts
    def remove_flanks(self,df):
        df['start'] +=1
        df['end'] -= 1
        df['exon_start'] +=1
        df['exon_end'] -= 1
        df['seq'] = df['seq'].apply(lambda x: x[1:-1]) #remove flanks 
        df['target_seq'] = df['target_seq'].apply(lambda x: x[1:-1]) #remove flanks 
        return df

        #need to keep track of poistion in gene vs introns 
    def extract_gene_seq_positions(self, gene_df):
        #add positions
       
        strand = gene_df.iloc[0]['strand']
        if strand == '+':
            gene_df = gene_df.sort_values(by='exon_start').reset_index(drop=True)

        else:
            gene_df = gene_df.sort_values(by='exon_end', ascending=False).reset_index(drop=True)

        seq = ''.join(gene_df['seq'].tolist())
        positions = list(chain.from_iterable(gene_df['positions']))
        subs = list(chain.from_iterable(gene_df['subs']))
        subs_alts = list(chain.from_iterable(gene_df['subs_alts']))
        return seq,  positions, subs, subs_alts

    #note for introns do not need to iterate by codon, as all muts are synon
    def generate_intron_mutations1(self, seq, start, end, strand):
        for i in range(0, len(seq) - 2):
            ref_base = seq[i]
            if ref_base not in self.bases:
                continue
            self.intron_poss_subs.append(3*[i+start] if strand == '+' else 3*[end-i])
        return
    
    
    #or work with numpy vectorised version
    def generate_intron_mutations(self, seq, start, end, strand):
        seq_arr = np.array(list(seq))  # convert string to array of characters

        # Boolean mask for valid bases
        valid_mask = np.isin(seq_arr, self.bases)
        idx = np.nonzero(valid_mask)[0]

        # Compute positions based on strand
        if strand == '+':
            positions = start + idx
        else:
            positions = end - idx

    # Each base has 3 possible mutations
        self.intron_poss_subs.extend(np.repeat(positions, 3).tolist())

        
    #overriding function in mutation matrix generator
    def generate_mutations(self, seq,positions):
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon not in self.codon_table:
                continue
            for m in range(3):
                seq_pos = i + m
                if seq_pos >= len(positions):
                    continue
                ref_base = seq[seq_pos]
                for alt in self.bases:
                    if alt == ref_base:
                        continue
                    alt_codon = list(codon)
                    alt_codon[m] = alt
                    alt_codon = ''.join(alt_codon)
                    mut_type = self.find_mutation_type(codon, alt_codon)
                    if mut_type == 0:
                        self.gene_records.append((positions[seq_pos],  alt))
                                
        return 

    
    def run(self):
        self.load_data()
        print('finished preprocessing exons')
        self.load_intron_data()
        print('finished preprocessing introns')
        #self.output_dir.mkdir(parents=True, exist_ok=True)
        l_s = 0
        l_s_start = 0
        l_s_end = 0
        bin_s_obs = 0
        s_obs_per_bin = []
        bin_size = 20000
        self.all_exons_df.sort_values(by =['chr','exon_start'], inplace =True)
        for (gene, strand, chr), gene_exon_df in self.all_exons_df.groupby(["gene","strand","chr"],sort=False):
            
            if gene not in self.intron_df_group.groups:
                continue
            self.gene_records = []
            seq, positions, gene_exon_subs, gene_exon_subs_alts = self.extract_gene_seq_positions(gene_exon_df)
            if len(seq) % 3 != 0:
                continue
            self.generate_mutations(seq, positions)

            #find synonymous subs
            gene_exon_syn_poss_dicts = {'pos': [r[0] for r in self.gene_records],
            'alt': [r[1] for r in self.gene_records]}
            gene_exon_poss_muts_df = pd.DataFrame(gene_exon_syn_poss_dicts)
            gene_exon_subs_dict = {'pos': gene_exon_subs, 'alt': gene_exon_subs_alts}
            gene_exon_subs_df = pd.DataFrame(gene_exon_subs_dict)
            gene_syn_exon_subs_df = gene_exon_subs_df.merge(gene_exon_poss_muts_df, on=['pos','alt'], how = 'inner')
            
            #extract intron all possible muts
            intron_gene_df = self.intron_df_group.get_group(gene)
            self.intron_poss_subs = []
            intron_gene_df.apply(lambda row : self.generate_intron_mutations(row['seq'], row['start'], row['end'], row['strand']), axis=1)
            
            #extra list of all pos syn muts in gene (from generate muts in intron and exon)
            gene_poss_muts = self.intron_poss_subs + list(gene_exon_poss_muts_df['pos']) #keep repeats-means 4 fold syn
            #print(f'possible syn muts in gene {len(gene_poss_muts)}')

            #extron intron subs (all 'syn' so do not need to merge with syn)
            gene_intron_subs = list(chain.from_iterable(intron_gene_df['subs']))
            
            #list of positions of all 'syn' subs that happened in gene 
            gene_syn_subs = gene_intron_subs + list(gene_syn_exon_subs_df['pos'])
            #print(f'syn subs in gene {len(gene_syn_subs)}')

            #sort order of poss muts
            '''
            if strand == '+':
                gene_poss_muts.sort()
            if strand == '-':
                gene_poss_muts.sort(reverse=True)
            '''
            gene_poss_muts.sort()

            #if start of bin is start of gene, then record start pos
            if l_s ==0:
                l_s_start = gene_poss_muts[0]

            #bin gene
            gene_l_s = len(gene_poss_muts)
            if l_s + gene_l_s < bin_size:
                l_s += gene_l_s
                bin_s_obs += len(gene_syn_subs)
            elif l_s + gene_l_s == bin_size:
                l_s += gene_l_s
                bin_s_obs += len(gene_syn_subs)
                l_s_end = max(gene_poss_muts)
                s_obs_per_bin.append([chr, l_s_start, l_s_end, bin_s_obs])
                l_s =0
                bin_s_obs =0
                
            else:
                while l_s + gene_l_s > bin_size:
                    #first bin from start of gene_df (or remaining gene_df)
                    bin1_gene_l_s = bin_size - l_s
                    l_s_end = gene_poss_muts[bin1_gene_l_s-1] #not sure about -1
                    #choose how to deal with mut at cut off positions (as 3 possibble muts not ordered)
                    #find s_obs for part of gene in bin1, start fog ene depends on strand!
                    '''
                    if strand == '+':
                        bin1_syn_subs = [x for x in gene_syn_subs if x<=l_s_end] 
                    elif strand == '-':
                        bin1_syn_subs = [x for x in gene_syn_subs if x>=l_s_end]
                    '''
                    bin1_syn_subs = [x for x in gene_syn_subs if x<l_s_end]
                    bin_s_obs +=  len(bin1_syn_subs)
                    s_obs_per_bin.append([chr, l_s_start, l_s_end, bin_s_obs]) 
                    
                    #reset subs and poss muts, removing binned section
                    gene_poss_muts = gene_poss_muts[bin1_gene_l_s:] #keep unbinned muts
                    gene_syn_subs = [x for x in gene_syn_subs if x not in bin1_syn_subs]
                    l_s_start = gene_poss_muts[0]
                    l_s = 0
                    gene_l_s = len(gene_poss_muts)
                    bin_s_obs =0
                 #keep making more bins in gene ie repeat loop above

                
                #now reset and calculate for last bin in gene
                bin_s_obs = len(gene_syn_subs) 
                l_s =  gene_l_s
                #CHECK WETHER NEED ANOTHER BIN  IN GENE, keep doing loop
            #print(s_obs_per_bin)
        df = pd.DataFrame({'chr':[x[0] for x in s_obs_per_bin],
                     'bin_start':[x[1] for x in s_obs_per_bin],
                     'bin_end':[x[2] for x in s_obs_per_bin],
                     's_obs':[x[3] for x in s_obs_per_bin]})
        print(df)
        df.to_csv(self.output_dir,index=False)
        
  
species = 'Anc4'
target_species = 'hg38'
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df2.bed'
intron_file = '/home/maria/synon_mut_rates/auxillary/introns_genes.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = '/home/maria/synon_mut_rates/auxillary/neutral_muts'

generator = NeutralGenerator(species, target_species, possible_species, exon_file, gene_annotations, intron_file,output_dir)
generator.run()


#test function 
'''
import pandas as pd
exon_df = pd.DataFrame({'exon_start':[1,1,1,10,10], 'exon_end': [8,8,8,18, 18],'start':[1,5, 6,10,15], 'end':[5,7,8,12, 17], 'hg38':['AAAA', 'GG','CC','TT','G'], 'other':[1,1,1,3,3], 'chr':[1,1,1,3,3]})
exon_df.to_csv('/home/maria/test.bed')
exon_annot_df = pd.DataFrame({'exon_start':[1,1,1,10,10, 20], 'exon_end': [8,8,8,18, 18, 25], 'gene':[1,1,1,3,3,3], 'chr':[1,1,1,3,3,3]})
exon_annot_df.to_csv('/home/maria/test_a.bed')

species = 'hg38'
possible_species = ['hg38']
exon_file = '/home/maria/test.bed'
gene_annotations = '/home/maria/test_a.bed'
output_dir = '/home/maria/non_syn_target_test'

generator = MutationMatrixGenerator(species, possible_species, exon_file, gene_annotations, output_dir)
generator.run()

'''

#cp /home/maria/cactus_target_size/scripts/MutationMatrixGenerator.py /home/maria/synon_mut_rates/scripts/MutationMatrixGenerator.py

#mut rate 10^-8 per gen
#generations = 6x10^6 / 30 = 2x10^5
#mut rate per base = 2x10^-3

#E[mut per bin] = 40