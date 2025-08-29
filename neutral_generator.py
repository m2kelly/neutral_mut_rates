
import pandas as pd
from pathlib import Path
from itertools import chain

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
target matrix of non syn mutations from trinuc to base, per gene
note never transform coordinates when reverse compliemnt within the intron 

Note - main difference is adding tracking of positions of all possible mutaions
#all added positions are with respect to the POSITIVE STRAND
# so when on - strand start counting fromm last positions  

PIPELINE-EXONS
1) load data for seq (ancestral) and target seq
2) merge segments of an exon if multiple exist in exon file (for seq &target_seq)
3) fill ends of exons, so correct lengths (with flanks still, for seq&target_seq)
4) reverse complement whole exons on negative strand (seq, target_seq)
5) add list of positions of all bases (so can order with introns)
Note - main difference is adding tracking of positions of all possible mutaions
#all added positions are with respect to the POSITIVE STRAND
# so when on - strand start counting fromm last positions  
6) for each exon record list of positions of each subsitions between seq & target seq
7) drop target seq-from now only work with seq
5) calculate trinucleotide contexts of full exons
6) remove flanks from exons
7) add any extra exons
8) reconstruct gene seqs (and gene trinucs and positions) by merging exons in genes 
per gene: 
9) iterate through codons and identify if all possible muts are syn or non-syn
10) for syn target extract correpsonding position and trinuc and add to gene table
 
#PIPELINE-INTRONS
add positions
rc seq, target seq and reverse order of positions for - strand
for each intron record list of positions of each subsitions between seq & target seq
discard target seq
per gene:
1) iterate through codons and identify all possible muts
2) save position and trinuc to df 

combined per gene:
calculate l_s = intron and exon syn muts per gene, 
if l_s < bin size: s_obs += subs and process next gene , 
else find pos where l_s = bin_size, and add subs before this to s_obs then reset counters
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
        self.reconstructed_df = self.reconstructed_df[:1000]

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

        
        #input segments of exons, output full exons and corresponding trinucs 
        exon_merged_df = self.apply_segment_merges()
        full_exons_df = self.fill_gaps_in_exons(exon_merged_df)
        full_exons_df['positions'] = full_exons_df.apply(lambda row : range(row['exon_start'], row['exon_end']), axis=1)
        #reverse comp seqs
        rc_df = self.apply_reverse_complement_to_exons(full_exons_df)
        rc_df['subs'] = rc_df.apply(lambda row : self.find_subs(row['seq'], row['target_seq']))
        rc_df.drop(columns='target')
        #here check for mutations and add column for target mut
        #then remove seq_target
        #reverse positions
        rc_df['positions'] = rc_df.apply(lambda row: row['positions'][::-1] if row['strand'] == '-' else row['positions'],axis=1)
        #subs columns contains lists of positions of subs
        #think no point in keeping alts? as just want to keep s obs
        #need to make function find_subs
        rc_df['subs'] = rc_df.apply(lambda row : self.find_subs(row['seq'], row['target_seq'], row['positions']))
        rc_df.drop(columns='target_seq', inplace=True)
        trinuc_df = self.add_trinucs_remove_flanks(rc_df) #also remove flanks from mutations
        self.all_exons_df = self.add_extra_exons(trinuc_df, cut_coord_df) #also fill mutations column
       
    def load_intron_data(self):
        header = ['chr', 'start', 'end'] + self.possible_species + ['gene', 'strand']
        intron_df = pd.read_csv(self.intron_file, sep="\t", names=header)
        for speci in self.possible_species:
            if speci == self.species:
                intron_df.rename(columns={speci : 'seq'}, inplace=True)
            elif speci == self.target_species:
                intron_df.rename(columns={speci : 'target_seq'}, inplace=True)
            else:
                intron_df.drop(columns=[speci], inplace=True)
        intron_df['positions'] = intron_df.apply(lambda row : range(row['start'], row['end']),axis=1)
        print(intron_df.columns)
        rc_intron_df = self.apply_reverse_complement_to_exons(intron_df)
        rc_intron_df['positions'] = rc_intron_df.apply(lambda row: row['positions'][::-1] if row['strand'] == '-' else row['positions'],axis=1)
        rc_intron_df['subs'] = rc_intron_df.apply(lambda row : self.find_subs(row['seq'], row['target_seq'], row['positions']))
        rc_intron_df.drop(columns='target_seq', inplace=True)
        self.intron_df_group = rc_intron_df.groupby('gene')


    def find_subs(self, seq, target_seq, positions):
        df = pd.DataFrame({'seq':list(seq), 'target_seq':list(target_seq), 'positions':positions})
        ok = df[['seq', 'target_seq']].isin(['A', 'C', 'G', 'T'])
        # rows where ALL species are valid
        valid_rows = ok.all(axis=1)
        df =df.loc[valid_rows, :]
        print(df)
        df = df[df['seq']!=df['target_seq']]
        print(df)
        if len(df)>1:
            return df['positions'].values
        elif len(df)==1:
            return [df['positions'][1]]
        else:
            return []

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
            missing['trinucs'] = missing.apply(lambda row: [-1] * (row['exon_end'] - row['exon_start']), axis=1)
            missing['positions'] = missing.apply(lambda row : range(row['exon_start'], row['exon_end']),axis=1)
            missing['subs'] = missing.apply(lambda row: [0] * (row['exon_end'] - row['exon_start']), axis=1)
            
            missing.drop(columns=['_merge'], inplace=True)
            return pd.concat([df, missing], ignore_index=True)

    #also editting to remove flanks from positions
    def add_trinucs_remove_flanks(self, df):
        df['trinucs'] = df['seq'].apply(lambda x: self.list_trinucs_of_flanked_seq(x))
        df['seq'] = df['seq'].apply(lambda x: x[1:-1]) #remove flanks 
        df['positions'] = df['positions'].apply(lambda x: x[1:-1])
        df['subs'] = df['subs'].apply(lambda x: x[1:-1])
        df['exon_start'] +=1
        df['exon_end'] -= 1
        return df


        #need to keep track of poistion in gene vs introns 
    def extract_gene_seq_trinucs(self, gene_df):
        #add positions
       
        strand = gene_df.iloc[0]['strand']
        if strand == '+':
            gene_df = gene_df.sort_values(by='exon_start').reset_index(drop=True)

        else:
            gene_df = gene_df.sort_values(by='exon_end', ascending=False).reset_index(drop=True)

        seq = ''.join(gene_df['seq'].tolist())
        trinucs = list(chain.from_iterable(gene_df['trinucs']))
        positions = list(chain.from_iterable(gene_df['positions']))
        subs = list(chain.from_iterable(gene_df['subs']))
        return seq, trinucs, positions, subs

    def generate_intron_mutations(self, seq, positions):
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            
            if codon not in self.codon_table:
                continue
            for m in range(3):
                seq_pos = i + m
                ref_base = seq[seq_pos]
                context = seq[seq_pos-1:seq_pos+2] 
                if context not in self.codon_table:
                    continue
        
                for alt in self.bases:
                    if alt == ref_base:
                        continue
                    alt_codon = list(codon)
                    alt_codon[m] = alt
                    alt_codon = ''.join(alt_codon)
                    self.gene_records.append((positions[seq_pos], context, alt))
                    
                       
        
    #overriding function in mutation matrix generator
    def generate_mutations(self, seq, trinucs, positions):
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon not in self.codon_table:
                continue
            for m in range(3):
                seq_pos = i + m
                if seq_pos >= len(trinucs):
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
                        context_index = trinucs[seq_pos]
                        if context_index != -1:
                            context = self.trinucleotides[context_index]
                            if context[1] == ref_base: #check that trinucs are correct
                                self.gene_records.append((positions[seq_pos], context, alt))
                                
                            else: 
                                print('mismatch', context, ref_base)
        return 

    
    def run(self):
        self.load_data()
        self.load_intron_data()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        l_s = 0
        bin_size = 20000
        for gene, strand, gene_exon_df in self.all_exons_df.groupby(["gene","strand"]):
            if gene not in intron_gene_df.groups:
                continue
            self.gene_records = []
            seq, trinucs, positions, exon_subs = self.extract_gene_seq_trinucs(gene_exon_df)
            if len(seq) % 3 == 0:
                self.generate_mutations(seq, trinucs, positions)
            intron_gene_df = self.intron_df_group.get_group(gene)
            intron_gene_df.apply(lambda row : self.generate_intron_mutations(row['seq'], row['positions']), axis=1)
            subs = exon_subs + intron_gene_df['subs'].values

            gene_dict = {
            'pos': [r[0] for r in self.gene_records],
            'trinuc': [r[1] for r in self.gene_records],
            'alt': [r[2] for r in self.gene_records],}

            gene_poss_muts_df = pd.DataFrame(gene_dict)
            if strand == '+':
                gene_poss_muts_df.sort_values(by = 'positions', inplace=True, ascending =True)
              
            elif strand == '-':
                gene_poss_muts_df.sort_values(by = 'positions', inplace=True, ascending =False)
                
            gene_l_s = len(gene_poss_muts_df)
            print(gene_poss_muts_df)
            if l_s + gene_l_s < bin_size:
                l_s = l_s + gene_l_s
                #use whole df
                #add all mutations from this region to mutation matrix !
                #must do intron allignemnt and exon allignments on both species at once !
                bin_s_obs += len(subs)
                #save then set to zero
                bin_s_obs = 0
            else:
                bin1_gene_l_s = bin_size - l_s
                bin1_gene_pos_cutoff = gene_poss_muts_df.loc[bin1_gene_l_s]['positions'][0]
                #choose how to dela with mut at cut off positions (as 3 possibble muts not ordered)
                #add muts also contained in anc here 
                if strand == '+':
                    bin1_subs = [x for x in subs if x < bin1_gene_pos_cutoff]
                elif strand == '-':
                    bin1_subs = [x for x in subs if x > bin1_gene_pos_cutoff]
                bin_s_obs +=  len(bin1_subs)
                #save 

                #now reset 
                bin2_subs = set(subs)- set(bin1_subs)
                bin_s_obs = len(bin2_subs)

                bin2_gene_l_s = gene_l_s - bin1_gene_l_s
                l_s = bin2_gene_l_s
           

species = 'Anc4'
target_species = 'hg38'
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df2.bed'
intron_file = '/home/maria/synon_mut_rates/auxillary/introns_genes.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = '/home/maria/synon_mut_rates/auxillary/neutral_matricies'

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
