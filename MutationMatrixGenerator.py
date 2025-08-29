#edit so rmeove flanks
#edit so dealing correctly with split up segments rom same exon


import pandas as pd
from pathlib import Path
from itertools import chain

#think should redo seq etraction, combining eon and exon end first-
# as causing contiuous extraction in extracted exon seqs (ie exon plus codon), 
# but in seqs not combining + adding headache to merging exons
# so preprocess gene coords, combining stop codons with exon 
# sort df and groupby gene, combine lastrow end = row start , apply to group, like merge func here
'''
DATA
exon_file - species data with reconstructed windows, all seq on + strand, 
bedtools intersected with full exon coords (without dup or cpg removed, but with overlapping gene segments removed)
all seq and coords contain flanks (added +- 1 to coords )
gene_annotations - contains full gene coords (without duplications removed) 
to add any missing exons in genes in the exon_file-without flanks in coords
OUTPUT
target matrix of non syn mutations from trinuc to base, per gene
PIPELINE
1) load data 
2) merge segments of an exon if multiple exist in exon file
3) fill ends of exons, so correct lengths (with flanks still)
4) reverse complement whole exons on negative strand
5) calculate trinucleotide contexts of full exons
6) remove flanks from exons
7) add any extra exons
8) reconstruct gene seqs (and gene trinucs) by merging exons in genes 
9) iterate through codons and identify if all possible muts are syn or non-syn
10) for non syn target extract correpsonding trinuc and add to matrix
11) save target matrix per gene  

'''
class MutationMatrixGenerator:
    def __init__(self, species, possible_species, exon_file, gene_annotations, output_dir):
        self.possible_species = possible_species
        self.species=species
        self.exon_file = exon_file
        self.gene_annotations = gene_annotations
        self.output_dir = Path(output_dir)
        self.bases = ['A', 'C', 'G', 'T']
        self.complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        self.trinucleotides = [a + b + c for a in self.bases for b in self.bases for c in self.bases]
        self.empty_mut_matrix = pd.DataFrame(0, index=self.trinucleotides, columns=self.bases)

        self.codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
        }

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
        rc_df = self.apply_reverse_complement_to_exons(full_exons_df)
        trinuc_df = self.add_trinucs_remove_flanks(rc_df)
        self.all_exons_df = self.add_extra_exons(trinuc_df, cut_coord_df)
       

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
                            new_row['start'] = last_row['start']
                            last_row = new_row
                         
                        else:
                            raise ValueError(f"Overlapping segments in group {exon_df.name}")
                    
            return pd.DataFrame([last_row])

    #try grouping, by exon start and chr, mostly will return fine, but if not need to merge
    def apply_segment_merges(self):
        # choose the grouping columns that uniquely define an exon
        group_cols = ['chr', 'exon_start', 'gene']  # consider ['chr','gene','strand','exon_start','exon_end'] if needed

        df = self.reconstructed_df.copy()
        df = df.sort_values(group_cols + ['start'])

        merged = (
            df.groupby(group_cols, sort=False, group_keys=False)
            .apply(self.merge_segments_of_exons)
            .reset_index(drop=True)
        )
        return merged



    def add_trinucs_remove_flanks(self, df):
        df['trinucs'] = df['seq'].apply(lambda x: self.list_trinucs_of_flanked_seq(x))
        df['seq'] = df['seq'].apply(lambda x: x[1:-1]) #remove flanks 
        df['exon_start'] +=1
        df['exon_end'] -= 1
        return df

    def list_trinucs_of_flanked_seq(self, seq):
        return [
            self.trinucleotides.index(seq[i-1:i+2]) if seq[i-1:i+2] in self.trinucleotides else -1
            for i in range(1, len(seq)-1)
        ]


#assuming all seqs on posiiv strand at this point
    def fill_gaps_in_exons(self, df):
        pad_start = (df['start'] - df['exon_start']).astype(int)
        pad_end = (df['exon_end'] - df['end']).astype(int)
        if (pad_start < 0).any() or (pad_end < 0).any():
            print(df[(df['start'] - df['exon_start'])<0])
            print(df[(df['exon_end'] - df['end'])<0])
            raise ValueError("padding problem: negative lengths found")
        
        df['seq'] = pad_start.apply(lambda x: 'K' * x) + df['seq'] + pad_end.apply(lambda x: 'K' * x)
        #df['trinucs'] = pad_start.apply(lambda x: [-1] * x) + df['trinucs'] + pad_end.apply(lambda x: [-1] * x)
        df.drop(columns=['start', 'end'], inplace=True)
        return df
    
    def reverse_complement(self, seq):
        bases = list(seq)
        comp_bases=[]
        for base in bases:
            if base in self.complement.keys():
                comp_base = self.complement[base]
            else:
                comp_base = base
            comp_bases.append(comp_base)
        comp_seq = ''.join(comp_bases)
        return comp_seq[::-1]
    
    def apply_reverse_complement_to_exons(self, df):
        df['seq'] = df.apply(
    lambda row: 
        self.reverse_complement(row['seq']) if row['strand'] == '-' else row['seq'],
    axis=1)
        return df
        
    #assuming exon_star and exon_coords do not contain flanks
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
            missing.drop(columns=['_merge'], inplace=True)
            return pd.concat([df, missing], ignore_index=True)
       

    def extract_gene_seq_trinucs(self, gene_df):
        strand = gene_df.iloc[0]['strand']
        if strand == '+':
            gene_df = gene_df.sort_values(by='exon_start').reset_index(drop=True)
        else:
            gene_df = gene_df.sort_values(by='exon_end', ascending=False).reset_index(drop=True)
        seq = ''.join(gene_df['seq'].tolist())
        trinucs = list(chain.from_iterable(gene_df['trinucs']))
        return seq, trinucs

    def find_mutation_type(self, ref, alt):
        if ref in self.codon_table and alt in self.codon_table:
            return 0 if self.codon_table[ref] == self.codon_table[alt] else 1
        return -1

    def generate_mutations(self, seq, trinucs, type):
        matrix = self.empty_mut_matrix.copy()
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
                    if mut_type == type:
                        context_index = trinucs[seq_pos]
                        if context_index != -1:
                            context = self.trinucleotides[context_index]
                            matrix.at[context, alt] += 1
        return matrix

    def run(self):
        self.load_data()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for gene, gene_df in self.all_exons_df.groupby("gene"):
            seq, trinucs = self.extract_gene_seq_trinucs(gene_df)
            if len(seq) % 3 == 0:
                matrix = self.generate_mutations(seq, trinucs, 1) #1 is for non syn
                matrix.to_csv(self.output_dir / f"{gene}")
              
               

'''
species = 'Anc4'
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df2_nocpg.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = '/home/maria/cactus_target_size/auxillary/non_syn_target_nocpg_anc'

generator = MutationMatrixGenerator(species, possible_species, exon_file, gene_annotations, output_dir)
exon_df = generator.load_data()
genes = set(exon_df['gene'])
print(len(genes))

'''
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
#tets
'''
species = 'hg38'
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df2.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = '/home/maria/cactus_target_size/auxillary/non_syn_target_hg38'

generator = MutationMatrixGenerator(species, possible_species, exon_file, gene_annotations, output_dir)
generator.run()
'''