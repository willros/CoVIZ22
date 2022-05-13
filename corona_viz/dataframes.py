import pandas as pd
from Bio import SeqIO
import sqlite3

class CoronaDataframe:
    '''Class that handles manipulation and generation of dataframes
    needed for plotting'''
    
    def __init__(self, query: str, ref: str, mutations: str) -> None:
        
        # Only include region between 256:29674 (same as pangolin)
        self.query = SeqIO.read(query, 'fasta').seq[256:29674]
        self.ref = SeqIO.read(ref, 'fasta').seq[256:29674]
        self.mutations = mutations
    
    
    def make_alignment_df(self) -> pd.DataFrame:
        '''Turns the sequence of the aligned fasta file into a df used to plot 
        the missing (N), deletions and mutations as a genome map'''
    
        reference = [x for x in self.ref]
        query = [x for x in self.query]
        converted_values = self._convert_alignment_value()

        alignment_df = pd.DataFrame({'nt': query, 'id': 'query', 
                                 'position': range(256, len(query) + 256), 
                                 'value': converted_values, 
                                 'ref': reference})
        return alignment_df
    
    
    def make_mutations_df(self) -> pd.DataFrame:
        '''Transforms the df of mutations to a suitable format'''
        
        mutations_df = pd.read_csv(self.mutations)
        mutations_df['unique'] = [1 if x == 1 else 0 for x in mutations_df.number_mutations]
        mutations_df['kind'] = [2 if '-' in str(x) else y for x,y in zip(mutations_df.mutation, 
                                                                         mutations_df.unique)]
        mutations_df['position'] = mutations_df['mutation'].str.extract('(\d+)').astype(int)

        return mutations_df
    
    def make_query_df(self) -> pd.DataFrame:
        '''Turns a list of mutations for the query into a df in the right format'''
        
        mutations_list = self._extract_snp()
        query_df = pd.DataFrame({'pango': 'QUERY', 'mutation': mutations_list})
        query_df['position'] = query_df['mutation'].str.extract('(\d+)').astype(int)
    
        return query_df
    
    def make_mutations_inspection_df(self) -> pd.DataFrame:
        '''Returns a dataframe with all mutations and lineages 
        for every given position in the query mutation set'''
        mutations = self.make_mutations_df()
        query_df = self.make_query_df()
        mutations = mutations[mutations['position'].isin(query_df.position)]
        query_df = query_df[query_df['position'].isin(mutations.position)]

        mutations_list = mutations.groupby('position')[['mutation', 'pango']].aggregate(list).reset_index()
        mutations_list['mutation_list'] = [list(zip(x, y)) for x,y in zip(mutations_list.pango, mutations_list.mutation)]
        mutations_list.drop(columns=['pango', 'mutation'], inplace=True)

        merged = query_df.merge(mutations_list)
        merged['pango_with_mutation'] = merged['mutation_list'].apply(len)
        merged['color'] = [1 if x.endswith('N') else 2 if x.endswith('-') else 3 for x in merged.mutation]

        return merged
    
    def _convert_alignment_value(self) -> list:
        '''Converts alignments to numbers depending on it is N, 
        deletion, substitution or the same'''
        values = []
        for q,r in zip(self.query, self.ref):
            if q != r:
                if q == 'N':
                    values.append(-1)
                elif q == '-':
                    values.append(-2)
                else:
                    values.append(1)
            else:
                values.append(0)  
        return values
    
    def _pull_out_pango_from_db(self, pango: str, connection: sqlite3.Connection) -> pd.DataFrame:
        '''Pulls out 10000 rows of data about a given pango from the database'''
        df = pd.read_sql_query(f'''
         SELECT *
         FROM mutations
         WHERE pango = "{pango}"
         ORDER BY date DESC 
         limit 10000;
         ''',
         connection)

        df['mutation'] = [set(x.split(',')) for x in df['mutation']]
        df['date'] = df['date'].astype('datetime64[ns]')

        return df
    
    def find_similar_sequences(self, pango: str, connection: sqlite3.Connection) -> pd.DataFrame:
        '''Searches for similar sequences in the database based on mutations (without N)'''
        mutation_list = self._extract_snp() 
        query = set([x for x in mutation_list if not x.endswith('N')])

        similarity_df = self._pull_out_pango_from_db(pango, connection)
        similarity_df['similarity'] = [query.intersection(x) for x in similarity_df['mutation']]
        similarity_df['number_similar'] = similarity_df['similarity'].apply(len)
        similarity_df['percent_similar'] = similarity_df['number_similar'] / len(query)
        similarity_df.sort_values(['percent_similar', 'date'], ascending=False, inplace=True)

        return similarity_df.head(50)

    
    # Adderat 1 till index 
    def _extract_snp(self) -> list:
        '''Helper function to extract information about mutation'''

        return [f'{a}{index + 1}{b}' for index,(a,b) in 
                enumerate(zip(self.ref, self.query), 256) if a != b]