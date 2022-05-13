import altair as alt
import pandas as pd
import numpy as np
import sqlite3

alt.data_transformers.disable_max_rows()

from .dataframes import CoronaDataframe

def gene_loci(df: pd.DataFrame) -> list:
    '''Helper function for adding gene information'''
    genes = ['E' if 26245 <= x <= 26472 else 
             'M' if 26523 <= x <= 27191 else
             'N' if 28274 <= x <= 29533 else
             'ORF1a' if 266 <= x <= 13468 else
             'ORF1b' if 13468 <= x <= 21555 else
             'ORF3a' if 25393 <= x <= 26220 else
             'ORF6' if 27202 <= x <= 27387 else
             'ORF7a' if 27394 <= x <= 27759 else
             'ORF7b' if 27756 <= x <= 27887 else
             'ORF8' if 27894 <= x <= 28259 else
             'ORF9b' if 28284 <= x <= 28577 else
             'S' if 21563 <= x <= 25384 else
             'None' for x in df.position]
    
    return genes

class CoronaPlot:
    '''Class for plotting corona data. Uses Altair.'''
    
    def __init__(self, query: str, ref: str, mutations: str) -> None:
        self.dataframes = CoronaDataframe(query, ref, mutations)
        
    def plot_alignment_map(self):
        ''' Returns chart and view which plots the alignment map of the query sequence'''
        ####### Chache with json file????? ######## 
        df = self.dataframes.make_alignment_df()
        df['gene'] = gene_loci(df)
        
        interval = alt.selection_interval(encodings=['x'])

        base = alt.Chart(df, title="Coverage and mutations").mark_rect().encode(
            x=alt.X('position:N', axis=alt.Axis(labels=False, ticks=False)),
            y='id:N',
            color=alt.Color('max(value):N', legend=None, scale=alt.Scale(domain=[0,-1, -2, 1], range=['white', 'red', 'black', 'green'])),
            tooltip=[alt.Tooltip('nt:N', title='nt'), 
                     alt.Tooltip('position:N', title='Position'), 
                     alt.Tooltip('ref:N', title='Reference'),
                     alt.Tooltip('gene:N', title='Gene')])

        chart = base.encode(x=alt.X('position:N', scale=alt.Scale(domain=interval.ref()), axis=alt.Axis(labels=False, ticks=False))
            ).properties(width=1200, height=300)

        view = base.add_selection(interval
        ).properties(width=1200, height=100)

        return chart, view
        
    def plot_mutation_heatmap(self):
        '''Plots a heatmap over all mutations in the mutation csv file'''
        mutations = self.dataframes.make_mutations_df()

        unique_mut = list(mutations.mutation.unique())
        order_level = sorted(unique_mut, key=lambda x: int("".join([i for i in x if i.isdigit()])))

        fig = alt.Chart(mutations, 
                        title='Mutations. Dark blue == deletions. Light blue == unique mutations').mark_rect(stroke='black').encode(
            x=alt.X('mutation', axis=alt.Axis(labels=False, ticks=False), sort=order_level),
            y='pango',
            color=alt.Color('kind', legend=None),
            tooltip=[
                alt.Tooltip('mutation', title='Mutation'),
                alt.Tooltip('pango', title='Pango')]).properties(width=1200)

        return fig
    
    def plot_compare(self, compare_to: str):
        '''Compares query mutatations against mutations for a given corona lineage.'''
        mutations_db = self.dataframes.make_mutations_df()
        mutations_db = mutations_db[['pango', 'mutation', 'kind', 'position']]
        lineage_df = mutations_db[mutations_db['pango'] == compare_to]

        query_df = self.dataframes.make_query_df()
        # Vad är bäst, jämföra mot bara mutationer i den man vill jämföra med eller ta alla mutationer som finns i query???
        # kommentera bort raderna nedan för att testa
        #query_df = query_df[query_df['position'].isin(mutations_db.position.to_list())]
        query_df = query_df[query_df['position'].isin(lineage_df.position.to_list())]

        query_df['kind'] = [3 if x in lineage_df.mutation.to_list() else 4 for x in query_df.mutation]

        concated = pd.concat([lineage_df, query_df])
        concated['gene'] = gene_loci(concated)
        
        # Order of x axis
        position_level = sorted(concated.position.unique())

        fig = alt.Chart(concated, title={'text': [f'Compared against {compare_to}'], 
                                         'subtitle': ['Green == same, pink == unique, blue == deletion, red == differs from comparison'], 
                                         'color': 'green'}).mark_rect(stroke='black').encode(
        x=alt.X('position:N', axis=alt.Axis(), sort=position_level),
        y='pango',
        color=alt.Color('kind', legend=None, 
                        scale=alt.Scale(domain=[0, 1, 2, 3, 4], range=['#e7fc98 ', '#f64fe7', '#98fcf3', '#e7fc98', '#fc9898'])),
        tooltip=[
            alt.Tooltip('mutation', title='Mutation'),
            alt.Tooltip('pango', title='Pango'),
            alt.Tooltip('gene', title='Gene')])#.properties(width=1000, height=50)

        return fig

    def plot_mutations_inspection(self):
        '''Plots how common each mutation is for every mutation in the query sequence'''
        df = self.dataframes.make_mutations_inspection_df()
        df['gene'] = gene_loci(df)
        
        fig = alt.Chart(df, title={'text': ['Mutations in'], 
                                         'subtitle': ['Green == regular, blue == deletion, red == N', 
                                                      'Size shows how common mutation is'], 
                                         'color': 'green'}).mark_point(stroke='black').encode(
            y=alt.Y('position:N', axis=alt.Axis()),
            x='pango',
            fill=alt.Fill('color:N', legend=None, scale=alt.Scale(domain=[1, 2, 3], range=['red', 'blue', 'green'])),
            size=alt.Size('pango_with_mutation', legend=None, scale=alt.Scale(range=[100, 1500])),
            tooltip=[
                alt.Tooltip('mutation', title='Mutation'),
                alt.Tooltip('gene', title='Gene'),
                alt.Tooltip('mutation_list', title='Pangos with mutation')])#.properties(width=1200, height=150)

        return fig
    
    
    def plot_gene_annotation(self):
        '''Plots gene annotation map'''
        mutations = pd.read_csv(self.dataframes.mutations)
        mutations['position'] = mutations['mutation'].str.extract('(\d+)').astype(int)
        mutations['gene'] = gene_loci(mutations)
        
        genes = mutations.groupby('gene').agg({'position': [np.min, np.max]}).reset_index()
        genes.columns = ['gene', 'start', 'end']
        genes = genes[lambda x: x.gene != 'Not annotated']
        genes = genes.sort_values('end')
        genes['len'] = genes['end'] - genes['start']
        genes = genes[lambda x: x.gene != 'None']

        fig = alt.Chart(genes).mark_bar().encode(
         alt.X('start:Q', axis=alt.Axis(title='Position')),
         alt.X2('end:Q'),
         alt.Fill('gene:N', sort=genes['gene'].to_list(), scale=alt.Scale(scheme='paired'),
                  legend=alt.Legend(orient='bottom', title=None)),
         alt.Text('gene'), 
         tooltip=[alt.Tooltip('gene', title='Gene'),
                  alt.Tooltip('start', title='Start'),
                  alt.Tooltip('end', title='End'),
                  alt.Tooltip('len', title='Length (nt)')]
        ).properties(width=1200)
        
        return fig
    
    def plot_geographical_similar(self, pango: str, connection: sqlite3.Connection):
        '''Plots the 50 most similar sequences to the query and the location of those'''
        df = self.dataframes.find_similar_sequences(pango, connection)
        df.drop(columns=['mutation', 'similarity'], inplace=True)
        
        fig = alt.Chart(df).mark_bar().encode(
         alt.X('date:T'),
         alt.Y('count()', axis=alt.Axis(title='Number of similar sequences')),
         alt.Color('percent_similar:Q', legend=alt.Legend(format='%')),
         tooltip=[alt.Tooltip('name:N'),
                  alt.Tooltip('place:N'),
                  alt.Tooltip('percent_similar:Q'),
                  alt.Tooltip('number_similar:Q'),
                  alt.Tooltip('date:T')]
        ).facet(
         'country:N', columns=1)
        
        return fig