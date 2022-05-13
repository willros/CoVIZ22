import streamlit as st
import altair as alt
from corona_viz.plotting import CoronaPlot
from corona_viz.nextclade import Nextclade

# Make the gene map above the smaller graph. Use position and a mark_rect with a heigth to plot the genes.

def app():
    nextclade = Nextclade()
    corona = CoronaPlot(nextclade.QUERY, nextclade.REF, 'mutations.csv')
    df = nextclade.make_meta_df()
    name = df.seqName.values[0]
    
    st.markdown(f'### Alignment of {name}')
    a, b = corona.plot_alignment_map()
    annotation_plot = corona.plot_gene_annotation()    
    st.altair_chart((a & b) & annotation_plot)
