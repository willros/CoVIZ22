import streamlit as st
import altair as alt
from corona_viz.plotting import CoronaPlot
from corona_viz.nextclade import Nextclade

# Cluster the heatmap in another way

def app():
    nextclade = Nextclade()
    corona = CoronaPlot(nextclade.QUERY, nextclade.REF, 'mutations.csv')
    
    st.markdown('### Heatmap of mutations in all lineages')
    heatmap = corona.plot_mutation_heatmap()
    st.altair_chart(heatmap)
    
