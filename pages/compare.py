import streamlit as st
import altair as alt
from corona_viz.plotting import CoronaPlot
from corona_viz.nextclade import Nextclade

def app():
    nextclade = Nextclade()
    corona = CoronaPlot(nextclade.QUERY, nextclade.REF, 'mutations.csv')
    df = nextclade.make_meta_df()
    pango_line = df.Nextclade_pango.values[0]
    clade = df.clade.values[0]
    name = df.seqName.values[0]

    st.markdown(f'##### Compare {name} to other lineages')
    st.markdown(f'**Predicted pangolin and clade:** {pango_line} | {clade}')
    compare_to = corona.dataframes.make_mutations_df()
    pango = st.selectbox('Select pango lineage to compare with', sorted(compare_to['pango'].unique()))
    compare_plot = corona.plot_compare(pango)
    st.altair_chart(compare_plot)
    
