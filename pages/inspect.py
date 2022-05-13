import streamlit as st
import altair as alt
from corona_viz.plotting import CoronaPlot
from corona_viz.nextclade import Nextclade

def app():
    nextclade = Nextclade()
    corona = CoronaPlot(nextclade.QUERY, nextclade.REF, 'mutations.csv')
    df = nextclade.make_meta_df()
    name = df.seqName.values[0]
    
    st.markdown(f'### Inspect mutations in {name}')
    inspect = corona.plot_mutations_inspection()
    st.altair_chart(inspect)
