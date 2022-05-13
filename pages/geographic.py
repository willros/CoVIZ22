import streamlit as st
import altair as alt
from corona_viz.plotting import CoronaPlot
from corona_viz.nextclade import Nextclade
import sqlite3

def app():
    nextclade = Nextclade()
    corona = CoronaPlot(nextclade.QUERY, nextclade.REF, 'mutations.csv')
    df = nextclade.make_meta_df()
    pango = df.Nextclade_pango.values[0]
    name = df.seqName.values[0]

    connection = sqlite3.connect('mutations_geographical.db')

    st.markdown(f'##### Similar sequences to {name} in different locations of the world')
    geo_plot = corona.plot_geographical_similar(pango, connection)
    st.altair_chart(geo_plot)
    