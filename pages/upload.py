import streamlit as st
from corona_viz.nextclade import Nextclade

def app():
    st.title('Upload file')
    fasta = st.file_uploader("FASTA file to run Nextclade on")
    if fasta is not None:
        st.markdown('#### Running Nextclade on the fasta file... Please wait a second.')
        with open('input.fasta', 'w+') as f:
            print(fasta.read().decode(), file=f)
            
        nextclade = Nextclade()
        nextclade.run_nextclade()
        st.markdown('##### Nextclade done!')
        #nextclade.make_meta_df()
        #st.write(nextclade.meta_df)
