import streamlit as st
from corona_viz.nextclade import Nextclade

def app():
    nextclade = Nextclade()
    df = nextclade.make_meta_df()
    pango = df.Nextclade_pango.values[0]
    clade = df.clade.values[0]
    overall = df['qc.overallStatus'].values[0].title()
    missing = df['qc.missingData.totalMissing'].values[0]
    substitutions = df['totalSubstitutions'].values[0]
    deletions = df['totalDeletions'].values[0]
    amino_sub = df['totalAminoacidSubstitutions'].values[0]
    amino_del = df['totalAminoacidDeletions'].values[0]
    reversion = str(df['privateNucMutations.reversionSubstitutions'].values[0]).split(',')
    labeled = str(df['privateNucMutations.labeledSubstitutions'].values[0]).split(',')
    unlabeled = str(df['privateNucMutations.unlabeledSubstitutions'].values[0]).split(',')

    st.markdown(f'### {df.seqName.values[0]}')
    st.markdown("""---""")

    col1, col2 = st.columns(2)
    
    col1.metric('Clade', f'{clade}')
    col2.metric('Pangolin', f'{pango}')
    
    col1.metric('Alignment score', f'{overall}')
    col2.metric('Missing nucleotides', f'{missing}')
    
    col1.metric('Total nt substitutions', f'{substitutions}')
    col2.metric('Total nt deletions', f'{deletions}')
    
    col1.metric('Total AA substitutions', f'{amino_sub}')
    col2.metric('Total AA deletions', f'{amino_del}')
    
    col1, col2, col3 = st.columns(3)
    col1.markdown('##### Reversion mutations')
    [col1.markdown(x) for x in reversion]
    col2.markdown('##### Labeled mutations')
    [col2.markdown(x) for x in labeled]
    col3.markdown('##### Novel mutations')
    [col3.markdown(x) for x in unlabeled]

    st.markdown("""---""")

