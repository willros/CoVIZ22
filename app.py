import os
import streamlit as st
import numpy as np
from PIL import  Image

from multipage import MultiPage
from pages import upload, compare, inspect, alignment, info, geographic

app = MultiPage()

# Main page
img = 'logo.png'
display = Image.open(img)
col1, col2, col3 = st.columns(3)
col2.image(display, width=300)

app.add_page('Upload Data', upload.app)
app.add_page('Compare against corona lineages', compare.app)
app.add_page('Inspect mutations', inspect.app)
app.add_page('Alignment information', alignment.app)
app.add_page('Sequence information', info.app)
app.add_page('Similar sequences', geographic.app)

app.run()


