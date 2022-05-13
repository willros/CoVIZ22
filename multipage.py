import streamlit as st
from PIL import  Image

class MultiPage: 
    def __init__(self) -> None:
        self.pages = []
    
    def add_page(self, title, func) -> None: 
        self.pages.append({"title": title, "function": func})

    def run(self):
        img = 'logo.png'
        display = Image.open(img)
        col1, col2, col3 = st.sidebar.columns(3)
        col2.image(display, width=100)
        page = st.sidebar.selectbox('Navigation', self.pages, format_func=lambda page: page['title'])
        page['function']()

        
