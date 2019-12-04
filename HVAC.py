import pandas as pd
import os


#perfdf = pd.read_excel('Perf.xlsx')
#perfdf.to_latex('perf.tex')
#statsdf = pd.read_excel('Stats.xlsx')
#statsdf.to_latex('Stats.tex')
#uncertdf = pd.read_excel('Uncert.xlsx')
#uncertdf.to_latex('Uncert.tex')

resdf = pd.read_excel('results.xlsx')
resdf.to_latex('res.tex')