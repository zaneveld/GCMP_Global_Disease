import pandas as pd
import matplotlib.pyplot as plt
import argparse
import altair as alt
from PIL import Image

parser = argparse.ArgumentParser(description="Summary of FRRP data")

parser.add_argument('-l', type=str, dest="lower_data", help="Enter path to lower table")
parser.add_argument('-u', type=str, dest="upper_data", help="Enter path to upper table")

# Handling the sourcing for the main chart 
opts = parser.parse_args()

cd_data= pd.read_csv('%s' % opts.lower_data)

cd_data = cd_data.melt(id_vars=['Species','Clade'])

domain = ['complex', 'robust', 'outgroup']
range_ = ['#f15362', "#4f7ac6", "#81c186"]

# sourcing for the sum chart

sum_data = pd.read_csv('%s' % opts.upper_data)

# ==== ENCODING ====

#Encoding X
x_encode = alt.X('Species:O')

#Encoding Y
y_encode = alt.Y('variable:O',\
	axis=alt.Axis(title = 'Disease'), \
	sort=alt.EncodingSortField(field = 'value', op = "average", order="descending")
)

#Encoding color
color_encode = alt.Color('Clade:O', scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(title="Clade"))

#Encoding Size
size_encode = alt.Size('value:Q', legend=alt.Legend(title="% Prevalence"))

# ==== CHART CREATION ====

main_chart = alt.Chart(cd_data).mark_circle().encode(
	x = x_encode,
	y = y_encode,
	color = color_encode,
	size = size_encode

	)

sum_chart = alt.Chart(sum_data).mark_bar().encode(
	x = 'Species:O',
	y = alt.Y('Total Disease:Q', title = '% with disease'),
	color = color_encode,
	).properties(height=100)

# combine the two charts into one image

full_chart = alt.vconcat(sum_chart, main_chart).configure(background="white")

# save graph

full_chart.save('../figures/summarydata.png')

#open image of chart

img = Image.open('../figures/summarydata.png')
img.show()


