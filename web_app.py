import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import sys
import os.path
import flask
import cowell as cw
import numpy as np
import math as mp
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go


app = dash.Dash()
app.title = 'Auth Propagation Model'

# # setup static folder
# STATIC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), 'static'))
#
# print(STATIC_PATH)
#
# @app.server.route('/static/<resource>')
# def serve_static(resource):
# 	return flask.send_from_directory(STATIC_PATH, resource)

# app.css.append_css({'external_url': '/static/style.css'})

app.layout = html.Div([

	html.H1('AUTH Propagation Model'),

	html.Div(children=['Satellite Orbit Determination using the propagation model AUTH'], style={'margin': '20px', 'font-size':'16pt'}),

	html.Div(children=[
		html.H3('Initial Parameters'),
		html.Label('State Vector (m - m/s)'),
		dcc.Input(value='4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03', type='text', style={'margin-left': '10px'}, id='state-vector'),
		html.Br(),

		html.Label('Initial Time (sec)'),
		dcc.Input(value='0', type='text', style={'margin-left': '40px'}, id='init-time'),
		html.Br(),

		html.Label('Final Time (sec)'),
		dcc.Input(value='10000', type='text', style={'margin-left': '46px'}, id='final-time'),
		html.Br(),

		html.Label('Step Size (sec)'),
		dcc.Input(value='200', type='text', style={'margin-left': '56px'}, id='step-size'),
	], style={'margin': '20px', }),

	html.Div(children=[
		html.H3('Force Model'),

		html.Label("Select Earth's Gravity Model", style={'margin-top': '10px'}),
		html.Br(),
		dcc.Dropdown(
			id="earth-model",
			options=[
				{'label': 'Keplerian', 'value': 1},
				{'label': 'Geodynamic Model', 'value': 2},
			],
			value=2
		),
		html.Br(),

		html.Label("Select other forces you want to include to the Model"),
		html.Br(),
		dcc.Checklist(
			id="other-forces",
			options=[
				{'label': "Sun's and Moon's Gravitational Force", 'value': 1},
				{'label': 'Solar Radiation', 'value': 2},
				{'label': 'Atmospheric Drag', 'value': 3},
			],
			values=[]
		),

	], style={'margin': '20px', 'width': '30%'}),


	html.Div(children=[

		html.Button(children="Determine Satellite Orbit", id="calculate-orbit", n_clicks=0)

	], style={'margin': '20px'}),

	html.Div(id='output-state', style={'margin': '40px'})


], style={'columnCount': 1})


@app.callback(Output('output-state', 'children'),
				[Input('calculate-orbit', 'n_clicks')], [
				State('state-vector', 'value'),
				State('init-time', 'value'),
				State('final-time', 'value'),
				State('step-size', 'value'),
				State('earth-model', 'value'),
				State('other-forces', 'values'),
				])
def update_output(n_clicks, input1, input2, input3, input4, input5, input6):

	y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
	index = 0
	parseInput1 = input1.split(",")
	for i in parseInput1:
		y[index] = float(i)
		index += 1

	if input2 == "" or input3 == "" or input4 == "":
		pass
	else:
		t0 = float(input2)
		tf = float(input3)
		step = float(input4)

		kepOrGeo = int(input5)
		other = list(input6)

		if kepOrGeo == 2:
			data = cw.read_geopotential_coeffs("restruct_EGM2008.gfc", False)
			C, S = cw.createNumpyArrayForCoeffs(data, 10, 10)
		else:
			C = 0
			S = 0

		sunAndMoon = False
		solar = False
		drag = False

		if 1 in other:
			sunAndMoon = True

		if 2 in other:
			solar = True

		if 3 in other:
			drag = True

		loopTf = step
		loopIndex = int((tf - t0) / step)
		final = np.zeros((loopIndex, 6))
		time = np.zeros(loopIndex)

		for i in range(0, loopIndex):
			final[i, :] = cw.rk4(y, t0, loopTf, 1, kepOrGeo, solar, sunAndMoon, drag, C, S)
			time[i] = t0
			t0 = loopTf
			loopTf = loopTf + step
			y = final[i, :]

		altitude = np.sqrt(final[:, 0]**2 + final[:, 1]**2 + final[:, 2]**2)

		returnTable = pd.DataFrame(data=final[:, :],
								   index=range(0, loopIndex),
								   columns=["rx(m)", "ry(m)", "rz(m)", "vx(m/s)", "vy(m/s)", "vz(m/s)"])

		return html.Table(
			# Header
			[html.Tr([html.Th(col) for col in returnTable.columns])] +

			# Body
			[html.Tr([
				html.Td(returnTable.iloc[i][col]) for col in returnTable.columns
			]) for i in range(min(len(returnTable), loopIndex))]
		), html.Br(), \
		dcc.Graph(
			id='alitutde-graph',
			figure={
				'data': [
					go.Scatter(
						x=time,
						y=altitude,
						opacity=0.7,
					)
				],
			}
		)



	# y = np.array([4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03])
	# t0, tf = 0, 100.00
	# final = np.zeros((100, 6))
	# step =
	#
	# # restruct_geopotential_file("EGM2008.gfc")


if __name__ == '__main__':
	app.run_server(debug=True)
