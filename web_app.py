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

app.layout = html.Div([

	html.H1('AUTH Propagation Model'),

	html.Div(children=['Satellite Orbit Determination using the propagation model AUTH'], style={'margin': '20px', 'font-size':'16pt'}),

	html.Div(children=[

		html.Div(children=[
			html.H3('Initial Parameters'),
			html.Label('State Vector (m - m/s)', className="main"),
			dcc.Textarea(value='4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03', id='state-vector', style={"margin-bottom": "5px", "min-height": "150px"}),

			html.Label('Initial Time (sec)', className="main"),
			dcc.Input(value='0', type='text', id='init-time', style={"margin-bottom": "5px"}),

			html.Label('Final Time (sec)', className="main"),
			dcc.Input(value='1000', type='text', id='final-time', style={"margin-bottom": "5px"}),

			html.Label('Step Size (sec)', className="main"),
			dcc.Input(value='100', type='text', id='step-size', style={"margin-bottom": "5px"}),
		], style={'margin': '20px', }),

		html.Div(children=[
			html.H3('Force Model'),

			html.Label("Select Earth's Gravity Model", style={'margin-top': '10px'}, className="main"),
			dcc.Dropdown(
				id="earth-model",
				options=[
					{'label': 'Keplerian', 'value': 1},
					{'label': 'Geodynamic Model', 'value': 2},
				],
				value=1
			),
			html.Br(),

			html.Div(id="geodynamic-parameters", style={'display': 'none'}, children=[
				dcc.Upload(
					id='upload-data',
					children=html.Div([
						html.A('Upload Geodynamic Model')
					]),
					style={
						'width': '100%',
						'height': '60px',
						'lineHeight': '60px',
						'borderWidth': '1px',
						'borderStyle': 'dashed',
						'borderRadius': '5px',
						'textAlign': 'center',
					},
					# Allow multiple files to be uploaded
					multiple=False
				),
				html.Label('Order (n)', className="main"),
				dcc.Input(value='20', type='text', id='order-geo', style={"margin-bottom": "5px"}),
				html.Label('Degree (m)', className="main"),
				dcc.Input(value='20', type='text', id='degree-geo', style={"margin-bottom": "5px"}),
			]),

			html.Label("Select other forces you want to include to the Model", className="main"),
			dcc.Checklist(
				id="other-forces",
				options=[
					{'label': "Sun's and Moon's Gravitational Force", 'value': 1},
					{'label': 'Solar Radiation', 'value': 2},
					{'label': 'Atmospheric Drag', 'value': 3},
				],
				values=[],
				style={"text-align": "left"}),
			html.Div(id="radiation-parameters", style={'display': 'none'}, children=[
				html.Label('Solar Radiation Reflectivity - Absorption (e)', className="main"),
				dcc.Input(value='0.8', type='text', id='solar-epsilon', style={"margin-bottom": "5px"}),
			]),
			html.Div(id="drag-parameters", style={'display': 'none'}, children=[
				html.Label('Atmospheric Drag Coefficient (Cd)', className="main"),
				dcc.Input(value='2.1', type='text', id='drag-coeff', style={"margin-bottom": "5px"}),
			]),

		], style={'margin': '20px'}),
	], style={'display': 'inline-flex'}),

	html.Div(children=[

		html.Button(children="Determine Satellite Orbit", id="calculate-orbit", n_clicks=0)

	], style={'margin': '20px'}),

	html.Div(id='output-state', style={'margin': '40px'})


], style={'columnCount': 1, "text-align": "center"})


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

	if n_clicks == 0:
		pass
	else:
		propagatorObj = cw.propagator()
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
				data = propagatorObj.read_geopotential_coeffs("restruct_EGM2008.gfc", False)
				C, S = propagatorObj.createNumpyArrayForCoeffs(data, 10, 10)
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
				final[i, :] = propagatorObj.rk4(y, t0, loopTf, step/50, kepOrGeo, solar, sunAndMoon, drag, C, S)
				time[i] = loopTf
				t0 = loopTf
				loopTf = loopTf + step
				y = final[i, :]
				print(t0, loopTf)

			altitude = np.sqrt(final[:, 0]**2 + final[:, 1]**2 + final[:, 2]**2)
			beforeDataFrame = np.zeros((loopIndex, 7))
			beforeDataFrame[:, 1:7] = final[:]
			beforeDataFrame[:, 0] = time[:]

			returnTable = pd.DataFrame(data=beforeDataFrame[:, :],
									   index=range(0, loopIndex),
									   columns=["Time(sec)", "rx(m)", "ry(m)", "rz(m)", "vx(m/s)", "vy(m/s)", "vz(m/s)"])


			return html.Div(children=[ html.Table(
				# Header
				[html.Tr([html.Th(col) for col in returnTable.columns])] +

				# Body
				[html.Tr([
					html.Td(returnTable.iloc[i][col]) for col in returnTable.columns
				]) for i in range(min(len(returnTable), loopIndex))]
			)],style={"padding-left": "20%"}), html.Br(), \
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


if __name__ == '__main__':
	app.run_server(debug=True)
