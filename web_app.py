import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import sys
import os.path
import flask


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
		html.Label('State Vector'),
		dcc.Input(value='', type='text', style={'margin-left': '10px'}, id='state-vector'),
		html.Br(),

		html.Label('Initial Time'),
		dcc.Input(value='', type='text', style={'margin-left': '13px'}, id='init-time'),
		html.Br(),

		html.Label('Final Time'),
		dcc.Input(value='', type='text', style={'margin-left': '18px'}, id='final-time'),
		html.Br(),

		html.Label('Step Size'),
		dcc.Input(value='', type='text', style={'margin-left': '28px'}, id='step-size'),
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
			  [Input('calculate-orbit', 'n_clicks')],
			  [State('state-vector', 'value'),
			   State('init-time', 'value'),
			   State('final-time', 'value'),
			   State('earth-model', 'value')])
def update_output(n_clicks, input1, input2, input3, input4):

	print(type(input1))
	print(type(input2))
	print(type(input3))
	print(input4)


if __name__ == '__main__':
	app.run_server(debug=True)
