# Authors: Gianvittorio Luria 
# License: BSD (3-clause)

import os
import socket
import tomli
import pandas as pd
import os.path as op
import plotly.express as px
import pingouin as pg

with open(op.join(os.getcwd(), "config.toml"), "rb") as f:
    config = tomli.load(f)
    
USER = config['user'][socket.gethostname()]

raw_data_folder = config['folders'][USER]['raw_data_folder']
print(f'Data will be loaded from {raw_data_folder}')
results_folder = config['folders'][USER]['results_folder']
print(f'Results will be saved in {results_folder}')

# Load files
db_drivers = pd.read_csv(op.join(raw_data_folder, 'drivers_summary.tsv'), sep='\t')
db_all = pd.read_csv(op.join(results_folder, 'csv', 'all_data.csv'), sep=';')

# Create save subdirectories
if not os.path.exists(op.join(results_folder, 'plots')):
    os.makedirs(op.join(results_folder, 'plots'))
if not os.path.exists(op.join(results_folder, 'plots', 'png')):
    os.makedirs(op.join(results_folder, 'plots', 'png'))
if not os.path.exists(op.join(results_folder, 'plots', 'html')):
    os.makedirs(op.join(results_folder, 'plots', 'html'))

db_list = list()
for index, row in db_drivers.iterrows():
	_disease_name = row['Outcome']
	_feature_type = row['Layer']
	_feature_name = row['FeatureName']
	print(f'*** {_disease_name} - {_feature_name} - {_feature_type} ***')
	db_aux = db_all[(db_all.feature_name==_feature_name)&(db_all.disease_name==_disease_name)].copy()
	db_aux.replace({'disease':{0:'No', 1:'Yes'}}, inplace=True)
	# MWU
	_all_times = [f'T{x}' for x in range(5)]
	_times_list = list()
	_outcome_list = list()
	_feature_type_list = list()
	_feature_name_list = list()
	_ann_dict_mwu = dict()
	_stat_list = list()
	_pval_list = list()
	_fc_list = list()
	for _t in _all_times:
		_times_list.append(_t)
		_outcome_list.append(_disease_name)
		_feature_type_list.append(_feature_type)
		_feature_name_list.append(_feature_name)
		_no = db_aux[(db_aux.timepoint==_t)&(db_aux.disease=='No')]['value'].to_numpy()
		_yes = db_aux[(db_aux.timepoint==_t)&(db_aux.disease=='Yes')]['value'].to_numpy()
		result_mwu = pg.mwu(_no,_yes).iloc[0]  
		_ann_dict_mwu[_t] = result_mwu['p-val']
		_stat_list.append(result_mwu['U-val'])
		_pval_list.append(result_mwu['p-val'])
		_fc_list.append(result_mwu['CLES'])
	db_partial_mwu = pd.DataFrame(data={'timepoint':_times_list, 'outcome':_outcome_list, 
										'feature':_feature_name_list,
										'feature_type':_feature_type_list, 'U':_stat_list,
										'pval':_pval_list, 'CLES':_fc_list})
	db_list.append(db_partial_mwu)
	
	# Plot
	_savefig = True
	_savename= f"BOXPLOT_{_feature_name.replace(' ', '_').replace(':','_').replace(';','_').replace('/','_')}_{_disease_name}"
	_annotation = _ann_dict_mwu#None

	fig = px.box(db_aux, x="timepoint", y="value", color="disease",
				hover_data=['sample'],
				category_orders={'disease':['No', 'Yes']}, 
				color_discrete_map={'Yes':'firebrick', 'No':'forestgreen'},
				labels={'disease':f"{_disease_name.replace('_', ' ').capitalize()}",
						'timepoint': 'Timepoint',
						'value':f' {_feature_name} log<sub>2</sub>(intensity)'}
				)

	fig.update_layout(
	height=300, width=600,
	boxgroupgap=0.3, boxgap=0.3,
	margin=dict(l=20, r=20, t=30, b=40),
	paper_bgcolor='aliceblue',
	)

	# statistical annotation
	if _annotation is not None:
		_savename += '_annotated'
		for _x, _pval in _annotation.items():
			if _pval >= 0.05:
				symbol = 'ns'
			elif _pval >= 0.01:
				symbol = '*'
			elif _pval >= 0.001:
				symbol = '**'
			else:
				symbol = '***'

			fig.add_annotation(
				yref='paper',x=_x, y=1.1,
				text=f"<b>{symbol}</b>", 
				font=dict(size=16),
				showarrow=False)


	## save
	if _savefig:
		# PNG
		fig.write_image(op.join(results_folder, 'plots', 'png', f"{_savename}.png"),
						scale=3)
		# HTML
		fig.write_html(op.join(results_folder, 'plots', 'html', f"{_savename}.html"))
		## show
		fig.show()

# Create statistical test results DataFrame
db_mwu = pd.concat(db_list)
db_mwu.reset_index(inplace=True, drop=True)

# Save DataFrame in csv format
db_mwu.to_csv(op.join(results_folder, 'csv', 'db_mwu.csv'), sep=';', index=False)
			
