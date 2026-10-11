# Authors: Gianvittorio Luria 
# License: BSD (3-clause)

import os
import socket
import tomli
import numpy as np
import pandas as pd
import os.path as op
import datetime as dt

with open(op.join(os.getcwd(), "config.toml"), "rb") as f:
    config = tomli.load(f)
    
USER = config['user'][socket.gethostname()]

raw_data_folder = config['folders'][USER]['raw_data_folder']
print(f'Data will be loaded from {raw_data_folder}')
results_folder = config['folders'][USER]['results_folder']
print(f'Results will be saved in {results_folder}')

# Load files
db_drivers = pd.read_csv(op.join(raw_data_folder, 'drivers_summary.tsv'), sep='\t')
db_metadata = pd.read_csv(op.join(raw_data_folder, 'metadata_filtered.tsv'), sep='\t')
db_proteins = pd.read_csv(op.join(raw_data_folder, 'proteins_filtered.tsv'), sep='\t')
db_lipids = pd.read_csv(op.join(raw_data_folder, 'lipids_filtered.tsv'), sep='\t')

# Create save subdirectories
if not os.path.exists(op.join(results_folder, 'npy')):
    os.makedirs(op.join(results_folder, 'npy'))
if not os.path.exists(op.join(results_folder, 'csv')):
    os.makedirs(op.join(results_folder, 'csv'))

_disease_name_list = list()
_feature_name_list = list()
_feature_type_list = list()
_timepoint_list = list()
_sample_list = list()
_value_list = list()
_sample_disease_list = list()

for index, row in db_drivers.iterrows():
    _disease = row['Outcome']
    _feature_type = row['Layer']
    _feature_name = row['FeatureName']
    print(f'*** {_disease} - {_feature_name} - {_feature_name} ***')
    
    for _timepoint in [f'T{x}' for x in range(5)]:
        print(_timepoint)
        db_aux1 = db_metadata[db_metadata.timepoint == _timepoint][['sample_id', _disease]].dropna()
        _yes_samples = db_aux1[db_aux1[_disease]==1]['sample_id'].to_list()
        _no_samples = db_aux1[db_aux1[_disease]==0]['sample_id'].to_list()
        
        match _feature_type:
            case 'Protein':
                print(f'protein {_feature_name}')
                db_aux2 = db_proteins[db_proteins.FeatureName==_feature_name][[_c for _c in db_proteins.columns if _c.endswith(_timepoint)]].T
            case 'Lipid':
                print(f'lipid {_feature_name}')
                db_aux2 = db_lipids[db_lipids.FeatureName==_feature_name][[_c for _c in db_lipids.columns if _c.endswith(_timepoint)]].T
        
        _yes_values = db_aux2[db_aux2.index.isin(_yes_samples)][db_aux2.columns[0]].to_numpy()
        _no_values = db_aux2[db_aux2.index.isin(_no_samples)][db_aux2.columns[0]].to_numpy()
        # save values into numpy arrays 
        np.save(op.join(results_folder, 'npy', f"{_timepoint}__{_disease}__{_feature_name.replace(' ', '_').replace(':','_').replace(';','_').replace('/','_')}__yes_values"), _yes_values)
        np.save(op.join(results_folder, 'npy', f"{_timepoint}__{_disease}__{_feature_name.replace(' ', '_').replace(':','_').replace(';','_').replace('/','_')}__no_values"), _no_values)
        
        for _s, _v in zip(_yes_samples, _yes_values):
            _disease_name_list.append(_disease)
            _feature_name_list.append(_feature_name)
            _feature_type_list.append(_feature_type)
            _timepoint_list.append(_timepoint)
            _sample_list.append(_s)
            _value_list.append(_v)
            _sample_disease_list.append(1)
            
        for _s, _v in zip(_no_samples, _no_values):
            _disease_name_list.append(_disease)
            _feature_name_list.append(_feature_name)
            _feature_type_list.append(_feature_type)
            _timepoint_list.append(_timepoint)
            _sample_list.append(_s)
            _value_list.append(_v)
            _sample_disease_list.append(0)

new_db = pd.DataFrame(data={
    'disease_name': _disease_name_list,
    'feature_name': _feature_name_list,
    'feature_type': _feature_type_list,
    'timepoint': _timepoint_list,
    'sample': _sample_list,
    'value': _value_list,
    'disease': _sample_disease_list    
})

# save new DataFrame in csv format
new_db.to_csv(op.join(results_folder, 'csv', 'all_data.csv'), sep=';', index=False)