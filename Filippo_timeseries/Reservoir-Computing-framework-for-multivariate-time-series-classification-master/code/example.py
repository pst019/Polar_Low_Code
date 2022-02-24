# General imports
import numpy as np
import scipy.io
from sklearn.preprocessing import OneHotEncoder

# Custom imports
from modules import RC_classifier

# ============ Network configuration and hyperparameter values ============
config = {}
config['dataset_name'] = 'JpVow'

config['seed'] = 1
np.random.seed(config['seed'])

# Parameters for the reservoir
config['n_internal_units'] = 450 # size of the reservoir
config['spectral_radius'] = 0.59 # largest eigenvalue of the reservoir
config['leak'] = 0.6 # amount of leakage in the reservoir state update (None or 1.0 --> no leakage)
config['connectivity'] = 0.25 # percentage of nonzero connections in the reservoir
config['input_scaling'] = 0.1 # scaling of the input weights
config['noise_level'] = 0.01 # noise in the reservoir state update
config['n_drop'] = 5 # transient states to be dropped
config['bidir'] = True # if True, use bidirectional reservoir
config['circ'] = False # use reservoir with circle topology

# Dimensionality reduction parameters
config['dimred_method'] ='tenpca' # options: {None (no dimensionality reduction), 'pca', 'tenpca'}
config['n_dim'] = 75 #  number of resulting dimensions after the dimensionality reduction procedure

# Type of MTS representation
config['mts_rep'] = 'reservoir' # MTS representation:  {'last', 'mean', 'output', 'reservoir'}
config['w_ridge_embedding'] = 10.0 # regularization parameter of the ridge regression

# Type of readout
config['readout_type'] = 'lin' # readout used for classification: {'lin', 'mlp', 'svm'}

# Linear readout parameters
config['w_ridge'] = 5.0 # regularization of the ridge regression readout

# SVM readout
config['svm_gamma'] = 0.005 # bandwith of the RBF kernel
config['svm_C'] = 5.0 # regularization for SVM hyperplane

# MLP readout parameters
config['mlp_layout'] = [10,10] # neurons in each MLP layer
config['num_epochs'] = 2000 # number of epochs 
config['p_drop'] = 0.1 # dropout probability (0=no drop)
config['w_l2'] = 0.001 # weight of the L2 regularization
config['nonlinearity'] = 'maxout' # type of activation function {'relu', 'tanh', 'maxout', 'kaf'}

print(config)

# ============ Load dataset ============
data = scipy.io.loadmat('../dataset/'+config['dataset_name']+'.mat')
X = data['X']  # shape is [N,T,V]
if len(X.shape) < 3:
    X = np.atleast_3d(X)
Y = data['Y']  # shape is [N,1]
Xte = data['Xte']
if len(Xte.shape) < 3:
    Xte = np.atleast_3d(Xte)
Yte = data['Yte']

print('Loading '+config['dataset_name']+' - Tr: '+ str(X.shape)+', Te: '+str(Xte.shape))

# One-hot encoding for labels
onehot_encoder = OneHotEncoder(sparse=False)
Y = onehot_encoder.fit_transform(Y)
Yte = onehot_encoder.transform(Yte)

# ============ Specify, train and evaluate model ============
classifier =  RC_classifier(
                          reservoir=None,     
                          n_internal_units=config['n_internal_units'],
                          spectral_radius=config['spectral_radius'],
                          leak=config['leak'],
                          connectivity=config['connectivity'],
                          input_scaling=config['input_scaling'],
                          noise_level=config['noise_level'],
                          circle=config['circ'],
                          n_drop=config['n_drop'],
                          bidir=config['bidir'],
                          dimred_method=config['dimred_method'], 
                          n_dim=config['n_dim'],
                          mts_rep=config['mts_rep'],
                          w_ridge_embedding=config['w_ridge_embedding'],
                          readout_type=config['readout_type'],            
                          w_ridge=config['w_ridge'],              
                          mlp_layout=config['mlp_layout'],
                          num_epochs=config['num_epochs'],
                          p_drop=config['p_drop'],
                          w_l2=config['w_l2'],
                          nonlinearity=config['nonlinearity'], 
                          svm_gamma=config['svm_gamma'],
                          svm_C=config['svm_C'],
                          seed=config['seed'])

tr_time = classifier.train(X,Y)
print('Training time = %.2f'%tr_time)

accuracy, f1 = classifier.test(Xte, Yte)
print('Accuracy = %.3f, F1 = %.3f'%(accuracy, f1))