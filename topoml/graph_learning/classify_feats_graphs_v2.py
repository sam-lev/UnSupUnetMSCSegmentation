import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report,f1_score
import torch.nn.functional as f
from .gbssl import CAMLP
from scipy import sparse

class graph_learners:
    def __init__(self):
        base_path = "./graph_learning/max_diadem/"
        print("working from base path: ", base_path)
        self.tests = "Deep walk, Confidence Label Prediction, and Random Forestswith ( Raw Features, Deep Walk, and Locally Filtered Features)"

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print('using', device)
        
        edges = pd.read_csv(base_path+'edges.csv')
        df = pd.read_csv(base_path+'arc_features.csv')
        labs = np.array(df.iloc[:,1])
        feats = np.array(df.iloc[:,2:])
        
        valid_ids = np.where(labs != 0)[0]
        unk_ids = np.where(labs == 0)[0]
        
        X = feats
        print('Features: {}'.format(X.shape))
        
        y = labs
        y[y == -1] = 0
        print('Labels: {}'.format(np.unique(y, return_counts=True)))
        
        A = np.array(edges)
        print('Edgelist: {}'.format(A.shape))
        print("edge list saved as 'edge.list'")
        np.savetxt(base_path+'edge.list', A.astype(int),fmt='%i', delimiter=' ')
        train_ids, test_ids = train_test_split(range(len(y)), test_size = 0.8, stratify=y)
        print(X[train_ids,:].shape)
        
        train_ids = np.setdiff1d(train_ids, unk_ids)
        test_ids = np.setdiff1d(test_ids, unk_ids)
        
        print('Purely Graph Based Approaches')
        # 1a. Semi-supervised Classifier
        G = sparse.lil_matrix((len(y),len(y)))
        for i in range(len(A)):
            G[A[i,0],A[i,1]] = 1
            G[A[i,1],A[i,0]] = 1
        G.tocsr()
        camlp = CAMLP(graph=G)#,beta=0, max_iter=50)
        camlp.fit(train_ids, y[train_ids])
        ytest_pred = camlp.predict(test_ids)
        print('Semi-Supervised Learning CAMLP')
        print('F1: ', f1_score(y[test_ids], ytest_pred))
        
        # 1b. Use deepwalk features for classification
        feats = np.loadtxt(base_path+'graph.embeddings',delimiter=' ')
        feats = feats[feats[:,0].argsort()]
        
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(feats[train_ids,:], y[train_ids])
        ytest_pred = clf.predict(feats[test_ids,:])
        print('Random Forests - Deepwalk embeddings')
        print('F1: ', f1_score(y[test_ids], ytest_pred))
        
        print('Feature Based Approaches')
        #
        # # 2a. Use the raw features for classiciation
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(X[train_ids,:], y[train_ids])
        ytest_pred = clf.predict(X[test_ids,:])
        print('Random Forests - Raw Features')
        print('F1: ', f1_score(y[test_ids], ytest_pred))
        #
        # 2b. Combine raw and deepwalk features
        X_new = np.concatenate((X,feats),axis=1)
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(X_new[train_ids,:], y[train_ids])
        ytest_pred = clf.predict(X_new[test_ids,:])
        print('RandomForests - Deep Walk + Features')
        print('F1: ', f1_score(y[test_ids], ytest_pred))
        
        # 2c. Locally filtered features
        nFeats = X.shape[1]
        W = np.eye(len(y),len(y))
        for i in range(len(A)):
            W[A[i,0],A[i,1]] = 1
            W[A[i,1],A[i,0]] = 1
        W =W/np.tile(np.sum(W, axis=1)[:,np.newaxis],[1,W.shape[1]])
        
        X_filt = W.dot(X)
        clf = RandomForestClassifier(n_estimators=50)
        clf.fit(X_filt[train_ids,:], y[train_ids])
        ytest_pred = clf.predict(X_filt[test_ids,:])
        print('RandomForests - Locally Filtered Features')
        print('F1: ', f1_score(y[test_ids], ytest_pred))
        
