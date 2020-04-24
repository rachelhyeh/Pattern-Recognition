# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 11:46:05 2018

@author: Rachel Yeh
"""

import numpy as np
import xlrd
import itertools
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import svm
from sklearn.neural_network import MLPClassifier

##########Data Analysis
def data (input_workbook, output_data):
    input_sheet = input_workbook.sheet_by_index(0)
    classifier = list()
    for row in range (1, input_sheet.nrows):
        output_data.append(input_sheet.row_values(row))
    tt = (np.array(output_data).T).tolist()
    tt.remove(tt[0])
    classifier.extend((np.array(tt[70]).T).tolist())
    tt.remove(tt[70])
    output_data = (np.array(tt).T).tolist()
    return output_data, classifier
##########Training Data
def training ():   
    training_data = list()
    training_workbook = xlrd.open_workbook("Training_Data.xlsm")
    return data(training_workbook, training_data)
##########Testing Data 
def testing ():
    testing_data = list()
    testing_workbook = xlrd.open_workbook("Testing_Data.xlsm")
    return data(testing_workbook, testing_data)
##########Exhaustive Search
def top2 (length, n):
    p = list(range(0,length))
    return np.asarray(list(itertools.combinations(p, n)))
##########Error Count
def error (orig, deter):
    error = 0;
    for i in range(len(orig.T)):
        if orig[i] != deter[i]:
            error = error + 1
    return error / len(deter.T)
    
##########Top 2 Genes/Echaustive Search
train = np.array(training()[0])
trn_lab = np.array(training()[1])
test = np.array(testing()[0])
tx_lab = np.array(testing()[1])
exh = top2(len(train[0]), 2)
appError_lda_exh = 1; appError_lsvm_exh = 1; appError_nsvm_exh = 1; appError_nn_exh = 1;

for i in range (len(exh.T[0])):
    trn_x = list(); tx_x = list();
    for j in range (2):
        trn_x.append(train.T[exh[i][j]])
        tx_x.append(test.T[exh[i][j]])
    trn_x = np.array(trn_x).T
    tx_x = np.array(tx_x).T
    ##########LDA, p=0.75  
    exh_lda = LinearDiscriminantAnalysis(priors = np.array([0.25, 0.75]))
    exh_lda.fit(trn_x,trn_lab)
    trn_y_lda = exh_lda.predict(trn_x)
    err_lda_exh = error(trn_lab, trn_y_lda)
    if err_lda_exh < appError_lda_exh:
        appError_lda_exh = err_lda_exh;
        top2_lda = exh[i]
        tx_y_lda_exh = exh_lda.predict(tx_x)
        txError_lda_exh = error(tx_lab, tx_y_lda_exh)
    ##########Linear SVM, C=1
    exh_lsvm = svm.SVC(kernel='linear')
    exh_lsvm.fit(trn_x,trn_lab)
    trn_y_lsvm = exh_lsvm.predict(trn_x)    
    err_lsvm_exh = error(trn_lab, trn_y_lsvm)
    if err_lsvm_exh < appError_lsvm_exh:
        appError_lsvm_exh = err_lsvm_exh;
        top2_lsvm = exh[i]
        tx_y_lsvm_exh = exh_lsvm.predict(tx_x)        
        txError_lsvm_exh = error(tx_lab, tx_y_lsvm_exh)
    ##########NonLinear SVM with Gaussian RBF kernel, C=1
    exh_nsvm = svm.SVC(kernel='rbf')
    exh_nsvm.fit(trn_x,trn_lab)
    trn_y_nsvm = exh_nsvm.predict(trn_x)    
    err_nsvm_exh = error(trn_lab, trn_y_nsvm)
    if err_nsvm_exh < appError_nsvm_exh:
        appError_nsvm_exh = err_nsvm_exh;
        top2_nsvm = exh[i]
        tx_y_nsvm_exh = exh_nsvm.predict(tx_x)        
        txError_nsvm_exh = error(tx_lab, tx_y_nsvm_exh)      
#    #########NN with 5 neurons in one hidden layer
#    exh_nn = MLPClassifier(solver='lbfgs', hidden_layer_sizes=(5.))
#    exh_nn.fit(trn_x,trn_lab)
#    trn_y_nn = exh_nn.predict(trn_x)    
#    err_nn_exh = error(trn_lab, trn_y_nn)
#    if err_nn_exh < appError_nn_exh:
#        appError_nn_exh = err_nn_exh;
#        top2_nn = exh[i]
#        tx_y_nn_exh = exh_nn.predict(tx_x)        
#        txError_nn_exh = error(tx_lab, tx_y_nn_exh)


##########Top 3-4 Genes/Sequencial Forward Search
lda = lsvm = nsvm = nn = np.array(range(len(train.T))).tolist()
appError_lda_sf = []; appError_lsvm_sf = []; appError_nsvm_sf = []; appError_nn_sf = [];
txError_lda_sf = []; txError_lsvm_sf = []; txError_nsvm_sf = []; txError_nn_sf = [];
top5_lda = top2_lda.tolist(); top5_lsvm = top2_lsvm.tolist(); top5_nsvm = top2_nsvm.tolist();
#top5_nn = top2_nn.tolist();
lda = list(set(lda)-set(top5_lda)); lsvm = list(set(lsvm)-set(top5_lsvm)); nsvm = list(set(nsvm)-set(top5_nsvm));
#nn = list(set(nn)-set(top5_nn));
for i_sf in range(3):
    ##########LDA, p=0.75    
    errorApp_lda = 1; num_lda = 0; 
    for j in range (len(lda)):
        trn_x_lda = list(); tx_x_lda = list();
        for k in range (len(top5_lda)):
            trn_x_lda.append(train.T[top5_lda[k]])
            tx_x_lda.append(test.T[top5_lda[k]])
        trn_x_lda.append(train.T[lda[j]])
        tx_x_lda.append(test.T[lda[j]]) 
        
        trn_x_lda = np.array(trn_x_lda).T
        tx_x_lda = np.array(tx_x_lda).T
        sf_lda = LinearDiscriminantAnalysis(priors = np.array([0.25, 0.75]))
        sf_lda.fit(trn_x_lda, trn_lab)
        trn_y_lda = sf_lda.predict(trn_x_lda)
        err_lda_sf = error(trn_lab, trn_y_lda)
        if err_lda_sf < errorApp_lda:
            errorApp_lda = err_lda_sf;
            num_lda = lda[j]
            tx_y_lda_sf = sf_lda.predict(tx_x_lda)
            errorTx_lda = error(tx_lab, tx_y_lda_sf)
    top5_lda.append(num_lda)
    lda = list(set(lda)-set(top5_lda))
    appError_lda_sf.append(errorApp_lda)
    txError_lda_sf.append(errorTx_lda)
           
    ##########Linear SVM, C=1       
    errorApp_lsvm = 1; num_lsvm = 0; 
    for j_lsvm in range (len(lsvm)):
        trn_x_lsvm = list(); tx_x_lsvm = list();
        for k_lsvm in range (len(top5_lsvm)):
            trn_x_lsvm.append(train.T[top5_lsvm[k_lsvm]])
            tx_x_lsvm.append(test.T[top5_lsvm[k_lsvm]])
        trn_x_lsvm.append(train.T[lsvm[j_lsvm]])
        tx_x_lsvm.append(test.T[lsvm[j_lsvm]]) 
        
        trn_x_lsvm = np.array(trn_x_lsvm).T
        tx_x_lsvm = np.array(tx_x_lsvm).T
        sf_lsvm = svm.SVC(kernel='linear')
        sf_lsvm.fit(trn_x_lsvm, trn_lab)
        trn_y_lsvm = sf_lsvm.predict(trn_x_lsvm)
        err_lsvm_sf = error(trn_lab, trn_y_lsvm)
        if err_lsvm_sf < errorApp_lsvm:
            errorApp_lsvm = err_lsvm_sf;
            num_lsvm = lsvm[j_lsvm];
            tx_y_lsvm_sf = sf_lsvm.predict(tx_x_lsvm)
            errorTx_lsvm = error(tx_lab, tx_y_lsvm_sf)
    top5_lsvm.append(num_lsvm)
    lsvm = list(set(lsvm)-set(top5_lsvm))
    appError_lsvm_sf.append(errorApp_lsvm)
    txError_lsvm_sf.append(errorTx_lsvm)        
    ##########NonLinear SVM with Gaussian RBF kernel, C=1            
    errorApp_nsvm = 1; num_nsvm = 0; 
    for j_nsvm in range (len(nsvm)):
        trn_x_nsvm = list(); tx_x_nsvm = list();
        for k_nsvm in range (len(top5_nsvm)):
            trn_x_nsvm.append(train.T[top5_nsvm[k_nsvm]])
            tx_x_nsvm.append(test.T[top5_nsvm[k_nsvm]])
        trn_x_nsvm.append(train.T[nsvm[j_nsvm]])
        tx_x_nsvm.append(test.T[nsvm[j_nsvm]]) 
        
        trn_x_nsvm = np.array(trn_x_nsvm).T
        tx_x_nsvm = np.array(tx_x_nsvm).T
        sf_nsvm = svm.SVC(kernel='rbf')
        sf_nsvm.fit(trn_x_nsvm, trn_lab)
        trn_y_nsvm = sf_nsvm.predict(trn_x_nsvm)
        err_nsvm_sf = error(trn_lab, trn_y_nsvm)
        if err_nsvm_sf < errorApp_nsvm:
            errorApp_nsvm = err_nsvm_sf;
            num_nsvm = nsvm[j_nsvm]
            tx_y_nsvm_sf = sf_nsvm.predict(tx_x_nsvm)
            errorTx_nsvm = error(tx_lab, tx_y_nsvm_sf)
    top5_nsvm.append(num_nsvm)
    nsvm = list(set(nsvm)-set(top5_nsvm))
    appError_nsvm_sf.append(errorApp_nsvm)
    txError_nsvm_sf.append(errorTx_nsvm)        
#    #########NN with 5 neurons in one hidden layer            
#    errorApp_nn = 1; num_nn = 0; 
#    for j_nn in range (len(nn)):
#        trn_x_nn = list(); tx_x_nn = list();
#        for k_nn in range (len(top5_nn)):
#            trn_x_nn.append(train.T[top5_nn[k_nn]])
#            tx_x_nn.append(test.T[top5_nn[k_nn]])
#        trn_x_nn.append(train.T[nn[j_nn]])
#        tx_x_nn.append(test.T[nn[j_nn]]) 
#        
#        trn_x_nn = np.array(trn_x_nn).T
#        tx_x_nn = np.array(tx_x_nn).T
#        sf_nn = MLPClassifier(solver='lbfgs', hidden_layer_sizes=(5.))
#        sf_nn.fit(trn_x_nn, trn_lab)
#        trn_y_nn = sf_nn.predict(trn_x_nn)    
#        err_nn_sf = error(trn_lab, trn_y_nn)
#        if err_nn_sf < errorApp_nn:
#            errorApp_nn = err_nn_sf;
#            num_nn = nn[j_nn]
#            tx_y_nn_sf = sf_nn.predict(tx_x_nn)
#            errorTx_nn = error(tx_lab, tx_y_nn_sf)
#    top5_nn.append(num_nn)
#    nn = list(set(nn)-set(top5_nn))
#    appError_nn_sf.append(errorApp_nn)
#    txError_nn_sf.append(errorTx_nn)   

##########All Genes/No Feature Selection
    ##########LDA, p=0.75
all_lda = LinearDiscriminantAnalysis(priors = np.array([0.25, 0.75]))
all_lda.fit(train,trn_lab)
trn_y_lda_all = all_lda.predict(train)  
appError_lda_all = error(trn_lab, trn_y_lda_all)
tx_y_lda_all = all_lda.predict(test)   
txError_lda_all = error(tx_lab, tx_y_lda_all)
    ##########Linear SVM, C=1 
all_lsvm = svm.SVC(kernel='linear')
all_lsvm.fit(train,trn_lab)
trn_y_lsvm_all = all_lsvm.predict(train)  
appError_lsvm_all = error(trn_lab, trn_y_lsvm_all)
tx_y_lsvm_all = all_lsvm.predict(test)  
txError_lsvm_all = error(tx_lab, tx_y_lsvm_all)
    ##########NonLinear SVM with Gaussian RBF kernel, C=1 
all_nsvm = svm.SVC(kernel='rbf')
all_nsvm.fit(train,trn_lab)
trn_y_nsvm_all = all_nsvm.predict(train)    
appError_nsvm_all = error(trn_lab, trn_y_nsvm_all)
tx_y_nsvm_all = all_nsvm.predict(test)    
txError_nsvm_all = error(tx_lab, tx_y_nsvm_all)
#    #########NN with 5 neurons in one hidden layer 
#all_nn = MLPClassifier(solver='lbfgs', hidden_layer_sizes=(5.))
#all_nn.fit(train,trn_lab)
#trn_y_nn_all = all_nn.predict(train)    
#appError_nn_all = error(trn_lab, trn_y_nn_all)
#tx_y_nn_all = all_nn.predict(test)    
#txError_nn_all = error(tx_lab, tx_y_nn_all)