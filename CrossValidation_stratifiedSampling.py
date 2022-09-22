#!/usr/bin/python3.7

### ROC cross-validation using sklearn ###

# import packages
#------------------

import numpy as np
import numpy
import pandas as pd
import matplotlib.pyplot as plt

from sklearn import svm
from sklearn.metrics import auc
from sklearn.pipeline import Pipeline
from sklearn.metrics import plot_roc_curve
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold # KFold stratified split 
from sklearn.model_selection import GridSearchCV # tune hyperparameters
from sklearn.model_selection import StratifiedShuffleSplit # stratified sampling
from sklearn.model_selection import train_test_split # split data to training and testset
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.preprocessing import LabelBinarizer
from sklearn.utils  import resample
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score

# read file
#------------
#FileName = '~/Documents/Mount_Sinai_Hospital_documents/Margaret/CRP_analysis_main/CRP_crossValidation/CRP_crossValidation_CRPMerge_alldata.csv' 
FileName = '~/Documents/Mount_Sinai_Hospital_documents/Margaret/CRP_analysis_main/CRP_crossValidation/CrossValidation_mainCrossVal2.csv' 

#############################################################################
# import files
# read datasets in csv fromat
CRP = pd.read_csv(FileName, index_col=0)
#print(CRP)

#Replace the species with 1,2 or 3 as appropriate
label_dict = dict()
label_dict['0'] = 'Male'
label_dict['1'] = 'Female'
CRP['sex'].replace(['Male', 'Female'], [0, 1], inplace=True)


# define the columns that I need
#columns = ["cxcl9", "antitnf_current"]
#columns = ["cxcl9"]
#columns = ["MeanCRP"]
columns = ["MeanCRP", "antitnf_current"]
#columns = ["cxcl9","MeanCRP","mmp1","il5","st1a1"]
#columns = ["cxcl9","MeanCRP","mmp1","il5","st1a1", "antitnf_current"]
#columns = ["cxcl9","MeanCRP","mmp1","il5","st1a1", "sample_age", "sex", "antitnf_current"]
#columns = ['il8','vegfa','mcp3','gdnf','cdcp1','cd244','il7','opg','laptgfbeta1','upa','il6','il17c','mcp1','il17a','cxcl11','axin1','trail','il20ra','cxcl9','cst5','il2rb','il1alpha','osm','il2','cxcl1','tslp','ccl4','cd6','scf','il18','slamf1','tgfalpha','mcp4','ccl11','tnfsf14','fgf23','il10ra','fgf5','mmp1','lifr','fgf21','ccl19','il15ra','il10rb','il22ra1','il18r1','pdl1','betangf','cxcl5','trance','hgf','il12b','il24','il13','artn','mmp10','il10','tnf','ccl23','cd5','ccl3','flt3l','cxcl6','cxcl10','ebp1','il20','sirt2','ccl28','dner','enrage','cd40','il33','ifngamma','fgf19','il4','lif','nrtn','mcp2','casp8','ccl25','cx3cl1','tnfrsf9','nt3','tweak','ccl20','st1a1','stampb','il5','ada','tnfb','csf1']

# read file as Dataframe taking the columns I need
x = pd.DataFrame(CRP, columns=columns)
X = x.to_numpy()
#print(x)

# extract y
columns = ['rutgeerts_score_bin']
y = pd.DataFrame(CRP, columns=columns)
y = y.to_numpy()
y = y.flatten() # transpose
#print(y)

# #############################################################################
#stratified sampling to split the data to training and test datasets no randomeState defined
# stratified sampling preserving the percentage of samples for each class
X_train, X_test, y_train, y_test = train_test_split( X, y, test_size=0.1, stratify=y)

NT=len(X_train)
#print(len(X_test))
#print(X_test)
#print(y_train)

# def average function
def Average(lst): 
    return sum(lst) / len(lst)


#################################################################################
# CrossValidation with Training dataset

# Run classifier with cross-validation and plot ROC curves
scaler = StandardScaler()
classifier = svm.SVC(kernel='linear', probability=True)
# stratified sampling for crossValidation


#pipe= Pipeline([('scaler', StandardScaler()), ('classifier', classifier)])
score = list()

# define how many time to run for a 10 flod crossvalidation
nRep = 5
for n in range(0,nRep,1):
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100) # for plotting gives equal random dots
	cv = StratifiedKFold(n_splits=5,shuffle=True)
	
	fig, ax = plt.subplots()
	for i, (train, test) in enumerate(cv.split(X_train, y_train)):
		#fit model
		classifier.fit(X_train[train], y_train[train])
		
		# predict y in the test dataset
		#y_pred = classifier.predict(X_test)
		#y_pred = classifier.predict(X_train[test])
		
		#yaccuracy = accuracy_score(y_test, y_pred)
		#YaccuracyMean=yaccuracy.mean()
		
		# append accuracy of the fit model for each n fold
		score.append(classifier.fit(X_train[train], y_train[train]).score(X_train[test], y_train[test]))
		average = Average(score)
		
		
		#print('Fold n ',n,i)
		#print(classifier.coef_)
		#print(confusion_matrix(y_test, y_pred, normalize='all'))
		#print('*'*50)
		
		# result visualization
		viz = plot_roc_curve(classifier, X_train[test], y_train[test],
	                         name='ROC fold{}'.format(i),
	                         alpha=0.3, lw=0.5, ax=ax)
		interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)

	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey', alpha=.8)
	
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	ax.plot(mean_fpr, mean_tpr, color='blue',
	        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
	        lw=2, alpha=.5)
	
	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
	                label=r'$\pm$ 1 std. dev.')
	                
	#plt.ylabel('Sensitivity')
	#plt.xlabel('1- Specifity')
	
	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
	       #title="ROC curve stratified sampling, sample number in training="+str(NT))
			title="ROC Curve")
	#ax.legend(loc="lower right", alpha=.2)
	plt.legend(loc="lower right", prop={'size':9})
	plt.show()
	
	print(score)
	print(average)
	
	#plt.savefig('CRPCV_allData'+str(n)+'.pdf')
	
	# Compute and print the confusion matrix and classification report
	#print('\n\nconfusion_matrix:')
	#print(confusion_matrix(y_test, y_pred, normalize='all'))
	#print('================================')
	
	#print('\n\nclassification_report:')
	#print(classification_report(y_test, y_pred))
	#print('================================')


	





