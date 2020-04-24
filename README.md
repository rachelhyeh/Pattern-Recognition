# Pattern-Recognition
Computer Project 2

### Data
The data come from the following the following cancer classification study: <br/>
**van de Vijver,M.J., He,Y.D., van't Veer,L.J., et al. (2002), "A gene-expression signature as a predictor of survival in breast cancer." New Eng. J. Med., 347, 1999-2009.** <br/>
This paper analyzes gene expression in breast tumor biopsies from 295 patients. The authors performed feature selection to obtain 70 genes; hence, the full data matrix is 70x295. This is a retrospective study, meaning that the patients were tracked over the years and their outcomes recorded. Using this clinical information, the authors labeled the tumor samples into two classes: the "good' prognosis" group were disease-free for at least five years after first treatment, whereas the "bad prognosis" group developed distant methatasis within the first five years. Of the 295 patients, 216 belong to the "good-prognosis" class, whereas the remaining 79 belong to the "poor-
prognosis" class.


### Goal
Search for gene feature sets and design classifiers that best discriminate the two prognosis classes on the training data, and use the testing data to determine their accuracy.
1. Classification rules:
   - LDA, p = 0.75
   - Linear SVM, C = 1
   - Nonlinear SVM with Caussian RBF kernel, C = 1
   - NN with 5 neurons in on hidden layer
2. Feature selection methods:
   - Top 2 genes (exhaustive search)
   - Top 3-5 genes (sequential forward search)
   - All genes (no feature selection)
   
### Files
- Description: Computer_Project_2
- Training data: Training_Data.txt
- Testing data: Testing_Data.txt
- Code: CP_2.py
