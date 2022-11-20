import sys
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import torch
from torch import nn
from torch.autograd import Variable
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter


from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import PandasTools

from sklearn.model_selection import train_test_split

from PyQt5 import QtWebEngineWidgets
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QUrl
from Model import *
from PandasModel import *
from DataBloc import *

writer = SummaryWriter()
class Linearregression(QtWidgets.QMainWindow):
    def __init__(self):
        super(Linearregression,self).__init__()
        
        self.ui= Ui_MainWindow()
        self.ui.setupUi(self)
        
        
        self.ui.DataButton.clicked.connect(self.tabData)
        self.ui.ModelButton.clicked.connect(self.tabModel)
        self.ui.EVALButton.clicked.connect(self.tabEval)
        
        
        self.ui.LoadButton.clicked.connect(self.LoadData)
        self.ui.radioButton.clicked.connect(self.DataSet)
        self.ui.SplitButton.clicked.connect(self.SplitData)
        self.ui.tableView_2.clicked.connect(self.DrawMolecule)
        self.ui.TrainButton.clicked.connect(self.Train)
        #self.ui.LossButton.clicked.connect(self.Loss)
        self.ui.ResultButton.clicked.connect(self.predictionplot)
        
        
    def tabData(self):
        self.ui.tabWidget.setCurrentIndex(0)
    def tabModel(self):
        self.ui.tabWidget.setCurrentIndex(1)
    def tabEval(self):
        self.ui.tabWidget.setCurrentIndex(2)
    
#Data Processing **********************************************************************************
    def LoadData(self):
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load file", "/home", "CSV (*.csv)")
        df = pd.read_csv(fileName)

        df['mol'] = df['SMILES sequence'].apply(lambda x: Chem.MolFromSmiles(x)) 
        df['mol'] = df['mol'].apply(lambda x: Chem.AddHs(x))
        df['num_of_atoms'] = df['mol'].apply(lambda x: x.GetNumAtoms())
        df['Num_of_bonds'] = df['mol'].apply(lambda x:x.GetNumBonds())
        df['num_of_heavy_atoms'] = df['mol'].apply(lambda x: x.GetNumHeavyAtoms())
        df['logP'] = df['mol'].apply(lambda x: Descriptors.MolLogP(x))
        df['ssr'] = df['mol'].apply(lambda x: Chem.GetSSSR(x))
        df['tpsa'] = df['mol'].apply(lambda x: Descriptors.TPSA(x))
        df['mol_exact_w'] = df['mol'].apply(lambda x: Descriptors.ExactMolWt(x))
        df['num_valence_electrons'] = df['mol'].apply(lambda x: Descriptors.NumValenceElectrons(x))
        df['num_heteroatoms'] = df['mol'].apply(lambda x: Descriptors.NumHeteroatoms(x))
        df['rotatable_bond'] = df['mol'].apply(lambda x: Descriptors.NumRotatableBonds(x))
        df['mol_mr'] = df['mol'].apply(lambda x: Descriptors.MolMR(x)) 
        df['mol_w'] = df['mol'].apply(lambda x: Descriptors.MolWt(x))
        df['NumAromaticCarbocycles'] = df['mol'].apply(lambda x: Descriptors.NumAromaticCarbocycles(x))
        df['NumAromaticHeterocycles'] = df['mol'].apply(lambda x: Descriptors.NumAromaticHeterocycles(x))
        df['NumAromaticRings'] = df['mol'].apply(lambda x: Descriptors.NumAromaticRings(x))
        df['NumHAcceptors'] = df['mol'].apply(lambda x: Descriptors.NumHAcceptors(x))
        df['NumHDonors'] = df['mol'].apply(lambda x: Descriptors.NumHDonors(x))
        df['NumSaturatedRings'] = df['mol'].apply(lambda x: Descriptors.NumSaturatedRings(x))
        df['NHOHCount'] = df['mol'].apply(lambda x: Descriptors.NHOHCount(x))
        df['NOCount'] = df['mol'].apply(lambda x: Descriptors.NOCount(x))
        df['HeavyAtomMolWt'] = df['mol'].apply(lambda x: Descriptors.HeavyAtomMolWt(x))
        df['FpDensityMorgan1'] = df['mol'].apply(lambda x: Descriptors.FpDensityMorgan1(x))
        df['FpDensityMorgan2'] = df['mol'].apply(lambda x: Descriptors.FpDensityMorgan2(x))
        df['FpDensityMorgan3'] = df['mol'].apply(lambda x: Descriptors.FpDensityMorgan3(x))
        df['fr_NH1'] = df['mol'].apply(lambda x: Descriptors.fr_NH1(x))
        df['fr_NH2'] = df['mol'].apply(lambda x: Descriptors.fr_NH2(x))
        df['fr_NH0'] = df['mol'].apply(lambda x: Descriptors.fr_NH0(x))
        df['fr_Al_COO'] = df['mol'].apply(lambda x: Descriptors.fr_Al_COO(x))
        df['fr_bicyclic'] = df['mol'].apply(lambda x: Descriptors.fr_bicyclic(x))
        df['MaxPartialCharge'] = df['mol'].apply(lambda x: Descriptors.MaxPartialCharge(x))
        df['HallKierAlpha'] = df['mol'].apply(lambda x: Descriptors.HallKierAlpha(x))
        df['RingCount'] = df['mol'].apply(lambda x: Descriptors.RingCount(x))
        df['BalabanJ'] = df['mol'].apply(lambda x: Descriptors.BalabanJ(x))

# Here we need a llop to cover all atoms of the molecule
        def number_of_atoms(atom_list, df):
            for i in atom_list:
               df['num_of_{}_atoms'.format(i)] = df['mol'].apply(lambda x: len(x.GetSubstructMatches(Chem.MolFromSmiles(i))))
        number_of_atoms(['C','O','N','S'], df)

#display in table
        
    def DataSet(self):
            V=df[['SMILES sequence']]
            model1=PandasModel(V)
            self.ui.tableView_2.setModel(model1)


        
            model=PandasModel(df)
            
            if self.ui.radioButton.isChecked():
               self.ui.tableView.setModel(model)
           
            
    def DrawMolecule(self):
        index=(self.ui.tableView_2.selectionModel().currentIndex())
        value=index.sibling(index.row(),index.column()).data()
        
        #print(value)
        molecule=Chem.MolFromSmiles(value)
        self.ui.label_4.setPixmap(Draw.MolToQPixmap(molecule))
    def SplitData(self):
            X = df.drop(columns=['SMILES sequence','mol', 'Binding Affinity'])
            y = df['Binding Affinity'].values
            
            percentage=(self.ui.spinBox.value())/100
            
            x_train,x_test,y_train,y_test = train_test_split(X, y, test_size=percentage, random_state=0)
            print(x_train)
            
            #transorm to numpy arrays*********************
            xn_train = np.array(x_train) 
            xn_test = np.array(x_test) 
            yn_train = np.array(y_train) 
            yn_test = np.array(y_test)
            xS_train=xn_train.reshape(-1,1).astype('float32')
            print(xS_train.shape)
            #print(xn_test.shape)
            
            #*********************************************
            
            #transfering numpy array to torch tensors for training later****************
            xt_train= torch.from_numpy(xn_train).float() 
            yt_train= torch.from_numpy(yn_train).float() 
            
            #print(xt_train.shape)
            xt_test= torch.from_numpy(xn_test).float() 
            yt_test= torch.from_numpy(yn_test).float()
            
            #print(xt_test.shape)
            self.x1=xt_train
            self.x2=xt_test
            
            self.y1=yt_train
            self.y2=yt_test
            #****************************************************************************
            self.ui.xTrainButton.clicked.connect(self.Xtrain)
    def Xtrain(self):
            model=PandasModel(x_train)
            if self.ui.xTrainButton.isChecked():
                self.ui.tableView.setModel(model)
            self.ui.xTestButton.clicked.connect(self.Xtest)
    def Xtest(self):
        model=PandasModel(x_test)
        if self.ui.xTestButton.isChecked():
            self.ui.tableView.setModel(model)
    
    #TRAINING MODULE*************************************************************************
    def Train(self):
        class linearRegression(nn.Module):
            def __init__(self):
                super(linearRegression, self).__init__()
                self.linear= nn.Linear(input_Dim, H)
                self.relu=torch.nn.ReLU()
                self.linear2=nn.Linear(H,output_Dim)
            
            def forward(self, xn_train):
                out=self.linear2(self.relu(self.linear(xn_train)))
                return out 
        H=2000
        #Input Size
        len(self.ui.lineEdit.text())
        input_Dim=int(self.ui.lineEdit.text())
        #Output Size
        len(self.ui.lineEdit_2.text())
        output_Dim=int(self.ui.lineEdit_2.text())
        #Learning Rate
        len(self.ui.lineEdit_4.text())
        learning_rate=float(self.ui.lineEdit_4.text())
        #Epochs
        len(self.ui.lineEdit_3.text())
        epochs=int(self.ui.lineEdit_3.text())
        
    
        model= linearRegression()
        self.lrmodel=model
    
        
        #loss functions**************************************************** 
        if str(self.ui.comboBox.currentText())=="MSE":
            criterion= nn.MSELoss()
        
        if str(self.ui.comboBox.currentText())=="MAE":
            criterion= nn.L1Loss()
        
        if str(self.ui.comboBox.currentText())=="Huber":
            criterion= nn.SmoothL1Loss() 
        
        #if str(self.ui.ComboBox.currentText())=="Log-Cosh":
            #criterion= nn. 
        
        #if str(self.ui.ComboBox.currentText())=="Quentile":
            #criterion= nn.
        
        #optimizers ******************************************************
        if str(self.ui.comboBox_2.currentText())=="SGD":
            optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
        
        if str(self.ui.comboBox_2.currentText())=="Adam":
            optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        
        if str(self.ui.comboBox_2.currentText())=="Momentum":
            optimizer = torch.optim.SGD(model.parameters(),momentum=0 , lr=learning_rate)
        
        if str(self.ui.comboBox_2.currentText())=="Adagrad":
            optimizer = torch.optim.Adagrad(model.parameters(), lr=learning_rate)
     
     #Start training
        num_epochs = epochs
        for epoch in range(num_epochs):
            inputs = self.x1
            target= self.y1
         
         #Forward
            out=model(inputs)
            loss=criterion(out,target)
            
            
            writer.add_scalar("Loss/train", loss, epoch)
        
         #Backward
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
         
            print(f'Epoch[{epoch+1}/{num_epochs}], loss: {loss.item():.6f}')
      
            
     #evalu
        model.eval()
        with torch.no_grad():
             predict = model(self.x1)
             y_pred=model(self.x2)
        predict = predict.data.numpy()
        y_pred=y_pred.data.numpy()
        self.pred=predict
        self.y_pred=y_pred
        
        writer.flush()
        self.ui.webEngineView.load(QUrl('http://localhost:6006'))
   
    
    def predictionplot(self):
     
        self.ui.PlotWidget.axes.plot(self.x1.numpy(), self.y1.numpy(),'ro')
        self.ui.PlotWidget.axes.plot(self.x1.numpy(),self.pred,'bo')
        self.ui.PlotWidget.axes.legend(('original data' ,'predicted'),loc='upper right')
        self.ui.PlotWidget.axes.set_title('Model predction')
        
        

        self.ui.PlotWidget.ax.scatter(self.y2,self.y_pred)
        self.ui.PlotWidget.ax.plot([self.y2.min(),self.y2.max()],[self.y2.min(),self.y2.max()],'k--',lw=4)
        self.ui.PlotWidget.draw()

        
        
        
if __name__=="__main__":
    app=QtWidgets.QApplication(sys.argv)
    MainWindow=QtWidgets.QMainWindow()
    MYsmn=Linearregression()
    MYsmn.show()
    sys.exit(app.exec_())
    