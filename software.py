import sys
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import torch
from torch import nn
from torch.autograd import Variable
from torch.utils.data import DataLoader,TensorDataset
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
from sklearn.metrics import r2_score

from PyQt5 import QtWebEngineWidgets
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QUrl
from plotcanvas import PlotCanvas
from Model import *
from PandasModel import *
from DataBloc import *

writer = SummaryWriter()
class APP(QtWidgets.QMainWindow):
    def __init__(self):
        super(APP,self).__init__()
        
        self.ui= Ui_MainWindow()
        self.ui.setupUi(self)
                
        #Navigation Buttons 
        self.ui.DataButton.clicked.connect(self.tabData)
        self.ui.ModelButton.clicked.connect(self.tabModel)
        self.ui.EVALButton.clicked.connect(self.tabEval)
        
        #Data related Buttons 
        self.ui.LoadButton.clicked.connect(self.LoadData)
        self.ui.radioButton.clicked.connect(self.DataSet)
        self.ui.tableView_2.clicked.connect(self.DrawMolecule)
        self.ui.SplitButton.clicked.connect(self.SplitData)
        self.ui.xTrainButton.clicked.connect(self.Xtrain)
        self.ui.xTestButton.clicked.connect(self.Xtest)
        
        
        #Model related Buttons
        self.ui.TrainButton.clicked.connect(self.Train)
        #self.ui.LossButton.clicked.connect(self.Loss)
        self.ui.ResultButton.clicked.connect(self.predictionplot)
        
#Move between Tabs******************************************************************      
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

        def number_of_atoms(atom_list, df):
            for i in atom_list:
               df['num_of_{}_atoms'.format(i)] = df['mol'].apply(lambda x: len(x.GetSubstructMatches(Chem.MolFromSmiles(i))))
        number_of_atoms(['C','O','N','S'], df)
        
        
        #Display all Data in table
    def DataSet(self):
        M=df.drop(columns=['mol'])
        model=PandasModel(M)
        if self.ui.radioButton.isChecked():
               self.ui.tableView.setModel(model)
               
               self.ui.lineEdit_5.setText(str(int(model.rowCount())))
               self.ui.lineEdit_6.setText(str(int(model.columnCount())))
        
        #In second table , display only SMILES Sequences 
        V=df[['SMILES sequence']]
        model2=PandasModel(V)
        self.ui.tableView_2.setModel(model2)
        
        #Drawing the molecules using Qpixmap
    def DrawMolecule(self):
        index=(self.ui.tableView_2.selectionModel().currentIndex())
        value=index.sibling(index.row(),index.column()).data()
        molecule=Chem.MolFromSmiles(value)
        
        self.ui.label_4.setPixmap(Draw.MolToQPixmap(molecule))

        #Splitting the data into two subsets (test and train)
    def SplitData(self):
        X = df.drop(columns=['SMILES sequence','mol', 'Binding Affinity'])
        y = df['Binding Affinity'].values
            
        percentage=(self.ui.spinBox.value())/100
            
        x_train,x_test,y_train,y_test = train_test_split(X, y, test_size=percentage, random_state=0)
        print(x_train)
        
        #X and Y train 
         #To numpy arrays 
        xn_train = np.array(x_train)
        
        yn_train=np.array(y_train)
        
        self.x_T=xn_train
        self.y_T=yn_train
        
        print(xn_train , yn_train)
        

         #To tensors 
        xt_train=torch.tensor(xn_train).float()
        yt_train=torch.tensor(yn_train).float()
        
        self.x1=xt_train
        self.y1=yt_train
        
         #DataLoader for Train 
        train_dataset = TensorDataset(xt_train, yt_train)
        train_loader = torch.utils.data.DataLoader(dataset=train_dataset, batch_size=256, shuffle=True)
        self.train_loader=train_loader
        
        #X and Y test 
         #To numpy arrays 
        xn_test = np.array(x_test)
        yn_test = np.array(y_test)
        
        self.x_TS=xn_test
        self.y_TS=xn_test
         #to tensors
        xt_test=torch.tensor(xn_test).float()
        yt_test=torch.tensor(yn_test).float()
        
        self.x2=xt_test
        self.y2=yt_test
        
        
         #DataLoader for Test
        test_dataset = TensorDataset(xt_test, yt_test)
        test_loader = torch.utils.data.DataLoader(dataset=test_dataset,batch_size=256, shuffle=False)
        self.test_loader=test_loader
        #Displaying training and testing sets in table
    def Xtrain(self):
        model=PandasModel(x_train)
        if self.ui.xTrainButton.isChecked():
           self.ui.tableView.setModel(model)
           
           self.ui.lineEdit_5.setText(str(int(model.rowCount())))
           self.ui.lineEdit_6.setText(str(int(model.columnCount())))
    def Xtest(self):
        model=PandasModel(x_test)
        if self.ui.xTestButton.isChecked():
           self.ui.tableView.setModel(model)
           
           self.ui.lineEdit_5.setText(str(int(model.rowCount())))
           self.ui.lineEdit_6.setText(str(int(model.columnCount())))
    
#Model Setup***************************************************************************
    def Train(self):
        class Network(nn.Module):
            def __init__(self,input_Dim,H,Dropout_rate,output_Dim):
                super(Network, self).__init__()
                self.linear= nn.Linear(input_Dim, H)
                self.h1=nn.Linear(H, H)
                self.h2=nn.Linear(H,H)
                self.linear_2=nn.Linear(H,output_Dim)
                
                #Layer Norm
                self.ln1=nn.LayerNorm(H,H)
                self.ln2=nn.LayerNorm(H,H)
                self.ln3=nn.LayerNorm(H,H)
                
                #Activation Function
                self.activation=nn.ReLU()
                
                #DroupOUT
                self.dropout=nn.Dropout(Dropout_rate)
            
            def forward(self, x):
                out= self.linear(x)
                out=self.ln1(out)
                out=self.activation(out)
                out=self.dropout(out)
                out=self.h1(out)
                out=self.ln2(out)
                out=self.activation(out)
                out=self.dropout(out)
                out=self.h2(out)
                out=self.ln3(out)
                out=self.activation(out)
                out=self.dropout(out)
                out=self.linear_2(out)
                return out 
    #parameters :    
        ##Hidden Size
        H=512
        Dropout_rate=0.05
        
        ##Input Size
        #len(self.ui.lineEdit.text())
        input_Dim=int(self.ui.lineEdit.text())
        ##Output Size
        len(self.ui.lineEdit_2.text())
        output_Dim=int(self.ui.lineEdit_2.text())
        ##Learning Rate
        len(self.ui.lineEdit_4.text())
        learning_rate=float(self.ui.lineEdit_4.text())
        ##Epochs
        len(self.ui.lineEdit_3.text())
        epochs=int(self.ui.lineEdit_3.text())
    #set model
        model= Network(input_Dim,H,Dropout_rate,output_Dim)
        self.lrmodel=model
        
    #Choice of Performance Function and optimizer 
        ##loss functions 
        if str(self.ui.comboBox.currentText())=="MSE":
            criterion= nn.MSELoss()
        
        if str(self.ui.comboBox.currentText())=="MAE":
            criterion= nn.L1Loss()
        
        if str(self.ui.comboBox.currentText())=="Huber":
            criterion= nn.SmoothL1Loss()
        
        ##optimizers 
        if str(self.ui.comboBox_2.currentText())=="SGD":
            optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
        
        if str(self.ui.comboBox_2.currentText())=="Adam":
            optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        
        if str(self.ui.comboBox_2.currentText())=="Momentum":
            optimizer = torch.optim.SGD(model.parameters(),momentum=0 , lr=learning_rate)
        
        if str(self.ui.comboBox_2.currentText())=="Adagrad":
            optimizer = torch.optim.Adagrad(model.parameters(), lr=learning_rate)
        
    #Training :
        model.train()
        num_epochs = epochs
        for epoch in range(num_epochs):
            running_loss=0
            for fps, Labels in self.train_loader:
                optimizer.zero_grad()

        #Forward
                output=model(fps)
                loss=criterion(output,Labels)
        #Backward
            
                loss.backward()
                optimizer.step()
                
                running_loss += loss.item()
                
                
        #///////////////////////////////////////////
                writer.add_scalar("Loss/train", running_loss, epoch)
        #///////////////////////////////////////////
            else:
                test_loss= torch.mean((self.y2 - model(self.x2))**2).item()
        
        
        #///////////////////////////////////////////
                writer.add_scalar("Loss/test", test_loss, epoch)
        #///////////////////////////////////////////

         #print Losses
                #print(f'Epoch[{epoch+1}/{num_epochs}], train loss: {(running_loss/len(self.train_loader)):.6f}, test loss : {test_loss:.6f}')
                print("Epoch: %3i Training loss: %0.2F Test loss: %0.2F"%(epoch,(running_loss/len(self.train_loader)), test_loss))
        
       #Evaluation
        model.eval()
        with torch.no_grad():
             ypred_train = model(self.x1)
             ypred_test = model(self.x2)
             
        ypred_train = ypred_train.detach().numpy()
        ypred_test = ypred_test.detach().numpy()
        
        self.pred_train=ypred_train
        self.pred_test= ypred_test
        #Display Loss in TensorBoard
        writer.flush()
        self.ui.webEngineView.load(QUrl('http://localhost:6006'))
        
        
    #Results :
    def predictionplot(self):
     
        self.ui.PlotWidget.axes.plot(self.y_T, self.ypred_train,'ro')
        self.ui.PlotWidget.axes.plot(self.y_TS,self.ypred_test,'bo')
        self.ui.PlotWidget.axes.plot([self.y_TS.min(), self.y_TS.max()], [self.y_TS.min(), self.y_TS.max()], 'k--', lw=3)
        self.ui.PlotWidget.axes.legend(('original data' ,'predicted'),loc='upper right')
        self.ui.PlotWidget.axes.set_title('Model predction')
        self.ui.PlotWidget.draw()
        
        
        
        #plot with RDKIT 
        self.ui.rdkitPlotwidget.axes.plot(self.pred_test[:200],"red",label="prediction",linewidth=2.0)
        self.ui.rdkitPlotwidget.axes.plot(self.pred_train[:200],"green",label="Measurment",linewidth=2.0)
        self.ui.rdkitPlotwidget.draw()
        

        
if __name__=="__main__":
    app=QtWidgets.QApplication(sys.argv)
    MainWindow=QtWidgets.QMainWindow()
    MYsmn=APP()
    MYsmn.show()
    sys.exit(app.exec_())