
import sys                                  #  A module assuring the interaction between Python script and the system. We have use it in sys.argv and sys.exit  
from PyQt5 import QtCore, QtGui, QtWidgets  # importing all the modules required to create a GUI into the current namespace. The QtWidgets module contains all the major widgets used. 
from fit import Ui_MainWindow               # importing our generated file. Using pyuic5 to convert the XML file (*.ui) into Python file (*.py) is better than using loadUi as you did, because compiling the XML file every time reduces the speed of the script 
import pandas as pd                         # importing Pandas to deal with loaded files
import numpy as np                          # importing the module numpy
from PandasModel import *                   # In the case of QTableView the data must be provided through a model.in the case of pandas there is no default model but we can create a custom one as shown in the file PandasModel.py 
                     #  Importing Scientific Graphics and GUI Library for Python and named pg
from DataBlocByMe import *
 
class Mytest(QtWidgets.QMainWindow):        # Creating a new Mytest class that inherits from the base class Mainwindow.
    def __init__(self):
        super(Mytest, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        

        self.ui.pushButton.clicked.connect(self.loadfile)     # a clicked() event connected to the method loadfile
        self.ui.pushButton_2.clicked.connect(self.plotPQG)    # a clicked() event connected to the method fitgraf
        self.ui.pushButton_4.clicked.connect(self.plotMPL)    # a clicked() event connected to the method fitmatp
        
        
        
    def loadfile(self, filename):
        #fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load file", "/home", "XLSX (*.xlsx)")  # Browsing the home drive (c:) and Loading an xlsx file
        
        #df = pd.read_excel(fileName)          # Loading the data from a xlsx file to a Pandas DataFrame named df
        #model = PandasModel(df)               # The custom Pandas model (see the file PandasModel.py) required to display the data stored into the DataFrame in a TableView
        #self.ui.tableView.setModel(model)     # Displaying the data in a tableView
        
        
        self.x1=df.x     # My test data loaded from the xlsx file, declared as "Global Variable" using "self" (self.x1) so we can use it into other functions or methods (in our case fitgraf )
        self.y1=df.y     # //
        

    def plotPQG(self):      # Plotting the data using pyqtgraph
       
        self.ui.graphwidget.plot(self.x1, self.y1, pen=(0,0,255))  
        

   
    def plotMPL(self):   # Plotting the data using MatLibPlot / it requires the file mplwidget.py  : 
                                  #  That  file is located into the actual directory and  must take  the same name as the widget (QWidget) that is promoted in Qt designer. 
        
        self.ui.MplWidget.canvas.axes.plot(self.x1,self.y1)
        self.ui.MplWidget.canvas.draw()
         
         
         
         
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)     #  an instance of QApplication, sys.argv is the list of command-line parameters that we can pass to the application when launching it through the shell  
    MainWindow = QtWidgets.QMainWindow()
    myapp = Mytest()                           #  An instance of the Mytest class is created with the name myapp
    myapp.show()                               #  Show our application's GUI
    sys.exit(app.exec_())                      # Passes the control over to Qt which will exit the application only when the user closes it from the GUI and then will release memory resources
    
    
 