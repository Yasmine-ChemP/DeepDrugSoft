from PyQt5 import QtWidgets, uic
from pyqtgraph import PlotWidget, plot
import numpy as np
import sys


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        #instead of converting .ui to .py , loading ui page directly (alternative)
        uic.loadUi('mainwindow.ui', self)
        x=np.linspace(-30,30,10)
        y=x**2
        self.plot(x,y)
    
    def plot(self, x, y):
         self.graphWidget.plot(x,y)



    
def main():
    app = QtWidgets.QApplication(sys.argv)
    grph = MainWindow()
    grph.show()
    sys.exit(app.exec_())

if __name__ == '__main__':      
    main()