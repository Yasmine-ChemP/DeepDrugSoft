#  Matplotlib consists of 3 different layers : "Backend" as the bottom layer , "Artist" and  "Scripting" as the top layers 
# Backend layer : FigureCanvas is the space (canvas) used to draw the graph
# Artist layer : used to display the title, the axis, the labels ...
# Script (pyplot) layer : handles the simplicity of the script and simplifies the usage of Matplotlib module. 

from PyQt5.QtWidgets import*                                                  # all components belonging to the QtWidgets class are added
from matplotlib.backends.backend_qt5agg import FigureCanvas    #   We add FigureCanvas  from Matplotlib module
from matplotlib.figure import Figure


# we define the class named MplWidget with the same name as the Widget created from the QWidget 
# and Promoted in Qt Designer

class MplWidget(QWidget):
    
    def __init__(self, parent = None):

        QWidget.__init__(self, parent)        
        self.canvas = FigureCanvas(Figure())
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)        
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.setLayout(vertical_layout)
