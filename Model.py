# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ModelApp.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(997, 862)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 130, 951, 631))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(9)
        font.setBold(True)
        font.setWeight(75)
        self.tabWidget.setFont(font)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.tabWidget_2 = QtWidgets.QTabWidget(self.tab)
        self.tabWidget_2.setGeometry(QtCore.QRect(20, 10, 891, 581))
        self.tabWidget_2.setObjectName("tabWidget_2")
        self.tab_4 = QtWidgets.QWidget()
        self.tab_4.setObjectName("tab_4")
        self.LoadButton = QtWidgets.QPushButton(self.tab_4)
        self.LoadButton.setGeometry(QtCore.QRect(70, 30, 151, 41))
        self.LoadButton.setObjectName("LoadButton")
        self.radioButton = QtWidgets.QRadioButton(self.tab_4)
        self.radioButton.setGeometry(QtCore.QRect(420, 20, 95, 20))
        self.radioButton.setObjectName("radioButton")
        self.xTrainButton = QtWidgets.QRadioButton(self.tab_4)
        self.xTrainButton.setGeometry(QtCore.QRect(550, 20, 95, 20))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(8)
        font.setBold(True)
        font.setWeight(75)
        self.xTrainButton.setFont(font)
        self.xTrainButton.setObjectName("xTrainButton")
        self.yTrainButton = QtWidgets.QRadioButton(self.tab_4)
        self.yTrainButton.setGeometry(QtCore.QRect(680, 20, 95, 20))
        font = QtGui.QFont()
        font.setPointSize(8)
        self.yTrainButton.setFont(font)
        self.yTrainButton.setObjectName("yTrainButton")
        self.xTestButton = QtWidgets.QRadioButton(self.tab_4)
        self.xTestButton.setGeometry(QtCore.QRect(550, 50, 95, 20))
        font = QtGui.QFont()
        font.setPointSize(8)
        self.xTestButton.setFont(font)
        self.xTestButton.setObjectName("xTestButton")
        self.yTestButton = QtWidgets.QRadioButton(self.tab_4)
        self.yTestButton.setGeometry(QtCore.QRect(680, 50, 95, 20))
        font = QtGui.QFont()
        font.setPointSize(8)
        self.yTestButton.setFont(font)
        self.yTestButton.setObjectName("yTestButton")
        self.tableView = QtWidgets.QTableView(self.tab_4)
        self.tableView.setGeometry(QtCore.QRect(40, 80, 761, 401))
        self.tableView.setObjectName("tableView")
        self.SplitButton = QtWidgets.QPushButton(self.tab_4)
        self.SplitButton.setGeometry(QtCore.QRect(40, 487, 151, 51))
        self.SplitButton.setObjectName("SplitButton")
        self.spinBox = QtWidgets.QSpinBox(self.tab_4)
        self.spinBox.setGeometry(QtCore.QRect(210, 500, 61, 31))
        self.spinBox.setMaximum(100)
        self.spinBox.setObjectName("spinBox")
        self.label = QtWidgets.QLabel(self.tab_4)
        self.label.setGeometry(QtCore.QRect(280, 510, 55, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.tab_4)
        self.label_2.setGeometry(QtCore.QRect(380, 510, 111, 21))
        self.label_2.setObjectName("label_2")
        self.lineEdit_5 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_5.setGeometry(QtCore.QRect(492, 510, 131, 22))
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.label_3 = QtWidgets.QLabel(self.tab_4)
        self.label_3.setGeometry(QtCore.QRect(630, 510, 16, 16))
        font = QtGui.QFont()
        font.setFamily("Courier New")
        font.setBold(False)
        font.setWeight(50)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.lineEdit_6 = QtWidgets.QLineEdit(self.tab_4)
        self.lineEdit_6.setGeometry(QtCore.QRect(650, 510, 131, 22))
        self.lineEdit_6.setObjectName("lineEdit_6")
        self.tabWidget_2.addTab(self.tab_4, "")
        self.tab_5 = QtWidgets.QWidget()
        self.tab_5.setObjectName("tab_5")
        self.tableView_2 = QtWidgets.QTableView(self.tab_5)
        self.tableView_2.setGeometry(QtCore.QRect(30, 30, 301, 491))
        self.tableView_2.setObjectName("tableView_2")
        self.label_4 = QtWidgets.QLabel(self.tab_5)
        self.label_4.setGeometry(QtCore.QRect(390, 100, 451, 381))
        self.label_4.setText("")
        self.label_4.setObjectName("label_4")
        self.tabWidget_2.addTab(self.tab_5, "")
        self.tab_6 = QtWidgets.QWidget()
        self.tab_6.setObjectName("tab_6")
        self.rdkitPlotwidget = PlotCanvas(self.tab_6)
        self.rdkitPlotwidget.setGeometry(QtCore.QRect(50, 30, 761, 481))
        self.rdkitPlotwidget.setObjectName("rdkitPlotwidget")
        self.tabWidget_2.addTab(self.tab_6, "")
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.groupBox = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox.setGeometry(QtCore.QRect(20, 30, 391, 451))
        self.groupBox.setObjectName("groupBox")
        self.lineEdit = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit.setGeometry(QtCore.QRect(190, 100, 131, 22))
        self.lineEdit.setObjectName("lineEdit")
        self.lineEdit_2 = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_2.setGeometry(QtCore.QRect(190, 170, 131, 22))
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.lineEdit_4 = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_4.setGeometry(QtCore.QRect(190, 310, 131, 22))
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.lineEdit_3 = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_3.setGeometry(QtCore.QRect(190, 380, 131, 22))
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.label_5 = QtWidgets.QLabel(self.groupBox)
        self.label_5.setGeometry(QtCore.QRect(50, 100, 111, 16))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.groupBox)
        self.label_6.setGeometry(QtCore.QRect(50, 170, 111, 16))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.groupBox)
        self.label_7.setGeometry(QtCore.QRect(10, 310, 161, 16))
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(self.groupBox)
        self.label_8.setGeometry(QtCore.QRect(10, 380, 181, 16))
        self.label_8.setObjectName("label_8")
        self.groupBox_2 = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox_2.setGeometry(QtCore.QRect(520, 30, 321, 331))
        self.groupBox_2.setObjectName("groupBox_2")
        self.comboBox = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox.setGeometry(QtCore.QRect(190, 80, 111, 31))
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox_2 = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox_2.setGeometry(QtCore.QRect(190, 210, 111, 31))
        self.comboBox_2.setObjectName("comboBox_2")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.label_9 = QtWidgets.QLabel(self.groupBox_2)
        self.label_9.setGeometry(QtCore.QRect(30, 90, 131, 16))
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(self.groupBox_2)
        self.label_10.setGeometry(QtCore.QRect(30, 220, 131, 16))
        self.label_10.setObjectName("label_10")
        self.TrainButton = QtWidgets.QPushButton(self.tab_2)
        self.TrainButton.setGeometry(QtCore.QRect(80, 510, 141, 51))
        self.TrainButton.setObjectName("TrainButton")
        self.LossButton = QtWidgets.QPushButton(self.tab_2)
        self.LossButton.setGeometry(QtCore.QRect(370, 510, 141, 51))
        self.LossButton.setObjectName("LossButton")
        self.ResultButton = QtWidgets.QPushButton(self.tab_2)
        self.ResultButton.setGeometry(QtCore.QRect(640, 510, 141, 51))
        self.ResultButton.setObjectName("ResultButton")
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.tabWidget_3 = QtWidgets.QTabWidget(self.tab_3)
        self.tabWidget_3.setGeometry(QtCore.QRect(20, 10, 891, 581))
        self.tabWidget_3.setObjectName("tabWidget_3")
        self.tab_7 = QtWidgets.QWidget()
        self.tab_7.setObjectName("tab_7")
        self.webEngineView = QtWebEngineWidgets.QWebEngineView(self.tab_7)
        self.webEngineView.setGeometry(QtCore.QRect(20, 20, 841, 511))
        self.webEngineView.setUrl(QtCore.QUrl("about:blank"))
        self.webEngineView.setObjectName("webEngineView")
        self.tabWidget_3.addTab(self.tab_7, "")
        self.tab_8 = QtWidgets.QWidget()
        self.tab_8.setObjectName("tab_8")
        self.PlotWidget = PlotCanvas(self.tab_8)
        self.PlotWidget.setGeometry(QtCore.QRect(20, 20, 841, 511))
        self.PlotWidget.setObjectName("PlotWidget")
        self.tabWidget_3.addTab(self.tab_8, "")
        self.tabWidget.addTab(self.tab_3, "")
        self.DataButton = QtWidgets.QPushButton(self.centralwidget)
        self.DataButton.setGeometry(QtCore.QRect(80, 10, 151, 111))
        self.DataButton.setText("")
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../../Pictures/Big-Data-blog1-16.9-1.jpg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.DataButton.setIcon(icon)
        self.DataButton.setIconSize(QtCore.QSize(350, 108))
        self.DataButton.setObjectName("DataButton")
        self.ModelButton = QtWidgets.QPushButton(self.centralwidget)
        self.ModelButton.setGeometry(QtCore.QRect(420, 10, 161, 111))
        self.ModelButton.setStyleSheet("")
        self.ModelButton.setText("")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("../../Pictures/nn.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.ModelButton.setIcon(icon1)
        self.ModelButton.setIconSize(QtCore.QSize(350, 108))
        self.ModelButton.setObjectName("ModelButton")
        self.EVALButton = QtWidgets.QPushButton(self.centralwidget)
        self.EVALButton.setGeometry(QtCore.QRect(780, 10, 151, 111))
        self.EVALButton.setText("")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("../../Pictures/eval.jpg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.EVALButton.setIcon(icon2)
        self.EVALButton.setIconSize(QtCore.QSize(350, 108))
        self.EVALButton.setObjectName("EVALButton")
        self.ExitButton = QtWidgets.QPushButton(self.centralwidget)
        self.ExitButton.setGeometry(QtCore.QRect(810, 770, 141, 51))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.ExitButton.setFont(font)
        self.ExitButton.setObjectName("ExitButton")
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(0)
        self.tabWidget_3.setCurrentIndex(1)
        self.ExitButton.clicked.connect(MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.LoadButton.setText(_translate("MainWindow", "Load Data"))
        self.radioButton.setText(_translate("MainWindow", "DataSet"))
        self.xTrainButton.setText(_translate("MainWindow", "X Train Data"))
        self.yTrainButton.setText(_translate("MainWindow", "Y Train Data"))
        self.xTestButton.setText(_translate("MainWindow", "X Test Data"))
        self.yTestButton.setText(_translate("MainWindow", "Y Test Data"))
        self.SplitButton.setText(_translate("MainWindow", "SPLIT DATA"))
        self.label.setText(_translate("MainWindow", "%"))
        self.label_2.setText(_translate("MainWindow", "Size Of Table :"))
        self.label_3.setText(_translate("MainWindow", "X"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_4), _translate("MainWindow", "Tabled Data"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_5), _translate("MainWindow", "Species Visualization"))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_6), _translate("MainWindow", "Plots with RdKit"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Data Processing"))
        self.groupBox.setTitle(_translate("MainWindow", "HYPERPARAMETERS :"))
        self.label_5.setText(_translate("MainWindow", "INPUT Size :"))
        self.label_6.setText(_translate("MainWindow", "OUTPUT Size :"))
        self.label_7.setText(_translate("MainWindow", "LEARNING RATE :"))
        self.label_8.setText(_translate("MainWindow", "NUMBER OF EPOCHS :"))
        self.groupBox_2.setTitle(_translate("MainWindow", "LEARNING FUNCTIONS :"))
        self.comboBox.setItemText(0, _translate("MainWindow", "MSE"))
        self.comboBox.setItemText(1, _translate("MainWindow", "MAE"))
        self.comboBox.setItemText(2, _translate("MainWindow", "Huber"))
        self.comboBox.setItemText(3, _translate("MainWindow", "Log-Cosh"))
        self.comboBox.setItemText(4, _translate("MainWindow", "Quantile"))
        self.comboBox_2.setItemText(0, _translate("MainWindow", "SGD"))
        self.comboBox_2.setItemText(1, _translate("MainWindow", "Adam"))
        self.comboBox_2.setItemText(2, _translate("MainWindow", "Momentum"))
        self.comboBox_2.setItemText(3, _translate("MainWindow", "Adagrad"))
        self.label_9.setText(_translate("MainWindow", "LOSS FUNCTION :"))
        self.label_10.setText(_translate("MainWindow", "OPTIMIZER :"))
        self.TrainButton.setText(_translate("MainWindow", "TRAIN"))
        self.LossButton.setText(_translate("MainWindow", "LOSS"))
        self.ResultButton.setText(_translate("MainWindow", "RESULTS"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Model Setup"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_7), _translate("MainWindow", "LOSS VISUALIZATION"))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_8), _translate("MainWindow", "PREDICTIONS"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindow", "Model Evaluation"))
        self.ExitButton.setText(_translate("MainWindow", "Exit"))
from PyQt5 import QtWebEngineWidgets
from plotcanvas import PlotCanvas


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
