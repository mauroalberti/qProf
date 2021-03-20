from qgis.PyQt.QtWidgets import *

from ..qt_utils.filesystem import *


class FigureExportDetailedDlg(QDialog):

    def __init__(self, plugin_name, export_params, parent=None):

        super(FigureExportDetailedDlg, self).__init__(parent)

        self.sPluginName = plugin_name
        self.dExportParams = export_params

        qlyMainLayout = QVBoxLayout()

        # main parameters groupbox

        qgbMainParams = QGroupBox("Main graphic parameters")

        qlyMainParams = QGridLayout()

        qlyMainParams.addWidget(QLabel(self.tr("Figure width (inches)")), 0, 0, 1, 1)
        self.qleFigWidthInch = QLineEdit(self.dExportParams["expfig_width_inch"])
        qlyMainParams.addWidget(self.qleFigWidthInch, 0, 1, 1, 1)

        qlyMainParams.addWidget(QLabel(self.tr("Resolution (dpi)")), 0, 2, 1, 1)
        self.qleFigResolutionDpi = QLineEdit(self.dExportParams["expfig_res_dpi"])
        qlyMainParams.addWidget(self.qleFigResolutionDpi, 0, 3, 1, 1)

        qlyMainParams.addWidget(QLabel(self.tr("Font size (pts)")), 0, 4, 1, 1)
        self.qleFigFontSizePts = QLineEdit(self.dExportParams["expfig_font_size_pts"])
        qlyMainParams.addWidget(self.qleFigFontSizePts, 0, 5, 1, 1)

        qgbMainParams.setLayout(qlyMainParams)

        qlyMainLayout.addWidget(qgbMainParams)

        # additional parameters groupbox

        qgbAddParams = QGroupBox(self.tr("Subplot configuration tools parameters"))

        qlyAddParams = QGridLayout()

        qlyAddParams.addWidget(QLabel("Top space"), 0, 2, 1, 1)
        self.qsbTopSpaceValue = QDoubleSpinBox()
        self.qsbTopSpaceValue.setRange(0.0, 1.0)
        self.qsbTopSpaceValue.setDecimals(2)
        self.qsbTopSpaceValue.setSingleStep(0.01)
        self.qsbTopSpaceValue.setValue(0.96)
        qlyAddParams.addWidget(self.qsbTopSpaceValue, 0, 3, 1, 1)

        qlyAddParams.addWidget(QLabel("Left space"), 1, 0, 1, 1)
        self.qsbLeftSpaceValue = QDoubleSpinBox()
        self.qsbLeftSpaceValue.setRange(0.0, 1.0)
        self.qsbLeftSpaceValue.setDecimals(2)
        self.qsbLeftSpaceValue.setSingleStep(0.01)
        self.qsbLeftSpaceValue.setValue(0.1)
        qlyAddParams.addWidget(self.qsbLeftSpaceValue, 1, 1, 1, 1)

        qlyAddParams.addWidget(QLabel("Right space"), 1, 4, 1, 1)
        self.qsbRightSpaceValue = QDoubleSpinBox()
        self.qsbRightSpaceValue.setRange(0.0, 1.0)
        self.qsbRightSpaceValue.setDecimals(2)
        self.qsbRightSpaceValue.setSingleStep(0.01)
        self.qsbRightSpaceValue.setValue(0.96)
        qlyAddParams.addWidget(self.qsbRightSpaceValue, 1, 5, 1, 1)

        qlyAddParams.addWidget(QLabel("Bottom space"), 2, 2, 1, 1)
        self.qsbBottomSpaceValue = QDoubleSpinBox()
        self.qsbBottomSpaceValue.setRange(0.0, 1.0)
        self.qsbBottomSpaceValue.setDecimals(2)
        self.qsbBottomSpaceValue.setSingleStep(0.01)
        self.qsbBottomSpaceValue.setValue(0.06)
        qlyAddParams.addWidget(self.qsbBottomSpaceValue, 2, 3, 1, 1)

        qlyAddParams.addWidget(QLabel("Blank width space between subplots"), 3, 0, 1, 2)
        self.qsbBlankWidthSpaceValue = QDoubleSpinBox()
        self.qsbBlankWidthSpaceValue.setRange(0.0, 1.0)
        self.qsbBlankWidthSpaceValue.setDecimals(2)
        self.qsbBlankWidthSpaceValue.setSingleStep(0.01)
        self.qsbBlankWidthSpaceValue.setValue(0.1)
        qlyAddParams.addWidget(self.qsbBlankWidthSpaceValue, 3, 2, 1, 1)

        qlyAddParams.addWidget(QLabel("Blank height space between subplots"), 3, 3, 1, 2)
        self.qsbBlankHeightSpaceValue = QDoubleSpinBox()
        self.qsbBlankHeightSpaceValue.setRange(0.0, 1.0)
        self.qsbBlankHeightSpaceValue.setDecimals(2)
        self.qsbBlankHeightSpaceValue.setSingleStep(0.01)
        self.qsbBlankHeightSpaceValue.setValue(0.1)
        qlyAddParams.addWidget(self.qsbBlankHeightSpaceValue, 3, 5, 1, 1)

        qlyAddParams.setRowMinimumHeight(3, 50)

        qgbAddParams.setLayout(qlyAddParams)

        qlyMainLayout.addWidget(qgbAddParams)

        # graphic parameters import and export

        qgbGraphicParamsIO = QGroupBox("Graphic parameters save/load")

        qlyGraphicParamsIO = QHBoxLayout()

        self.qpbGraphicParamsSave = QPushButton("Save")
        self.qpbGraphicParamsSave.clicked.connect(self.output_graphic_params_save)
        qlyGraphicParamsIO.addWidget(self.qpbGraphicParamsSave)

        self.qpbGraphicParamsLoad = QPushButton("Load")
        self.qpbGraphicParamsLoad.clicked.connect(self.output_graphic_params_load)
        qlyGraphicParamsIO.addWidget(self.qpbGraphicParamsLoad)

        qgbGraphicParamsIO.setLayout(qlyGraphicParamsIO)

        qlyMainLayout.addWidget(qgbGraphicParamsIO)

        # output file parameters

        qgbOutputFile = QGroupBox(self.tr("Output file - available formats: tif, pdf, svg"))

        qlyOutputFile = QGridLayout()

        self.qleFigureOutPath = QLineEdit()
        qlyOutputFile.addWidget(self.qleFigureOutPath, 3, 0, 1, 1)

        self.qpbFigureOutPath = QPushButton(self.tr("Choose"))
        self.qpbFigureOutPath.clicked.connect(self.define_figure_outpath)
        qlyOutputFile.addWidget(self.qpbFigureOutPath, 3, 1, 1, 1)

        qgbOutputFile.setLayout(qlyOutputFile)

        qlyMainLayout.addWidget(qgbOutputFile)

        # execution buttons

        qwdgDecide = QWidget()

        qhblButtons = QHBoxLayout()
        qhblButtons.addStretch()

        qpbOk = QPushButton("&OK")
        qpbCancel = QPushButton("Cancel")

        qhblButtons.addWidget(qpbOk)
        qhblButtons.addWidget(qpbCancel)

        qwdgDecide.setLayout(qhblButtons)

        qlyMainLayout.addWidget(qwdgDecide)

        self.setLayout(qlyMainLayout)

        self.connect(qpbOk, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(qpbCancel, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Export figure")

    def output_graphic_params_save(self):

        output_file_path = new_file_path(self, "Define output configuration file", "*.txt", "txt")

        if not output_file_path:
            return

        out_configuration_string = """figure width = %f
resolution (dpi) = %d
font size (pts) = %f
top space = %f
left space = %f        
right space = %f        
bottom space = %f  
blank width space = %f
blank height space = %f""" % (float(self.qleFigWidthInch.text()),
                              int  (self.qleFigResolutionDpi.text()),
                              float(self.qleFigFontSizePts.text()),
                              float(self.qsbTopSpaceValue.value()),
                              float(self.qsbLeftSpaceValue.value()),
                              float(self.qsbRightSpaceValue.value()),
                              float(self.qsbBottomSpaceValue.value()),
                              float(self.qsbBlankWidthSpaceValue.value()),
                              float(self.qsbBlankHeightSpaceValue.value()))

        with open(output_file_path, "w") as ofile:
            ofile.write(out_configuration_string)

        self.info("Graphic parameters saved")

    def output_graphic_params_load(self):

        input_file_path = old_file_path(self, "Choose input configuration file", "*.txt", "txt")

        if not input_file_path:
            return

        with open(input_file_path, "r") as ifile:
            config_lines = ifile.readlines()

        try:
            figure_width_inches = float(config_lines[0].split("=")[1])
            figure_resolution_dpi = int(config_lines[1].split("=")[1])
            figure_fontsize_pts = float(config_lines[2].split("=")[1])
            top_space_value = float(config_lines[3].split("=")[1])
            left_space_value = float(config_lines[4].split("=")[1])
            right_space_value = float(config_lines[5].split("=")[1])
            bottom_space_value = float(config_lines[6].split("=")[1])
            blank_width_space = float(config_lines[7].split("=")[1])
            blank_height_space = float(config_lines[8].split("=")[1])
        except:
            self.warn("Error in configuration file")
            return

        self.qleFigWidthInch.setText(str(figure_width_inches))
        self.qleFigResolutionDpi.setText(str(figure_resolution_dpi))
        self.qleFigFontSizePts.setText(str(figure_fontsize_pts))
        self.qsbTopSpaceValue.setValue(top_space_value)
        self.qsbLeftSpaceValue.setValue(left_space_value)
        self.qsbRightSpaceValue.setValue(right_space_value)
        self.qsbBottomSpaceValue.setValue(bottom_space_value)
        self.qsbBlankWidthSpaceValue.setValue(blank_width_space)
        self.qsbBlankHeightSpaceValue.setValue(blank_height_space)

    def define_figure_outpath(self):

        outfile_path = new_file_path(self, "Create", "", "Images (*.svg *.pdf *.tif)")

        self.qleFigureOutPath.setText(outfile_path)

    def info(self, msg):

        QMessageBox.information(self, self.sPluginName, msg)

    def warn(self, msg):

        QMessageBox.warning(self, self.sPluginName, msg)


class FigureExportDlg(QDialog):

    def __init__(self, plugin_name, export_params, parent=None):

        super(FigureExportDlg, self).__init__(parent)

        self.sPluginName = plugin_name
        self.dExportParams = export_params

        qlyMainLayout = QVBoxLayout()

        # main parameters groupbox

        qgbMainParams = QGroupBox("Main graphic parameters")

        qlyMainParams = QGridLayout()

        qlyMainParams.addWidget(QLabel(self.tr("Resolution (dpi)")), 0, 2, 1, 1)
        self.qleFigResolutionDpi = QLineEdit(self.dExportParams["expfig_res_dpi"])
        qlyMainParams.addWidget(self.qleFigResolutionDpi, 0, 3, 1, 1)

        qgbMainParams.setLayout(qlyMainParams)

        qlyMainLayout.addWidget(qgbMainParams)

        # output file parameters

        qgbOutputFile = QGroupBox(self.tr("Output file - suggested formats: tif, pdf, svg"))

        qlyOutputFile = QGridLayout()

        self.qleFigureOutPath = QLineEdit()
        qlyOutputFile.addWidget(self.qleFigureOutPath, 3, 0, 1, 1)

        self.qpbFigureOutPath = QPushButton(self.tr("Choose"))
        self.qpbFigureOutPath.clicked.connect(self.define_figure_outpath)
        qlyOutputFile.addWidget(self.qpbFigureOutPath, 3, 1, 1, 1)

        qgbOutputFile.setLayout(qlyOutputFile)

        qlyMainLayout.addWidget(qgbOutputFile)

        # execution buttons

        qwdgDecide = QWidget()

        qhblButtons = QHBoxLayout()
        qhblButtons.addStretch()

        qpbOk = QPushButton("&OK")
        qpbCancel = QPushButton("Cancel")

        qhblButtons.addWidget(qpbOk)
        qhblButtons.addWidget(qpbCancel)

        qwdgDecide.setLayout(qhblButtons)

        qlyMainLayout.addWidget(qwdgDecide)

        self.setLayout(qlyMainLayout)

        self.connect(qpbOk, SIGNAL("clicked()"),
                     self, SLOT("accept()"))
        self.connect(qpbCancel, SIGNAL("clicked()"),
                     self, SLOT("reject()"))

        self.setWindowTitle("Export figure")

    def define_figure_outpath(self):

        outfile_path = new_file_path(self, "Create", "", "Images (*.svg *.pdf *.tif)")

        self.qleFigureOutPath.setText(outfile_path)

    def info(self, msg):

        QMessageBox.information(self, self.sPluginName, msg)

    def warn(self, msg):

        QMessageBox.warning(self, self.sPluginName, msg)


