from abaqusGui import getAFXApp, Activator, AFXMode, afxCreatePNGIcon
from abaqusConstants import ALL
import os
thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerGuiMenuButton(
    buttonText='Elastic Properties Calculator (EasyPBC)', 
    object=Activator(os.path.join(thisDir, 'easyPBCDB.py')),
    kernelInitString='import easypbc',
    messageId=AFXMode.ID_ACTIVATE,
    icon=afxCreatePNGIcon(os.path.join(thisDir,'easypbc.png')),
    applicableModules=ALL,
    version='1.3',
    author='Sadik Omairey',
    description='N/A',
    helpUrl='N/A'
)
