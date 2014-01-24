from distutils.core import setup
import os
VERSION = '2.1.1'
SKJAR_FILE = 'SpatialKappa-v' + VERSION + '.jar'
SKJAR_FILE_PATH = os.path.join('..', 'mlm', 'SpatialKappa', SKJAR_FILE)
ANTLRJAR_FILE = 'ant-antlr-3.2.jar'
ANTLRJAR_FILE_PATH = os.path.join('..', 'mlm', 'SpatialKappa', ANTLRJAR_FILE)


setup(name='SpatialKappa',
      version=VERSION,
      py_modules=['SpatialKappa'],
      data_files=[('share/SpatialKappa',[SKJAR_FILE_PATH, ANTLRJAR_FILE_PATH])],)
