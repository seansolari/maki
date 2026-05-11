
import os
import ctypes

lib_path = os.path.join(os.path.dirname(__file__), "libs", "_maki.cpython-313-x86_64-linux-gnu.so")
ctypes.CDLL(lib_path)

from _maki import *
