import ctypes as C
import sys
sys.path.insert(0, './MapMatching')

nav = C.CDLL('./MapMatching/libnav.so', mode=C.RTLD_GLOBAL)

class sotest:
    def __init__(self):
        self.lib = nav

    def test(self):
        self.lib.main()