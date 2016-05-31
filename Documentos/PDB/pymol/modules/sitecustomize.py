# NOTE: This only works for module-based PyMOL builds.
# Embedded versions of PyMOL call PyUnicode_SetDefaultEncoding at startup
import sys
sys.setdefaultencoding("utf-8")

# add "modules" directory as site dir
try:
    import site, os.path
    site.addsitedir(os.path.dirname(__file__))
except ImportError:
    pass
