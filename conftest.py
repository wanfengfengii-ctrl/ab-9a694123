# 使 pytest 无论从何处启动都能导入仓库根目录下的 app 包。
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
