import os

from sys import platform
import shutil

def build_evoef():

    os.mkdir('perm')
    shutil.move('DynaBench', 'perm')

    import DynaBench

    cp = os.getcwd()
    path = DynaBench.__file__[:-12]

    os.chdir(path)


    if platform == "linux" or platform == "linux2":
        os.system('git clone https://github.com/tommyhuangthu/EvoEF.git')
        os.chdir("./EvoEF/")
        os.system("g++ -O3 --fast-math -o EvoEF src/*.cpp")
        os.chdir("../")

    elif platform == "darwin":
        os.system('git clone https://github.com/tommyhuangthu/EvoEF.git')
        os.chdir("./EvoEF/")
        os.system("g++ -O3 -ffast-math -o EvoEF src/*.cpp")
        os.chdir("../")

    elif platform == 'win32':
        os.system('wsl git clone https://github.com/tommyhuangthu/EvoEF.git')
        os.chdir("EvoEF/")
        os.system("wsl g++ -O3 --fast-math -o EvoEF src/*.cpp")
        os.chdir("../")

    os.chdir(cp)
    
    shutil.move('./perm/DynaBench', cp)
    os.rmdir('perm')

build_evoef()

