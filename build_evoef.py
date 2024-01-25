import os

from sys import platform
import shutil

import argparse

parser = argparse.ArgumentParser()

def build_evoef():
    
    cp = os.getcwd()

    os.mkdir('perm')
    shutil.move('DynaBench', 'perm')

    import DynaBench

    path = DynaBench.__file__[:-12]

    try:

        if platform == "win32":
            def install_evoef(path):
                
                dynabench_destination = DynaBench.__file__[:-12]
                print(dynabench_destination)
                shutil.move(path, dynabench_destination)


            parser.add_argument("--evoef_path", type=str)
            args = parser.parse_args()

            install_evoef(args.evoef_path)
            
        else:

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

            os.chdir(cp)


        shutil.move('./perm/DynaBench', cp)
        os.rmdir('perm')
    finally:
        if os.path.exists("./perm/DynaBench"):
            shutil.move('./perm/DynaBench', cp)
            os.rmdir('perm')  

build_evoef()

