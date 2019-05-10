import os

for root, dirs, files in os.walk("."):  
    for filename in files:
        if ('overlay' in filename) or ('bnb' in filename):
            os.system("project.py --xml "+filename+" --stage nuecc --clean")
            os.system("project.py --xml "+filename+" --stage nuecc --submit")
