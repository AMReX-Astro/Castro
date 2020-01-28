import re
import sys

def process_analysis(filename):
    outtxt = ""
    with open(filename, 'r') as f:
        r = re.compile(r'^(.\.\/\.\.\/\.\.\/Source\/.*?(?:warning|note).*?)(?=\n\S)', flags=re.M|re.S)

        n_warnings = 0

        for m in r.finditer(f.read()):
            outtxt += m.group(1) + '\n'
            n_warnings += 1

        if outtxt == "":
            # no warnings were found in Castro
            print("Static analysis found no warnings in Castro")
        else:
            print('Static analysis warnings from Castro:\n-------------------------------------\n')
            print(outtxt)
            print(f'{n_warnings} bugs found')
            return False 

    return True



if __name__=="__main__":
    process_analysis(sys.argv[1])