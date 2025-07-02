import re
import sys

def process_analysis(filename):

    with open(filename, 'r') as f:
        r = re.compile(r'^(.\.\/\.\.\/\.\.\/Source\/[\w/]+\.cpp.*?(?:warning|note).*?)(?=\n\S)', flags=re.M|re.S)

        bugs = set()

        for m in r.finditer(f.read()):
            if 'ignoring #pragma gpu box' not in m.group(1):
                bugs.add(m.group(1))

        n_warnings = len(bugs)

        if n_warnings == 0:
            # no warnings were found in Castro
            print("Static analysis found no warnings in Castro")
        else:
            print('Static analysis warnings from Castro:\n-------------------------------------\n')
            print('\n'.join(bugs))
            print(f'{n_warnings} bugs found')
            sys.exit(1)

    return True

if __name__=="__main__":
    process_analysis(sys.argv[1])