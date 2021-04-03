import os
import subprocess
import sys
import warnings
# import pybedtools as bedtools


# Biopython
def check_for_biopython():
    try:
        import Bio
        print("Biopython was successfully found")

        if float(Bio.__version__) >= 1.78:
            print("Biopython version is valid")
            return True

        warnings.warn('Biopython seems to be an illegal version! Make sure it is version = 1.78 or above.')
        return False

    except ModuleNotFoundError:
        warnings.warn("Biopython was not found")
        exit()


def check_for_pybedtools():
    try:
        import pybedtools
        print("pybedtools was successfully found")

        if pybedtools.__version__ >= '0.8.2':
            print("pybedtools version is valid")
            return True

        warnings.warn('pybedtools seems to be an illegal version! Make sure it is version = 0.8.2 or above.')
        return False

    except ModuleNotFoundError:
        warnings.warn("pybedtools was not found")
        exit()

def check_for_bedtools():
    try:
        cmd_return = subprocess.run(['bedtools', '-h'], capture_output=True)

        if cmd_return.returncode == 0:
            print("Bedtools was successfully found")

            cmd_return = subprocess.run(['bedtools', '-version'], capture_output=True)
            if cmd_return.stdout.decode().rsplit('v', 1)[-1] >= '2.29.2':
                print("Bedtools version is valid")
                return True

            else:
                warnings.warn('Bedtools seems to be an illegal version! Make sure it is version = 2.29.2 or above.')
                return False
        else:
            warnings.warn('The command "bedtools -h" does not seem to give return code = 0')
            exit()

    except FileNotFoundError:
        warnings.warn('Bedtools could not be found!')
        exit()


# Samtools
def check_for_samtools():
    try:
        cmd_return = subprocess.run(['samtools', '--help'], capture_output=True)

        if cmd_return.returncode == 0:
            print("Samtools was successfully found")

            cmd_return = subprocess.run(['samtools', '--version'], capture_output=True)
            # Isolate version tag
            version = cmd_return.stdout.decode().split(' ', 1)[-1]
            version = version.split('\n', 1)[0]

            if version >= '1.19' or version == '1.12':
                print("Samtools version is valid")
                return True

            else:
                warnings.warn('Samtools seems to be an illegal version! Make sure it is version = 1.19, 1.12 or above.')
                return False
        else:
            warnings.warn('The command "samtools -h" does not seem to give return code = 0')
            exit()

    except FileNotFoundError:
        warnings.warn('Samtools could not be found!')
        exit()


if __name__ == '__main__':
    #TODO - add name of tool!
    print('\n------ Checking dependencies for Magphi ------')
    biopython_presence = check_for_biopython()

    pybedtools_presence = check_for_pybedtools()

    bedtools_presence = check_for_bedtools()

    samtools_presence = check_for_samtools()

    print('\n------ Summary of dependency check ------')

    if all([biopython_presence, pybedtools_presence, bedtools_presence, samtools_presence]):
        print("All dependencies were found to be okay! We are ready to go!")
    else:
        print("!!! Some dependency does not seem be the correct version, please check the warnings above !!!")
