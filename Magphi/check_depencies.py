import subprocess
import warnings
import sys

try:
    from Magphi.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error
EXIT_DEPENDENCY_ERROR = 4


def check_for_biopython(verbose):
    """
    Function to test presence and version of BioPython
    :param verbose: To print or not to print
    :return: Version of BioPython if present else False
    """
    try:
        import Bio
        if verbose:
            print("Biopython was successfully found")

        if float(Bio.__version__) >= 1.78:
            if verbose:
                print("Biopython version is valid")
            return float(Bio.__version__)

        warnings.warn('Biopython seems to be an untested version! Make sure it is version = 1.78 or above.')
        return False

    except ModuleNotFoundError as exception:
        warnings.warn("Biopython was not found")
        exit_with_error(str(exception), EXIT_DEPENDENCY_ERROR)


def check_for_pybedtools(verbose):
    """
    Function to test presence and version of PyBedtools
    :param verbose: To print or not to print
    :return: Version of PyBedtools if present else False
    """
    try:
        import pybedtools
        if verbose:
            print("pybedtools was successfully found")

        if pybedtools.__version__ >= '0.8.2':
            if verbose:
                print("pybedtools version is valid")
            return pybedtools.__version__

        warnings.warn('pybedtools seems to be an untested version! Make sure it is version = 0.8.2 or above.')
        return False

    except ModuleNotFoundError as exception:
        warnings.warn("pybedtools was not found")
        exit_with_error(str(exception), EXIT_DEPENDENCY_ERROR)


def check_for_bedtools(verbose):
    """
    Function to test presence and version of Bedtools
    :param verbose: To print or not to print
    :return: Version of Bedtools if present else False
    """
    try:
        cmd_return = subprocess.run(['bedtools', '-h'], capture_output=True)

        if cmd_return.returncode == 0:
            if verbose:
                print("Bedtools was successfully found")

            cmd_return = subprocess.run(['bedtools', '-version'], capture_output=True)
            version = cmd_return.stdout.decode().rsplit('v', 1)[-1]
            if version >= '2.29.2':
                if verbose:
                    print("Bedtools version is valid")
                return version

            else:
                warnings.warn('Bedtools seems to be an untested version! Make sure it is version = 2.29.2 or above.')
                return False
        else:
            warnings.warn('The command "bedtools -h" does not seem to give return code = 0')
            sys.exit(EXIT_DEPENDENCY_ERROR)

    except FileNotFoundError as exception:
        warnings.warn('Bedtools could not be found!')
        exit_with_error(str(exception), EXIT_DEPENDENCY_ERROR)


# Samtools
def check_for_samtools(verbose):
    """
    Function to test presence and version of samtools
    :param verbose: To print or not to print
    :return: Version of samtools if present else False
    """
    ''' Check for the presence of samtools and if the version satisfies the required version 1.12, 1.19 or above,
            if this can not be satisfied then warn the use of a missing dependency and exit.
            '''
    try:
        cmd_return = subprocess.run(['samtools', '--help'], capture_output=True)

        if cmd_return.returncode == 0:
            if verbose:
                print("Samtools was successfully found")

            cmd_return = subprocess.run(['samtools', '--version'], capture_output=True)
            # Isolate version tag
            version = cmd_return.stdout.decode().split(' ', 1)[-1]
            version = version.split('\n', 1)[0]

            if version >= '1.19' or version == '1.12':
                if verbose:
                    print("Samtools version is valid")
                return version

            else:
                warnings.warn('Samtools seems to be an untested version! Make sure it is version = 1.19, 1.12 or above.')
                return False
        else:
            warnings.warn('The command "samtools -h" does not seem to give return code = 0')
            sys.exit(EXIT_DEPENDENCY_ERROR)

    except FileNotFoundError as exception:
        warnings.warn('Samtools could not be found!')
        exit_with_error(str(exception), EXIT_DEPENDENCY_ERROR)


def check_for_blast_plus(verbose=False):
    """
    Function to test presence and version of Blast+
    :param verbose: To print or not to print
    :return: Version of Blast+ if present else False
    """
    ''' Check for the presence of Blast+ and if the version satisfies the required version ,
                if this can not be satisfied then warn the use of a missing dependency and exit. '''

    try:
        cmd_return = subprocess.run(['makeblastdb', '-version'], capture_output=True)

        if cmd_return.returncode == 0:
            if verbose:
                print("Blast+ was successfully found")

            # Isolate version tag
            version = cmd_return.stdout.decode().split('\n')[1]
            version = version.split(' ')[3]
            version = version.rsplit(',', 1)[0]

            if version == '2.6.0':
                if verbose:
                    print("Blast+ version is valid")
                return version

            else:
                warnings.warn('Blast+ seems to be an untested version! Make sure it is version = 1.19, 1.12 or above.')
                return False
        else:
            warnings.warn('The command "makeblastdb -version" does not seem to give return code = 0')
            sys.exit(EXIT_DEPENDENCY_ERROR)

    except FileNotFoundError as exception:
        warnings.warn('Blast+ could not be found!')
        exit_with_error(str(exception), EXIT_DEPENDENCY_ERROR)


def check_dependencies_for_main(verbose=False):
    """
    Function to be called from main of Magphi to check the dependecies before running
    :param verbose: To print or not to print
    :return: Versions of programs of False if one of them are not present.
    """
    ''' Function to check dependencies when Magphi is executed. If all programs are found their versions are returned,
    if they are not found False is returned.'''
    biopython_present = check_for_biopython(verbose)

    pybedtools_present = check_for_pybedtools(verbose)

    bedtools_present = check_for_bedtools(verbose)

    samtools_present = check_for_samtools(verbose)

    if all([biopython_present, pybedtools_present, bedtools_present, samtools_present]):
        return [biopython_present, pybedtools_present, bedtools_present, samtools_present]
    else:
        return False


def check_dependencies_only():
    print('\n------ Checking dependencies for Magphi ------')
    biopython_presence = check_for_biopython(verbose=True)

    pybedtools_presence = check_for_pybedtools(verbose=True)

    bedtools_presence = check_for_bedtools(verbose=True)

    samtools_presence = check_for_samtools(verbose=True)

    check_for_blast_plus(verbose=True)

    print('\n------ Summary of dependency check ------')

    if all([biopython_presence, pybedtools_presence, bedtools_presence, samtools_presence]):
        print("All dependencies were found to be okay! We are ready to go! \n")
    else:
        print("!!! Some dependency does not seem be the correct version, please check the warnings above !!! \n")

    sys.exit(0)

if __name__ == '__main__':
    pass
