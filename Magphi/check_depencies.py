import subprocess
import warnings
import sys
# from Magphi.exit_with_error import exit_with_error
from exit_with_error import exit_with_error
EXIT_DEPENDENCY_ERROR = 4


def check_for_biopython(verbose):
    '''Check for the presence of Biopython and if the version satisfies the required version 1.78 or above,
    if this can not be satisfied then warn the use of a missing dependency and exit.
    '''
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
    '''Check for the presence of pybedtools and if the version satisfies the required version 0.8.2 or above,
        if this can not be satisfied then warn the use of a missing dependency and exit.
        '''
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
    '''Check for the presence of bedtools and if the version satisfies the required version 2.29.2 or above,
            if this can not be satisfied then warn the use of a missing dependency and exit.
            '''
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
    '''Check for the presence of samtools and if the version satisfies the required version 1.12, 1.19 or above,
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


def check_dependencies_for_main(verbose = False):
    ''' Function to check dependencies when Magphi is executed. If all programs are found their versions are returned,
    if they are not found False is returned.'''
    # TODO - log the version of tools used for the run. - possibly do this in main when TRUE is returned. or return the versions instead of TRUE
    biopython_presence = check_for_biopython(verbose)

    pybedtools_presence = check_for_pybedtools(verbose)

    bedtools_presence = check_for_bedtools(verbose)

    samtools_presence = check_for_samtools(verbose)

    # if all([biopython_presence, pybedtools_presence, bedtools_presence, samtools_presence]):
    if all([biopython_presence, pybedtools_presence]):
        return [biopython_presence, pybedtools_presence]
        # return [biopython_presence, pybedtools_presence, bedtools_presence, samtools_presence]
    else:
        return False


if __name__ == '__main__':
    #TODO - add name of tool!
    print('\n------ Checking dependencies for Magphi ------')
    biopython_presence = check_for_biopython(verbose=True)

    pybedtools_presence = check_for_pybedtools(verbose=True)

    bedtools_presence = check_for_bedtools(verbose=True)

    samtools_presence = check_for_samtools(verbose=True)

    print('\n------ Summary of dependency check ------')

    if all([biopython_presence, pybedtools_presence, bedtools_presence, samtools_presence]):
        print("All dependencies were found to be okay! We are ready to go!")
    else:
        print("!!! Some dependency does not seem be the correct version, please check the warnings above !!!")
