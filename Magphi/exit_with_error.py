import sys
import logging
import os


def exit_with_error(message, exit_status, tmp_folder=None):
    """
    Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.
    :param message: Message to give the user upon exit
    :param exit_status: Status returned as exit status
    :param tmp_folder: Temporary folder for Magphi to be deleted under some circumstances.
    :return: None
    """

    # Delete tmp files and folder
    try:
        if tmp_folder is not None:
            tmp_files = os.listdir(tmp_folder)
            for file in tmp_files:
                os.remove(os.path.join(tmp_folder, file))
            os.rmdir(tmp_folder)
        else:
            pass
    except FileNotFoundError:
        pass

    logging.error(message)
    print(f"Magphi ERROR: {message}, exiting", file=sys.stderr)
    sys.exit(exit_status)
