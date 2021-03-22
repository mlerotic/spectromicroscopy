import sys


def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    # Do argument parsing here (eg. with argparse) and anything else
    # you want your project to do. Return values are exit codes.
    
    from . import mantis_qt
    mantis_qt.main()  # Open the GUI


if __name__ == "__main__":
    sys.exit(main())
