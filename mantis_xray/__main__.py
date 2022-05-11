import sys, getopt


def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    # Do argument parsing here (eg. with argparse) and anything else
    # you want your project to do. Return values are exit codes.
    try:
        options, extraParams = getopt.getopt(args, '', ['batch', 'nnma'])
    except:
        print('Error - wrong command line option used. Available options are --batch and --nnma')
        return
    
    batch_mode = False
    for opt, arg in options:
        if opt in '--batch':
            batch_mode = True
    
    if batch_mode:
        from . import mantis
        mantis.main()
    else:
        from . import mantis_qt
        mantis_qt.main()  # Open the GUI

if __name__ == "__main__":
    sys.exit(main())
