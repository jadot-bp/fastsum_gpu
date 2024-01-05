import lyncs_io as lio
import sys


def main(argv):
    """
    Converts from openqcd to 'C' format
    if filename ends in .oqcd, replace that with .C
    else append .C
    """
    if len(argv) > 1:
        # Get the file name from the first argument
        file_name = argv[1]
    else:
        # Print an error message if no file name was given
        print("Please provide a full file path as an argument.")
    # setup output filename
    if '.oqcd' in file_name:
        C_file_name = file_name.replace('.oqcd', '.C')
    else:
        C_file_name = file_name + '.C'
    # Load the openqcd gaugefield
    U = lio.load(file_name, format='openqcd')
    # Dump out to .C format
    U.tofile(C_file_name)
    # Print the output filename to screen
    print(C_file_name)

if __name__ == '__main__':
    main(sys.argv)
