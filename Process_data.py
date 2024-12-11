import argparse
import glob
from loadcorrprocess import *

def main(data_folder, fig_folder):
    # Load the data
    print('data_folder:', data_folder)
    fileL = sorted(glob.glob(data_folder + '*' + '.hdf5'))
    
    if not fileL:
        print("No HDF5 files found in the specified data folder.")
        return

    filename = fileL[0]
    print('Loading %s' % filename)
    if 'Corr' in filename:
        print(f"{filename} is a correlator file")
    
    realname = filename.split('/')[-1]  # Remove path
    reallabel = realname.split('.')[0]  # Remove suffix
    reallabel = reallabel.split(' ')[0]  # Remove Gary's name, adjust if needed
    print('reallabel:', reallabel)

    if 'SouthWest' in filename:
        print('SW')
        mode = 'SA_SW'
    elif 'NorthEast' in filename:
        print('NE')
        mode = 'SA_NE'
    else:
        print('Correlator mode')
        mode = 'Correlator'

    if 'Corr' in filename:
        aveX, f_GHz, tUtcBytes, az_degX, el_degX, mjd, \
        satTles, aveY, az_degY, el_degY, mjd, \
        XYave, hdrD = getCorrData(filename)
        f_corr = f_GHz
        f_GHz = np.linspace(10.7, 10.7 + 2, len(f_GHz))
        mjd_corr = mjd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process dynamic spectra numpy arrays from Onsala data.')
    parser.add_argument('--data_folder', type=str, required=True, help='Folder where the data files are located')
    parser.add_argument('--output_folder', type=str, default='npdata/',required=True, help='Folder where the output will be saved')

    args = parser.parse_args()
    main(args.data_folder, args.fig_folder)