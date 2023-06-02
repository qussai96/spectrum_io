from spectrum_io.raw.msraw import MSRaw
import pandas as pd
from pathlib import Path

class SpectraFilter:
    def __init__(self, mzml_dir, msms_file):
        self.mzml_dir = mzml_dir
        self.msms_file = msms_file

    def filter_spectra(self):
        unmatched_spectra = self.get_unmatched_spectra()
        self.process_mzml_directory(unmatched_spectra)

    def get_unmatched_spectra(self):
        """
        get a list of spectra that didn't match to any peptide sequence
        """
        raw_df = MSRaw.read_mzml(self.mzml_dir)
        msms_df = pd.read_csv(self.msms_file, sep='\t', index_col=False, header=0)
        msms_df.rename(columns={'Scan number': 'SCAN_NUMBER'}, inplace=True)
        merged_df = pd.merge(raw_df, msms_df, on='SCAN_NUMBER', how='inner')
        unmatched_df = raw_df[~raw_df['SCAN_NUMBER'].isin(merged_df['SCAN_NUMBER'])]
        unmatched_spectra = unmatched_df['SCAN_NUMBER'].tolist()
        print(f'msms_file has: {len(msms_df)} scans')
        print(f'number of matched spectra: {len(raw_df)} scans')
        print(f'number of unmatched spectra: {len(unmatched_df)} scans')
        return unmatched_spectra

    def process_mzml_directory(self, unmatched_spectra):
        """
        iterate over all mzml files in the given directory
        """
        path = Path(self.mzml_dir)
        mzml_files = path.glob("*.mzML")
        for mzml_file in mzml_files:
            unmatched_file = mzml_file.with_name(f"{mzml_file.stem}_unmatched_spectra.mzML")
            self.remove_spectra(mzml_file, unmatched_file, unmatched_spectra)

    @staticmethod
    def remove_spectra(mzml_file, unmatched_file, unmatched_spectra):
    """
    take mzml file and return a new file with the unmatched spectra
    """
        with open(mzml_file, 'r') as file:
            lines = file.readlines()

        output_lines = []
        skip_next_lines = False

        for line in lines:
            if skip_next_lines:
                if line.strip() == '</spectrum>':
                    skip_next_lines = False
                continue

            if '<spectrum id=' in line:
                scan_number = int(line.split('scan=')[1].split()[0].strip('"'))
                if scan_number not in unmatched_spectra:
                    skip_next_lines = True
                    continue

            output_lines.append(line)

        with open(unmatched_file, 'w') as file:
            file.writelines(output_lines)
