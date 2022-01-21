import pandas as pd
import numpy as np
from .spectral_library import SpectralLibrary
from fundamentals.mod_string import internal_without_mods, internal_to_spectronaut


class Spectronaut(SpectralLibrary):
    # Check spectronaut folder for output format.

    def write(self):
        """
        Writing method; splits intermediate dataframe into chunks of 10k lines, explodes them, and filters for relevant
        peaks before storing it in self.out_path.
        @return: None
        """
        out = open(self.out_path, "w")
        out.truncate()
        out.close()

        n = 7000  # split df into chunks of size n
        initial = True
        for group, segment in self.spectra_output.groupby(np.arange(len(self.spectra_output)) // n):
            segment = segment.explode(['intensities', 'fragment_mz', 'fragment_types', 'fragment_numbers', 'fragment_charges'])
            segment = segment[segment['intensities'] > 0]  # set to >= if 0 should be kept
            segment.rename(columns={'intensities': 'RelativeIntensity', 'fragment_mz': 'FragmentMz', 'fragment_types': 'FragmentType', 'fragment_numbers':'FragmentNumber', 'fragment_charges':'FragmentCharge'},inplace=True)
            segment['FragmentLossType'] = 'noloss'
            segment = segment[['RelativeIntensity','FragmentMz','ModifiedPeptide','LabeledPeptide','StrippedPeptide','PrecursorCharge','PrecursorMz','iRT','proteotypicity','FragmentNumber','FragmentType','FragmentCharge','FragmentLossType']]
            segment.to_csv(self.out_path, mode='a', header=initial,index=False)
            if initial:
                initial = False

    def prepare_spectrum(self):
        """
        Converts grpc output and metadata dataframe into spectronaut format
        @return: intermediate dataframe
        """
        intensities = self.grpc_output[list(self.grpc_output)[0]]['intensity']
        fragment_mz = self.grpc_output[list(self.grpc_output)[0]]['fragmentmz']
        annotation = self.grpc_output[list(self.grpc_output)[0]]['annotation']
        fragment_types = annotation['type']
        fragment_numbers = annotation['number']
        fragment_charges = annotation['charge']
        irt = self.grpc_output[list(self.grpc_output)[1]]
        irt = irt.flatten()
        proteotypicity = self.grpc_output[list(self.grpc_output)[2]]
        proteotypicity = proteotypicity.flatten()
        modified_sequences_spec = internal_to_spectronaut(self.spectra_input['MODIFIED_SEQUENCE'].apply(lambda x: '_' + x + '_'))
        modified_sequences = self.spectra_input['MODIFIED_SEQUENCE']

        labelled_sequences = internal_without_mods(modified_sequences)
        stripped_peptide = internal_without_mods(modified_sequences)
        charges = self.spectra_input['PRECURSOR_CHARGE']
        precursor_masses = self.spectra_input['MASS']
        precursor_mz = (precursor_masses + charges) / charges

        inter_df = pd.DataFrame(data={'ModifiedPeptide': modified_sequences_spec, 'LabeledPeptide': labelled_sequences,
                                      'StrippedPeptide': stripped_peptide, 'PrecursorCharge': charges,
                                      'PrecursorMz': precursor_mz})
        inter_df['iRT'], inter_df['proteotypicity'] = irt.tolist(), proteotypicity.tolist()
        inter_df['intensities'], inter_df['fragment_mz'] = intensities.tolist(), fragment_mz.tolist()
        inter_df['fragment_types'] = fragment_types.tolist()
        inter_df['fragment_numbers'] = fragment_numbers.tolist()
        inter_df['fragment_charges'] = fragment_charges.tolist()

        self.spectra_output = inter_df
