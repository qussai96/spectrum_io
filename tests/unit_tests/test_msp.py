import tempfile

import numpy as np
import pandas as pd
import pytest

import spectrum_io.spectral_library.msp as msp


class TestMspPrepareSpectrum:
    """Class to test msp."""

    def test_prepare_spectrum(self, spectra_input, grpc_dict):
        """Test preparation of spectrum."""
        output_path = ""
        msp_lib = msp.MSP(spectra_input, grpc_dict, output_path)
        msp_lib.prepare_spectrum()

    def test_write(self, spectra_input, grpc_dict):
        """Test write to file."""
        out_file = tempfile.NamedTemporaryFile(delete=False)
        msp_lib = msp.MSP(spectra_input, grpc_dict, out_file.name)
        msp_lib.prepare_spectrum()
        msp_lib.write()
        file_content = out_file.read().decode()
        file_content = file_content.replace("\r", "")  # explicitly remove to work for windows
        anticipated_content = (
            "Name: AAACCCCKR/1\n"
            "MW: 124.407276467\n"
            "Comment: Parent=124.407276467 Collision_energy=10.0 Mods=2/3,C,Carbamidomethyl/5,C,Carbamidomethyl "
            "ModString=AAACCCCKR//Carbamidomethyl@C3; Carbamidomethyl@C5/1 iRT=982.12 proteotypicity=123.1\n"
            "Num peaks: 2\n"
            '0.9	0.1	"b1/0.0ppm"\n'
            '0.8	0.2	"y1^2/0.0ppm"\n'
            "Name: AAACILKKR/2\n"
            "MW: 1617.057276467\n"
            "Comment: Parent=1617.057276467 Collision_energy=20.0 Mods=0 ModString=AAACILKKR///2 iRT=382.12 proteotypicity=234.2\n"
            "Num peaks: 2\n"
            '0.6	0.4	"b1^2/0.0ppm"\n'
            '0.5	0.5	"y3^3/0.0ppm"\n'
        )
        assert file_content == anticipated_content


@pytest.fixture
def spectra_input():
    """Test spectra input."""
    spectra_input = pd.DataFrame()
    spectra_input["MODIFIED_SEQUENCE_SPEC"] = ["AAACCCC", "AAACILK"]
    spectra_input["MODIFIED_SEQUENCE"] = ["AAAC[UNIMOD:4]CC[UNIMOD:4]CKR", "AAACILKKR"]
    spectra_input["MASS"] = [123.4, 3232.1]
    spectra_input["COLLISION_ENERGY"] = [10.0, 20.0]
    spectra_input["PRECURSOR_CHARGE"] = [1, 2]
    return spectra_input


@pytest.fixture
def grpc_dict():
    """Creates grpc dictionary."""
    grpc_dict = {
        "model": {
            "intensity": np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
            "fragmentmz": np.array([[0.9, 0.8, 0.7], [0.6, 0.5, 0.4]]),
            "annotation": {
                "charge": np.array([[1, 2, 3], [2, 3, 1]]),
                "number": np.array([[1, 1, 2], [1, 3, 5]]),
                "type": np.array([["b", "y", "N"], ["b", "y", "N"]]),
            },
        },
        "model_irt": np.array([[982.12], [382.12]]),
        "model_proteotypicity": np.array([[123.1], [234.2]]),
    }
    return grpc_dict
