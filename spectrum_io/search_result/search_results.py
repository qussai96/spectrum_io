import logging
from abc import abstractmethod
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from spectrum_io.file import csv

logger = logging.getLogger(__name__)


class SearchResults:
    """Handle search results from different software."""

    orig_res: pd.DataFrame
    fake_msms: pd.DataFrame

    def __init__(self, path: Union[str, Path]):
        """
        Init Searchresults object.

        :param path: path to file
        """
        if isinstance(path, str):
            path = Path(path)
        self.path = path

    @abstractmethod
    def read_result(self, path: Union[str, Path], tmt_labeled: str):
        """Read result."""
        raise NotImplementedError

    def generate_internal(self, tmt_labeled: str, out_path: Optional[Union[str, Path]] = None) -> str:
        """
        Generate df and save to out_path.

        :param out_path: path to output
        :param tmt_labeled: tmt label as str
        :return: path to output file
        """
        if out_path is None:
            out_path = self.path.with_suffix(".prosit")
        if isinstance(out_path, str):
            out_path = Path(out_path)

        if out_path.is_file():
            logger.info(f"Found search results in internal format at {out_path}, skipping conversion")
            return out_path

        df = self.read_result(self.path, tmt_labeled)
        csv.write_file(df, out_path)

        return out_path

    def read_internal(self, path: Union[str, Path]) -> pd.DataFrame:
        """
        Read file from path.

        :param path: path to file
        :return: dataframe after reading the file
        """
        return csv.read_file(path)
