"""manifest.py: Code to handle the file manifest to locate raw data.

author(s): @christinehc
"""

# =========================================================
# Imports
# =========================================================
from typing import Optional, Union

import pandas as pd


# =========================================================
# Classes
# =========================================================
class DataManifest:
    """Class to manage and retrieve data from master manifest."""

    def __init__(self, filepath: str):
        """Initialize the DataManifest.

        Parameters
        ----------
        filepath : str
            Path to manifest file (.CSV)
        """
        self.filepath = filepath
        self.manifest = self._load_manifest()
        self._multiple_input = {
            "data_type": False,
            "sample_type": False,
            "name": False,
            "version": False,
        }

    def _load_manifest(self) -> pd.DataFrame:
        """Load manifest file (.CSV)

        Returns
        -------
        pd.DataFrame
            Dataframe containing manifest contents

        Raises
        ------
        ValueError
            If `filepath` is invalid or else manifest fails to load.
        """
        try:
            # Load from external URL
            if self.filepath.startswith(("http://", "https://")):
                manifest = pd.read_csv(self.filepath)
            # Otherwise, load from local file
            else:
                manifest = pd.read_csv(self.filepath)
            return manifest
        except Exception as e:
            raise ValueError(f"Failed to load manifest from {self.filepath}: {e}")

    def _filter_manifest(
        self,
        data_type: Union[str, list[str], None] = None,
        sample_type: Union[str, list[str], None] = None,
        name: Union[str, list[str], None] = None,
        version: Union[str, list[str], None] = None,
    ) -> pd.DataFrame:
        """Filter manifest based on input criteria.

        Parameters
        ----------
        data_type : Union[str, list[str], None], optional
            Data type (value for `data_type` col in manifest),
                by default None
        sample_type : Union[str, list[str], None], optional
            Sample type (value for `sample_type` col in manifest),
                by default None
        name : Union[str, list[str], None], optional
            Name (value for `name` col in manifest),
                by default None
        version : Union[str, list[str], None], optional
            Version (value for `version` col in manifest),
                by default None

        Returns
        -------
        pd.DataFrame
            Dataframe filtered to match input parameters.

        Raises
        ------
        ValueError
            If all input parameters are None.
        """
        # Require at least 1 filtering criterion
        if all(param is None for param in [data_type, sample_type, name, version]):
            raise ValueError("All input parameters cannot be None.")

        # Apply filters
        df = self.manifest
        for col, param in {
            "data_type": data_type,
            "sample_type": sample_type,
            "name": name,
            "version": version,
        }.items():
            if param is not None:
                if isinstance(param, list):
                    df = df[df[col].isin(param)]
                    self._multiple_input[col] = True
                else:
                    df = df[df[col] == param]
                    self._multiple_input[col] = False
        return df

    def get(
        self,
        data_type: Union[str, list[str], None] = None,
        sample_type: Union[str, list[str], None] = None,
        name: Union[str, list[str], None] = None,
        version: Union[str, list[str], None] = None,
        return_first: bool = True,
        output_format: str = "str",
    ) -> Optional[str]:
        """Get file location/URL based on input parameters.

        Parameters
        ----------
        data_type : Union[str, list[str], None], optional
            Data type (value for `data_type` col in manifest),
                by default None
        sample_type : Union[str, list[str], None], optional
            Sample type (value for `sample_type` col in manifest),
                by default None
        name : Union[str, list[str], None], optional
            Name (value for `name` col in manifest),
                by default None
        version : Union[str, list[str], None], optional
            Version (value for `version` col in manifest),
                by default None
        return_first : bool, optional
            If True, returns only first result.
            If False, returns comma-separated list of results.
                by default True
        output_format : str, optional
            Output format if `return_first` = False. Ignored if
                `return_first` = True, by default "str"
            Options: ["str", "list"]
                "str" : Comma-separated string list
                    e.g. "filepath1.csv,filepath2.csv"
                "list" : List of results
                    e.g. ["filepath1.csv", "filepath2.csv"]

        Returns
        -------
        Optional[str]
            Returns one of the following:
                1. URL string
                2. Comma-separated list of URL strings,
                    if `return_first` = False and only one
                    `data_type` selected
                3. Grouped list of lists where each data type is
                    grouped by name across data types.
        """
        filtered = self._filter_manifest(data_type, sample_type, name, version)

        if filtered.empty:
            return None
        elif len(filtered) == 1:
            return filtered.iloc[0]["location"]
        else:
            if return_first:
                # If returning one file, return either highest (newest)
                # version or first match
                if "version" in filtered.columns:
                    return filtered.nlargest(1, "version").iloc[0]["location"]
                else:
                    return filtered.iloc[0]["location"]
            elif self._multiple_input["data_type"] is False:
                return ",".join(list(filtered["location"]))

            # NOTE: currently assumes any multiple without return_first are grouped
            else:
                pivoted = filtered.pivot(
                    index="name", columns="data_type", values="location"
                )
                grouped = [
                    [row[dt] for dt in data_type] for _, row in pivoted.iterrows()
                ]
                return grouped

    def summary(self) -> pd.DataFrame:
        """Summarize available data in manifest by data/sample type.

        Returns
        -------
        pd.DataFrame
            Table of manifest grouped by `data_type` and `sample_type`
        """
        return (
            self.manifest.groupby(["data_type", "sample_type"])
            .agg({"name": "first", "version": "max", "location": "count"})
            .rename(columns={"location": "count"})
            .reset_index()
        )
