"""data.py: Handle data from Figshare.

author(s): @christinehc

Note: Written with assistance from Claude Sonnet 4.0
    - FigshareDataDownloader:
        - Caching methods written using Claude
        - Downloading methods originally written using Claude
            and modified by @christinehc
    - FigshareDataLoader: Written entirely by @christinehc
"""

# =========================================================
# Imports
# =========================================================
import json
import os
from pathlib import Path
from typing import Optional, Union

import pandas as pd
import requests


# =========================================================
# Class Definitions
# =========================================================
class FigshareDownloader:
    """Class for downloading data from Figshare"""

    def __init__(
        self,
        cache_dir: str = "./.figshare_cache",
        api_token: Optional[str] = None,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ):
        # Cache setup
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.downloaded_files = {}

        # Authentication setup
        self.api_token = api_token or os.getenv("FIGSHARE_API_TOKEN")
        self.username = username or os.getenv("FIGSHARE_USERNAME")
        self.password = password or os.getenv("FIGSHARE_PASSWORD")

        # Session for authenticated requests
        self.session = requests.Session()
        self._setup_authentication()
        self._load_cache_manifest()

    def _setup_authentication(self):
        """Setup authentication headers and session"""
        if self.api_token:
            # Token-based authentication (preferred)
            self.session.headers.update({"Authorization": f"token {self.api_token}"})
        elif self.username and self.password:
            # Basic authentication fallback
            from requests.auth import HTTPBasicAuth

            self.session.auth = HTTPBasicAuth(self.username, self.password)
        else:
            print("No authentication configured - only public files accessible")

    def _load_cache_manifest(self):
        """Load manifest of previously downloaded files"""
        manifest_path = self.cache_dir / "manifest.json"
        if manifest_path.exists():
            with open(manifest_path, "r") as f:
                self.downloaded_files = json.load(f)

    def _save_cache_manifest(self):
        """Save manifest of downloaded files"""
        manifest_path = self.cache_dir / "manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(self.downloaded_files, f, indent=2)

    def download(
        self,
        file_id: Union[str, int],
        filename: Optional[str] = None,
        force_redownload: bool = False,
    ) -> Path:
        """
        Download a file from Figshare

        Args:
            file_id: Figshare file ID
            filename: Optional custom filename
            force_redownload: Re-download even if file exists

        Returns:
            Path to downloaded file
        """
        file_id = str(file_id)

        # Check if already downloaded
        if file_id in self.downloaded_files and not force_redownload:
            file_path = Path(self.downloaded_files[file_id]["path"])
            if file_path.exists():
                print(f"File {file_id} already cached at {file_path}")
                return file_path
        print(f"Downloading file {file_id}...")

        # Get file metadata for download URL
        url = f"https://ndownloader.figshare.com/files/{file_id}"

        response = self.session.get(url, stream=True)
        response.raise_for_status()

        # Get filename
        if filename is None:
            content_disposition = response.headers.get("content-disposition", "")
            if "filename=" in content_disposition:
                filename = content_disposition.split("filename=")[-1].strip('"')
            else:
                filename = f"figshare_file_{file_id}"

        file_path = self.cache_dir / filename

        return self._save_file(response, file_id, filename, url)

    def _save_file(
        self,
        response: requests.Response,
        file_id: str,
        filename: str,
        url: str,
    ) -> Path:
        """Save downloaded file and update manifest"""
        file_path = self.cache_dir / filename

        with open(file_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

            # Update manifest
            self.downloaded_files[file_id] = {
                "path": str(file_path),
                "filename": filename,
                "url": url,
                "file_size_mb": file_path.stat().st_size / 1000,
            }
            self._save_cache_manifest()

            print(f"Downloaded: {filename} -> {file_path}")
            return file_path

    def get_file_path(self, file_id: Union[str, int]) -> Optional[Path]:
        """Get path to cached file if it exists"""
        file_id = str(file_id)
        if file_id in self.downloaded_files:
            path = Path(self.downloaded_files[file_id]["path"])
            return path if path.exists() else None
        return None

    def list_cached_files(self) -> dict[str, Union[str, int]]:
        """List all cached files"""
        return self.downloaded_files

    def clear_cache(self):
        """Clear all cached files"""
        for file_info in self.downloaded_files.values():
            file_path = Path(file_info["path"])
            if file_path.exists():
                file_path.unlink()
        self.downloaded_files.clear()
        self._save_cache_manifest()


class FigshareDataLoader(FigshareDownloader):
    """Data loader class for Figshare data."""

    def __init__(
        self,
        cache_dir: str = "./.figshare_cache",
        api_token: Optional[str] = None,
        username: Optional[str] = None,
        password: Optional[str] = None,
    ):
        super().__init__(cache_dir, api_token, username, password)
        self.data_cache = {}

    def load_data(
        self,
        file_id: Union[str, int],
        force_reload: bool = False,
    ) -> pd.DataFrame:
        """Download and load data from Figshare.

        Parameters
        ----------
        file_id : Union[str, int]
            Figshare file ID
        force_reload : bool, optional
            If True, forces reload even if data already cached,
                by default False

        Returns
        -------
        pd.DataFrame
            Dataframe of loaded data
        """
        file_id = str(file_id)

        # Check cache for existing data
        if file_id in self.data_cache and not force_reload:
            print(f"Data for {file_id} already loaded (cached)")
            return self.data_cache[file_id]

        # If not cached, download and load data
        file_path = self.download(file_id)
        data = self._load_file(file_path)

        # Save data to cache
        self.data_cache[file_id] = data

        return data

    def _load_file(self, file_path: str) -> pd.DataFrame:
        """Load data from file.

        Parameters
        ----------
        file_path : str
            Path to data file
                Note: Accepts .csv, .tsv, and .xlsx files.

        Returns
        -------
        pd.DataFrame
            Dataframe of loaded data

        Raises
        ------
        ValueError
            If file extension is not one of .csv, .tsv, or .xlsx
        """
        file_type = os.path.splitext(file_path)[1]

        if file_type == ".csv":
            return pd.read_csv(file_path)

        elif file_type == ".tsv":
            return pd.read_csv(file_path, sep="\t")

        elif file_type == ".xlsx":
            return pd.read_excel(file_path)

        else:
            raise ValueError(
                f"Unexpected file extension: {file_type}."
                "Valid extensions are .csv, .tsv, or .xlsx."
            )
