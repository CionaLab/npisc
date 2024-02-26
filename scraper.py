from typing import List, Optional
from urllib.parse import urlparse
import hashlib
import os
import io
import mimetypes

import requests
from lxml import etree


def http_get(url: str, verify_ssl: bool = True) -> requests.Response:
    """
    Send a GET request to the specified URL and return the response.

    Args:
    - url (str): The URL to send the GET request to.
    - verify_ssl (bool): Whether to verify SSL certificates.

    Returns:
    - requests.Response: The response from the GET request.
    """
    try:
        response = requests.get(url, verify=verify_ssl)
        response.raise_for_status()
    except requests.HTTPError as e:
        raise requests.HTTPError(f"{e} when trying to download {url}") from None
    else:
        if response.status_code != 200:
            raise requests.HTTPError(f"Unexpected error when trying to download {url}")

    return response


def fetch_parse(url: str, verify_ssl: bool = True) -> etree._Element:
    """
    Fetch the content of the URL and parse it into an HTML tree.

    Args:
    - url (str): The URL to fetch.
    - verify_ssl (bool): Whether to verify SSL certificates.

    Returns:
    - etree._Element: The parsed HTML tree.
    """
    response = http_get(url, verify_ssl)
    parser = etree.HTMLParser()
    tree = etree.fromstring(response.content, parser)

    return tree


def get_base_url(url: str) -> str:
    """
    Extract the base URL from the given URL.

    Args:
    - url (str): The full URL from which to extract the base URL.

    Returns:
    - str: The base URL.
    """
    parsed_url = urlparse(url)
    return f"{parsed_url.scheme}://{parsed_url.hostname}"


def extract_info(tree: etree._ElementTree, xpath: str) -> Optional[str]:
    """
    Helper function to extract information using a given XPath.

    Args:
    - tree (etree._ElementTree): The parsed HTML tree.
    - xpath (str): The XPath to extract information.

    Returns:
    - Optional[str]: Extracted information or None if not found.
    """
    info = tree.xpath(xpath)
    return info[0] if info else None


def download_file(
    url: str, path_prefix: str = "images", verify_ssl: bool = True
) -> str:
    """
    Download a file from a URL, save it to a StringIO buffer, compute its sha512 value and use it as the filename.

    Args:
    - url (str): The URL of the file to download.
    - path_prefix (str): The path prefix for the downloaded file.
    - verify_ssl (bool): Whether to verify SSL certificates.

    Returns:
    - str: The filename of the downloaded file.
    """
    response = http_get(url, verify_ssl)

    # Save the file to a StringIO buffer
    buffer = io.BytesIO(response.content)

    # Compute the sha512 value of the file
    sha512 = hashlib.sha512()
    sha512.update(buffer.getvalue())
    filename = sha512.hexdigest()

    # Get the extension name from the MIME type
    mime = response.headers.get("content-type")
    extension = mimetypes.guess_extension(mime)

    # Use the sha512 value as the filename and the extension from the MIME type
    filename_with_extension = f"{filename}{extension}"

    # Add the path prefix to the filename
    filename_with_path = (
        os.path.join(path_prefix, filename_with_extension)
        if path_prefix
        else filename_with_extension
    )

    # Create directories in the path prefix if they do not exist
    os.makedirs(os.path.dirname(filename_with_path), exist_ok=True)

    # Write the file to disk
    with open(filename_with_path, "wb") as f:
        f.write(buffer.getvalue())

    return filename_with_path
