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

    :param url: The URL to send the GET request to.
    :type url: str
    :param verify_ssl: Whether to verify SSL certificates.
    :type verify_ssl: bool
    :return: The response from the GET request.
    :rtype: requests.Response
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

    :param url: The URL to fetch.
    :type url: str
    :param verify_ssl: Whether to verify SSL certificates.
    :type verify_ssl: bool
    :return: The parsed HTML tree.
    :rtype: etree._Element
    """
    response = http_get(url, verify_ssl)
    parser = etree.HTMLParser()
    tree = etree.fromstring(response.content, parser)

    return tree


def get_base_url(url: str) -> str:
    """
    Extract the base URL from the given URL.

    :param url: The full URL from which to extract the base URL.
    :type url: str
    :return: The base URL.
    :rtype: str
    """
    parsed_url = urlparse(url)
    return f"{parsed_url.scheme}://{parsed_url.hostname}"


def extract_info(tree: etree._ElementTree, xpath: str) -> Optional[str]:
    """
    Helper function to extract information using a given XPath.

    :param tree: The parsed HTML tree.
    :type tree: etree._ElementTree
    :param xpath: The XPath to extract information.
    :type xpath: str
    :return: Extracted information or None if not found.
    :rtype: Optional[str]
    """
    info = tree.xpath(xpath)
    return info[0] if info else None


def download_file(
    url: str, path_prefix: str = "images", verify_ssl: bool = True
) -> str:
    """
    Download a file from a URL, save it to a StringIO buffer, compute its sha512
    value and use it as the filename.

    :param url: The URL of the file to download.
    :type url: str
    :param path_prefix: The path prefix for the downloaded file.
    :type path_prefix: str
    :param verify_ssl: Whether to verify SSL certificates.
    :type verify_ssl: bool
    :return: The filename of the downloaded file.
    :rtype: str
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
