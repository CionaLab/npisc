from typing import List, Optional
from urllib.parse import urlparse

import requests
from lxml import etree


def fetch_parse(url: str, verify_ssl: bool = True) -> etree._Element:
    """
    Fetch the content of the URL and parse it into an HTML tree.

    Args:
    - url (str): The URL to fetch.
    - verify_ssl (bool): Whether to verify SSL certificates.

    Returns:
    - etree._Element: The parsed HTML tree.
    """
    response = requests.get(url, verify=verify_ssl)
    if response.status_code != 200:
        raise Exception(
            f"Failed to fetch URL. HTTP Status Code: {response.status_code}"
        )

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
