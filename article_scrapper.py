import argparse
import sys
import time
from typing import List, Tuple, Optional, Generator
from urllib.parse import urljoin

from scraper import get_base_url, fetch_parse


def scrape_page(url: str, verify_ssl: bool = True) -> Tuple[List[str], Optional[str]]:
    """
    Scrape the provided webpage and retrieve article links.

    Args:
    - url (str): The URL of the webpage to scrape.
    - verify_ssl (bool): Whether to verify SSL certificates.

    Returns:
    - Tuple[List[str], Optional[str]]: A tuple containing a list of article
    links and the URL of the next page (if available).
    """
    tree = fetch_parse(url, verify_ssl)

    base_url = get_base_url(url)

    # Extract links inside articles and convert them to absolute URLs if they're relative
    article_links = tree.xpath('//*[@id="informations"]/article/header/p/a/@href')
    article_links = [urljoin(base_url, link) for link in article_links]

    # Get the next page link and convert it to an absolute URL if it's relative
    next_page_links = tree.xpath(
        '//*[@id="informations"]/ul[@class="pagelinks"]/li[@class="pagelinks-next"]/a/@href'
    )
    next_page_url = urljoin(base_url, next_page_links[0]) if next_page_links else None

    return article_links, next_page_url


def scrape_website(
    start_url: str, sleep_duration: float, verify_ssl: bool = True
) -> Generator[str, None, None]:
    """
    Generator that yields article links from the website starting from the provided URL,
    navigating through paginated content.

    Args:
    - start_url (str): The starting URL for the scraping process.
    - sleep_duration (float): Duration to wait between requests in seconds.
    - verify_ssl (bool): Whether to verify SSL certificates.

    Yields:
    - str: Article link found on the page.
    """
    current_url = start_url

    while current_url:
        article_links, next_page_url = scrape_page(current_url, verify_ssl)
        for link in article_links:
            yield link
        current_url = next_page_url
        time.sleep(sleep_duration)


def main():
    parser = argparse.ArgumentParser(description="Web Scraper for paginated content")
    parser.add_argument(
        "start_url", type=str, help="The starting URL for the scraping process"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file where the results should be written. Default is stdout.",
    )
    parser.add_argument(
        "--insecure",
        "-n",
        action="store_true",
        help="Ignore SSL certificate verification.",
    )
    parser.add_argument(
        "--sleep",
        "-s",
        type=float,
        default=0.5,
        help="Duration to wait between requests in seconds. Default is 0.5 seconds.",
    )

    args = parser.parse_args()

    with args.output as output:
        for link in scrape_website(args.start_url, args.sleep, not args.insecure):
            output.write(f"{link}\n")


if __name__ == "__main__":
    main()
