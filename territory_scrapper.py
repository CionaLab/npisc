import argparse
import csv
import sys
import time
from typing import List, Tuple, Optional

import requests
from lxml import etree


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


def get_expression_info(
    url: str, verify_ssl: bool = True
) -> Tuple[Optional[str], Optional[str], List[str]]:
    """
    Retrieve the "Stage", "Expression Territory", and "Predicted Gene" from the specified URL.

    Args:
    - url (str): The URL containing the desired information.
    - verify_ssl (bool): Whether to verify SSL certificates.

    Returns:
    - Tuple[Optional[str], Optional[str], List[str]]: A tuple containing the stage, the predicted gene,
      and a list of territories  or None if not reported.
    """
    response = requests.get(url, verify=verify_ssl)
    if response.status_code != 200:
        raise Exception(
            f"Failed to fetch URL. HTTP Status Code: {response.status_code}"
        )

    parser = etree.HTMLParser()
    tree = etree.fromstring(response.content, parser)

    # Extracting the information using the refactored function
    stage = extract_info(tree, '//*[@id="informations"]/div[3]/div/div[2]/div[2]/p/text()')
    territories = tree.xpath(
        '//*[@id="informations"]/div[3]/div/div[3]/div/div/table/tr/td[1]/a/text()'
    )
    predicted_gene = extract_info(
        tree, '//*[@id="informations"]/div[3]/div/div[2]/div[3]/p/a/text()'
    )

    return stage, predicted_gene, territories


def main():
    parser = argparse.ArgumentParser(
        description="Web Scraper for Expression Territory and Predicted Gene information"
    )
    parser.add_argument(
        "--input",
        "-i",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input file containing URLs to scrape. Default is stdin.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file to write the results in TSV format. Default is stdout.",
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

    with args.input as input_file, args.output as output:
        writer = csv.writer(output, delimiter="\t")
        # Writing the header
        writer.writerow(["URL", "Stage", "Gene", "Territory"])

        for line in input_file:
            url = line.strip()
            stage, predicted_gene, territories = get_expression_info(
                url, not args.insecure
            )

            for territory in territories:
                writer.writerow(
                    [
                        url,
                        stage if stage else "None",
                        predicted_gene if predicted_gene else "None",
                        territory,
                    ]
                )

            time.sleep(args.sleep)


if __name__ == "__main__":
    main()
