import argparse
import csv
import sys
import time
from typing import List, Tuple, Optional

from scraper import extract_info, fetch_parse


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
    tree = fetch_parse(url, verify_ssl)

    # Extracting the information using the refactored function
    stage = extract_info(
        tree,
        '//section[@id="informations"]/div[@class="content_title"]/div[@class="content_yellow"]/div[@class="content_frame mod"]/div[2]/p/text()',
    )
    territories = (
        t
        if (
            t := tree.xpath(
                '//section[@id="informations"]/div[@class="content_title"]/div[@class="content_yellow"]/div[@class="table"]/div[@class="results_yellow"]//div[@class="results_content"]/table/tr/td[1]/a/text()',
            )
        )
        else ["None"]
    )
    predicted_gene = extract_info(
        tree,
        '//section[@id="informations"]/div[@class="content_title"]/div[@class="content_yellow"]/div[@class="content_frame mod"]/div[3]/p/a/text()',
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
