import argparse
import csv
import sys
import time

from urllib.parse import urljoin

from scraper import fetch_parse, download_file, get_base_url

XPATH_IMAGE = '//div[contains(@class, "content_yellow")]/div/div[@id="picture_description"]/div[contains(@class, "mini_picture")]/a/@href'


def main():
    parser = argparse.ArgumentParser(
        description="Web Scraper for in situ images table."
    )
    # TODO: remove argparse.FileType("r") and argparse.FileType("w") and use str
    # as they are not safe.
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
    parser.add_argument(
        "--prefix",
        "-p",
        type=str,
        default="images",
        help="Prefix for the images path. Default is 'images'.",
    )

    args = parser.parse_args()

    with args.input as input_file, args.output as output:
        writer = csv.writer(output, delimiter="\t")
        # Writing the header
        writer.writerow(["URL", "Image_URL", "Image"])

        for line in input_file:
            url = line.strip()

            tree = fetch_parse(url, not args.insecure)

            for i in tree.xpath(XPATH_IMAGE):
                url_image = urljoin(get_base_url(url), i)
                file_image = download_file(url_image, args.prefix, not args.insecure)
                writer.writerow(
                    [
                        url,
                        url_image,
                        file_image,
                    ]
                )

            time.sleep(args.sleep)


if __name__ == "__main__":
    main()
