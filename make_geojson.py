import argparse
import json
import sys
from typing import Dict


def generate_grid(rows: int, cols: int) -> Dict:
    """
    Generate a grid in geoJSON format.

    Args:
    - rows (int): Number of rows for the grid.
    - cols (int): Number of columns for the grid.

    Returns:
    - Dict: A dictionary representing the geoJSON format of the grid.
    """

    features = []

    for i in range(rows):
        for j in range(cols):
            cell = {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [[j, i], [j, i + 1], [j + 1, i + 1], [j + 1, i], [j, i]]
                    ],
                },
                "properties": {"name": f"({i}, {j})"},
            }
            features.append(cell)

    grid = {"type": "FeatureCollection", "features": features}

    return grid


def main():
    parser = argparse.ArgumentParser(
        description="Generate a geoJSON grid based on provided parameters."
    )
    parser.add_argument(
        "--rows", "-r", type=int, default=1, help="Number of rows. Default is 1."
    )
    parser.add_argument(
        "--cols", "-c", type=int, default=1, help="Number of columns. Default is 1."
    )
    parser.add_argument(
        "--output",
        "-o",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file to write the geoJSON grid. Default is stdout.",
    )

    args = parser.parse_args()

    grid_geojson = generate_grid(args.rows, args.cols)
    json.dump(grid_geojson, args.output, indent=4)


if __name__ == "__main__":
    main()
