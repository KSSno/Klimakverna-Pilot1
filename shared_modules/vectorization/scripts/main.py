import argparse
import json
import logging
import sys
from logging.handlers import TimedRotatingFileHandler
from pathlib import Path

from .config import DEBUG_FOLDER, LOG_DAYS, LOG_FOLDER, LOG_LEVEL, DefultContourConfig
from .models import ContourGeneratorConfig
from .utils import (
    contour_generator,
    get_contour_definition,
    get_feature_layer,
    list_feature_layers,
)

logger = logging.getLogger()
logger.setLevel(LOG_LEVEL)
logger.handlers = []

formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# Set up timed rotating file handler
if LOG_FOLDER:
    LOG_FOLDER.mkdir(parents=True, exist_ok=True)
    log_file = LOG_FOLDER / "klimakart.log"
    file_handler = TimedRotatingFileHandler(
        log_file, when="midnight", interval=1, backupCount=LOG_DAYS
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)


def main():
    # Use argparse to parse command line arguments

    parser = argparse.ArgumentParser(
        description="Generate contour lines from raster datasets."
    )
    # Load dataset configuration from JSON file
    parser.add_argument(
        "-d",
        "--datasets",
        type=Path,
        default=None,
        help="Path to the datasets configuration JSON file.",
    )
    # Load contour levels from JSON file
    parser.add_argument(
        "-c",
        "--contour_levels",
        type=Path,
        default=None,
        help="Path to the contour levels configuration JSON file.",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode for detailed logging.",
    )

    # Add option to list available feature layers
    parser.add_argument(
        "--list-layers",
        action="store_true",
        help="List all available feature layers in the package.",
    )

    # Add option to list available contour configurations
    parser.add_argument(
        "--list-contours",
        action="store_true",
        help="List all available default contour configurations.",
    )

    # Add option to generate default datasets configuration
    parser.add_argument(
        "--default-datasets-config",
        nargs="?",
        const="datasets.json",
        metavar="FILE",
        help="Generate default datasets configuration to a JSON file (default: datasets.json)",
    )

    # Add option to generate default contour levels configuration
    parser.add_argument(
        "--default-contours-config",
        nargs="?",
        const="contour_levels.json",
        metavar="FILE",
        help="Generate default contour levels configuration to a JSON file (default: contour_levels.json)",
    )

    # Get the arguments
    args = parser.parse_args()

    def fix_datasets(datasets):
        """
        Fix the datasets by converting Path objects to strings.
        """
        return [
            {
                key: value.as_posix() if isinstance(value, Path) else value
                for key, value in dataset.items()
            }
            for dataset in datasets
        ]

    # Handle the generate options before regular processing
    generate_files = False

    if args.default_datasets_config:
        output_file = Path(args.default_datasets_config)
        logger.info(f"Generating default datasets configuration to {output_file}")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(fix_datasets(DefultContourConfig.datasets), f, indent=4)
        print(f"Default datasets configuration written to {output_file}")
        generate_files = True

    if args.default_contours_config:
        output_file = Path(args.default_contours_config)
        logger.info(f"Generating default contour levels configuration to {output_file}")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(DefultContourConfig.contour_levels, f, indent=4)
        print(f"Default contour levels configuration written to {output_file}")
        generate_files = True

    # Exit if files were generated
    if generate_files:
        logger.info("Exiting after generating configuration files")
        sys.exit(0)

    # Handle the list options before regular processing
    if args.list_layers:
        layers = list_feature_layers()
        print("Available feature layers:")
        for layer in layers:
            print(f"  - {layer}")
        sys.exit(0)  # Exit after listing layers

    # Handle the list-contours option
    if args.list_contours:
        contour_configs = DefultContourConfig.contour_levels
        print("Available contour configurations:")
        for name, details in contour_configs.items():
            print(f"  - {name}")
            print(f"    Levels: {details['levels']}")
            print(f"    Colors: {details['colormap']}")
            print()  # Add a blank line for readability
        sys.exit(0)  # Exit after listing contours

    # Load the datasets and contour levels from the provided JSON files
    if args.datasets:
        with args.datasets.open("r", encoding="utf-8") as f:
            datasets = json.load(f)
            logger.info(f"Loaded datasets from {args.datasets}")
    else:
        datasets = DefultContourConfig.datasets
        logger.info("Using default datasets configuration")

    if args.contour_levels:
        with args.contour_levels.open("r", encoding="utf-8") as f:
            contour_levels = json.load(f)
            logger.info(f"Loaded contour levels from {args.contour_levels}")
    else:
        contour_levels = DefultContourConfig.contour_levels
        logger.info("Using default contour levels configuration")

    debug = args.debug
    if args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug mode enabled")
        if DEBUG_FOLDER:
            DEBUG_FOLDER.mkdir(parents=True, exist_ok=True)
            with open(DEBUG_FOLDER / ".gitignore", "w") as f:
                f.write("*\n")

    logger.info("Starting contour generation process")

    for i, ds in enumerate(datasets):
        logger.info(f"Processing contour configuration {i + 1}/{len(datasets)}")

        if isinstance(ds["clip_layer"], str):
            # If clip_layer is a string, get the definition from the dictionary
            ds["clip_layer"] = get_feature_layer(ds["clip_layer"])

        if isinstance(ds["contour_levels"], str):
            # If contour_levels is a string, get the definition from the dictionary
            ds["contour_levels"] = get_contour_definition(
                ds["contour_levels"], contour_dict=contour_levels
            )

        if "nodata_feature_layer" in ds:
            if isinstance(ds["nodata_feature_layer"], str):
                # If nodata_feature_layer is a string, get the definition from the dictionary
                ds["nodata_feature_layer"] = get_feature_layer(
                    ds["nodata_feature_layer"]
                )

        contour_generator(config=ContourGeneratorConfig(**ds), debug=debug)

    logger.info("Processing complete!")


if __name__ == "__main__":
    main()
