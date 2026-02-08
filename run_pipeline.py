import argparse
import yaml
import os

from utils.logger import setup_logger
from integration.normalize_dlbc import normalize_dlbc
from integration.rank_histones import run_ranking


logger = setup_logger()


def load_config():

    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)


def run_autism_pipeline(config):

    logger.info("Running autism expression pipeline")

    expression_file = config["datasets"]["autism"]["expression_file"]

    if not os.path.exists(expression_file):

        logger.error(f"Expression file not found: {expression_file}")
        return False

    logger.info("Autism expression data verified")

    return True


def run_dlbc_pipeline(config):

    logger.info("Running DLBC genomic alteration pipeline")

    alteration_file = config["datasets"]["dlbc"]["alteration_file"]

    if not os.path.exists(alteration_file):

        logger.error(f"DLBC alteration file not found: {alteration_file}")
        return False

    logger.info("Normalizing DLBC alteration data")

    normalize_dlbc()

    logger.info("DLBC normalization complete")

    return True


def run_integration_pipeline(config):

    logger.info("Running multi-omics integration pipeline")

    run_ranking()

    logger.info("Integration and ranking complete")

    return True


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--disease",
        type=str,
        required=True,
        help="Disease to analyze (autism, dlbc, or all)",
    )

    args = parser.parse_args()

    config = load_config()

    logger.info("Pipeline started")

    if args.disease == "autism":

        run_autism_pipeline(config)

    elif args.disease == "dlbc":

        run_dlbc_pipeline(config)

    elif args.disease == "all":

        run_autism_pipeline(config)
        run_dlbc_pipeline(config)
        run_integration_pipeline(config)

    else:

        logger.error("Invalid disease option")

    logger.info("Pipeline finished")


if __name__ == "__main__":

    main()
