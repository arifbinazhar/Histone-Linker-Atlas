# import argparse
# import yaml
# import os

# from utils.logger import setup_logger
# from integration.normalize_dlbc import normalize_dlbc
# from integration.rank_histones import run_ranking


# logger = setup_logger()


# def load_config():

#     with open("config.yaml", "r") as f:
#         return yaml.safe_load(f)


# def run_autism_pipeline(config):

#     logger.info("Running autism expression pipeline")

#     expression_file = config.get("datasets", {}).get("autism", {}).get("expression_file")


#     if not os.path.exists(expression_file):

#         logger.error(f"Expression file not found: {expression_file}")
#         return False

#     logger.info("Autism expression data verified")

#     return True


# def run_dlbc_pipeline(config):

#     logger.info("Running DLBC genomic alteration pipeline")

#     alteration_file = config["datasets"]["dlbc"]["alteration_file"]

#     if not os.path.exists(alteration_file):

#         logger.error(f"DLBC alteration file not found: {alteration_file}")
#         return False

#     logger.info("Normalizing DLBC alteration data")

#     normalize_dlbc()

#     logger.info("DLBC normalization complete")

#     return True


# def run_integration_pipeline(config):

#     logger.info("Running multi-omics integration pipeline")

#     run_ranking()

#     logger.info("Integration and ranking complete")

#     return True


# def main():

#     parser = argparse.ArgumentParser()

#     parser.add_argument(
#         "--disease",
#         type=str,
#         required=True,
#         help="Disease to analyze (autism, dlbc, or all)",
#     )

#     args = parser.parse_args()

#     config = load_config()

#     logger.info("Pipeline started")

#     if args.disease == "autism":

#         run_autism_pipeline(config)

#     elif args.disease == "dlbc":

#         run_dlbc_pipeline(config)

#     elif args.disease == "all":

#         run_autism_pipeline(config)
#         run_dlbc_pipeline(config)
#         run_integration_pipeline(config)

#     else:

#         logger.error("Invalid disease option")

#     logger.info("Pipeline finished")


# if __name__ == "__main__":

#     main()


"""-------------------------------------------------------------------------------"""


#!/usr/bin/env python

import argparse
import yaml
import os

from utils.logger import setup_logger

# Import fetchers
from data_fetchers.geo_fetch import run_autism_pipeline
from data_fetchers.cbio_fetch import run_dlbc_pipeline

# Import processing modules
from integration.normalize_dlbc import normalize_dlbc
from integration.rank_histones import run_ranking


logger = setup_logger()


def load_config():

    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)


# -----------------------------
# GEO FETCH PIPELINE
# -----------------------------

def run_autism_pipeline_2(config):

    logger.info("Starting GEO autism data fetch pipeline")

    try:

        run_autism_pipeline()

        logger.info("GEO autism pipeline completed successfully")

    except Exception as e:

        logger.error(f"GEO pipeline failed: {e}")
        raise


# -----------------------------
# DLBC FETCH PIPELINE
# -----------------------------

def run_dlbc_fetch_pipeline(config):

    logger.info("Starting DLBC cBioPortal fetch pipeline")

    try:

        run_dlbc_pipeline()

        logger.info("DLBC fetch pipeline completed successfully")

    except Exception as e:

        logger.error(f"DLBC fetch failed: {e}")
        raise


# -----------------------------
# NORMALIZATION PIPELINE
# -----------------------------

def run_dlbc_normalization(config):

    logger.info("Starting DLBC normalization")

    try:

        normalize_dlbc()

        logger.info("DLBC normalization completed")

    except Exception as e:

        logger.error(f"DLBC normalization failed: {e}")
        raise


# -----------------------------
# INTEGRATION PIPELINE
# -----------------------------

def run_integration(config):

    logger.info("Starting multi-omics integration")

    try:

        run_ranking()

        logger.info("Integration and ranking completed")

    except Exception as e:

        logger.error(f"Integration failed: {e}")
        raise


# -----------------------------
# MASTER CONTROLLER
# -----------------------------

def main():

    parser = argparse.ArgumentParser(
        description="HistoneLinker-Atlas Pipeline"
    )

    parser.add_argument(
        "--disease",
        type=str,
        required=True,
        choices=["autism", "dlbc", "integration", "all"],
        help="Pipeline stage to run",
    )

    args = parser.parse_args()

    config = load_config()

    logger.info("Pipeline started")

    if args.disease == "autism":

        run_autism_pipeline_2(config)

    elif args.disease == "dlbc":

        run_dlbc_fetch_pipeline(config)
        run_dlbc_normalization(config)

    elif args.disease == "integration":

        run_integration(config)

    elif args.disease == "all":

        run_autism_pipeline_2(config)
        run_dlbc_fetch_pipeline(config)
        run_dlbc_normalization(config)
        run_integration(config)

    logger.info("Pipeline finished successfully")


if __name__ == "__main__":

    main()
