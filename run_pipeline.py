import yaml
import logging
from fetcher.geo_fetch import run_autism_pipeline

logging.basicConfig(level=logging.INFO)

def main():
    with open("config.yaml") as f:
        config = yaml.safe_load(f)

    if "autism" in config["diseases"]:
        logging.info("Running Autism GEO pipeline...")
        run_autism_pipeline()

if __name__ == "__main__":
    main()
