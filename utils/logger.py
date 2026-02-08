import logging
import os


def setup_logger(name="pipeline_logger"):

    os.makedirs("logs", exist_ok=True)

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    if not logger.handlers:

        file_handler = logging.FileHandler("logs/pipeline.log")
        console_handler = logging.StreamHandler()

        formatter = logging.Formatter(
            "%(asctime)s — %(levelname)s — %(message)s"
        )

        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

    return logger
