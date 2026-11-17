"""Replay the extended frozen sample using the original catalog selection."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'repeat-catalogs'))
from compare import compare, main, project
from fetch_catalogs import product_url


if __name__ == '__main__':
    main(HERE)
