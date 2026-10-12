import random
import shutil

import pytest

from tests.synthetic import make_reference_dir


@pytest.fixture(scope="session")
def refs(tmp_path_factory):
    if shutil.which("blastn") is None:
        pytest.skip("BLAST+ not installed")
    d = tmp_path_factory.mktemp("references")
    genes = make_reference_dir(d)
    return d, genes


@pytest.fixture(scope="session")
def engine(refs):
    from septatyper.prediction_engine import SeptaTyperEngine
    return SeptaTyperEngine(refs[0], threads=2)


@pytest.fixture()
def rng():
    return random.Random(11)
