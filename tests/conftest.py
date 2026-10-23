import pytest

from wellqc.pipeline import run_audit
from wellqc.synthetic import make_demo_dataset


@pytest.fixture(scope="session")
def demo_projects(tmp_path_factory):
    return make_demo_dataset(tmp_path_factory.mktemp("demo_data"), seed=7)


@pytest.fixture(scope="session")
def audit(demo_projects, tmp_path_factory):
    return run_audit(demo_projects, tmp_path_factory.mktemp("out"), figures=False)


def messages(audit, category=None, well=None):
    return [i.message for i in audit.issues
            if (category is None or i.category == category) and (well is None or i.well == well)]
