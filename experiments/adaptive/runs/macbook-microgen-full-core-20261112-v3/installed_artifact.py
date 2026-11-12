from pathlib import Path
import sys
import microgen

def pytest_collection_modifyitems(items):
    root = Path(microgen.__file__).parent
    assert root.is_relative_to(Path(sys.prefix))
    assert len(list(root.rglob('*.py'))) == 51
    for item in items:
        if item.module.__name__.endswith('test_no_top_level_ocp_imports'):
            item.module.MICROGEN_ROOT = root
