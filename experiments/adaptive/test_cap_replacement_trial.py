"""Never apply overlapping plans or silently overwrite unrelated cells."""

import numpy as np
import pytest

from cap_replacement_trial import replace_plans


def test_replacement_preserves_unselected_cell_and_boundary():
    cells = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
    surface = np.array([[0, 1, 2], [1, 2, 4]])
    plan = {"incident_elements": [0], "alternative_cells": [[0, 1, 2, 5]],
        "old_cap_faces": [[0, 1, 2]], "new_cap_faces": [[0, 1, 5]]}
    replaced, boundary = replace_plans(cells, surface, [plan])
    np.testing.assert_array_equal(replaced, [[1, 2, 3, 4], [0, 1, 2, 5]])
    np.testing.assert_array_equal(boundary, [[1, 2, 4], [0, 1, 5]])
    with pytest.raises(ValueError, match="overlap"):
        replace_plans(cells, surface, [plan, plan])


def test_wrong_source_face_is_rejected():
    plan = {"incident_elements": [0], "alternative_cells": [[0, 1, 2, 3]],
        "old_cap_faces": [[0, 1, 4]], "new_cap_faces": [[0, 1, 2]]}
    with pytest.raises(ValueError, match="match source"):
        replace_plans(np.array([[0, 1, 2, 3]]), np.array([[0, 1, 2]]), [plan])
