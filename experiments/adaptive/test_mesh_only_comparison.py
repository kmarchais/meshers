"""Benchmark controls preserve the measured specimens and element ordering."""

import json
import numpy as np
import pytest

import microgen_workflow as w
from mesh_only_comparison import configs,spin_field


@pytest.mark.parametrize('case,source',[
    ('spinodoid','g3-spinodoid-partial-settle-curved24-tol003-solve-v1'),
    ('smooth-lattice','g3-lattice-adapt26-tol0028-preflight-v1')])
def test_known_curved_recipe_matches_frozen_measurement(case,source):
    expected=json.loads((w.HERE/'runs'/source/'result.json').read_text())
    actual=next(c for c in configs(case,True) if c['name']=='meshers-curved-known')
    for current,old in [('cells','cells'),('order','order'),('adapt','adaptive'),
        ('tolerance','native_tolerance'),('adapt_tolerance','adapt_tolerance'),('min_length','min_length')]:
        assert actual[current]==expected[old]


def test_compiled_spinodoid_preserves_frozen_field():
    w.imports()
    import meshers
    case,shape,_=w.specimen('spinodoid')
    compiled=spin_field(shape,meshers,np)
    points=np.random.default_rng(20261018).uniform(-.5,.5,(4096,3))
    np.testing.assert_allclose(compiled.values(points),case.function(*points.T),atol=1e-12,rtol=1e-12)


def test_comparison_spectral_expression_matches_original_compiler_inputs():
    import comparison_fixtures
    import meshers
    case,shape,_=w.specimen('spinodoid')
    original=spin_field(shape,meshers,np)
    current=meshers.compile_field(comparison_fixtures.spectral_field(
        comparison_fixtures.coefficients_path('spinodoid',w.HERE)))
    points=np.random.default_rng(20261026).uniform(-.5,.5,(4096,3))
    np.testing.assert_array_equal(current.values(points),original.values(points))


@pytest.mark.parametrize('case,source',[
    ('spinodoid','g3-spinodoid-compiled32-v1'),
    ('smooth-lattice','g3-lattice-void-clamp32-v1')])
def test_historical_linear_route_is_retained(case,source):
    expected=json.loads((w.HERE/'runs'/source/'result.json').read_text())
    actual=next(c for c in configs(case,True) if c['name']=='meshers-tet4-historical32')
    assert actual['cells']==expected['cells']
    assert actual['min_length']==expected['min_length']
    assert actual['tolerance']==expected.get('native_tolerance',.002)
    assert actual['order']==1 and actual['adapt'] and actual['adapt_tolerance']==.002
    if case=='smooth-lattice':
        assert actual['lattice_field']=='original'
