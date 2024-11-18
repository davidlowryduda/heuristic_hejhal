from sage.misc.cython import cython

FNAMES = [
    "lpkbessel.spyx",
    "maass_levelone_computations.spyx",
    "maass_levelone_computations_tests.spyx",
    "maass_sqfreelevel_computations.spyx",
    "maass_sqfreelevel_tests.spyx",
]

def recompile_all():
    for fname in FNAMES:
        cython(fname, compile_message=True, use_cache=True, create_local_so_file=True)


if __name__ == "__main__":
    recompile_all()
