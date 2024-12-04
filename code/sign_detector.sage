from maass_sqfreelevel_computations import *
import sys


VERBOSE = True


def setup_level_data(level, Rmax = 26):
    gd = group_data(level)
    msd = MaassSpaceData(Rmax, gd)
    return gd, msd


def single_diff(cand, n):
    """
    Return the minimum distance from `cand` to +- 1/sqrt(n)
    """
    return min(abs(cand - 1/sqrt(n)), abs(cand + 1/sqrt(n)))


def total_diff(cands, ns):
    assert len(cands) == len(ns)
    return sum(single_diff(c, n) for (c, n) in zip(cands, ns))


def hecke_diff(coeffs):
    """
    Check a(2)*a(3) = a(6) and a(2)*a(5) = a(10) and a(3)*a(5) = a(15)
    """
    a2 = coeffs[int(2-2)]
    a3 = coeffs[int(3-2)]
    a5 = coeffs[int(5-2)]
    a6 = coeffs[int(6-2)]
    a10 = coeffs[int(10-2)]
    a15 = coeffs[int(15-2)]
    return abs(a2 * a3 - a6) + abs(a2 * a5 - a10) + abs(a3 * a5 - a15)


def find_correct_signs(level, gd, msd, Y, R, zdata, symmetry):
    checks = [d for d in divisors(level) if d!= 1 and d != level]
    checks = [p for p in checks if is_prime(p)]
    limit = sqrt(msd.Q)
    check_ns = []
    for p in checks:
        j = 1
        while p^j < limit:
            check_ns.append(p^j)
            j += 1
    record = 1e10
    record_signs = -1
    with open(f"{level}.{R}.verbose.data", "w", encoding="utf8") as vfile:
        for s in short_all_signs(level):
            coeffs = maass_form_coeffs(Y, R, zdata, s, symmetry=symmetry)
            check_cands = [coeffs[int(n - 2)] for n in check_ns]
            cand_diff = total_diff(check_cands, check_ns) + hecke_diff(coeffs)
            if VERBOSE: vfile.write(f"{s}: {cand_diff}: {coeffs[:20]}\n")
            print(s, cand_diff, coeffs[:8])
            if cand_diff < record:
                record = cand_diff
                record_signs = s
    return record, record_signs


def main(level):
    fname = f"missing.{level}.record"
    print(f"Reading missing forms from {fname}")
    with open(fname, "r", encoding="utf8") as mfile:
        lines = mfile.readlines()
    checklevel, maxR = list(map(int, lines[0].split(';')))
    assert level == checklevel
    lines = lines[1:]
    print(f"Found {len(lines)} missing forms to complete.")

    gd, msd = setup_level_data(level, maxR)
    Y = msd.Y
    zdata = msd.zdata

    outfname = f"guess.{level}.data"
    with open(outfname, "w", encoding="utf8") as outfile:
        for line in lines:
            _, label, R, sym = line.split(";")
            sym = sym.strip()
            # Annoying difference: 0/1 vs 1/-1 for even/odd.
            if sym == '0':
                sym = 1
            else:
                sym = -1
            R = float(R)
            record, signs = find_correct_signs(level, gd, msd, Y, R, zdata, sym)
            outfile.write(f"{label}; {signs[0]}; {record}\n")


if __name__ == "__main__":
    level = int(sys.argv[1])
    main(level)
