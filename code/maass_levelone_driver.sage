from maass_levelone_computations import find_single_ev_linearized, find_evs
import sys
import time

def main():
    if len(sys.argv) < 4:
        print("  Usage:   sage progname R radius symtype")
        print("  example: sage progname 9.5 0.5 -1")
        sys.exit()
    R = float(sys.argv[1])
    radius = float(sys.argv[2])
    symtype = int(sys.argv[3])
    print("Searching for an eigenvalue in B({}, {}) with symtype {}".format(R, 2*radius, symtype))
    val = find_single_ev_linearized(R, radius, symtype, verbosity=2)
    if val != None:
        print("{} is an eigenvalue.".format(val[0]))
    else:
        print("None found.")


def other_main():
    if len(sys.argv) < 4:
        print("  Usage:   sage progname R1 R2 symtype")
        print("  example: sage progname 9 11 -1")
        sys.exit()
    R1 = float(sys.argv[1])
    R2 = float(sys.argv[2])
    symtype = int(sys.argv[3])

    start = time.time()
    evs = find_evs(R1, R2, symmetry=symtype, verbosity=2)
    print("Finished in {} seconds".format(time.time() - start))
    if evs:
        prev = evs[0]
        print('----- ev -----')
        print(prev[0])
        for val in prev[1]:
            print("  {}".format(val))
        print()

        for ev in evs[1:]:
            if abs(ev[0] - prev[0]) < 1e-4:
                continue
            print('----- ev -----')
            print(ev[0])
            for val in ev[1]:
                print("  {}".format(val))
            print()
            prev = ev

if __name__ == "__main__":
    #other_main()
    main()
